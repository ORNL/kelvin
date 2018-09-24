/**----------------------------------------------------------------------------
 Copyright  2018-, UT-Battelle, LLC
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 * Neither the name of fern nor the names of its
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 Author(s): Jay Jay Billings (billingsjj <at> ornl <dot> gov)
 -----------------------------------------------------------------------------*/
#include <Grid.h>
#include <memory>

using namespace mfem;
using namespace std;

namespace Kelvin {

Grid::Grid(MeshContainer & meshContainer) : _meshContainer(meshContainer){
	// TODO Auto-generated constructor stub

}

Grid::~Grid() {
	// TODO Auto-generated destructor stub
}

void Grid::assemble(const std::vector<Kelvin::MaterialPoint> & particles) {

	// Get the nodal coordinates from the mesh
	auto & mesh = _meshContainer.getMesh();
	int numVerts = mesh.GetNV();
	double * mcCoords;
	int dim = _meshContainer.dimension();
	for (int i = 0; i < numVerts; i++) {
		mcCoords = mesh.GetVertex(i);
		Point point;
		// Create the node. Only the position needs to be assigned.
		for (int j = 0; j < dim; j++) {
			point.pos[j] = mcCoords[j];
		}
		_nodes.push_back(point);
	}

	// The shape matrix is very sparse, so this computation exploits that by only
	// adding shape values for nodes that exist for the given particle.
	int numParticles = particles.size();
	_shapeMatrix = make_unique<SparseMatrix>(numParticles,numVerts);
	for (int i = 0; i < numParticles; i++) {
		// Get the shape
		auto nodeIds = _meshContainer.getSurroundingNodeIds(particles[i].pos);
		auto shape = _meshContainer.getNodalShapes(particles[i].pos);
		for (int j = 0; j < shape.size(); j++) {
			_shapeMatrix->Add(i,nodeIds[j],shape[j]);
			_nodeSet.insert(nodeIds[j]);
		}
		// Get the gradients associated with the particle
		auto gradients = _meshContainer.getNodalGradients(particles[i].pos);
		_gradientMap[i] = gradients;
	}
	_shapeMatrix->Finalize();
	_shapeMatrix->SortColumnIndices();

	// Construct the mass matrix associated with the grid nodes
	_massMatrix = make_unique<MassMatrix>(particles);

	// Resize the internal and external force vectors
	ForceVector emptyForceWithCorrectDimensionality(dim);
	int numForces = _nodeSet.size();
	_internalForces.resize(numForces,emptyForceWithCorrectDimensionality);
	_externalForces.resize(numForces,emptyForceWithCorrectDimensionality);
	// Set the force node ids
	set<int>::iterator it;
	int index = 0;
	for (it = _nodeSet.begin(); it != _nodeSet.end(); it++) {
		_internalForces[index].nodeId = *it;
		_externalForces[index].nodeId = *it;
		index++;
	}

	return;
}

void Grid::updateShapeMatrix() {

}

void Grid::updateMassMatrix() {

}

void Grid::update() {

	// Get the diagonalized form of the mass matrix
	auto diagonalMassMatrix = _massMatrix->lump();

}

const std::map<int,std::vector<Gradient>> & Grid::gradients() const {
	return _gradientMap;
}

const std::vector<ForceVector> & Grid::internalForces(
		const std::vector<Kelvin::MaterialPoint> & particles) {

	/**
	 * This operation computes the internal forces according to Sulsky's 1994
	 * paper:
	 * f = -\sum_p M_p * G_ip^T * stress_p
	 */

	// Compute the forces
	int numForces = _internalForces.size();
	int dim = _meshContainer.dimension();
	int numParticles = particles.size();
	for (int i = 0; i < numForces; i++) {
		auto & forceVector = _internalForces[i];
		forceVector.clear();
		int forceNodeId = forceVector.nodeId;
		// Computer the component of the force due to each particle. The
		// gradient map contains the gradients for each particle for each
		// surrounding node, so we need to look to see if the node id
		// of any of those gradients matches the node id of force vector.
		for (int j = 0; j < numParticles; j++) {
			auto & mPoint = particles[j];
			auto & gradients = _gradientMap[j];
			auto numGrads = gradients.size();
			// Compute the component due to each element of the gradient...
			for (int k = 0; k < numGrads; k++) {
				auto & grad = gradients[k];
				int nodeId = grad.nodeId;
				// ... iff its node id matches
				if (nodeId == forceNodeId) {
					// Need to compute it over all dimensions
					for (int l = 0; l < dim; l++) {
						forceVector.values[l] -= grad.values[l]*
								mPoint.mass*mPoint.stress[l];
					}
				}
			}
		}
	}

	return _internalForces;
}

const std::vector<ForceVector> & Grid::externalForces(
		const std::vector<Kelvin::MaterialPoint> & particles) {

	/**
	 * This operation computes the internal forces according to Sulsky's 1994
	 * paper:
	 * f = \sum_p M_p * S_ip^T * bodyForces_p
	 */

	Vector rowI;
	Array<int> colsI;

	// Compute the forces
	int numForces = _externalForces.size();
	int dim = _meshContainer.dimension();
	int numParticles = particles.size();
	for (int i = 0; i < numForces; i++) {
		auto & forceVector = _externalForces[i];
		forceVector.clear();
		int forceNodeId = forceVector.nodeId;
		// Computer the component of the force due to each particle. The shape
		// matrix contains all the shapes and the nodeset describes which nodes
		// are massive.
		for (int j = 0; j < numParticles; j++) {
			auto & mPoint = particles[j];
			_shapeMatrix->GetRow(j, colsI, rowI);
			int colI = colsI.Find(i);
			if (colI >= 0) {
				// Need to compute it over all dimensions
				for (int l = 0; l < dim; l++) {
					forceVector.values[l] += rowI[colI]*
							mPoint.mass * mPoint.bodyForce[l];
				}
			}
		}
	}

	// TESTING GRID FUNCTIONS - NOW GETTING SLEEPY AND WILL MESS SOMETHING UP. GOING TO LAY DOWN.

	cout << "GRID FUNCTION TEST JUNK!" << endl;

	auto & space = _meshContainer.getSpace();
	GridFunction gf(&space);
	gf = 1.0;
	Vector nVals(dim);

	int numNodes = 6;
	gf.GetNodalValues(nVals);
	nVals.Print();

	// Find the element that contains the point
	mfem::IntegrationPoint intPoint;
	intPoint.Set2(0.5,0.5);
	std::vector<double> point = {0.5,0.5};

	cout << gf.GetValue(0,intPoint) << endl;
	auto shape = _meshContainer.getNodalShapes(point);
	for (int i = 0; i < shape.size(); i++) {
		cout << shape[i] << " ";
	}
	cout << endl;

	gf[1] = 2.0;
	gf[4] = 2.0;
	cout << gf.GetValue(0,intPoint) << endl;

	// FYI - IT WORKED AS EXPECTED.

	return _externalForces;
}

MassMatrix & Grid::massMatrix(
		const std::vector<Kelvin::MaterialPoint> & particles) {
	_massMatrix->assemble(*_shapeMatrix,_nodeSet);
	return *_massMatrix;
}

const std::vector<Point> & Grid::nodes() const {
	return _nodes;
}

void Grid::updateNodalAccelerations(const double & timeStep,
		const std::vector<Kelvin::MaterialPoint> & particles) {

	int dim = _meshContainer.dimension();
	// Compute the mass matrix (Sulsky steps 8 & 10a)
	auto & massMat = massMatrix(particles);
	// Lump the mass matrix (Sulsky step 10b)
	auto lumpedMassMat = massMat.lump();
	// Compute the forces (Sulsky step 9)
	auto & intForces = internalForces(particles);
	auto & exForces = externalForces(particles);
	// Making an assumption here that the node sets for the internal and
	// external forces match. That's reasonable given the implementation of the
	// base class, but may not be in the case of a subclass.

	// Only proceed if the sizes of these systems match. If the do not, this
	// is a critical error, so it should be checked regularly.
	bool sizesMatch = (intForces.size() == exForces.size())
			&& (lumpedMassMat.size() == intForces.size());
	if (sizesMatch) {
		int numNodes = lumpedMassMat.size();
		// a_i = (f^int_i + f^ex_i)/m_i
		for (int i = 0; i < numNodes; i++) {
			auto & intForceVec = intForces[i];
			auto & exForceVec = exForces[i];
			// Under the nodeset assumption above, just get the node id from
			// the internal force vector
			int nodeId = intForceVec.nodeId;
			for (int j = 0; j < dim; j++) {
				// Note that the mass is non-dimensional (thus i, not j)
				_nodes[nodeId].acc[j] = (intForceVec.values[j]
						+ exForceVec.values[j])/lumpedMassMat[i];
			}
		}
	} else {
		throw "Mass matrix diagonal and force vector size mismatch (unequal).";
	}
	// Compute the acceleration and update the grid (Sulsky step 1)
	// ALL GRID POINTS?

	return;
}

void Grid::updateNodalVelocities(const double & timeStep,
		const std::vector<Kelvin::MaterialPoint> & particles) {
	return;
}

} /* namespace Kelvin */
