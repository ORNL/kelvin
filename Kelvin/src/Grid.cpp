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

void Grid::assemble(const std::vector<Kelvin::Point> & particles) {

	// Get the nodal coordinates from he mesh
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
			nodeSet.insert(nodeIds[j]);
		}
		// Get the gradients associated with the particle
		auto gradients = _meshContainer.getNodalGradients(particles[i].pos);
		_gradientMap[i] = gradients;
	}
	_shapeMatrix->Finalize();
	_shapeMatrix->SortColumnIndices();

	// Create the gradient matrix. Use a std::map<int,std::vector<Gradient>> to
	// store the gradient set. Use a std::map<int,std::map<int,int>> to handle
	// a fast-indexed node set. Make sure all map accesses are by ref!


	// Construct the mass matrix associated with the grid nodes
	_massMatrix = make_unique<MassMatrix>(particles);
	double particleMass = 1.0; // FIXME! Needs to be something real and from input.
	_massMatrix->assemble(*_shapeMatrix,nodeSet);

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

const MassMatrix & Grid::massMatrix() const {
	return *_massMatrix;
}

const std::vector<Point> Grid::nodes() const {
	return _nodes;
}

} /* namespace Kelvin */
