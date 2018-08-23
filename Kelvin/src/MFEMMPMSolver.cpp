/**----------------------------------------------------------------------------
 Copyright (c) 2018-, UT-Battelle, LLC
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
#include <MFEMMPMSolver.h>
#include <set>
#include <MassMatrix.h>

using namespace std;
using namespace mfem;

namespace Kelvin {

MFEMMPMSolver::MFEMMPMSolver() {
	// TODO Auto-generated constructor stub

}

MFEMMPMSolver::~MFEMMPMSolver() {
	// TODO Auto-generated destructor stub
}

void MFEMMPMSolver::solve(MFEMMPMData & data) {

	// Set of all nodes currently near particles
	set<int> nodeSet;
	set<int>::iterator outerIt;
	set<int>::iterator innerIt;

	// Compute initial mapping matrices
	auto & meshCont = data.meshContainer();
	auto & mesh = meshCont.getMesh();
	int numElems = mesh.GetNE();
	int numVerts = mesh.GetNV();
	auto & particles = data.particles();
	int numParticles = particles.size();
	// The shape matrix is very sparse, so this computation exploits that by only
	// adding shape values for nodes that exist for the given particle.
	SparseMatrix shapeMatrix(numParticles,numVerts);
	for (int i = 0; i < numParticles; i++) {
		auto nodeIds = meshCont.getSurroundingNodeIds(particles[i].coords);
		auto shape = meshCont.getNodalShapes(particles[i].coords);
		for (int j = 0; j < shape.size(); j++) {
			shapeMatrix.Add(i,nodeIds[j],shape[j]);
			nodeSet.insert(nodeIds[j]);
		}
	}
	shapeMatrix.Finalize();
	shapeMatrix.SortColumnIndices();

	shapeMatrix.Print();

	cout << shapeMatrix.Size() << endl;

	// Construct the mass matrix associated with the grid nodes
	double particleMass = 1.0; // FIXME! Needs to be something real and from input.
	MassMatrix massMatrix;
	massMatrix.setParticles(particles);
	massMatrix.assemble(shapeMatrix,nodeSet);
	// Get the diagonalized form of the mass matrix
	auto diagonalMassMatrix = massMatrix.lump();

	// Compute the acceleration at the grid nodes


	// Compute the velocity through explicit integration

	// Compute grad(v) at the material points

	// Use the velocity gradient to compute the strain rate at material points

	// Compute/update the stress at material points using the constitutive
	// equation

	// Use mapping functions to compute the velocity and acceleration at the
	// material points

	// Compute updates to the material point positions and velocity using
	// explicit integration

	// Update mapping matrices

	// Update gradients

	// Update the consistent mass matrix on the grid

	// Update the diagonal mass matrix

	// Use conservation of momentum to solve for the new grid velocities

	// Loop


}


} /* namespace Kelvin */
