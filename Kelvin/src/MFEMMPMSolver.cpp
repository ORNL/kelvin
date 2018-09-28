
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

 * Neither the name of Kelvin nor the names of its
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

	// Assemble the grid
	auto & particles = data.particles();
	auto & grid = data.grid();
	grid.assemble(particles);

	// Set basic start time parameters - FIXME! Will read from input
	double tFinal = 10.0, dt = 1.0;

	// Compute the acceleration at the grid nodes
	grid.updateNodalAccelerations(dt,particles);
	// Compute the initial velocity from the momenta
	grid.updateNodalVelocitiesFromMomenta(particles);
	// Compute the velocity update
	grid.updateNodalVelocities(dt,particles);

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
