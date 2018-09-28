
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
#include <BasicMFEMGridMapper.h>

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
	// Setup a mapper for mapping from the grid to the material points
	BasicMFEMGridMapper mapper(data.meshContainer().getMesh());
	// Create a storage vector for velocities from the mapper
	int dim = data.meshContainer().dimension();
	int numParticles = particles.size();
	std::vector<double> velUpdate(numParticles*dim);

	// Set basic start time parameters - FIXME! Will read from input
	double tInit = 0.0, tFinal = 10.0, dt = 1.0; // dtOverstep = tFinal % dt;
	int numTimeSteps = (int) tFinal/dt, printStepFrequency = 5;

	// FIXME! time stepping issues
	// 1) Handle final time overstepping
	// 2) Follow time step limitations described in the mpm-libre article

	// Integrate over time. At the moment this will not integrate exactly to
	// tFinal. See time stepping issues above - just a place holder for now.
	for (int ts = 0; ts < numTimeSteps; ts++) {

		// Compute the acceleration at the grid nodes
		grid.updateNodalAccelerations(dt, particles);
		// Compute the initial velocity from the momenta
		grid.updateNodalVelocitiesFromMomenta(particles);
		// Compute the velocity update
		grid.updateNodalVelocities(dt, particles);

		// USE A CONSTITUTIVE EQUATION CLASS INSTEAD! PASS GRID AND PARTICLES
		// Compute grad(v) at the material points - Do I need this?
		// Use the velocity gradient to compute the strain rate at material
		// points
		// Compute/update the stress at material points using the constitutive
		// equation

		// Use mapping functions to compute the velocity and acceleration at
		// the material points
		mapper.updateParticleAccelerations(grid, particles);
		mapper.updateParticleVelocities(grid, particles, velUpdate);

		// Compute updates to the material point positions and velocity using
		// explicit integration. This is just a simple explicit Euler update.
		for (int i = 0; i < numParticles; i++) {
			auto & mPoint = particles[i];
			for (int j = 0; j < dim; j++) {
				mPoint.pos[j] += dt * velUpdate[i * dim + j];
				mPoint.vel[j] += dt * mPoint.acc[j];
			}
		}

		// Print stepping information
		if (!(ts % printStepFrequency)) {
			cout << "dt = " << dt << ", ts = " << ts << endl;
		}

	}

	// Loop


}


} /* namespace Kelvin */
