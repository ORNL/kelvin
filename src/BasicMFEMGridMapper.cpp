/**----------------------------------------------------------------------------
 Copyright (c) 2017-, UT-Battelle, LLC
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 * Neither the name of the copyright holder nor the names of its
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
#include <BasicMFEMGridMapper.h>
#include <Grid.h>
#include <cmath>

using namespace Kelvin;
using namespace std;
using namespace mfem;

namespace Kelvin {

BasicMFEMGridMapper::BasicMFEMGridMapper(mfem::Mesh & mesh) : _mesh(mesh) {
	// TODO Auto-generated constructor stub

}

BasicMFEMGridMapper::~BasicMFEMGridMapper() {
	// TODO Auto-generated destructor stub
}

void BasicMFEMGridMapper::updateParticleAccelerations(const Kelvin::Grid & grid,
		std::vector<Kelvin::MaterialPoint> & particles) const {

	// Create the H1 field
	int dim = _mesh.Dimension();
	H1_FECollection accCol(1,dim);
	FiniteElementSpace accSpace(&_mesh,&accCol,dim,Ordering::byVDIM);
	// Create and fill the grid function
	GridFunction accGf(&accSpace);
	auto & nodes = grid.nodes();
	for (int i = 0; i < nodes.size(); i++) {
		auto & nodeAcc = nodes[i].acc;
		for (int j = 0; j < dim; j++) {
			accGf[i*dim+j] = nodeAcc[j];
			cout << i << " " << nodeAcc[j] << " ";
		}
		cout << endl;
	}

    // Find the point quickly since the background mesh is known to be a cube.
    int numNodes = _mesh.GetNV();
    int nodesPerSide = (dim == 2) ? sqrt(numNodes) - 1 : cbrt(numNodes) - 1;
    double * node1 = _mesh.GetVertex(0);
    double * node2 = _mesh.GetVertex(1);
    double sideLength = node2[0] - node1[0];
    int id = 0;

	// Map the grid accelerations to the particles
	mfem::IntegrationPoint intPoint;
	Vector acc(dim);
	for (int i = 0; i < particles.size(); i++) {
		auto & mPoint = particles[i];

		// Find the point - FIXME! Put it on the Point class
		if (dim == 2) {
			id = ((int) (mPoint.pos[0] / sideLength))
					+ nodesPerSide * ((int) (mPoint.pos[1] / sideLength));
		} else if (dim == 3) {
			id = ((int) (mPoint.pos[0] / sideLength))
					+ nodesPerSide * ((int) (mPoint.pos[1] / sideLength))
					+ nodesPerSide * nodesPerSide
							* ((int) (mPoint.pos[2] / sideLength));
		}

		intPoint.Set(mPoint.pos.data(),dim);
		accGf.GetVectorValue(id,intPoint,acc);
		// Write the acceleration data back to the point
		for (int j = 0; j < dim; j++) {
			mPoint.acc[j] = acc[j];
		}
	}

	return;
}

void BasicMFEMGridMapper::updateParticleVelocities(const Kelvin::Grid & grid,
		std::vector<Kelvin::MaterialPoint> & particles) const {

	// Create the H1 field
	int dim = _mesh.Dimension();
	H1_FECollection velCol(1,dim);
	FiniteElementSpace velSpace(&_mesh,&velCol,dim,Ordering::byVDIM);
	// Create and fill the grid function
	GridFunction velGf(&velSpace);
	auto & nodes = grid.nodes();
	for (int i = 0; i < nodes.size(); i++) {
		auto & nodeVel = nodes[i].vel;
		for (int j = 0; j < dim; j++) {
			velGf[i*dim+j] = nodeVel[j];
		}
	}

    // Find the point quickly since the background mesh is known to be a cube.
    int numNodes = _mesh.GetNV();
    int nodesPerSide = (dim == 2) ? sqrt(numNodes) - 1 : cbrt(numNodes) - 1;
    double * node1 = _mesh.GetVertex(0);
    double * node2 = _mesh.GetVertex(1);
    double sideLength = node2[0] - node1[0];
    int id = 0;

	// Map the grid velocities to the particles
	mfem::IntegrationPoint intPoint;
	Vector vel(dim);
	for (int i = 0; i < particles.size(); i++) {
		auto & mPoint = particles[i];

		// Find the point - FIXME! Put it on the Point class
		if (dim == 2) {
			id = ((int) (mPoint.pos[0] / sideLength))
					+ nodesPerSide * ((int) (mPoint.pos[1] / sideLength));
		} else if (dim == 3) {
			id = ((int) (mPoint.pos[0] / sideLength))
					+ nodesPerSide * ((int) (mPoint.pos[1] / sideLength))
					+ nodesPerSide * nodesPerSide
							* ((int) (mPoint.pos[2] / sideLength));
		}

		intPoint.Set(mPoint.pos.data(),dim);
		velGf.GetVectorValue(id,intPoint,vel);
		// Write the velocity data back to the point
		for (int j = 0; j < dim; j++) {
			mPoint.vel[j] = vel[j];
		}
	}

	return;

}

void BasicMFEMGridMapper::updateParticleVelocities(const Kelvin::Grid & grid,
		const std::vector<Kelvin::MaterialPoint> & particles,
		std::vector<double> & velocities) const {

	// Create the H1 field
	int dim = _mesh.Dimension();
	H1_FECollection velCol(1, dim);
	FiniteElementSpace velSpace(&_mesh, &velCol, dim, Ordering::byVDIM);
	// Create and fill the grid function
	GridFunction velGf(&velSpace);
	auto & nodes = grid.nodes();
	for (int i = 0; i < nodes.size(); i++) {
		auto & nodeVel = nodes[i].vel;
		for (int j = 0; j < dim; j++) {
			velGf[i * dim + j] = nodeVel[j];
		}
	}

    // Find the point quickly since the background mesh is known to be a cube.
    int numNodes = _mesh.GetNV();
    int nodesPerSide = (dim == 2) ? sqrt(numNodes) - 1 : cbrt(numNodes) - 1;
    double * node1 = _mesh.GetVertex(0);
    double * node2 = _mesh.GetVertex(1);
    double sideLength = node2[0] - node1[0];
    int id = 0;

	// Map the grid velocities to the storage vector
	mfem::IntegrationPoint intPoint;
	Vector vel(dim);
	for (int i = 0; i < particles.size(); i++) {
		auto & mPoint = particles[i];

		// Find the point - FIXME! Put it on the Point class
		if (dim == 2) {
			id = ((int) (mPoint.pos[0] / sideLength))
					+ nodesPerSide * ((int) (mPoint.pos[1] / sideLength));
		} else if (dim == 3) {
			id = ((int) (mPoint.pos[0] / sideLength))
					+ nodesPerSide * ((int) (mPoint.pos[1] / sideLength))
					+ nodesPerSide * nodesPerSide
							* ((int) (mPoint.pos[2] / sideLength));
		}

		intPoint.Set(mPoint.pos.data(), dim);
		velGf.GetVectorValue(id, intPoint, vel);
		// Write the velocity data back to the storage vector
		for (int j = 0; j < dim; j++) {
			velocities[i*dim+j] = vel[j];
		}
	}

	return;

}

} /* namespace Kelvin */
