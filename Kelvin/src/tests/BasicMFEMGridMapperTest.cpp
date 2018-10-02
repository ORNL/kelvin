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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE kelvin

#include <boost/test/included/unit_test.hpp>
#include <mfem.hpp>
#include <vector>
#include <Grid.h>
#include <MeshContainer.h>
#include <H1FESpaceFactory.h>
#include <INIPropertyParser.h>
#include <BasicMFEMGridMapper.h>

using namespace std;
using namespace mfem;
using namespace Kelvin;
using namespace fire;

// Test file names
static std::string inputFile = "2SquaresInput-smallerMesh.ini";

/**
 * This operation checks the basic functionality of the Grid class.
 */
BOOST_AUTO_TEST_CASE(checkGridMapper) {

	// Create the space factory
	H1FESpaceFactory spaceFactory;

	// Load the input file
	INIPropertyParser propertyParser;
	propertyParser.setSource(inputFile);
    propertyParser.parse();

	// Load the mesh
    MeshContainer mc(propertyParser.getPropertyBlock("mesh"),spaceFactory);

    // Create the grid
    Grid grid(mc);
    // Just use the quadrature points as the particles
    auto points = mc.getQuadraturePoints();
    // FIXME! Convert to material points
    std::vector<MaterialPoint> mPoints;
    for (int i = 0; i < points.size(); i++) {
    	MaterialPoint point(points[i]);
    	mPoints.push_back(point);
    }

    // Assemble it to create the arrays of nodes, shapes, masses, etc.
    grid.assemble(mPoints);

    // Set the particle properties to 1 to get a value for the internal forces
    // that is equal to the gradient/sum of gradients.
    for (int i = 0; i < mPoints.size(); i++) {
    	mPoints[i].stress[0] = {1.0,1.0};
    	mPoints[i].stress[1] = {1.0,1.0};
    	mPoints[i].mass = 1.0;
    	mPoints[i].bodyForce[0] = 1.0;
    	mPoints[i].bodyForce[1] = 1.0;
    }

	// Compute the accelerations at the grid points
	grid.updateNodalAccelerations(1.0,mPoints);

	// Set the initial particle velocities and update the grid
    for (int i = 0; i < mPoints.size(); i++) {
    	mPoints[i].vel[0] = 1.0;
    	mPoints[i].vel[1] = 1.0;
    }
	// Compute the velocities at the grid points from the momenta
	grid.updateNodalVelocitiesFromMomenta(mPoints);

	// Compute the velocity at the grid point for time dt = 1.0
	grid.updateNodalVelocities(1.0,mPoints);

	// Create the grid mapper
	BasicMFEMGridMapper mapper(mc.getMesh());

	// Map grid node accelerations to the particle accelerations
	cout << "----- Updated Accelerations" << endl;
	mapper.updateParticleAccelerations(grid,mPoints);
    for (int i = 0; i < mPoints.size(); i++) {
    	cout << mPoints[i].acc[0] << " " << mPoints[i].acc[1] << endl;
    }
    // p1
    BOOST_REQUIRE_CLOSE(3.0,mPoints[0].acc[0],1.0e-15);
    BOOST_REQUIRE_CLOSE(3.0,mPoints[0].acc[1],1.0e-15);
    // p2
    BOOST_REQUIRE_CLOSE(3.0,mPoints[1].acc[0],0.0);
    BOOST_REQUIRE_CLOSE(3.0,mPoints[1].acc[1],1.0e-15);

    // Map grid node velocities back to particles, but store in a vector for
    // further processing
    cout << "----- Updated Velocities (by vector)" << endl;
    std::vector<double> vel(mPoints.size()*2);
    mapper.updateParticleVelocities(grid,mPoints,vel);
    for (int i = 0; i < mPoints.size(); i++) {
    	for (int j = 0; j < 2; j++) {
    		cout << vel[i*2+j] << " ";
    	}
    	cout << endl;
    }
    // p1
    BOOST_REQUIRE_CLOSE(4.0,vel[0],1.0e-15);
    BOOST_REQUIRE_CLOSE(4.0,vel[1],1.0e-15);
    // p2
    BOOST_REQUIRE_CLOSE(4.0,vel[2],1.0e-15);
    BOOST_REQUIRE_CLOSE(4.0,vel[3],1.0e-15);

	// Map grid node velocities to the particle velocities
	cout << "----- Updated Velocities" << endl;
	mapper.updateParticleVelocities(grid,mPoints);
    for (int i = 0; i < mPoints.size(); i++) {
    	cout << mPoints[i].vel[0] << " " << mPoints[i].vel[1] << endl;
    }
    // p1
    BOOST_REQUIRE_CLOSE(4.0,mPoints[0].vel[0],1.0e-15);
    BOOST_REQUIRE_CLOSE(4.0,mPoints[0].vel[1],1.0e-15);
    // p2
    BOOST_REQUIRE_CLOSE(4.0,mPoints[1].vel[0],1.0e-15);
    BOOST_REQUIRE_CLOSE(4.0,mPoints[1].vel[1],1.0e-15);

	return;
}
