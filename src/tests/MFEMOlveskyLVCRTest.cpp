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
#include <MFEMOlevskyLVCR.h>

using namespace std;
using namespace mfem;
using namespace Kelvin;
using namespace fire;

// Test file names
static std::string inputFile = "2SquaresInput-smallerMesh.ini";

/**
 * This operation checks that the strain rate can be computed correctly from
 * the MFEMOlevskyLVCR class.
 */
BOOST_AUTO_TEST_CASE(checkStrain) {

	// Create the space factory
	H1FESpaceFactory spaceFactory;

	// Create the data holder
	MFEMData data;
	data.load(inputFile);

	// Load the input file
	INIPropertyParser propertyParser;
	propertyParser.setSource(inputFile);
    propertyParser.parse();

	// Load the mesh
    MeshContainer mc(propertyParser.getPropertyBlock("mesh"),spaceFactory);
    int dim = mc.dimension();

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

    // Set the particle properties to 1 to get a value for the internal forces
    // that is equal to the gradient/sum of gradients.
    for (int i = 0; i < mPoints.size(); i++) {
    	mPoints[i].stress[0] = {1.0,1.0};
    	mPoints[i].stress[1] = {1.0,1.0};
    	mPoints[i].mass = 1.0;
    	mPoints[i].bodyForce[0] = 1.0;
    	mPoints[i].bodyForce[1] = 1.0;
    	mPoints[i].vel[0] = 0.0;
    	mPoints[i].vel[1] = 0.0;
    }

    // Assemble it to create the arrays of nodes, shapes, masses, etc.
    grid.assemble(mPoints);

	// Compute the accelerations at the grid points
	grid.updateNodalAccelerations(1.0,mPoints);

	// Compute the velocities at the grid points from the momenta
	grid.updateNodalVelocitiesFromMomenta(mPoints);

	// Compute the velocity at the grid point for time dt = 1.0
	grid.updateNodalVelocities(1.0,mPoints);

	// Create the constitutive equation
	MFEMOlevskyLVCR conRel(data);

	// Compute the strain
	for (int i = 0; i < mPoints.size(); i++) {
		conRel.updateStrainRate(grid,mPoints[i]);
	}
	// Setup reference values
	std::vector<double> p1Strain{0,-1,-1,-2};
	std::vector<double> p2Strain{-4,-3,-3,-2};
	cout << "----- Strains" << endl;
	for (int i = 0; i < mPoints.size(); i++) {
		for (int j = 0; j < dim; j++) {
			for (int k = 0; k < dim; k++) {
				cout << mPoints[i].strain[j][k] << " ";
			}
			cout << endl;
		}
		cout << "-----" << endl;
	}
	// Check point 1
	BOOST_REQUIRE_CLOSE(p1Strain[0],mPoints[0].strain[0][0],1.0e-15);
	BOOST_REQUIRE_CLOSE(p1Strain[1],mPoints[0].strain[0][1],1.0e-15);
	BOOST_REQUIRE_CLOSE(p1Strain[2],mPoints[0].strain[1][0],1.0e-15);
	BOOST_REQUIRE_CLOSE(p1Strain[3],mPoints[0].strain[1][1],1.0e-15);
	// Check point 2
	BOOST_REQUIRE_CLOSE(p2Strain[0],mPoints[1].strain[0][0],1.0e-15);
	BOOST_REQUIRE_CLOSE(p2Strain[1],mPoints[1].strain[0][1],1.0e-15);
	BOOST_REQUIRE_CLOSE(p2Strain[2],mPoints[1].strain[1][0],1.0e-15);
	BOOST_REQUIRE_CLOSE(p2Strain[3],mPoints[1].strain[1][1],1.0e-15);

	return;
}

/**
 * This operation checks that the stress can be computed correctly from the
 * MFEMOlevskyLVCR class.
 */
BOOST_AUTO_TEST_CASE(checkStress) {

	// Create the space factory
	H1FESpaceFactory spaceFactory;

	// Load the input file
	INIPropertyParser propertyParser;
	propertyParser.setSource(inputFile);
    propertyParser.parse();

	// Create the data holder
	MFEMData data;
	data.load(inputFile);

	// Load the mesh
    MeshContainer mc(propertyParser.getPropertyBlock("mesh"),spaceFactory);
    int dim = mc.dimension();

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

    // Set the particle properties to 1 to get a value for the internal forces
    // that is equal to the gradient/sum of gradients.
    for (int i = 0; i < mPoints.size(); i++) {
    	mPoints[i].stress[0] = {1.0,1.0};
    	mPoints[i].stress[1] = {1.0,1.0};
    	mPoints[i].mass = 1.0;
    	mPoints[i].bodyForce[0] = 1.0;
    	mPoints[i].bodyForce[1] = 1.0;
    	mPoints[i].vel[0] = 1.0;
    	mPoints[i].vel[1] = 1.0;
    }

    // Assemble it to create the arrays of nodes, shapes, masses, etc.
    grid.assemble(mPoints);

	// Compute the accelerations at the grid points
	grid.updateNodalAccelerations(1.0,mPoints);

	// Compute the velocities at the grid points from the momenta
	grid.updateNodalVelocitiesFromMomenta(mPoints);

	// Compute the velocity at the grid point for time dt = 1.0
	grid.updateNodalVelocities(1.0,mPoints);

	// Create the constitutive equation
	MFEMOlevskyLVCR conRel(data);

	for (int i = 0; i < mPoints.size(); i++) {
		// Compute the strain
		conRel.updateStrainRate(grid,mPoints[i]);
		// Compute the stress
		conRel.updateStress(grid,mPoints[i]);
	}

	// Setup reference values
	std::vector<double> p1Stress{-169444,-508333,-508333,-1.18611e6};
	std::vector<double> p2Stress{-2.54167e6,-1.525e6,-1.525e6,-1.525e6};
//  I'm guessing that these will change really soon... but for now:
//	----- Stress
//	-169444 -508333
//	-508333 -1.18611e+06
//	-----
//	-2.54167e+06 -1.525e+06
//	-1.525e+06 -1.525e+06
//	-----


	cout << "----- Stress" << endl;
	for (int i = 0; i < mPoints.size(); i++) {
		for (int j = 0; j < dim; j++) {
			for (int k = 0; k < dim; k++) {
				cout << mPoints[i].stress[j][k] << " ";
			}
			cout << endl;
		}
		cout << "-----" << endl;
	}

	// Check point 1
	BOOST_REQUIRE_CLOSE(p1Stress[0],mPoints[0].stress[0][0],1.0e-3);
	BOOST_REQUIRE_CLOSE(p1Stress[1],mPoints[0].stress[0][1],1.0e-3);
	BOOST_REQUIRE_CLOSE(p1Stress[2],mPoints[0].stress[1][0],1.0e-3);
	BOOST_REQUIRE_CLOSE(p1Stress[3],mPoints[0].stress[1][1],1.0e-3);
	// Check point 2
	BOOST_REQUIRE_CLOSE(p2Stress[0],mPoints[1].stress[0][0],1.0e-3);
	BOOST_REQUIRE_CLOSE(p2Stress[1],mPoints[1].stress[0][1],1.0e-3);
	BOOST_REQUIRE_CLOSE(p2Stress[2],mPoints[1].stress[1][0],1.0e-3);
	BOOST_REQUIRE_CLOSE(p2Stress[3],mPoints[1].stress[1][1],1.0e-3);

	return;
}
