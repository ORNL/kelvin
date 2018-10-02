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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE kelvin

#include <boost/test/included/unit_test.hpp>
#include <mfem.hpp>
#include <vector>
#include <Grid.h>
#include <MeshContainer.h>
#include <H1FESpaceFactory.h>
#include <INIPropertyParser.h>

using namespace std;
using namespace mfem;
using namespace Kelvin;
using namespace fire;

// Test file names
static std::string inputFile = "2SquaresInput-smallerMesh.ini";

/**
 * This operation checks the basic functionality of the Grid class.
 */
BOOST_AUTO_TEST_CASE(checkGrid) {

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

    // Check the dimension
    BOOST_REQUIRE_EQUAL(mc.dimension(),grid.dimension());

    // Assemble it to create the arrays of nodes, shapes, masses, etc.
    grid.assemble(mPoints);

    // Check the nodal positions
    auto & nodes = grid.nodes();
    cout << "Num nodes = " << nodes.size() << endl;
    BOOST_REQUIRE_EQUAL(6,nodes.size());
    // (0.0,0.0) - Note that I am copying the value of the point here. That's
    // OK for a small test, but don't do it in production code! Loop instead.
    auto point = nodes[0].pos;
    BOOST_REQUIRE_CLOSE(0.0,point[0],0.0);
    BOOST_REQUIRE_CLOSE(0.0,point[1],0.0);
    // (1.0,0.0)
    point = nodes[1].pos;
    BOOST_REQUIRE_CLOSE(1.0,point[0],1.0e-15);
    BOOST_REQUIRE_CLOSE(0.0,point[1],0.0);
    // (2.0,0.0)
    point = nodes[2].pos;
    BOOST_REQUIRE_CLOSE(2.0,point[0],1.0e-15);
    BOOST_REQUIRE_CLOSE(0.0,point[1],0.0);
    // (0.0,1.0)
    point = nodes[3].pos;
    BOOST_REQUIRE_CLOSE(0.0,point[0],0.0);
    BOOST_REQUIRE_CLOSE(1.0,point[1],1.0e-15);
    // (1.0,1.0)
    point = nodes[4].pos;
    BOOST_REQUIRE_CLOSE(1.0,point[0],1.0e-15);
    BOOST_REQUIRE_CLOSE(1.0,point[1],1.0e-15);
    // (2.0,1.0)
    point = nodes[5].pos;
    BOOST_REQUIRE_CLOSE(2.0,point[0],1.0e-15);
    BOOST_REQUIRE_CLOSE(1.0,point[1],1.0e-15);

    // Check the gradient matrix. A basic sanity check is enough since the
    // gradients are thoroughly tested by the mesh container.
    auto & gradients = grid.gradients();
    BOOST_REQUIRE_EQUAL(2,gradients.size());
    // Note that .at() must be used because the return value is const &
    // and *this cannot be accepted as in put argument to [].
    auto & gradList1 = gradients.at(0);
    // 4 gradients, one per shape functiond: N1/dx, dN1/dy, dN2/dx, dN2/dy
    BOOST_REQUIRE_EQUAL(4,gradList1.size());
    // Just make sure that grad3 has the correct node id
    vector<int> ids = {0,1,4,3};
    for (int i = 0; i < gradList1.size(); i++) {
    	auto & grad = gradList1[i];
    	BOOST_REQUIRE_EQUAL(ids[i],grad.nodeId);
    }

    // Handy loop for debugging the gradients that I don't want to rewrite!
    cout << "----- Gradients -----" << endl;
    for (int i = 0; i < gradients.size(); i++) {
    	auto & elemGrads = gradients.at(i);
    	for (int j = 0; j < elemGrads.size(); j++) {
    		cout << elemGrads[j].values[0] << " "
    				<< elemGrads[j].values[1]
					<< " | " << elemGrads[j].nodeId << endl;
    	}
    	cout << endl;
    }
    cout << "-----" << endl;

    // This still doesn't deal with my fast sort problem though. Do I need that
    // here?
    // Why did I decide that I didn't need it in the grid class?

    // Check the mass matrix. A detailed check is warranted since the Grid
    // creates it.
    std::vector<double> refMasses = {0.0625, 0.0625, 0.0, 0.0625, 0.0625, 0.0,
    		0.0625, 0.125, 0.0625, 0.0625, 0.125, 0.0625, 0.0, 0.0625, 0.0625,
			0.0, 0.0625, 0.0625, 0.0625, 0.0625, 0.0, 0.0625, 0.0625, 0.0,
			0.0625,	0.125, 0.0625, 0.0625, 0.125, 0.0625, 0.0, 0.0625, 0.0625,
			0.0, 0.0625, 0.0625};
    auto & massMatrix = grid.massMatrix(mPoints);
    int k = 0;
    for (int i = 0; i < 6; i++) {
    	for (int j = 0; j < 6; j++) {
    		auto mij = massMatrix(i,j);
    		BOOST_REQUIRE_CLOSE(refMasses[k],mij,1.0e-15);
    		k++;
    	}
    }
    // Check the other version of the massMatrix() call that returns the
    // previously computed mass matrix.
    auto & storedMassMatrix = grid.massMatrix();
    k = 0;
    for (int i = 0; i < 6; i++) {
    	for (int j = 0; j < 6; j++) {
    		auto mij = storedMassMatrix(i,j);
    		BOOST_REQUIRE_CLOSE(refMasses[k],mij,1.0e-15);
    		k++;
    	}
    }

    // Not checking lumping - that is done in the mass matrix test.

    // Set the particle properties to 1 to get a value for the internal forces
    // that is equal to the gradient/sum of gradients.
    for (int i = 0; i < mPoints.size(); i++) {
    	mPoints[i].stress[0] = {1.0,1.0};
    	mPoints[i].stress[1] = {1.0,1.0};
    	mPoints[i].mass = 1.0;
    	mPoints[i].bodyForce[0] = 1.0;
    	mPoints[i].bodyForce[1] = 1.0;
    }

    // Check the computation of the internal forces. Should be equal to the
    // gradients given the configuration of the particles above.
    auto & internalForces = grid.internalForces(mPoints);
    // For this mesh, this should be equal to the number of nodes.
    BOOST_REQUIRE_EQUAL(nodes.size(),internalForces.size());
    for (int i = 0; i < internalForces.size(); i++) {
    	auto & forceVector = internalForces[i];
    	BOOST_REQUIRE_EQUAL(2,forceVector.dimension());
    }
    // n1
	BOOST_REQUIRE_CLOSE(1.0,internalForces[0].values[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,internalForces[0].values[1],1.0e-15);
	// n2
	BOOST_REQUIRE_CLOSE(1.0,internalForces[1].values[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,internalForces[1].values[1],1.0e-15);
	// n3
	BOOST_REQUIRE_CLOSE(0.0,internalForces[2].values[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(0.0,internalForces[2].values[1],1.0e-15);
	// n4
	BOOST_REQUIRE_CLOSE(0.0,internalForces[3].values[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(0.0,internalForces[3].values[1], 1.0e-15);
	// n5
	BOOST_REQUIRE_CLOSE(1.0,internalForces[4].values[0], 1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,internalForces[4].values[1], 1.0e-15);
	// n6
	BOOST_REQUIRE_CLOSE(-1.0,internalForces[5].values[0], 1.0e-15);
	BOOST_REQUIRE_CLOSE(-1.0,internalForces[5].values[1], 1.0e-15);

	// Check the computation of the external forces. Does not current include
	// traction, so it should be equal to the shapes given the particle
	// configuration above.
	auto & externalForces = grid.externalForces(mPoints);
	BOOST_REQUIRE_EQUAL(nodes.size(),externalForces.size());
    for (int i = 0; i < externalForces.size(); i++) {
    	auto & forceVector = externalForces[i];
    	BOOST_REQUIRE_EQUAL(2,forceVector.dimension());
    }
    // n1
	BOOST_REQUIRE_CLOSE(0.25,externalForces[0].values[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(0.25,externalForces[0].values[1],1.0e-15);
	// n2
	BOOST_REQUIRE_CLOSE(0.5,externalForces[1].values[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(0.5,externalForces[1].values[1],1.0e-15);
	// n3
	BOOST_REQUIRE_CLOSE(0.25,externalForces[2].values[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(0.25,externalForces[2].values[1],1.0e-15);
	// n4
	BOOST_REQUIRE_CLOSE(0.25,externalForces[3].values[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(0.25,externalForces[3].values[1], 1.0e-15);
	// n5
	BOOST_REQUIRE_CLOSE(0.5,externalForces[4].values[0], 1.0e-15);
	BOOST_REQUIRE_CLOSE(0.5,externalForces[4].values[1], 1.0e-15);
	// n6
	BOOST_REQUIRE_CLOSE(0.25,externalForces[5].values[0], 1.0e-15);
	BOOST_REQUIRE_CLOSE(0.25,externalForces[5].values[1], 1.0e-15);

	// Compute the accelerations at the grid points
	grid.updateNodalAccelerations(1.0,mPoints);
	auto lumpedMasses = massMatrix.lump();
	// Check them. Should be a = (f_int + f_ex)/m.
	cout << "----- Accelerations" << endl;
	for (int i = 0; i < nodes.size(); i++) {
		cout << i << " ";
		for (int j = 0; j < 2; j++) {
			cout << nodes[i].acc[j] << " ";
		}
		cout << "| " << lumpedMasses[i] << endl;
	}

	// FIXME! Get your notebook out and compare against the newly computed f_int values

	// n1
	BOOST_REQUIRE_CLOSE(5.0,nodes[0].acc[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(5.0,nodes[0].acc[1],1.0e-15);
	// n2
	BOOST_REQUIRE_CLOSE(3.0,nodes[1].acc[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(3.0,nodes[1].acc[1],1.0e-15);
	// n3
	BOOST_REQUIRE_CLOSE(1.0,nodes[2].acc[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,nodes[2].acc[1],1.0e-15);
	// n4
	BOOST_REQUIRE_CLOSE(1.0,nodes[3].acc[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,nodes[3].acc[1],1.0e-15);
	// n5
	BOOST_REQUIRE_CLOSE(3.0,nodes[4].acc[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(3.0,nodes[4].acc[1],1.0e-15);
	// n6
	BOOST_REQUIRE_CLOSE(-3.0,nodes[5].acc[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(-3.0,nodes[5].acc[1],1.0e-15);

	// Set the initial particle velocities and update the grid
    for (int i = 0; i < mPoints.size(); i++) {
    	mPoints[i].vel[0] = 1.0;
    	mPoints[i].vel[1] = 1.0;
    }
	// Compute the velocities at the grid points from the momenta
	grid.updateNodalVelocitiesFromMomenta(mPoints);
	// Check them.
	cout << "----- Velocities from Momenta" << endl;
	for (int i = 0; i < nodes.size(); i++) {
		cout << i << " ";
		for (int j = 0; j < 2; j++) {
			cout << nodes[i].vel[j] << " ";
		}
		cout << "| " << lumpedMasses[i] << endl;
	}
	// n1
	BOOST_REQUIRE_CLOSE(1.0,nodes[0].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,nodes[0].vel[1],1.0e-15);
	// n2
	BOOST_REQUIRE_CLOSE(1.0,nodes[1].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,nodes[1].vel[1],1.0e-15);
	// n3
	BOOST_REQUIRE_CLOSE(1.0,nodes[2].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,nodes[2].vel[1],1.0e-15);
	// n4
	BOOST_REQUIRE_CLOSE(1.0,nodes[3].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,nodes[3].vel[1],1.0e-15);
	// n5
	BOOST_REQUIRE_CLOSE(1.0,nodes[4].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,nodes[4].vel[1],1.0e-15);
	// n6
	BOOST_REQUIRE_CLOSE(1.0,nodes[5].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(1.0,nodes[5].vel[1],1.0e-15);

	// Compute the velocity at the grid point for time dt = 1.0
	grid.updateNodalVelocities(1.0,mPoints);
	// Check them.
	cout << "----- Velocities with acceleration update" << endl;
	for (int i = 0; i < nodes.size(); i++) {
		cout << i << " ";
		for (int j = 0; j < 2; j++) {
			cout << nodes[i].vel[j] << " ";
		}
		cout << "| " << lumpedMasses[i] << endl;
	}
	// n1
	BOOST_REQUIRE_CLOSE(6.0,nodes[0].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(6.0,nodes[0].vel[1],1.0e-15);
	// n2
	BOOST_REQUIRE_CLOSE(4.0,nodes[1].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(4.0,nodes[1].vel[1],1.0e-15);
	// n3
	BOOST_REQUIRE_CLOSE(2.0,nodes[2].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(2.0,nodes[2].vel[1],1.0e-15);
	// n4
	BOOST_REQUIRE_CLOSE(2.0,nodes[3].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(2.0,nodes[3].vel[1],1.0e-15);
	// n5
	BOOST_REQUIRE_CLOSE(4.0,nodes[4].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(4.0,nodes[4].vel[1],1.0e-15);
	// n6
	BOOST_REQUIRE_CLOSE(-2.0,nodes[5].vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(-2.0,nodes[5].vel[1],1.0e-15);

	// Get the massive node set and make sure it contains [0,5].
	auto & nodeSet = grid.massiveNodeSet();
	for (int i = 0; i < 6; i++) {
		auto result = nodeSet.find(i);
		BOOST_REQUIRE(result != nodeSet.end());
	}

	return;
}
