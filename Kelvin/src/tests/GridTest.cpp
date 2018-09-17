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
    // Assemble it to create the arrays of nodes, shapes, masses, etc.
    grid.assemble(points);

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

    // This still doesn't deal with my fast sort problem though. Do I need that here?
    // Why did I decide that I didn't need it in the grid class?

    // Check the mass matrix. A detailed check is warranted since the Grid
    // creates it.
    std::vector<double> refMasses = {0.0625, 0.0625, 0.0, 0.0625, 0.0625, 0.0,
    		0.0625, 0.125, 0.0625, 0.0625, 0.125, 0.0625, 0.0, 0.0625, 0.0625,
			0.0, 0.0625, 0.0625, 0.0625, 0.0625, 0.0, 0.0625, 0.0625, 0.0,
			0.0625,	0.125, 0.0625, 0.0625, 0.125, 0.0625, 0.0, 0.0625, 0.0625,
			0.0, 0.0625, 0.0625};
    auto & massMatrix = grid.massMatrix();
    int k = 0;
    for (int i = 0; i < 6; i++) {
    	for (int j = 0; j < 6; j++) {
    		auto mij = massMatrix(i,j);
    		BOOST_REQUIRE_CLOSE(refMasses[k],mij,1.0e-15);
    		k++;
    	}
    }

    // Not checking lumping - that is done in the mass matrix test.

    // Check the acceleration at the nodes
//    auto & acc = grid.acc();

    // FIXME! Tests!

    // Check the velocity at the nodes
//    auto & vel = grid.vel();

    // FIXME! Tests!

    BOOST_FAIL("Not yet implemented.");

	return;

}
