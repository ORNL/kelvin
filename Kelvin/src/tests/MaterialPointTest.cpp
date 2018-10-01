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
#include <MaterialPoint.h>

using namespace std;
using namespace Kelvin;

/**
 * This operation insures that the point can be constructed and read correctly.
 */
BOOST_AUTO_TEST_CASE(checkConstruction) {

	// Check the basic size of the point
	int dim = 3;
	MaterialPoint point;
	BOOST_REQUIRE_EQUAL(dim,point.dimension());

	// Check component dimensionality
	BOOST_REQUIRE_EQUAL(dim,point.pos.size());
	BOOST_REQUIRE_EQUAL(dim,point.vel.size());
	BOOST_REQUIRE_EQUAL(dim,point.acc.size());
	BOOST_REQUIRE_EQUAL(dim,point.stress.size());
	BOOST_REQUIRE_EQUAL(dim,point.strain.size());
	BOOST_REQUIRE_EQUAL(dim,point.bodyForce.size());
	// And mass
	BOOST_REQUIRE_EQUAL(0.0,point.mass);
	// Check the sizes of the rows of the stress and strain tensors
	for (int i = 0; i < dim; i++) {
		BOOST_REQUIRE_EQUAL(dim,point.stress[i].size());
		BOOST_REQUIRE_EQUAL(dim,point.strain[i].size());
	}

	// Loading the points to test copy construction
	for (int i = 0; i < dim; i++) {
		point.pos[i] = (double) i;
		point.vel[i] = (double) i;
		point.acc[i] = (double) i;
		point.stress[i] = {(double) i,(double) i, (double) i};
		point.strain[i] = {(double) i,(double) i, (double) i};
		point.bodyForce[i] = (double) i;
	}
	point.mass = 1.0;
	point.materialId = 5.0;

	// Check copy construction
	MaterialPoint point2(point);
	BOOST_REQUIRE_EQUAL(point.dimension(),point2.dimension());
	for (int i = 0; i < dim; i++) {
		BOOST_REQUIRE_CLOSE(point.pos[i],point2.pos[i],1.0e-15);
		BOOST_REQUIRE_CLOSE(point.vel[i],point2.vel[i],1.0e-15);
		BOOST_REQUIRE_CLOSE(point.acc[i],point2.acc[i],1.0e-15);
		for (int j = 0; j < dim; j++) {
			BOOST_REQUIRE_CLOSE(point.stress[i][j],point2.stress[i][j],1.0e-15);
			BOOST_REQUIRE_CLOSE(point.strain[i][j],point2.strain[i][j],1.0e-15);
		}
		BOOST_REQUIRE_CLOSE(point.bodyForce[i],point2.bodyForce[i],1.0e-15);
	}
	BOOST_REQUIRE_CLOSE(point.mass,point2.mass,1.0e-15);
	BOOST_REQUIRE_EQUAL(point.materialId,point2.materialId);

	// Check copy construction from a point
	Point point3;
	point3.pos[0] = 1.0;
	point3.vel[0] = 1.0;
	point3.acc[0] = 1.0;
	MaterialPoint point4(point3);
	BOOST_REQUIRE_CLOSE(point3.pos[0],point4.pos[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(point3.vel[0],point4.vel[0],1.0e-15);
	BOOST_REQUIRE_CLOSE(point3.acc[0],point4.acc[0],1.0e-15);
	// Other attributes should be zero
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			BOOST_REQUIRE_CLOSE(0.0,point4.stress[i][j],1.0e-15);
			BOOST_REQUIRE_CLOSE(0.0,point4.strain[i][j],1.0e-15);
		}
	}
	BOOST_REQUIRE_CLOSE(0.0,point4.mass,0.0);
	BOOST_REQUIRE_EQUAL(0,point4.materialId);

	std::string foo("bar");
	cout << sizeof(foo) << " " << foo.length() << endl;
	foo.resize(foo.length());
	cout << sizeof(foo) << " " << foo.length() << endl;
	cout << sizeof(int) << endl;

	return;
}


