
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
#include <KelvinBaseTypes.h>

using namespace std;
using namespace Kelvin;

/**
 * This operation insures that the NodalValueVector can be constructed and read
 * correctly.
 */
BOOST_AUTO_TEST_CASE(checkConstruction) {

	// Check the basic size of the NodalValueVector
	NodalValueVector grad;
	BOOST_REQUIRE_EQUAL(3,grad.dimension());
	BOOST_REQUIRE_EQUAL(-1,grad.nodeId);
	NodalValueVector twoDGrad(2);
	BOOST_REQUIRE_EQUAL(2,twoDGrad.dimension());
	BOOST_REQUIRE_EQUAL(-1,twoDGrad.nodeId);

	// Check component dimensionality
	BOOST_REQUIRE_EQUAL(3,grad.values.size());
	BOOST_REQUIRE_EQUAL(2,twoDGrad.values.size());

	// Loading the NodalValueVectors to test copy construction
	for (int i = 0; i < 2; i++) {
		twoDGrad.values[i] = (double) i;
	}

	// Check copy construction
	twoDGrad.nodeId = 2;
	NodalValueVector grad2(twoDGrad);
	BOOST_REQUIRE_EQUAL(twoDGrad.dimension(),grad2.dimension());
	BOOST_REQUIRE_EQUAL(twoDGrad.nodeId,grad2.nodeId);
	for (int i = 0; i < 2; i++) {
		twoDGrad.values[i] = (double) i;
		BOOST_REQUIRE_CLOSE(twoDGrad.values[i],grad2.values[i],1.0e-15);
	}

	return;
}

/**
 * This operation insures that the NodalValueVector can be cleared correctly.
 */
BOOST_AUTO_TEST_CASE(checkClearing) {
\
	NodalValueVector twoDGrad(2);
	BOOST_REQUIRE_EQUAL(2,twoDGrad.dimension());
	BOOST_REQUIRE_EQUAL(-1,twoDGrad.nodeId);
	BOOST_REQUIRE_EQUAL(2,twoDGrad.values.size());

	// Loading the NodalValueVectors to test copy construction
	for (int i = 0; i < 2; i++) {
		twoDGrad.values[i] = (double) i;
	}

	// Clear the vector and check it
	twoDGrad.clear();
	for (int i = 0; i < 2; i++) {
		BOOST_REQUIRE_CLOSE(0.0,twoDGrad.values[i],0.0);
	}

	return;
}


