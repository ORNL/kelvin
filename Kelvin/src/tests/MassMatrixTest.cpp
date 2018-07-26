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
#define BOOST_TEST_MODULE Kelvin

#include <boost/test/included/unit_test.hpp>
#include <MeshContainer.h>
#include <mfem.hpp>
#include <vector>
#include <INIPropertyParser.h>
#include <IFESpaceFactory.h>
#include <H1FESpaceFactory.h>

using namespace std;
using namespace mfem;
using namespace Kelvin;
using namespace fire;

// Test file names
static std::string meshFileName = "2Squares.mesh";
static std::string inputFile = "2SquaresInput.ini";

/**
 *
 */
BOOST_AUTO_TEST_CASE(checkConstruction) {


	return;
}

/**
 */
BOOST_AUTO_TEST_CASE(checkMesh) {

	// Create the space factory
	H1FESpaceFactory spaceFactory;

	// Load the input file
	INIPropertyParser propertyParser;
	propertyParser.setSource(inputFile);
    propertyParser.parse();

	// Load the mesh
    MeshContainer mc(propertyParser.getPropertyBlock("mesh"),spaceFactory);

    // Check some mesh properties
    auto & mesh = mc.getMesh();
    BOOST_REQUIRE_EQUAL(2,mesh.Dimension());
    BOOST_REQUIRE_EQUAL(1,mesh.GetElementTransformation(0)->Order());

	return;

}

/**
 */
BOOST_AUTO_TEST_CASE(checkShape) {



	return;

}
