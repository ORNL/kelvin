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

#include <ConstitutiveRelationshipService.h>
#include <ConstitutiveRelationship.h>
#include <Grid.h>
#include <MaterialPoint.h>
#include <boost/test/included/unit_test.hpp>
#include <memory>
#include <mfem.hpp>
#include <vector>
#include <MeshContainer.h>
#include <H1FESpaceFactory.h>
#include <INIPropertyParser.h>

using namespace std;
using namespace mfem;
using namespace Kelvin;
using namespace fire;

class TestConstitutiveRelationship : public ConstitutiveRelationship {
public:

	/**
	 * This operation updates the strain at the material points.
	 * @param grid the computational grid on which nodal quantities are defined.
	 * @param the material points at which the strains should be updated.
	 */
	virtual void updateStrainRate(const Kelvin::Grid & grid,
			std::vector<Kelvin::MaterialPoint> & matPoints) {

		// Just set the strains to 5.0.
		int dim = matPoints[0].dimension();
		for (int i = 0; i < matPoints.size(); i++) {
			auto & point = matPoints[i];
			for (int j = 0; j < dim; j++) {
				for (int k = 0; k < dim; k++) {
					point.strain[j][k] = 5.0;
				}
			}
		}

	}

	/**
	 * This operation updates the stress at the material points.
	 * @param grid the computational grid on which nodal quantities are defined.
	 * @param the material points at which the strains should be updated.
	 */
	virtual void updateStress(const Kelvin::Grid & grid,
			std::vector<Kelvin::MaterialPoint> & matPoints) {

	}

	virtual ~TestConstitutiveRelationship() {};
};

// Test file names
static std::string inputFile = "2SquaresInput-smallerMesh.ini";

/**
 * This operation checks the basic functionality of the constitutive
 * relationship service.
 */
BOOST_AUTO_TEST_CASE(checkServiceRegistration) {

	// Create the test pointer
	unique_ptr<ConstitutiveRelationship> testRelPtr =
			make_unique<TestConstitutiveRelationship>();
	// Add it to the service
	ConstitutiveRelationshipService::add(0,std::move(testRelPtr));

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
    std::vector<MaterialPoint> mPoints;
    for (int i = 0; i < points.size(); i++) {
    	MaterialPoint point(points[i]);
    	mPoints.push_back(point);
    }

	// Get the test relationship
	auto & relationship = ConstitutiveRelationshipService::get(0);

	// Update the strain using the relationship
	relationship.updateStrainRate(grid,mPoints);

	// Check them
	int dim = mPoints[0].dimension();
	for (int i = 0; i < mPoints.size(); i++) {
		auto & point = mPoints[i];
		for (int j = 0; j < dim; j++) {
			for (int k = 0; k < dim; k++) {
				BOOST_REQUIRE_CLOSE(5.0,point.strain[j][k],1.0e-5);
			}
		}
	}

	return;
}
