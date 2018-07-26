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
#include <limits>
#include <INIPropertyParser.h>
#include <IFESpaceFactory.h>
#include <H1FESpaceFactory.h>
#include <memory>

using namespace std;
using namespace mfem;
using namespace Kelvin;
using namespace fire;

// Test file names
static std::string meshFileName = "2Squares.mesh";
static std::string inputFile = "2SquaresInput.ini";

/* This simple mesh looks like the following:
 *
 * d ---- e ---- f
 * |	  |		 |
 * |	  |		 |
 * a ---- b ---- c
 *
 * with nodal coordinates in real space
 * a = 0,0
 * b = 1,0
 * c = 2,0
 * d = 0,1
 * e = 1,1
 * f = 2,1
 *
 * The following test points are used in the tests below and
 * the value of the shape function for each node is checked.
 * g = 1/2,1/2
 * h = 3/2,1/2
 *
 */

/**
 * This operation checks the basic construction of the mesh container.
 */
BOOST_AUTO_TEST_CASE(checkConstruction) {

	// Create the space factory
	H1FESpaceFactory spaceFactory;

	// Load the input file
	INIPropertyParser propertyParser;
	propertyParser.setSource(inputFile);
    propertyParser.parse();

	// Load the mesh
    MeshContainer mc(propertyParser.getPropertyBlock("mesh"),spaceFactory);

    // Check some basic mesh properties
    BOOST_REQUIRE_EQUAL("2Squares",mc.name());
    BOOST_REQUIRE_EQUAL(1,mc.order());
    BOOST_REQUIRE_EQUAL(2,mc.dimension());

	return;
}

/**
 * This operation checks the mesh to insure that it is properly loaded by the
 * MeshContainer. It includes some basic tests that are similar to those in
 * checkConstruction, but work directly on the underlying mfem mesh to insure
 * that it was properly created.
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

    // Check the surrounding ids for a point in the second element
    std::vector<double> point = {1.5,0.5};
    auto nodeIds = mc.getSurroundingNodeIds(point);
    // There should be four nodes
    BOOST_REQUIRE_EQUAL(4,nodeIds.size());
    // And the ids should match those of the second element
    BOOST_REQUIRE_EQUAL(1,nodeIds[0]);
    BOOST_REQUIRE_EQUAL(2,nodeIds[1]);
    BOOST_REQUIRE_EQUAL(5,nodeIds[2]);
    BOOST_REQUIRE_EQUAL(4,nodeIds[3]);

    // Get and check the quadrature points of the elements.
    cout << "----- Checking quadrature points -----" << endl;
    auto points = mc.getQuadraturePoints();
    // There should be two points for this mesh - one quadrature point per
    // element.
    BOOST_REQUIRE_EQUAL(2,points.size());
    // Check the first point
    auto & firstPoint = points[0];
    BOOST_REQUIRE_CLOSE(0.5,firstPoint.coords[0],1.0e-15);
    BOOST_REQUIRE_CLOSE(0.5,firstPoint.coords[1],1.0e-15);
    // Should be zero, but the only good way to check it is to make sure that
    // it is very small since mfem returns ~6.9e-310.
    BOOST_REQUIRE(firstPoint.coords[2] < std::numeric_limits<double>::min());
    // Check the second point
    auto & secondPoint = points[1];
    BOOST_REQUIRE_CLOSE(1.5,secondPoint.coords[0],1.0e-15);
    BOOST_REQUIRE_CLOSE(0.5,secondPoint.coords[1],1.0e-15);
    BOOST_REQUIRE(secondPoint.coords[2] < std::numeric_limits<double>::min());

    cout << "--------------------------------------" << endl;

//    	std::cout << "Element " << i << ":"<< std::endl;
//    			std::cout << "Local | Global" << endl;
//    			std::cout << intPoint.x << " " << intPoint.y << " " << intPoint.z << " | ";
//    			std::cout << vPoint(0) << " " << vPoint(1) << " " << vPoint(2) << std::endl;

   cout << "End checkMesh" << endl;

	return;

}

/**
 * This operation checks the mesh container to make sure that it can properly
 * pull the quadrature points off the mesh.
 */
BOOST_AUTO_TEST_CASE(checkQuadraturePoints) {

	// Create the space factory
	H1FESpaceFactory spaceFactory;

	// Load the input file
	INIPropertyParser propertyParser;
	propertyParser.setSource(inputFile);
    propertyParser.parse();

	// Load the mesh
    MeshContainer mc(propertyParser.getPropertyBlock("mesh"),spaceFactory);

}

/**
 * This operation checks the mesh container to insure that it can correctly
 * report the value of the nodal shape functions for a given element.
 */
BOOST_AUTO_TEST_CASE(checkShape) {

	// Create the space factory
	H1FESpaceFactory spaceFactory;

	// Load the input file
	INIPropertyParser propertyParser;
	propertyParser.setSource(inputFile);
    propertyParser.parse();

	// Load the mesh
    MeshContainer mc(propertyParser.getPropertyBlock("mesh"),spaceFactory);

    // This mesh has two elements, so grab the transformations for both from the
    // mesh.
    auto & mesh = mc.getMesh();
    auto * e1Transform = mesh.GetElementTransformation(0);
    auto * e2Transform = mesh.GetElementTransformation(0);

    // Get the finite element types (they should be the same)
    auto type1 = e1Transform->GetGeometryType();
    auto type2 = e2Transform->GetGeometryType();
    // Get the finite elements for each
    auto * feCollection = mc.getSpace().FEColl();
    auto * e1 = feCollection->FiniteElementForGeometry(type1);
    auto * e2 = feCollection->FiniteElementForGeometry(type2);

    /* The shapes are computed by transforming the physical coordinates into the
     * coordinates in reference space and then calculating the shape on the
     * canonical element. This requires:
     * 1) Creating a dense matrix with each point as a column. Thus the matrix
     * is n-points by d-dimensions.
     * 2) Locating the elements in which the points reside and the values of the
     * reference coordinates (integration points) in those elements.
     * 3) Calculating the shape on the finite element.
     */

    // Create the matrix of points
    mfem::DenseMatrix points(2);
    // row 0 contains x coordinates
    points(0,0) = 0.5; // x1
    points(0,1) = 1.5; // x2
    // row 1 contains y coordinates
    points(1,0) = 0.5; // y1
    points(1,1) = 0.5; // y2

    // Create return containers for the point data from the mesh
    mfem::Array<int> elementIds(2);
    mfem::Array<mfem::IntegrationPoint> intPoints(2);
    // Find the points in the mesh
    mesh.FindPoints(points,elementIds,intPoints);

    // Print for debugging
    auto & intPt1 = intPoints[0];
    auto & intPt2 = intPoints[1];
    BOOST_TEST_MESSAGE("----- Test Points -----");
    points.Print();
    BOOST_TEST_MESSAGE("intPt1 in e1 = (" << intPt1.x << ", " << intPt1.y << ")");
    BOOST_TEST_MESSAGE("intPt2 in e2 = (" << intPt2.x << ", " << intPt2.y << ")");
    BOOST_TEST_MESSAGE("-----------------------");

    // Calculate the shape on the finite element
    mfem::Vector shape1(e1->GetDof());
    e1->CalcShape(intPt1,shape1);
    BOOST_TEST_MESSAGE("----- Reference Shape Values -----");
    shape1.Print();
    mfem::Vector shape2(e2->GetDof());
    e2->CalcShape(intPt2,shape2);
    shape2.Print();
    BOOST_TEST_MESSAGE("----------------------------------");

    /**
     * Now get the shapes through the MeshContainer API and compare them.
     */

    // Create a vector for each point
    vector<double> pt1Vec = {points(0,0),points(1,0)};
    vector<double> pt2Vec = {points(0,1),points(1,1)};
    // Get the shapes through the mesh container API
    auto pt1Shape = mc.getNodalShapes(pt1Vec);
    auto pt2Shape = mc.getNodalShapes(pt2Vec);

    // Check the shapes - first, they should be the same size as each other
    BOOST_REQUIRE_EQUAL(pt1Shape.size(),pt2Shape.size());
    // They should be the same size as MFEM's vectors.
    BOOST_REQUIRE_EQUAL(pt1Shape.size(),shape1.Size());
    BOOST_REQUIRE_EQUAL(pt2Shape.size(),shape2.Size());
    // The values should match. e1...
    for (int i = 0; i < pt1Shape.size(); i++) {
    	BOOST_REQUIRE_CLOSE(pt1Shape[i],shape1[i],1.0e-15);
    }
    // e2
    for (int i = 0; i < pt2Shape.size(); i++) {
       	BOOST_REQUIRE_CLOSE(pt2Shape[i],shape2[i],1.0e-15);
    }

	return;

}
