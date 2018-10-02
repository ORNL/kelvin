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
#include <StringCaster.h>
#include <mfem.hpp>
#include <MeshContainer.h>

using namespace fire;
using namespace mfem;

namespace Kelvin {

MeshContainer::MeshContainer(
		const std::map<std::string, std::string> & meshProps,
		IFESpaceFactory & spaceFactory) :
				meshFilename(meshProps.at("file")), mesh(meshFilename.c_str()),
			   _order(StringCaster<int>::cast(meshProps.at("order"))),
			    dim(mesh.Dimension()), _name(meshProps.at("name")),
				feCollection(spaceFactory.getCollection(_order,dim)),
				space(spaceFactory.getFESpace(mesh,feCollection)) {

	// Helpful diagnostic information.
	cout << "Loaded mesh " << meshFilename << ". Mesh dimension = " << dim
			<< " with " << mesh.GetNE() << " elements and " << mesh.GetNV()
			<< " vertices." << endl;

}

MeshContainer::MeshContainer(const char * meshFile, const int & order,
		IFESpaceFactory & spaceFactory) :
		meshFilename(meshFile), mesh(meshFile), _order(order),
		dim(mesh.Dimension()), _name("mesh"),
		feCollection(spaceFactory.getCollection(_order,dim)),
		space(spaceFactory.getFESpace(mesh,feCollection)) {

}

Mesh & MeshContainer::getMesh() {
	return mesh;
}

FiniteElementSpace & MeshContainer::getSpace() {
	return space;
}

int MeshContainer::order() {
	return _order;
}

int MeshContainer::dimension() {
	return dim;
}

const std::string & MeshContainer::name() {
	return _name;
}

void MeshContainer::setDirichletBoundaryCondition(double value) {
	dbConditions.emplace_back(mesh,space,0,value);
}

void MeshContainer::setDirichletBoundaryCondition(int meshSide, double value) {
	dbConditions.emplace_back(mesh,space,meshSide,value);
}

std::vector<DirichletBoundaryCondition> &
MeshContainer::getDirichletBoundaryConditions() {
	return dbConditions;
}

mfem::DenseMatrix MeshContainer::convertPointToMatrix(const std::vector<double> & point) {

	// Pack the point
	mfem::DenseMatrix pointMatrix(point.size(),1);
	for (int i = 0; i < point.size(); i++) {
		pointMatrix(i,0) = point[i];
	}

	return pointMatrix;
}

std::vector<int> MeshContainer::getSurroundingNodeIds(
		const std::vector<double> & point) {
	std::vector<int> ids;

	// Find the element that contains the point
	mfem::Array<int> elementId(1);
	mfem::Array<mfem::IntegrationPoint> intPoint(1);
	auto pointMatrix = convertPointToMatrix(point);
    mesh.FindPoints(pointMatrix,elementId,intPoint);

    // Only proceed if an element containing the point was found
    if (elementId[0] > -1) {
    	// Get the element
        auto * element = mesh.GetElement(elementId[0]);
    	// Get the vertex ids
        auto numVertices = element->GetNVertices();
        auto * vertexIds = element->GetVertices();
    	// Repack the ids
        for (int i = 0; i < numVertices; i++) {
        	ids.push_back(vertexIds[i]);
        }
    }

	return ids;
}

std::vector<double> MeshContainer::getNodalShapes(const std::vector<double> & point) {
	std::vector<double> shapes;

	// Find the element that contains the point
	mfem::Array<int> elementId(1);
	mfem::Array<mfem::IntegrationPoint> intPoint(1);
	auto pointMatrix = convertPointToMatrix(point);
    mesh.FindPoints(pointMatrix,elementId,intPoint);

	// Only proceed if the element was found.
    if (elementId[0] > -1) {

    	// Get the element transform, type and the finite element itself.
    	auto * elemTransform = mesh.GetElementTransformation(elementId[0]);
    	auto type = elemTransform->GetGeometryType();
    	auto * feCollection = space.FEColl();
    	auto * fElement= feCollection->FiniteElementForGeometry(type);

    	// Compute the shape
    	int numShapes = fElement->GetDof();
    	mfem::Vector shapeVec(numShapes);
    	fElement->CalcShape(intPoint[0],shapeVec);

    	// Repack the shapes to return them
    	shapes.resize(numShapes);
    	for (int i = 0; i < numShapes; i++) {
    		shapes[i] = shapeVec[i];
    	}
    }

	return shapes;
}

std::vector<Gradient> MeshContainer::getNodalGradients(const std::vector<double> & point) {
	std::vector<Gradient> gradients;

	// Find the element that contains the point
	mfem::Array<int> elementId(1);
	mfem::Array<mfem::IntegrationPoint> intPoint(1);
	auto pointMatrix = convertPointToMatrix(point);
    mesh.FindPoints(pointMatrix,elementId,intPoint);

	// Only proceed if the element was found.
    if (elementId[0] > -1) {

    	// Get the element transform, type and the finite element itself.
    	auto * elemTransform = mesh.GetElementTransformation(elementId[0]);
    	auto type = elemTransform->GetGeometryType();
    	auto * feCollection = space.FEColl();
    	auto * fElement= feCollection->FiniteElementForGeometry(type);

    	// Compute the gradient
    	int numDof = fElement->GetDof();
    	mfem::DenseMatrix gradientMatrix(numDof,dim);
    	fElement->CalcDShape(intPoint[0],gradientMatrix);

    	cout << "Gradient MC" << endl;
    	gradientMatrix.Print();
    	cout << "-----" << endl;

    	// Repack the gradients to return them. We only need numDof entries in
    	// the vector, so shrink it.
    	//
    	// Repacking this data, especially like this, will turn out to be a
    	// performance bottleneck. So... FIXME!
    	gradients.resize(numDof);
    	// Get the vertex ids
    	// Get the element
        auto * element = mesh.GetElement(elementId[0]);
        auto * vertexIds = element->GetVertices();
        // Assuming numVerts = numDof
    	for (int i = 0; i < numDof; i++) {
    		Gradient grad(dim);
    		grad.nodeId = vertexIds[i];
    		for (int j = 0; j < dim; j++) {
    			grad.values[j] = gradientMatrix(i,j);
    		}
    		gradients[i] = grad;
    	}
    }

	return gradients;
}

std::vector<Point> MeshContainer::getQuadraturePoints() {
	std::vector<Point> points;

	// Get the number of elements and pull the quadrature points for each.
	int numElements = space.GetNE();
	for (int i = 0; i < numElements; i++) {
		// Get the element
		auto * element = space.GetFE(i);
		// Get the coordinate transformation for local->global
		auto * transform = space.GetElementTransformation(i);
		Vector vPoint;
		// Get the quadrature rue
		auto & intRule = IntRules.Get(element->GetGeomType(),_order);
		// Loop over all the quadrature points in the element, transform them,
		// and put them into the list.
		int numIntPoints = intRule.GetNPoints();
		for (int j = 0; j < numIntPoints; j++) {
			// Get the quadrature point
			auto & intPoint = intRule.IntPoint(j);
			// Transform it to global coordinates
			transform->Transform(intPoint,vPoint);
			// Load the point
			Point point(dim);
			for (int j = 0; j < dim; j++) {
				point.pos[j] = vPoint(j);
			}
			points.push_back(point);
		}
	}

	return points;
}

} /* namespace Kelvin */
