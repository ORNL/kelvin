/**----------------------------------------------------------------------------
 Copyright (c) 2015-, UT-Battelle, LLC
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

 Author(s): Jay Jay Billings (jayjaybillings <at> gmail <dot> com)
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
			<< " with " << mesh.GetNE() << " elements." << endl;

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

std::vector<DirichletBoundaryCondition> & MeshContainer::getDirichletBoundaryConditions() {
	return dbConditions;
}


} /* namespace Kelvin */
