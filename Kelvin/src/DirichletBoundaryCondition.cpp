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
#include <DirichletBoundaryCondition.h>

using namespace mfem;

namespace Kelvin {

DirichletBoundaryCondition::DirichletBoundaryCondition(
		Mesh & _mesh, FiniteElementSpace & _feSpace,
		int _sideId, double _value) :
			mesh(_mesh), feSpace(_feSpace), sideId(_sideId),
			coefficient(_value) {

	// Apply essential (Dirichlet) boundary conditions.
	//
	// How does this work? Boundary attributes should be numbered sequentially
	// from side 1 such that bdr_attributes.Max() returns the maximum boundary
	// attribute value, which may well be equal to the size and to the number
	// of sides. For simplicity, boundary attributes can be set to zero and
	// then only those that have essential constraints can be reset to a
	// non-zero value, with one being most commonly used in the MFEM examples.
	// Below, side 6 in the mesh has essential boundary conditions imposed,
	// which puts it at index 5 in the array and that value is thus set to 1.
	// In the case where multiple sides had essential boundary conditions, each
	// side would have a non-zero value set in the array and all other values
	// would be zero, (see MFEM example 17 for an example of this case).
    boundaryAttributes.SetSize(mesh.bdr_attributes.Max());
    // Check the side id. Per the definition, it if is zero then the condition
    // is applied to all sides.
    if (sideId == 0) {
    	boundaryAttributes = 1;
    } else {
    	// Otherwise the boundation condition is only applied to the specified
    	// side.
        boundaryAttributes = 0;
        boundaryAttributes[sideId-1] = 1;
    }

    return;
}

DirichletBoundaryCondition::DirichletBoundaryCondition(
		const DirichletBoundaryCondition & otherCond) : mesh(otherCond.mesh),
				feSpace(otherCond.feSpace), sideId(otherCond.sideId),
				coefficient(otherCond.coefficient) {

}

Array<int> & DirichletBoundaryCondition::getElements() {
	feSpace.GetEssentialTrueDofs(boundaryAttributes,elements);
	return elements;
}

Array<int> & DirichletBoundaryCondition::getBoundaryAttributes() {
	return boundaryAttributes;
}

Coefficient & DirichletBoundaryCondition::getCoefficient() {
    return coefficient;
}

} /* namespace Kelvin */
