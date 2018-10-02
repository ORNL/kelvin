/**----------------------------------------------------------------------------
 Copyright (c) 2018-, UT-Battelle, LLC
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
#ifndef DIRICHLETBOUNDARYCONDITION_H_
#define DIRICHLETBOUNDARYCONDITION_H_

#include <mfem.hpp>

namespace Kelvin {

/**
 * This class represents a Dirichlet (or "essential") Boundary Condition.
 */
class DirichletBoundaryCondition {
private:

	/**
	 * The mesh
	 */
	mfem::Mesh & mesh;

	/**
	 * The finite element space
	 */
	mfem::FiniteElementSpace & feSpace;

	/**
	 * The side id where the condition is set
	 */
	int sideId;

	/**
	 * The elements that are on the Dirichlet Boundary
	 */
	mfem::Array<int> elements;

	/**
	 * An array of boundary attributes that describes which elements of the
	 * mesh are on the boundary.
	 */
	mfem::Array<int> boundaryAttributes;

	/**
	 * A constant coefficient that describes the value of the boundary
	 * condition.
	 */
	mfem::ConstantCoefficient coefficient;

public:

	/**
	 * Construction
	 * @param _mesh the mesh
	 * @param _feSpace the finite element space
	 * @param _sideId 0 for all sides, otherwise an integer greater than 0 for
	 * each side on which the condition exists.
	 * @param _value the value of the condition on the given side
	 */
	DirichletBoundaryCondition(mfem::Mesh & _mesh, mfem::FiniteElementSpace & _space,
			int _sideId, double _value);

	/**
	 * Copy constructor.
	 * @param otherCond the condition to copy
	 */
	DirichletBoundaryCondition(const DirichletBoundaryCondition & otherCond);

	/**
	 * This operation returns the elements that are covered by this boundary
	 * condition.
	 * @return the elements
	 */
	mfem::Array<int> & getElements();

	/**
	 * This operation returns the constant coefficient value of this condition.
	 * @return the constant coeffficient with the value of the condition
	 */
	mfem::Coefficient & getCoefficient();

	/**
	 * This operation returns an MFEM-friendly array of boundary attributes
	 * with this boundary marked accordingly in the array.
	 * @return the MFEM boundary attributes array.
	 */
	mfem::Array<int> & getBoundaryAttributes();

};
} /* namespace Kelvin */

#endif /* DIRICHLETBOUNDARYCONDITION_H_ */
