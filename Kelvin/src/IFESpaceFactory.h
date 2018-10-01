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
#ifndef IFESPACEFACTORY_H_
#define IFESPACEFACTORY_H_

#include <mfem.hpp>

namespace Kelvin {

/**
 * This class is a factory interface for creating MFEM finite element
 * collections and spaces. It can be implemented and passed to other classes
 * as a mechanism of dynamically selecting the collection and finite element
 * space types.
 */
class IFESpaceFactory {
public:

	/**
	 * Constructor
	 */
	IFESpaceFactory() {};

	/**
	 * Destructor
	 */
	virtual ~IFESpaceFactory() {};

	/**
	 * This operation creates a finite element collection with the given order
	 * and dimensionality. The exact type of the finite element collection
	 * depends on the implementation of this interface.
	 * @param order the order of the mesh elements
	 * @param dim the dimension of the mesh
	 * @return the finite element collection
	 */
	virtual mfem::FiniteElementCollection & getCollection(int order, int dim) = 0;

	/**
	 * This operation creates a finite element space from the mesh and
	 * collection of elements. The exact type of the space is left up to the
	 * implementation.
	 * @param mesh the mesh
	 * @param collection the collection of elements
	 * @return the finite element space
	 */
	virtual mfem::FiniteElementSpace & getFESpace(
			mfem::Mesh & mesh, mfem::FiniteElementCollection & collection) = 0;
};
} /* namespace Kelvin */
#endif /* IFESPACEFACTORY_H_ */

