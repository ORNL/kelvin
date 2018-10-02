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
#ifndef SRC_KELVINBASETYPES_H_
#define SRC_KELVINBASETYPES_H_

#include <vector>

namespace Kelvin {

/**
 * This class represents a basic NodalValueVector. Its dimension is set on
 * creation, with a default value of n=3. NodalValueVectors should be used in
 * collections such that the NodalValueVectors of many points are a vector of
 * NodalValueVector objects, with each NodalValueVector storing one vector and
 * node Id.
 *
 * This is a basic data class, so access to some member variables is
 * unrestricted. The size of the members is equal to the dimension of the
 * NodalValueVector unless otherwise noted.
 *
 * NodalValueVectors allocate their space on construction, so they should only
 * be created when they are needed.
 */
class NodalValueVector {
protected:

	/**
	 * The number of dimensions of the NodalValueVector.
	 */
	int nDim;

public:

	/**
	 * The ID of the node for the shape function that was used to compute the
	 * NodalValueVector. The default value is -1, which should be an unrealistic
	 * value in most schemes.
	 */
	int nodeId;

	/**
	 * The vector holding the values of the NodalValueVector.
	 */
	std::vector<double> values;

	/**
	 * Constructor
	 */
	NodalValueVector(int dim = 3);

	/**
	 * Copy constructor
	 */
	NodalValueVector(const NodalValueVector & otherNodalValueVector);

	/**
	 * Destructor
	 */
	virtual ~NodalValueVector() {};

	/**
	 * This operation returns the dimension of the NodalValueVector, (i.e. - 1,
	 * 2 or 3D).
	 * @return the dimension
	 */
	int dimension() const;

	/**
	 * This operation clears the contents of the values, but leaves the node id
	 * unchanged.
	 */
	void clear();

};

/**
 * This class represents a basic Gradient. Its dimension is set on
 * creation, with a default value of n=3. Gradients should be used in
 * collections such that the gradients of many points are a vector of
 * Gradient objects, with each Gradient storing one value of dN_i/dU_j.
 */
using Gradient = NodalValueVector;

/**
 * This is a basic vector that stores components of forces.
 */
using ForceVector = NodalValueVector;

} /* namespace Kelvin */

#endif /* SRC_KELVINBASETYPES_H_ */
