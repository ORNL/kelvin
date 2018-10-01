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
#ifndef SRC_POINT_H_
#define SRC_POINT_H_

#include <vector>

namespace Kelvin {

/**
 * This class represent a basic point. Its dimension is set on creation, with
 * a default value of n=3.
 *
 * This is a basic data class, so access to some member variables is
 * unrestricted. The size of the members is equal to the dimension of the point
 * unless otherwise noted.
 *
 * Points allocate their space on construction, so they should only be created
 * when they are needed.
 */
class Point {
protected:

	/**
	 * The number of dimensions of the point.
	 */
	int nDim;

public:

	/**
	 * The position vector of the point.
	 */
	std::vector<double> pos;

	/**
	 * The velocity vector of the point.
	 */
	std::vector<double> vel;

	/**
	 * The acceleration vector of the point.
	 */
	std::vector<double> acc;

	/**
	 * Constructor
	 */
	Point(int dim = 3);

	/**
	 * Copy constructor
	 */
	Point(const Point & otherPoint);

	/**
	 * Destructor
	 */
	virtual ~Point() {};

	/**
	 * This operation returns the dimension of the point, (i.e. - 1, 2 or 3D).
	 * @return the dimension
	 */
	int dimension() const;
};

} /* namespace Kelvin */

#endif /* SRC_POINT_H_ */
