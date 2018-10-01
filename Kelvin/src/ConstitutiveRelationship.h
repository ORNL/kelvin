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
#ifndef SRC_CONSTITUTIVERELATIONSHIP_H_
#define SRC_CONSTITUTIVERELATIONSHIP_H_

#include <MaterialPoint.h>
#include <Grid.h>
#include <map>

namespace Kelvin {

/**
 * This is the base class for materials deformations provided through
 * constitutive equations in Kelvin. It is a simple interface for stress and
 * strain updates.
 */
class ConstitutiveRelationship {

public:

	/**
	 * Destructor
	 */
	virtual ~ConstitutiveRelationship() {};

	/**
	 * This operation updates the strain at the material points.
	 * @param grid the computational grid on which nodal quantities are defined.
	 * @param the material points at which the strains should be updated.
	 */
	virtual void updateStrain(const Kelvin::Grid & grid,
			std::vector<Kelvin::MaterialPoint> & matPoints) = 0;

	/**
	 * This operation updates the stress at the material points.
	 * @param grid the computational grid on which nodal quantities are defined.
	 * @param the material points at which the strains should be updated.
	 */
	virtual void updateStress(const Kelvin::Grid & grid,
			std::vector<Kelvin::MaterialPoint> & matPoints) = 0;

};

} /* namespace Kelvin */

#endif /* SRC_CONSTITUTIVERELATIONSHIP_H_ */
