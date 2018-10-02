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
#ifndef SRC_CONSTITUTIVERELATIONSHIPSERVICE_H_
#define SRC_CONSTITUTIVERELATIONSHIPSERVICE_H_

#include <ConstitutiveRelationship.h>
#include <map>
#include <memory>

namespace Kelvin {

/**
 * This is a global service for registering and retrieving constitutive
 * relationships to enable multi-material support. Clients register
 * constitutive equations with the register() operaiton and retrieve them with
 * get(). All constitutive relationships are assigned a unique id.
 *
 * The easiest way for Clients to access this service is to call it in or near
 * main().
 */
class ConstitutiveRelationshipService {
private:

	/**
	 * The map that indexes constitutive equations against integer ids.
	 */
	static std::map<int,std::unique_ptr<ConstitutiveRelationship>>
		_relationships;

public:

	/**
	 * This operation registers a new constitutive relationship with the
	 * service that can be used by clients.
	 * @param relationship the constitutive relationship to register
	 * @param id the material id for the constitutive relationship
	 */
	static void add(const int & id,
			std::unique_ptr<ConstitutiveRelationship> relationship);

	/**
	 * This operation retrieves a previously registered constitutive
	 * relationship. It will throw an exception if the relationship was not
	 * found in the registry.
	 * @paramid the material id for the constitutive relationship
	 * @return the constitutive relationship
	 */
	static ConstitutiveRelationship & get(const int & id);

	/**
	 * Destructor
	 */
	virtual ~ConstitutiveRelationshipService();
};

} /* namespace Kelvin */

#endif /* SRC_CONSTITUTIVERELATIONSHIPSERVICE_H_ */