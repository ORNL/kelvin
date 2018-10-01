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
#include <ConstitutiveRelationshipService.h>
#include <stdexcept>
#include <iostream>

using namespace std;

namespace Kelvin {

// Initialize the map
std::map<int,std::unique_ptr<ConstitutiveRelationship>>
	ConstitutiveRelationshipService::_relationships;

ConstitutiveRelationshipService::~ConstitutiveRelationshipService() {
	// TODO Auto-generated destructor stub
}

void ConstitutiveRelationshipService::add(const int & id,
		std::unique_ptr<ConstitutiveRelationship> relationship) {

	if (_relationships.count(id) == 0) {
		_relationships[id] = std::move(relationship);
	} else {
		std::stringstream error;
		error << "Constitutive relationship with id "
				<< id << " already exists." << std::endl;
		throw std::runtime_error(error.str());
	}

	return;
}

ConstitutiveRelationship &
	ConstitutiveRelationshipService::get(const int & id) {
	return *_relationships[0];
}

} /* namespace Kelvin */
