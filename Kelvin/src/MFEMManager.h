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

 Author(s): Jay Jay Billings (billingsjj <at> ornl <dot> gov)
 -----------------------------------------------------------------------------*/
#ifndef SRC_MFEMMANAGER_H_
#define SRC_MFEMMANAGER_H_

#include <MFEMData.h>
#include <MFEMSolver.h>

namespace Kelvin {

/**
 * This class is a simple manager that provides functions to initialize and
 * solve problems using MFEM.
 */
class MFEMManager {

	/**
	 * The data used in the problem solved by MFEM.
	 */
	MFEMData data;

	/**
	 * The solver used with MFEM to compute the solution.
	 */
	MFEMSolver solver;

public:

	/**
	 * Constructor
	 */
	MFEMManager();

	/**
	 * Destructor
	 */
	virtual ~MFEMManager();

	/**
	 * This operation sets up the MFEM problem in the manager by configuring
	 * input options and data.
	 */
	void setup(const std::string & inputFile, const int argc, char * argv[]);

	/**
	 * Run the solve. This will delegate to the pre-configured solver.
	 */
	void solve();
};

} /* namespace Kelvin */

#endif /* SRC_MFEMMANAGER_H_ */
