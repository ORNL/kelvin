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
#include <TimeIntegrator.h>
#include <stdio.h>
#include <stdlib.h>
#include <StringCaster.h>

using namespace mfem;
using namespace fire;

namespace Kelvin {

TimeIntegrator::TimeIntegrator(ThermalOperator & timeEvOp,
		const std::map<std::string, std::string> & solverProps,
		DataCollection & dc) :
		timeOperator(timeEvOp), dataColl(dc) {

	t = StringCaster<double>::cast(solverProps.at("startTime"));
	tFinal = StringCaster<double>::cast(solverProps.at("finalTime"));
	dt = StringCaster<double>::cast(solverProps.at("initialTimeStep"));
	outputStep = StringCaster<int>::cast(solverProps.at("outputStepFrequency"));
	done = false;
	timeOperator.SetTime(t);
	solver.Init(timeOperator);

}

void TimeIntegrator::integrate() {

	// The time was already set on the operator and the solver initialized in
	// the constructor. No need to repeat it. Just get the solution vector.
	auto & x = timeOperator.solution();

	// Do the time integration
	for (int ti = 1; !done; ti++) {
		// Check the final stepping condition.
		if (t + dt >= tFinal - dt / 2) {
			done = true;
		}
		// Do the step
		solver.Step(x, t, dt);

		// Recover the FEM solution
		timeOperator.recoverSolution();

		// Update the state and re-impose the essential boundary
		// conditions.
		timeOperator.update();

		// Output the result at the current step
		if (done || (ti % outputStep) == 0) {
			cout << "time step: " << ti << ", time: " << t << endl;
			dataColl.SetCycle(ti);
			dataColl.SetTime(t);
			dataColl.Save();
		}
	}

	return;
}

} /* namespace Kelvin */
