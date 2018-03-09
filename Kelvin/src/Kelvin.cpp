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
#include <stdio.h>
#include <stdlib.h>
#include <memory>
#include <INIPropertyParser.h>
#include <StringCaster.h>
#include <mfem.hpp>
#include <IFESpaceFactory.h>
#include <H1FESpaceFactory.h>
#include <MeshContainer.h>
#include <DirichletBoundaryCondition.h>
#include <TimeIntegrator.h>
#include <ThermalOperator.h>

using namespace std;
using namespace mfem;
using namespace fire;
using namespace Kelvin;

class MFEMData {
	fire::INIPropertyParser propertyParser;

	std::unique_ptr<MeshContainer> mc;

	H1FESpaceFactory spaceFactory;

	std::unique_ptr<mfem::DataCollection> dc;

public:

	MFEMData() {

	}

	void load(const std::string & inputFile) {
		// Load the input file
		propertyParser = build<INIPropertyParser, const string &>(
				inputFile);

		// Load the mesh
		mc = make_unique<MeshContainer>(
				propertyParser.getPropertyBlock("mesh"),spaceFactory);

		// Load the data container
		dc = make_unique<VisItDataCollection>(mc->name().c_str(),
				&mc->getMesh());
	}

	fire::INIPropertyParser & properties() {
		return propertyParser;
	}

	MeshContainer & meshContainer() {
		return *mc;
	}

	mfem::DataCollection & collection() {
		return *dc;
	}

};

class MFEMSolver {
public:
	void solve() {

	}
};

class MFEMManager {

	MFEMData data;
	MFEMSolver solver;



public:

	MFEMManager(const string & inputFile, const int argc, char * argv[]) {

		// Create the default command line arguments
		mfem::OptionsParser args(argc, argv);
		const char * inputFilePtr = inputFile.c_str();
		args.AddOption(&inputFilePtr, "-i", "--input", "Input file to use.");

		// Load the data -- FIXME! - This should pull the actual value from args.
		data.load(inputFile);

	}

	/**
	 * Run the thermal solve.
	 */
	void solve() {
		// Get the heat transfer coefficient for the material. Assumed to be
		// constant and provided as a property at the moment.
		auto & thermalProps = data.properties().getPropertyBlock("thermal");

	    // Create the thermal solver
		ThermalOperator thermalOperator(data.meshContainer(),thermalProps,
				data.collection());

		// Do the time integration. Get solver properties first.
		auto & solverProps = data.properties().getPropertyBlock("solver");
		TimeIntegrator integrator(thermalOperator,solverProps,
				data.collection());
		integrator.integrate();
	}

};

/**
 * Main program
 * @param argc the number of input arguments
 * @param argv the input arguments array of argc elements
 * @return EXIT_SUCCESS if successful, otherwise another value.
 */
int main(int argc, char * argv[]) {

	// Input file name - default is input.ini in the present directory.
	string inputFile("input.ini");

	// Create the MFEM problem manager
	MFEMManager manager(inputFile,argc,argv);

	// Do the thermal solve
	manager.solve();

	return EXIT_SUCCESS;
}
