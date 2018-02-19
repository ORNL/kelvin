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

/**
 * Run the thermal solve.
 * @param meshContainer the container holding the mesh and FE space
 * @param parser the property parser for thermal and solver properties
 * @param dc the data collection to store output data
 */
void solve(MeshContainer & meshContainer,
		INIPropertyParser & parser,
		DataCollection & dc) {
	// Get the heat transfer coefficient for the material. Assumed to be
	// constant and provided as a property at the moment.
	auto & thermalProps = parser.getPropertyBlock("thermal");

    // Create the thermal solver
	ThermalOperator thermalOperator(meshContainer,thermalProps,dc);

	// Do the time integration. Get solver properties first.
	auto & solverProps = parser.getPropertyBlock("solver");
	TimeIntegrator integrator(thermalOperator,solverProps,dc);
	integrator.integrate();
}

/**
 * This operation configures a data container to store simulation output. By
 * default it configures a VisItDataCollection and stores output in the VisIt
 * native format.
 * @param meshContainer the mesh container that holds the FE space and mesh.
 * @return the data collection
 */
DataCollection setupDataCollection(MeshContainer & meshContainer) {
	// Create data collection for solution output in the VisIt format
	VisItDataCollection dc(meshContainer.name().c_str(),
			&meshContainer.getMesh());
	return dc;
}

/**
 * This operation creates the space factory.
 * @return the space factory
 */
H1FESpaceFactory createSpaceFactory() {
	H1FESpaceFactory spaceFactory;
	return spaceFactory;
}

/**
 * This operation loads the mesh and finite element space into a mesh
 * container.
 * @param parser the property parser
 * @param spaceFactory the space factory that the mesh container should use
 * @return the mesh container
 */
MeshContainer loadMesh(INIPropertyParser & parser,
		IFESpaceFactory & spaceFactory) {
	MeshContainer meshContainer(parser.getPropertyBlock("mesh"),
			spaceFactory);
	return meshContainer;
}

/**
 * Main program
 * @param argc the number of input arguments
 * @param argv the input arguments array of argc elements
 * @return EXIT_SUCCESS if successful, otherwise another value.
 */
int main(int argc, char * argv[]) {

	// Input file name - default is input.ini in the present directory.
	const char *input_file = "input.ini";

	// Create the default command line arguments
	OptionsParser args(argc, argv);
	args.AddOption(&input_file, "-i", "--input", "Input file to use.");

	// Load the input file
	INIPropertyParser propertyParser;
	propertyParser = build<INIPropertyParser, const string &>(
			string(input_file));

	// Create the space factory
	H1FESpaceFactory spaceFactory = createSpaceFactory();
	// Load the mesh
	auto meshContainer = loadMesh(propertyParser, spaceFactory);

	// Load the data container
	auto dc = setupDataCollection(meshContainer);

	// Do the thermal solve
	solve(meshContainer,propertyParser,dc);

	return EXIT_SUCCESS;
}
