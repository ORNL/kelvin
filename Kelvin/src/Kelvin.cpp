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

INIPropertyParser propertyParser;

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
	propertyParser = build<INIPropertyParser, const string &>(
			string(input_file));

	H1FESpaceFactory spaceFactory;
	MeshContainer meshContainer(propertyParser.getPropertyBlock("mesh"),
			spaceFactory);
	auto & mesh = meshContainer.getMesh();

	// Create data collection for solution output in the VisIt format
	VisItDataCollection dc(meshContainer.name().c_str(), &mesh);

	// Get the heat transfer coefficient for the material. Assumed to be
	// constant and provided as a property at the moment.
	auto & thermalProps = propertyParser.getPropertyBlock("thermal");

//	// Get the furnace parameters
//	auto & furnaceProps = propertyParser.getPropertyBlock("furnace");
//	double elementTemp = StringCaster<double>::cast(
//			furnaceProps.at("bottomTemperature"));
//	int elementSide = StringCaster<double>::cast(
//			furnaceProps.at("bottomSide"));
//	// Create the furnace element boundary condition
//	meshContainer.setDirichletBoundaryCondition(elementSide,elementTemp);

    // Create the thermal solver
	ThermalOperator thermalOperator(meshContainer,thermalProps,dc);

	// Get the temperature, which will currently reflect the initial
	// conditions.
	auto & temp = thermalOperator.solution();
	// Save the mesh and output the initial conditions
	int precision = 8;
	auto meshName = meshContainer.name();
	string meshOutputName = meshName + ".mesh";
    ofstream omesh(meshOutputName.c_str());
	omesh.precision(precision);
	mesh.Print(omesh);
	ofstream osol("thermal_initial.gf");
	osol.precision(precision);
	temp.Save(osol);

	// Do the time integration. Get solver properties first.
	auto & solverProps = propertyParser.getPropertyBlock("solver");
	TimeIntegrator integrator(thermalOperator,solverProps,dc);
	integrator.integrate();

	// Print a vtk file of the mesh
	ofstream vtkStream("mesh.vtk");
	mesh.PrintVTK(vtkStream,0,1);

	// Save the final result to a serialized grid function file.
    ofstream final_osol("thermal_final.gf");
    final_osol.precision(precision);
    temp.Save(final_osol);

	return EXIT_SUCCESS;
}
