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

using namespace std;
using namespace mfem;
using namespace fire;

/** A time-dependent operator for the right-hand side of the ODE. The DG weak
 form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
 and advection matrices, and b describes the flow on the boundary. This can
 be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
 used to evaluate the right-hand side. */
class TimeEvolutionOperator: public TimeDependentOperator {
private:
	SparseMatrix &M, &K;
	const Vector &b;
	DSmoother M_prec;
	CGSolver M_solver;

	mutable Vector z;

public:
	TimeEvolutionOperator(SparseMatrix &_M, SparseMatrix &_K, const Vector &_b);

	virtual void Mult(const Vector &x, Vector &y) const;

	virtual ~TimeEvolutionOperator() {
	}
};

// Implementation of the time evolution operator
TimeEvolutionOperator::TimeEvolutionOperator(SparseMatrix &_M, SparseMatrix &_K,
		const Vector &_b) :
		TimeDependentOperator(_M.Size()), M(_M), K(_K), b(_b), z(_M.Size()) {
	M_solver.SetPreconditioner(M_prec);
	M_solver.SetOperator(M);

	M_solver.iterative_mode = false;
	M_solver.SetRelTol(1e-9);
	M_solver.SetAbsTol(0.0);
	M_solver.SetMaxIter(100);
	M_solver.SetPrintLevel(0);
}

void TimeEvolutionOperator::Mult(const Vector &x, Vector &y) const {
	// y = M^{-1} (K x + b)
	K.Mult(x, z);
	z += b;
	M_solver.Mult(z, y);
}

/**
 * This is the initial condition function that returns the value of u0 at the
 * position x. Righ now this is just a simple function for testing and should
 * be reconfigured as needed.
 * @param x the position vector where the temperature should be initially
 * configured
 * @return roughly room temperature in K
 */
double initialConditions(const Vector &x) {
	return 25.0;
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
	INIPropertyParser propertyParser = build<INIPropertyParser, const string &>(
			string(input_file));

	// Load the mesh - first get the mesh name from the properties
	auto meshProps = propertyParser.getPropertyBlock("mesh");
	auto meshFilename = meshProps.at("file");
	int order = StringCaster<int>::cast(meshProps.at("order"));
	Mesh *mesh = new Mesh(meshFilename.c_str());
	int dim = mesh->Dimension();

	// Helpful diagnostic information.
	cout << "Loaded mesh " << meshFilename << ". Mesh dimension = " << dim
			<< " with " << mesh->GetNE() << " elements." << endl;

	// Create the time integrator - Ideally get this from input. Just RK6 for
	// now.
	ODESolver *odeSolver = new RK6Solver;

	// Define the finite element space
	H1_FECollection feCollection(order, dim);
	FiniteElementSpace feSpace(mesh, &feCollection);

	// Define the bilinear form for heat transfer
	BilinearForm massMatrix(&feSpace), stiffnessMatrix(&feSpace);
	// Get the heat transfer coefficient for the material. Assumed to be
	// constant and provided as a property at the moment.
	auto thermalProps = propertyParser.getPropertyBlock("thermal");
	double alpha = StringCaster<double>::cast(
			thermalProps.at("constantThermalDiffusivity"));
	// Define the mfem constant coefficient
	ConstantCoefficient thermalCoeff(alpha);
	// Add appropriate domain integrators
	massMatrix.AddDomainIntegrator(new MassIntegrator);
	stiffnessMatrix.AddDomainIntegrator(new DiffusionIntegrator(thermalCoeff));

	// Both matrices can be assembled, but not yet finalized
	int skip_zeros = 0;
	stiffnessMatrix.Assemble(skip_zeros);
	massMatrix.Assemble();

	// Create an empty linear form vector since that portion is zero.
	LinearForm forcingVector(&feSpace);
	ConstantCoefficient zero(0.0);
	forcingVector.AddDomainIntegrator(new DomainLFIntegrator(zero));
	forcingVector.Assemble();

	// Define the initial conditions (temp0) and solution vector (temp). Point
	// to function that takes a vector of coordinates x. "temp" defines the
	// function across the grid that describes the value of the temperature as
	// a function of the vector x. The temp0 function coefficient describes the
	// initial temperature at time 0.
	FunctionCoefficient temp0(initialConditions);
	GridFunction temp(&feSpace);
	temp.ProjectCoefficient(temp0);

	// Get the furnace parameters
	auto furnaceProps = propertyParser.getPropertyBlock("furnace");
	double elementTemp = StringCaster<double>::cast(
			furnaceProps.at("elementTemperature"));
	int elementSide = StringCaster<double>::cast(
			furnaceProps.at("elementSide"));

	// Apply essential (Dirichlet) boundary conditions. In this case, the
	// boundary condition with attribute 99 is the top of the system.
	// FIXME! - will do something more general later.
	//
	// How does this work? Boundary attributes should be numbered sequentially
	// from side 1 such that bdr_attributes.Max() returns the maximum boundary
	// attribute value, which may well be equal to the size and to the number
	// of sides. For simplicity, boundary attributes can be set to zero and
	// then only those that have essential constraints can be reset to a
	// non-zero value, with one being most commonly used in the MFEM examples.
	// Below, side 6 in the mesh has essential boundary conditions imposed,
	// which puts it at index 5 in the array and that value is thus set to 1.
	// In the case where multiple sides had essential boundary conditions, each
	// side would have a non-zero value set in the array and all other values
	// would be zero, (see MFEM example 17 for an example of this case).
	Array<int> dirichletBoundaries(mesh->bdr_attributes.Max()), essBoundsTrueDoF;
	dirichletBoundaries = 0;
	dirichletBoundaries[elementSide-1] = 1;
	feSpace.GetEssentialTrueDofs(dirichletBoundaries,essBoundsTrueDoF);

	// Set the furnace element temperature as a Dirichlet/essential boundary
	// condition
	ConstantCoefficient boundaryTemp(elementTemp);
	temp.ProjectBdrCoefficient(boundaryTemp,dirichletBoundaries);

	// Form the linear system and remove essential degrees of freedom.
	SparseMatrix A;
	Vector B,X;
	stiffnessMatrix.FormLinearSystem(essBoundsTrueDoF,temp,forcingVector,A,X,B);

	// Finalize the bilinear form matrices.
	massMatrix.Finalize();
	stiffnessMatrix.Finalize(skip_zeros);

	int precision = 8;
	// Create data collection for solution output in the VisIt format
	DataCollection *dc = NULL;
	dc = new VisItDataCollection("Example9", mesh);
	dc->SetPrecision(precision);
	dc->RegisterField("solution", &temp);
	dc->SetCycle(0);
	dc->SetTime(0.0);
	dc->Save();

	// Save the mesh and output the initial conditions
	auto meshName = meshProps.at("name");
	string meshOutputName = meshName + ".mesh";
    ofstream omesh(meshOutputName.c_str());
	omesh.precision(precision);
	mesh->Print(omesh);
	ofstream osol("thermal_initial.gf");
	osol.precision(precision);
	temp.Save(osol);

	// Define the time-dependent operator. Point to mass and stiffness matrices
	// from the bilinear forms. Set initial time & initialize.
	TimeEvolutionOperator temperatureOverTime(massMatrix.SpMat(),
			stiffnessMatrix.SpMat(), forcingVector);

	// Do the time integration. Get solver properties first.
	auto solverProps = propertyParser.getPropertyBlock("solver");
	double t = StringCaster<double>::cast(solverProps.at("startTime"));
	double tFinal = StringCaster<double>::cast(solverProps.at("finalTime"));
	double dt = StringCaster<double>::cast(solverProps.at("initialTimeStep"));
	int outputStep = StringCaster<int>::cast(
			solverProps.at("outputStepFrequency"));
	// Set the initial time
	temperatureOverTime.SetTime(t);
	odeSolver->Init(temperatureOverTime);
	// Create loop and manually increment steps.
	bool done = false;
	for (int ti = 0; !done;) {
		// Compute the real time step to avoid overstepping
		double dt_real = min(dt, tFinal - t);
		// Do the step
		odeSolver->Step(temp, t, dt_real);
		// Increment the counter
		ti++;
		// Check the final condition
		done = (t >= tFinal - 1e-8 * dt);
		// Output the result
		if (done || ti % outputStep == 0) {
			cout << "time step: " << ti << ", time: " << t << endl;
			dc->SetCycle(ti);
			dc->SetTime(t);
			dc->Save();
		}
	}

	// Save the final result to a serialized grid function file.
    ofstream final_osol("thermal_final.gf");
    final_osol.precision(precision);
    temp.Save(final_osol);

	// Clean up
	delete odeSolver;

	return EXIT_SUCCESS;
}
