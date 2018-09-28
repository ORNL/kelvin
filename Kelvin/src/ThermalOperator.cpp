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
#include <StringCaster.h>
#include <ThermalOperator.h>

using namespace std;
using namespace fire;
using namespace mfem;

namespace Kelvin {

ThermalOperator::ThermalOperator(MeshContainer & _meshContainer,
		const std::map<std::string, std::string> & thermalProps,
		DataCollection & _dataColl) :
		TimeDependentOperator(_meshContainer.getSpace().GetTrueVSize()),
		meshContainer(_meshContainer),
		feSpace(_meshContainer.getSpace()),
		dataColl(_dataColl), forcingVector(&_meshContainer.getSpace()),
		temperature(&_meshContainer.getSpace()), zero(0.0), mPreconditioner(),
		mSolver(), tempSparseMat(), tempPreconditioner(), tempSolver(),
		z(forcingVector.Size()),
		conductionCoeff(_meshContainer.dimension(),temperature),
		sparseMassMatrix(), K(), b(), x() {

	// Set the thermal properties
	alpha = StringCaster<double>::cast(
							thermalProps.at("constantThermalDiffusivity"));
	density = StringCaster<double>::cast(
			thermalProps.at("density"));
	specificHeat = StringCaster<double>::cast(
			thermalProps.at("constantSpecificHeatCapacity"));
	conductivity = StringCaster<double>::cast(
			thermalProps.at("constantConductivity"));
	initialSurfaceTemperature = StringCaster<double>::cast(
			thermalProps.at("surfaceTemperature"));
	initialInteriorTemperature = StringCaster<double>::cast(
			thermalProps.at("initialTemperature"));

	// Define the initial conditions (temp0) and solution vector (temp). Point
	// to function that takes a vector of coordinates x. "temp" defines the
	// function across the grid that describes the value of the temperature as
	// a function of the vector x. The temp0 function coefficient describes the
	// initial temperature at time 0.
	ConstantCoefficient temp0(initialInteriorTemperature);
	temperature.ProjectCoefficient(temp0);

	// Set the surface temperature - a Dirichlet Boundary Condition.
	setSurfaceBoundaryConditions();

	// Get the boundary conditions - this will only apply a single surface
	// temperature across the entire boundary for now.
	auto & dbcs = meshContainer.getDirichletBoundaryConditions();
	auto & dbc = dbcs[0]; // FIXME! - Only considering one DBC right now!

	// Set the boundary conditions due to an incoming heat flux on the surface.
	// setFluxes();

	// I will need to think about how I combine the heat flux and surface
	// temperature conditions.


	// Ideally I would only attempt to eliminate the boundary conditions if
	// dbcs.Size() > 0. This could be done by removing the EliminateEssentialBC
	// calls from the statements below to a separate function, I think, that
	// would be called before finalization, which could happen at the end of
	// of this constructor instead of its current position. It may be that
	// the bilinear forms need to be finalized before being passed to the
	// so, OK.
	//
	// Those calls in question could presumably be moved to
	// setSurfaceBoundaryConditions, right?

	// Assemble the forcing functions

	ConstantCoefficient zero_coeff(0.0);
    forcingVector.AddDomainIntegrator(new DomainLFIntegrator(zero_coeff));
	forcingVector.Assemble();

	// Create the mass matrix as a bilinear form with an integrator..
	massMatrix = make_unique<BilinearForm>(&feSpace);
	massMatrix->AddDomainIntegrator(new MassIntegrator);
	massMatrix->Assemble(skip_zeros);
	massMatrix->FormSystemMatrix(dbc.getElements(), sparseMassMatrix);
//	massMatrix->Finalize(skip_zeros);

	// Configure the mass matrix solver.
	mSolver.SetPreconditioner(mPreconditioner);
	mSolver.SetOperator(massMatrix->SpMat());
	mSolver.iterative_mode = false;
	mSolver.SetRelTol(1.0e-9);
	mSolver.SetAbsTol(0.0);
	mSolver.SetMaxIter(100);
	mSolver.SetPrintLevel(0);

	// Setup the stiffness matrix. The stiffness matrix uses a diffusion
	// integrator because the thermal kernel is the diffusion/laplace/poisson
	// operator.
	ConstantCoefficient alphaCoeff(alpha);
	stiffnessMatrix = make_unique<BilinearForm>(&feSpace);
	stiffnessMatrix->AddDomainIntegrator(new DiffusionIntegrator(alphaCoeff));
	stiffnessMatrix->Assemble(skip_zeros);
//	stiffnessMatrix->Finalize(skip_zeros);

	// Setup the temporary solver
	tempSolver.iterative_mode = false;
	tempSolver.SetRelTol(1.0e-8);
	tempSolver.SetAbsTol(0.0);
	tempSolver.SetMaxIter(100);
	tempSolver.SetPrintLevel(0);
	tempSolver.SetPreconditioner(tempPreconditioner);

	// Update the stiffness matrix (initially configure it) and the temperature
	// array, including applying Dirichlet (essential) boundary conditions.
	update();

	// FIXME! Do I need that call to update above?

	// Setup initial properties for the data collection before the solve.
	// Note: Make this a private function in the solver.
	int precision = 8;
	dataColl.SetPrecision(precision);
	dataColl.RegisterField("temperature", &temperature);
	dataColl.SetCycle(0);
	dataColl.SetTime(0.0);
	dataColl.Save();

	return;
}

//void ThermalOperator::setFluxes() {
//
//	// Create a heat flux coefficient that is restricted to flow only over the
//	// part of the boundary where the user has indicated the flux exists. This
//	// requires creating a vector coefficient first. In this case the vector
//	// coefficient is a constant coefficient because the heat flux is constant
//	// at the boundary.
//	//
//	// In the heat equation, heat flux is scaled by the inverse volumetric heat
//	// capacity, 1/(density*specificHeat). It is also negative, which indicates
//	// that the heat is flowing from hot to cold.
//	double volHeatCap = -1.0/(density*specificHeat);
//	Array<int> bdrs(7);
//	bdrs = 0;
//	bdrs[5] = 1;
//	Vector heatFlux(3);
//	heatFlux(0) = 0.0;
//	heatFlux(1) = 0.0;
//	heatFlux(2) = 5000.0*volHeatCap;
//	VectorConstantCoefficient HFVectorCoeff(heatFlux);
//	// Restrict it to the boundary defined by the flux boundary condition.
//	VectorRestrictedCoefficient HFVectorRestrictedCoeff(HFVectorCoeff,bdrs);
//
//	// Create the linear form for the heat flux. This linear form is given by
//	// dT/dt \dot n = q \dot n where q is the vector flux at the boundary and n
//	// is the normal vector and \dot indicates the dot product.
//	forcingVector.AddBoundaryIntegrator(
//			new BoundaryNormalLFIntegrator(HFVectorRestrictedCoeff));
//
//	return;
//}

void ThermalOperator::setSurfaceBoundaryConditions() {
	// Set the temperature on all boundaries to the surface temperature
	meshContainer.setDirichletBoundaryCondition(initialSurfaceTemperature);
	// Get the surface boundary condition and map it to the solution vector
	auto & dbcs = meshContainer.getDirichletBoundaryConditions();
	auto & dbc = dbcs[0]; // FIXME! - Only considering one BC right now!
	temperature.ProjectBdrCoefficient(dbc.getCoefficient(),
			dbc.getBoundaryAttributes());
	// Set the heat conduction flux that results from a constant temperature
	// on the boundary, q = k*grad(u) where k is the thermal conductivity and
	// u is the temperature.
//	conductionCoeff.setConductivity(conductivity);
//	forcingVector.AddBoundaryIntegrator(
//			new BoundaryNormalLFIntegrator(conductionCoeff));
}

void ThermalOperator::Mult(const Vector &u, Vector &dudt) const {

	static int counter;

//	stiffnessMatrix->SpMat().Print();

	counter++;
	cout << "----- RK4 Step: " << counter << " -----" << endl;
	cout << "u = " << endl;
	u.Print();

	// y = M^{-1} (K x + b)
	K.Mult(u, z);
	//z must be negated because the bilinear form is actually on the right hand
	//side.
	z.Neg();
	// Add in b
	z += forcingVector;

	mSolver.Mult(z, dudt);

	cout << "dudt = " << endl;
	dudt.Print();

	return;
}

void ThermalOperator::ImplicitSolve(const double dt, const Vector &u, Vector &duDt) {

	// Allocate the temporary matrix to solve
	// du/dt = M^{-1}*(K(u+dt*du/dt) + b)
	if (!tempSparseMat) {
		tempSparseMat.reset(Add(1.0, massMatrix->SpMat(), dt,
				K));
		currentDt = dt;
		tempSolver.SetOperator(*tempSparseMat);
	}

	//Make sure the dt is correct to satisfy the SDIRK method
	MFEM_VERIFY(dt == currentDt, "");

	// Finish the solve
	K.Mult(u,z);
	z.Neg();
	z += b;
	tempSolver.Mult(z,duDt);

	return;
}

void ThermalOperator::update() {

	// Get the boundary conditions
	auto & dbcs = meshContainer.getDirichletBoundaryConditions();
	auto & dbc = dbcs[0]; // FIXME! - Only considering one BC right now!
	// Project the dirichlet boundary values
	temperature.ProjectBdrCoefficient(dbc.getCoefficient(),
			dbc.getBoundaryAttributes());

	// Form the linear system
	b = Vector();
	stiffnessMatrix->FormLinearSystem(dbc.getElements(),temperature,
			forcingVector,K,x,b,true);

	// Clear the temporary storage matrix so it is computed on the next step.
	tempSparseMat.reset(nullptr);

	return;
}

void ThermalOperator::recoverSolution() {
	stiffnessMatrix->RecoverFEMSolution(x,b,temperature);
}

ThermalOperator::~ThermalOperator() {

}

GridFunction & ThermalOperator::solution() {
	return temperature;
}

} /* namespace Kelvin */
