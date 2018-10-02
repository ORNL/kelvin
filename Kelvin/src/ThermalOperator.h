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
#ifndef THERMALOPERATOR_H_
#define THERMALOPERATOR_H_

#include <mfem.hpp>
#include <MeshContainer.h>
#include <memory>

namespace Kelvin {

/**
 * This coefficient computes the heat conduction flux, q = -k*grad(T) on the
 * boundary.
 */
class ConductionFluxCoefficient : public mfem::VectorCoefficient {

	mfem::GridFunction & u;

	double k;

public:
	ConductionFluxCoefficient(int _dim, mfem::GridFunction _u) :
		mfem::VectorCoefficient::VectorCoefficient(_dim), u(_u), k(1.0) {

	}

	virtual ~ConductionFluxCoefficient() {};

	virtual void setConductivity(double _k) {
		k = _k;
	}

	virtual void Eval(mfem::Vector & v, mfem::ElementTransformation & T,
			const mfem::IntegrationPoint & ip)  {
		// Set the size of the vector to match the dimensionality
	    v.SetSize(vdim);
	    // v = grad(u)
	    u.GetGradient(T,v);
	    // v = k*v = k*grad(u)
	    v *= -k;
	}
};

/**
 * This class represents the main operator for the thermal field. It is
 * time-dependent MFEM operator that otherwise behaves as all MFEM
 * TimeDependentOperator subclasses.
 */
class ThermalOperator: public mfem::TimeDependentOperator {
private:

	MeshContainer & meshContainer;
	mfem::FiniteElementSpace & feSpace;
	mfem::DataCollection & dataColl;

	std::unique_ptr<mfem::BilinearForm> massMatrix;
	std::unique_ptr<mfem::BilinearForm> stiffnessMatrix;
	mfem::LinearForm forcingVector;
	mfem::GridFunction temperature;
	mfem::ConstantCoefficient zero;

	/**
	 * Preconditioner for the mass matrix solve.
	 */
	mfem::DSmoother mPreconditioner;

	/**
	 * Solver for the mass matrix.
	 */
	mfem::CGSolver mSolver;

	/**
	 * A temporary matrix used for storing T = M + dt*K
	 */
	std::unique_ptr<mfem::SparseMatrix> tempSparseMat;

	/**
	 * Preconditioner for the temp matrix solve.
	 */
	mfem::DSmoother tempPreconditioner;

	/**
	 * Solver for the temp matrix.
	 */
	mfem::CGSolver tempSolver;

	/**
	 * Temporary storage vector.
	 */
	mutable mfem::Vector z;

	int skip_zeros = 0;
	double currentDt = 0.0;

	/**
	 * The thermal diffusivity
	 */
	double alpha;

	/**
	 * The density of the material
	 */
	double density;

	/**
	 * The specific heat of the material
	 */
	double specificHeat;

	/**
	 * The thermal conductivity of the material
	 */
	double conductivity;

	/**
	 * The initial temperature of the surface of the system.
	 */
	double initialSurfaceTemperature;

	/**
	 * The initial temperature of the interior of the system.
	 */
	double initialInteriorTemperature;

	/**
	 * A coefficient that computes the heat flux due to conduction at the
	 * surface of the material, q = -k*grad(T).
	 */
	ConductionFluxCoefficient conductionCoeff;

	mfem::SparseMatrix sparseMassMatrix;

	/**
	 * The matrix in the linear system formed from the stiffness matrix, Kx=b.
	 */
	mfem::SparseMatrix K;

	/**
	 * The RHS vector in the linear system formed from the stiffness matrix,
	 * Kx=b.
	 */
	mfem::Vector b;

	/**
	 * The unknown vector in the linear system formed from the stiffness
	 * matrix, Kx=b.
	 */
	mfem::Vector x;

	/**
	 * This operation configures the boundary conditions related to the surface
	 * temperature of the system. It maps the new boundary values to the final
	 * solution vector as well.
	 */
	void setSurfaceBoundaryConditions();

public:

	/**
	 * Constructor
	 */
	ThermalOperator(MeshContainer & _meshContainer,
			const std::map<std::string, std::string> & thermalProps,
			mfem::DataCollection & _dataColl);

	/**
	 * Destructor
	 */
	virtual ~ThermalOperator();

	void update();

	void recoverSolution();

	/** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
	 This is the only requirement for high-order SDIRK implicit integration.*/
	virtual void ImplicitSolve(const double dt, const mfem::Vector &u, mfem::Vector &k);

	mfem::GridFunction & solution();

	/**
	 * This operation applies the operator to x and returns the result in y,
	 * i.e. - y=Ax.
	 * @param x the vector to which the operator should be applied
	 * @param y the vector where the result should be stored.
	 */
	virtual void Mult(const mfem::Vector &x, mfem::Vector &y) const;

};

} /* namespace Kelvin */

#endif /* THERMALOPERATOR_H_ */
