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
#ifndef TIMEEVOLUTIONOPERATOR_H_
#define TIMEEVOLUTIONOPERATOR_H_

#include <mfem.hpp>

namespace Kelvin {

/**
 * This is a simple time-dependent operator for solving the ODE:
 * du/dt = (M^{-1})(Ku + b)
 * where M is mass matrix, K is the stiffness matrix and b is
 * the forcing vector.
 */
class TimeEvolutionOperator : public mfem::TimeDependentOperator {
private:

	/**
	 * Sparse mass matrix
	 */
	mfem::SparseMatrix & M;

	/**
	 * Spare stiffness matrix
	 */
	mfem::SparseMatrix & K;

	/**
	 * Forcing vector
	 */
	const mfem::Vector &b;

	/**
	 * Preconditioner for the mass matrix solve.
	 */
	mfem::DSmoother mPreconditioner;

	/**
	 * Solver for the mass matrix.
	 */
	mfem::CGSolver mSolver;

	/**
	 * Temperator storage vector.
	 */
	mutable mfem::Vector z;

public:

	/**
	 * Constructor
	 * @param _M The mass matrix in sparse matrix form
	 * @param _K The stiffness matrix in sparse matrix form
	 * @param _b the forcing vector (normally denoted "b")
	 */
	TimeEvolutionOperator(mfem::SparseMatrix &_M, mfem::SparseMatrix &_K,
			const mfem::Vector &_b);

	/**
	 * This operation applies the operator to x and returns the result in y,
	 * i.e. - y=Ax.
	 * @param x the vector to which the operator should be applied
	 * @param y the vector where the result should be stored.
	 */
	virtual void Mult(const mfem::Vector &x, mfem::Vector &y) const;

	virtual ~TimeEvolutionOperator();
};

} /* namespace Kelvin */

#endif /* TIMEEVOLUTIONOPERATOR_H_ */
