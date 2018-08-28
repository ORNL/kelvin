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
#ifndef SRC_MASSMATRIX_H_
#define SRC_MASSMATRIX_H_

#include <mfem.hpp>
#include <vector>
#include <set>
#include <Point.h>
#include <functional>

namespace Kelvin {

/**
 * The mass matrix represents the mass of the system discretized on the
 * computational grid. A mass matrix element m_{ij} represents the mass shared
 * between nodes i and j.
 *
 * The mass can be retrieved as matrix elements m_{ij} or as a list of diagonal
 * entries on mass-lumped diagonalized matrix, m_{D,ij}. For descriptions of
 * both see Sulsky's 1994 paper "A particle method for history-dependent
 * materials."
 *
 * This is a sparse structure over all massive nodes where particles exist,
 * instead of all nodes regardless of whether or not they have particles
 * nearby. Taking advantage of the sparse matrix greatly increases the
 * performance.
 *
 */
class MassMatrix {

private:

	/**
	 * A private dummy to get around injecting references in the constructor.
	 */
	std::set<int> nodesDummy;

	/**
	 * A private dummy to get around injecting references in the constructor.
	 */
	mfem::SparseMatrix shapesDummy;

protected:

	/**
	 * The list of nodes with particles near by, i.e. - "massive nodes."
	 */
	std::set<int> & nodes;

	/**
	 * The full set of particles that give rise to mass in the grid.
	 */
	const std::vector<Point> & particles;

	/**
	 * The shape matrix that defines the relationship between nodes and
	 * particles. That is, the map of particles to grid nodes/cells.
	 *
	 * This is stored as a pointer because of a bug in MFEM that adversely
	 * affects storing it as a reference. This will be fully investigated
	 * in the future, but for now the pointer is fine. Note that this class
	 * does not free the pointer.
	 */
	mfem::SparseMatrix * shapes;

public:

	/**
	 * Constructor. The mass matrix must have a reference to the master list of
	 * particles for the system.
	 * @param particlesList the particles in the system
	 */
	MassMatrix(const std::vector<Point> & particleList);

	/**
	 * Destructor
	 */
	virtual ~MassMatrix();

	/**
	 * This operator assembles the mass matrix from a shape matrix and the list
	 * of nodes that have mass in the background mesh. This matrix is sparse
	 * and uses a sparse form of the shape matrix.
	 * @param shapeMatrix A sparse matrix that contains the shape at nodes in
	 * the background mesh
	 * @param nodeSet the list of nodes in the background mesh that actually have
	 * mass
	 */
	void assemble(mfem::SparseMatrix & shapeMatrix,
			std::set<int> & nodeSet);

	/**
	 * This operator accesses the element in the matrix at the i-th row and the
	 * j-th column.
	 *
	 * Since this matrix is sparse, matrix elements are computed on the fly and
	 * not stored. Elements that are outside the nodeset always have
	 * m_ij = 0.0.
	 *
	 * This method will indiscriminately compute m_ij, even if i and j are not
	 * in the node set, and does not exploit the sparse nature of the matrix.
	 *
	 * @param i row number
	 * @param j column number
	 * @return the mass element at row i and column j if it is non-zero,
	 * otherwise 0.0 exactly;
	 */
	double operator()(int i, int j);

	/**
	 * This operation creates a diagonalized form of the mass matrix by summing
	 * across the rows of the original mass matrix, which is called
	 * "mass lumping."
	 * @return the non-zero diagonal/lumped mass matrix entries, one for
	 * entry in the massive node set. The result is moved using C++11 move
	 * semantics. (It would be nice to have a by-input-reference version for
	 * high performance concerns about repeated allocations, etc.)
	 */
	std::vector<double> lump();

};

} /* namespace Kelvin */

#endif /* SRC_MASSMATRIX_H_ */
