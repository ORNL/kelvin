/**----------------------------------------------------------------------------
 Copyright  2018-, UT-Battelle, LLC
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
#ifndef SRC_GRID_H_
#define SRC_GRID_H_

#include <Point.h>
#include <mfem.hpp>
#include <set>
#include <MassMatrix.h>
#include <MeshContainer.h>

namespace Kelvin {

/**
* This is the background Eulerian grid. It computes and stores the raw mesh
* (through the meshContainer), velocity, acceleration, mass, stress,
* strain, gradients of these quantities, etc.
*/
class Grid {

protected:

	/**
	 * Nodal positions
	 */
	std::vector<Point> _pos;

	/**
	 * Nodal velocities
	 */
	std::vector<Point> _vel;

	/**
	 * Nodal acceleration
	 */
	std::vector<Point> _acc;

	/**
	 * The list of massive nodes
	 */
	std::set<int> nodeSet;

	/**
	 * The mesh container that holds the original finite element mesh.
	 */
	MeshContainer & _meshContainer;

	/**
	 * The shape/mapping matrix
	 */
	std::unique_ptr<mfem::SparseMatrix> _shapeMatrix;

	/**
	 * The mass matrix that shows the amount of mass shared between nodes due
	 * to the particles.
	 */
	std::unique_ptr<MassMatrix> _massMatrix;

	/**
	 * This private operation updates the state of the shape matrix.
	 */
	void updateShapeMatrix();

	/**
	 * This private operation updates the state of the mass matrix.
	 */
	void updateMassMatrix();

public:

	/**
	 * Constructor. Requires a mesh container as input to extract nodal and
	 * gradient information.
	 * @param meshContainer the finite element mesh used to create the grid
	 */
	Grid(MeshContainer & meshContainer);

	/**
	 * Destructor
	 */
	virtual ~Grid();

	/**
	 * Assemble all matrices. Rock and roll.
	 * @param particles the particle list that represent a mass distributed
	 * across the grid. Used to create shape and mass matrices.
	 */
	void assemble(const std::vector<Kelvin::Point> & particles);

	/// Update kinematics, fields, etc.
	void update();

	/**
	 * Get the mass matrix associated with the grid.
	 * @return the mass matrix
	 */
	const MassMatrix & massMatrix() const;

	/**
	 * This operation returns the present position of the grid nodes
	 * @return the position of the nodes, one for each massive node on the
	 * grid.
	 */
	const std::vector<Point> pos() const;

	/**
	 * This operation returns the present velocity at the grid nodes
	 * @return the velocity at the nodes, one for each massive node on the
	 * grid.
	 */
	const std::vector<Point> vel() const;

	/**
	 * This operation returns the present acceleration at the grid nodes
	 * @return the acceleration at the nodes, one for each massive node on the
	 * grid.
	 */
	const std::vector<Point> acc() const;

};

} /* namespace Kelvin */

#endif /* SRC_GRID_H_ */
