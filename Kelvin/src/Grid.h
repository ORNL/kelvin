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

#include <MaterialPoint.h>
#include <mfem.hpp>
#include <set>
#include <MassMatrix.h>
#include <MeshContainer.h>
#include <KelvinBaseTypes.h>
#include <functional>
#include <map>

namespace Kelvin {

/**
* This is the background Eulerian grid. It computes and stores the raw mesh
* (through the meshContainer), velocity, acceleration, mass, stress,
* strain, gradients of these quantities, etc.
*
* Nodes are stored as Points, which have position, velocity, and
* acceleration. Accessing these quantities should be done by pulling and
* indexing the list of nodes through the nodes() operation.
* @code
* auto & nodes = grid.nodes();
* auto & pos5 = nodes[4].pos;
* auto & vel5 = nodes[4].vel;
* auto & acc5 = nodes[4].acc;
* ...
* @endcode
*/
class Grid {

protected:

	/**
	 * Nodal positions
	 */
	std::vector<Point> _nodes;

	/**
	 * The list of massive nodes
	 */
	std::set<int> _nodeSet;

	/**
	 * The mesh container that holds the original finite element mesh.
	 */
	MeshContainer & _meshContainer;

	/**
	 * The shape/mapping matrix
	 */
	std::unique_ptr<mfem::SparseMatrix> _shapeMatrix;

	/**
	 * A map representing a sparse matrix holding the shape function gradients.
	 */
	std::map<int,std::vector<Gradient>> _gradientMap;

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

	/**
	 * The values of the internal forces applied at the grid nodes.
	 */
    std::vector<ForceVector> _internalForces;

	/**
	 * The values of the external forces applied at the grid nodes.
	 */
    std::vector<ForceVector> _externalForces;

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
	void assemble(const std::vector<Kelvin::MaterialPoint> & particles);

	/// Update kinematics, fields, etc.
	void update();

	/**
	 * Get the mass matrix associated with the grid. This operation computes
	 * the mass matrix from the particle distribution.
	 * @return the mass matrix
	 */
	MassMatrix & massMatrix(
			const std::vector<Kelvin::MaterialPoint> & particles);

	/**
	 * Get the mass matrix associated with the grid. This returns the mass
	 * matrix computed previously, or the default, empty matrix.
	 * @return the mass matrix
	 */
	MassMatrix & massMatrix();

	/**
	 * This operation returns the gradients of the nodal shape functions for
	 * the nodes that are near particles.
	 * @return a map representing the matrix of gradients (technically a rank
	 * 3 tensor). The map key represents the row id which is equal to the
	 * particle id in the particles list provided to assemble(), and the value
	 * set contains the Gradients and represents the columns assigned to that
	 * row.
	 */
	const std::map<int,std::vector<Gradient>> & gradients() const;

	/**
	 * This operation computes and returns the internal forces at the massive
	 * grid nodes. Any node that does not have mass assigned to it is excluded.
	 * @param particles the list of particles. This should be the same list as
	 * that passed to assemble().
	 * @return the internal forces
	 */
	const std::vector<ForceVector> & internalForces(
			const std::vector<Kelvin::MaterialPoint> & particles);

	/**
	 * This operation computes and returns the external forces at the massive
	 * grid nodes. Any node that does not have mass assigned to it is excluded.
	 * @param particles the list of particles. THis should be the same list as
	 * that passed to assemble().
	 */
	const std::vector<ForceVector> & externalForces(
			const std::vector<Kelvin::MaterialPoint> & particles);

	/**
	 * This operation returns the present kinematic information at the nodes
	 * as Points with positions, velocities, and accelerations. This includes
	 * the whole grid, not just the massive nodes.
	 * @return the kinematic information of the nodes
	 */
	const std::vector<Point> & nodes() const;

	/**
	 * This operation computes the acceleration at the nodes using the
	 * internal and external forces, and the mass matrix. The result is stored
	 * directly on the nodes() array. The base class implementation does not
	 * use the time step since the acceleration is computed directly from the
	 * equation of motion.
	 * @param timeStep the amount of time that has passed between the present
	 * and new values (which will be computed)
	 * @param particles the present particle configuration on the grid
	 */
	void updateNodalAccelerations(const double & timeStep,
			const std::vector<Kelvin::MaterialPoint> & particles);

	/**
	 * This operation updates the nodal velocities at the present time using
	 * momenta calculated from the particle configuration. The result is stored
	 * directly on the nodes() array. Only velocities at the massive nodes are
	 * updated. The mass matrix is not recomputed; the present mass matrix is
	 * used.
	 */
	void updateNodalVelocitiesFromMomenta(
			const std::vector<Kelvin::MaterialPoint> & particles);

	/**
	 * This operation computes the acceleration at the nodes using the
	 * internal and external forces, and the mass matrix. The result is stored
	 * directly on the nodes() array. Only velocities at the massive nodes are
	 * updated.
	 * @param timeStep the amount of time that has passed between the present
	 * and new values (which will be computed)
	 * @param particles the present particle configuration on the grid
	 */
	void updateNodalVelocities(const double & timeStep,
			const std::vector<Kelvin::MaterialPoint> & particles);

	/*
	 * This operation returns the set of node ids representing the nodes on the
	 * grid that have mass
	 * @return the set of node ids
	 */
	const std::set<int> & massiveNodeSet();

};

} /* namespace Kelvin */

#endif /* SRC_GRID_H_ */
