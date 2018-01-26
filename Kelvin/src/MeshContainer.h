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
#ifndef MESHCONTAINER_H_
#define MESHCONTAINER_H_

#include <IFESpaceFactory.h>
#include <DirichletBoundaryCondition.h>
#include <mfem.hpp>
#include <string>
#include <vector>

namespace Kelvin {

/**
 * This class is a simple container for MFEM meshes, finite element
 * collections, and finite element spaces. It makes it easier to create and
 * modify these classes without writing separate functions.
 *
 * This class also ties these elements to a property map that is injected up
 * construction. The valid key-value pairs are:
 * name=meshName
 * file=meshFileName
 * order=meshOrder
 */
class MeshContainer {
private:

	/**
	 * The name of the mesh file
	 */
	const std::string & meshFilename;

	/**
	 * The mesh
	 */
	mfem::Mesh mesh;

	/**
	 * The order of the mesh
	 */
	int _order;

	/**
	 * The dimension of the mesh
	 */
	int dim;

	/**
	 * The name of the mesh, which is a property that is different than its
	 * file name.
	 */
	std::string _name;

	/**
	 * The finite element collection for this mesh which is created by the
	 * space factory.
	 */
	mfem::FiniteElementCollection & feCollection;

	/**
	 * The finte element space created by the space factory.
	 */
	mfem::FiniteElementSpace & space;

	/**
	 * The Dirichlet (essential) boundary conditions defined on the mesh.
	 */
	std::vector<DirichletBoundaryCondition> dbConditions;

public:

	/**
	 * Constructor
	 * @param meshProps a map of the mesh properties
	 * @param spaceFactory a factory that creates a particular type of finite
	 * element collection and space for this mesh.
	 */
	MeshContainer(const std::map<std::string, std::string> & meshProps,
			IFESpaceFactory & spaceFactory);

	/**
	 * This operation returns the MFEM mesh object.
	 * @return mesh the mesh
	 */
	mfem::Mesh & getMesh();

	/**
	 * This operation returns the MFEM finite element space.
	 * @return space the space
	 */
	mfem::FiniteElementSpace & getSpace();

	/**
	 * This operation returns the order of the mesh
	 * @return the order
	 */
	int order();

	/**
	 * This operation returns the dimension of the mesh
	 * @return the dimension
	 */
	int dimension();

	/**
	 * This operation returns the name of the mesh
	 * @return name the name of the mesh
	 */
	const std::string & name();

	/**
	 * This operation sets a Dirichlet (essential) boundary condition
	 * on all sides - faces - of the mesh.
	 * @param value the constant value of the Dirichlet condition
	 */
	void setDirichletBoundaryCondition(double value);

	/**
	 * This operation sets a Dirichlet (essential) boundary condition
	 * on the mesh.
	 * @param meshSide the side of the mesh, which should be consistent with
	 * boundary attributes defined on the mesh
	 * @param value the constant value of the Dirichlet condition
	 */
	void setDirichletBoundaryCondition(int meshSide, double value);

	/**
	 * This operation return the set of Dirichlet (essential) boundary
	 * conditions defined on the mesh.
	 * @return the boundary conditions
	 */
	std::vector<DirichletBoundaryCondition> & getDirichletBoundaryConditions();
};

} /* namespace Kelvin */

#endif /* MESHCONTAINER_H_ */
