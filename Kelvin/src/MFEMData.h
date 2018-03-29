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
#ifndef SRC_MFEMDATA_H_
#define SRC_MFEMDATA_H_

#include <INIPropertyParser.h>
#include <mfem.hpp>
#include <memory>
#include <IFESpaceFactory.h>
#include <H1FESpaceFactory.h>
#include <MeshContainer.h>

namespace Kelvin {

/**
 * This is a simple data holder that collects all of the data used by MFEM and
 * the problem under consideration.
 */
class MFEMData {

	/**
	 * Parser for INI properties.
	 */
	fire::INIPropertyParser propertyParser;

	/**
	 * The mesh container that holds the mesh and finite element space.
	 */
	std::unique_ptr<MeshContainer> mc;

	/**
	 * The space factory that generates the finite element space.
	 */
	H1FESpaceFactory spaceFactory;

	/**
	 * The data collection that holds the final output.
	 */
	std::unique_ptr<mfem::DataCollection> dc;

public:

	/**
	 * Constructor
	 */
	MFEMData();

	/**
	 * Destructor
	 */
	virtual ~MFEMData();

	/**
	 * This operation loads the data based on the input file and its contents.
	 * @param inputFile the input file that points to other data and provides
	 * initil values.
	 */
	void load(const std::string & inputFile);

	/**
	 * This operation returns a reference to the INI property parser to load
	 * properties.
	 * @return the INI property parser
	 */
	fire::INIPropertyParser & properties();

	/**
	 * This operation returns the mesh container loaded in load().
	 * @return the mesh container
	 */
	MeshContainer & meshContainer();

	/**
	 * This operation returns the output data collection associated with this
	 * data set.
	 * @return the data collection
	 */
	mfem::DataCollection & collection();

};

} /* namespace Kelvin */

#endif /* SRC_MFEMDATA_H_ */
