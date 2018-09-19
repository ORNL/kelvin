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
#ifndef SRC_MFEMMPMDATA_H_
#define SRC_MFEMMPMDATA_H_

#include <MFEMData.h>
#include <vector>
#include <memory>
#include <Grid.h>

namespace Kelvin {

/**
 * This subclass of MFEMData adds support for particle data for the Material
 * Point Method.
 *
 * This classes adds a method to retrieve a list of particles. The particles
 * are read from a file that is passed in as input in a block of the simulation
 * input file. The format of the file must be csv with x,y,z coordinates.
 */
class MFEMMPMData: public MFEMData {
private:

	/**
	 * This is the background Eulerian grid and it includes the raw mesh
	 * (through the meshContainer), velocity, acceleration, mass, stress,
	 * strain, gradients of these quantities, etc.
	 */
	std::unique_ptr<Grid> _grid;

	/**
	 * This vector contains the point particles read from input and associated
	 * with the mesh stored in the mesh container of this Data instance.
	 */
	std::vector<MaterialPoint> _particles;

public:
	/**
	 * Constructor
	 */
	MFEMMPMData();

	/**
	 * Destructor
	 */
	virtual ~MFEMMPMData();

	/**
	 * This operation returns a set of points that represent Lagrangian
	 * material points.
	 * @return the set of Lagrangian Material Points.
	 */
	std::vector<MaterialPoint> & particles();

	/**
	 * This operation returns the computational grid that defines the Eulerian
	 * components of the MPM method.
	 * @return the grid
	 */
	Grid & grid();

	virtual void load(const std::string & inputFile);

};

} /* namespace Kelvin */

#endif /* SRC_MFEMMPMDATA_H_ */
