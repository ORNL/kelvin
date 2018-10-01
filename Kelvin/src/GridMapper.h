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
#ifndef SRC_GRIDMAPPER_H_
#define SRC_GRIDMAPPER_H_

#include <Grid.h>

namespace Kelvin {

/**
 * This is a base class for a utility that maps grid fields to particles.
 * Separating these functions from the grid allows for greater variability in
 * implementation without modifying or subclassing the grid.
 */
class GridMapper {
public:

	/**
	 * Destructor
	 */
	virtual ~GridMapper() {};

	/**
	 * This function maps the acceleration on the grid nodes to particle
	 * positions. This operation updates the particle arrays in-place.
	 * @param grid the grid where nodal accelerations are stored
	 * @param particles the material points for which accelerations should be
	 * calculated
	 */
	virtual void updateParticleAccelerations(const Kelvin::Grid & grid,
			std::vector<Kelvin::MaterialPoint> & particles) const = 0;

	/**
	 * This function maps the velocity on the grid nodes to particle
	 * positions. This operation updates the particle arrays in-place.
	 * @param grid the grid where nodal velocities are stored
	 * @param particles the material points for which accelerations should be
	 * calculated
	 */
	virtual void updateParticleVelocities(const Kelvin::Grid & grid,
			std::vector<Kelvin::MaterialPoint> & particles) const = 0;

	/**
	 * This function maps the velocity on the grid nodes to particle
	 * positions. This operation provides the updates to the particle arrays
	 * in an out of place vector.
	 * @param grid the grid where nodal velocities are stored
	 * @param particles the material points for which accelerations should be
	 * calculated
	 * @param velocities a vector in which the velocities should be stored with
	 * total length equal to the number of particles by the number of
	 * dimensions
	 */
	virtual void updateParticleVelocities(const Kelvin::Grid & grid,
			const std::vector<Kelvin::MaterialPoint> & particles,
			std::vector<double> & velocities) const = 0;
};

} /* namespace Kelvin */

#endif /* SRC_GRIDMAPPER_H_ */
