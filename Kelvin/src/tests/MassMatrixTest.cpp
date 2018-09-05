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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE kelvin

#include <boost/test/included/unit_test.hpp>
#include <mfem.hpp>
#include <vector>
#include <MassMatrix.h>

using namespace std;
using namespace mfem;
using namespace Kelvin;

/**
 * This operation checks the basic functionality of the MassMatrix class.
 */
BOOST_AUTO_TEST_CASE(checkMass) {

	// Set of all nodes currently near particles
	set<int> nodeSet;
	set<int>::iterator outerIt;
	set<int>::iterator innerIt;

	// This test uses a particle mesh with four particles along the diagonal
	// at half integer steps. The full 4x4 mesh would have 16 elements and 25
	// vertices, but only the elements that have particles are created here.

	// Create some particles
	int numParticles = 4;
	vector<Kelvin::Point> particles(numParticles);
	for (int i = 0; i < numParticles; i++) {
		particles[i].pos[0] = 0.5 + ((double) i);
		particles[i].pos[1] = 0.5 + ((double) i);
		particles[i].pos[2] = 0.5 + ((double) i);
		std::cout << particles[i].pos[0] << " " << particles[i].pos[1]
				<< " " << particles[i].pos[2] << endl;
	}

	// Setup the shape matrix
	int numVerts = 25;
	SparseMatrix shapeMatrix(numParticles,numVerts);
	// Particle 1
	shapeMatrix.Add(0,0,0.25);
	nodeSet.insert(0);
	shapeMatrix.Add(0,1,0.25);
	nodeSet.insert(1);
	shapeMatrix.Add(0,5,0.25);
	nodeSet.insert(5);
	shapeMatrix.Add(0,6,0.25);
	nodeSet.insert(6);
	// Particle 2
	shapeMatrix.Add(1,6,0.25);
	shapeMatrix.Add(1,7,0.25);
	nodeSet.insert(7);
	shapeMatrix.Add(1,11,0.25);
	nodeSet.insert(11);
	shapeMatrix.Add(1,12,0.25);
	nodeSet.insert(12);
	// Particle 3
	shapeMatrix.Add(2,12,0.25);
	shapeMatrix.Add(2,13,0.25);
	nodeSet.insert(13);
	shapeMatrix.Add(2,17,0.25);
	nodeSet.insert(17);
	shapeMatrix.Add(2,18,0.25);
	nodeSet.insert(18);
	// Particle 4
	shapeMatrix.Add(3,18,0.25);
	shapeMatrix.Add(3,19,0.25);
	nodeSet.insert(19);
	shapeMatrix.Add(3,23,0.25);
	nodeSet.insert(23);
	shapeMatrix.Add(3,24,0.25);
	nodeSet.insert(24);

	// Finalize the matrix and put its columns in order from least to greatest.
	shapeMatrix.Finalize();
	shapeMatrix.SortColumnIndices();

	// Debug
	shapeMatrix.Print();

	// Create the mass matrix
	MassMatrix massMatrix(particles);
	massMatrix.assemble(shapeMatrix,nodeSet);

	// Construct the mass matrix associated with the grid nodes
	double particleMass = 1.0; // FIXME! Needs to be something real and from input.

	// Test some of the expected values of the matrix. Not worth it to test all.
	// The important thing to note here is that expectedValue1 is the expected
	// mass if one particle is near by, and expectedValued2 is is the expected
	// value if two particles are nearby.
	double expectedValue1 = 0.0625, expectedValue2 = 0.125, percentEPS = 1.0e-15;
	// 0 - One nearby particles
	BOOST_REQUIRE_CLOSE(expectedValue1,massMatrix(0,0),percentEPS);
	BOOST_REQUIRE_CLOSE(expectedValue1,massMatrix(0,1),percentEPS);
	BOOST_REQUIRE_CLOSE(expectedValue1,massMatrix(0,5),percentEPS);
	BOOST_REQUIRE_CLOSE(expectedValue1,massMatrix(0,6),percentEPS);
	// 12 - Two nearby particles nearby for 12 coupled to itself, but only one
	// for the others.
	BOOST_REQUIRE_CLOSE(expectedValue2,massMatrix(12,12),percentEPS);
	BOOST_REQUIRE_CLOSE(expectedValue1,massMatrix(12,13),percentEPS);
	BOOST_REQUIRE_CLOSE(expectedValue1,massMatrix(12,17),percentEPS);
	BOOST_REQUIRE_CLOSE(expectedValue1,massMatrix(12,18),percentEPS);
	// 20 - No nearby particles, not in the massive node set, m_ij = 0
	BOOST_REQUIRE_CLOSE(0.0,massMatrix(20,20),percentEPS);
	BOOST_REQUIRE_CLOSE(0.0,massMatrix(20,21),percentEPS);
	BOOST_REQUIRE_CLOSE(0.0,massMatrix(20,15),percentEPS);
	BOOST_REQUIRE_CLOSE(0.0,massMatrix(20,16),percentEPS);

	// Print the values at the massive nodes for debugging.
	for (outerIt = nodeSet.begin(); outerIt != nodeSet.end(); outerIt++) {
		for(innerIt = nodeSet.begin(); innerIt != nodeSet.end(); innerIt++) {
			cout << "(" << *outerIt << ", " << *innerIt << ") "
					<< massMatrix(*outerIt,*innerIt) << ", ";
		}
		cout << endl;
	}

	// Check the diagonalized/lumped form
	auto lumpedMassMatrix = massMatrix.lump();
	outerIt = nodeSet.begin();
	for (int i = 0; i < lumpedMassMatrix.size(); i++) {
		double expectedValue = 0.0;
		if (*outerIt != 6 && *outerIt != 12 && *outerIt != 18) {
			expectedValue = 4.0*expectedValue1;
		} else {
			expectedValue = 4.0*expectedValue2;
		}
		BOOST_REQUIRE_CLOSE(expectedValue,lumpedMassMatrix[i],percentEPS);
		cout << *outerIt << " " << lumpedMassMatrix[i] << endl;
		outerIt++;
	}

	return;

}
