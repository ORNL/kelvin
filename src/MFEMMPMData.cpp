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
#include <MFEMMPMData.h>
#include <memory>
#include <vector>
#include <DelimitedTextParser.h>
#include <iostream>

using namespace std;
using namespace fire;

namespace Kelvin {

MFEMMPMData::MFEMMPMData() {
	// TODO Auto-generated constructor stub

}

MFEMMPMData::~MFEMMPMData() {
	// TODO Auto-generated destructor stub
}

void MFEMMPMData::load(const std::string & inputFile) {

	// Load the rest of the input data - mesh, etc. - first and then pull the
	// particle data
	MFEMData::load(inputFile);

	// Get the particles file
	auto & block = propertyParser.getPropertyBlock("particles");
	auto & particlesFile = block.at("file");
	// Load the particles
	DelimitedTextParser<vector<vector<double>>,double> parser(",","#");
	parser.setSource(particlesFile);
	parser.parse();
	shared_ptr<vector<vector<double>>> data = parser.getData();

	cout << "Loaded " << data->size() << " particles from "
			<< particlesFile << endl;

	// Convert to material points and pack the particles vector
	for (int i = 0; i < data->size(); i++) {
		auto & rawCoords = data->at(i);
		int numCoords = rawCoords.size();
		MaterialPoint point(rawCoords.size());
		for (int j = 0; j < numCoords; j++) {
			point.pos[j] = rawCoords[j];
		}
		_particles.push_back(point);
	}

	// Configure the data needed by the grid
	_grid = make_unique<Grid>(*mc);

	return;
}

Grid & MFEMMPMData::grid() {
	if (!loaded) throw "Data not loaded!";
	return *_grid;
}

std::vector<MaterialPoint> & MFEMMPMData::particles() {
	if (!loaded) throw "Data not loaded!";
	return _particles;
}

} /* namespace Kelvin */
