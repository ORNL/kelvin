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
#include <MFEMOlevskyLVCR.h>
#include <mfem.hpp>
#include <MeshContainer.h>
#include <StringCaster.h>

using namespace mfem;
using namespace std;

namespace Kelvin {

MFEMOlevskyLVCR::MFEMOlevskyLVCR(MFEMData & data) : _data(data), porosity(0.0),
		shearModulus(0.0), phi(0.0), psi(0.0), density(1.0) {

	// Get the properties
	auto & propertiesParser = data.properties();
	auto & properties = propertiesParser.getPropertyBlock("material");
	// Store the porosity and shear modulus
	porosity = fire::StringCaster<double>::cast(properties.at("porosity"));
	shearModulus = fire::StringCaster<double>::cast(
			properties.at("shearModulus"));
	density = fire::StringCaster<double>::cast(properties.at("density"));
	// Compute phi and psi for the stress updates
	phi = (1.0-porosity)*(1.0-porosity);
	psi = (2.0/3.0)*phi*(1-porosity)/porosity;

	return;
}

MFEMOlevskyLVCR::~MFEMOlevskyLVCR() {
	// TODO Auto-generated destructor stub
}

void MFEMOlevskyLVCR::updateStrainRate(const Kelvin::Grid & grid,
		std::vector<Kelvin::MaterialPoint> & matPoints) {

	// Create the H1 field
	int dim = _data.meshContainer().dimension();
	auto & meshContainer = _data.meshContainer();
	auto & _mesh = meshContainer.getMesh();
	H1_FECollection velCol(1,dim);
	FiniteElementSpace velSpace(&_mesh,&velCol,dim,Ordering::byVDIM);
	// Create and fill the grid function
	GridFunction velGf(&velSpace);
	auto & nodes = grid.nodes();
	for (int i = 0; i < nodes.size(); i++) {
		auto & nodeVel = nodes[i].vel;
		for (int j = 0; j < dim; j++) {
			velGf[i*dim+j] = nodeVel[j];
		}
	}

	// Compute the strain rate and strain. Is this correct? I get the
	// impression that the share gradients are off.
	int numParticles = matPoints.size();
	for (int i = 0; i < numParticles; i++) {
		auto & point = matPoints[i];
		// Find the element that contains the point
		mfem::Array<int> elementId(1);
		mfem::Array<mfem::IntegrationPoint> intPoints(1);
		auto pointMatrix = meshContainer.convertPointToMatrix(point.pos);
	    _mesh.FindPoints(pointMatrix,elementId,intPoints);
	    // Get the element transformation
		auto * elemTrans = _mesh.GetElementTransformation(elementId[0]);
		// Set the material point position in reference coordinates
		elemTrans->SetIntPoint(&intPoints[0]);
		// Get the gradient of the velocity at the material point
		DenseMatrix gradVel(dim,dim);
		velGf.GetVectorGradient(*elemTrans,gradVel);
		// Compute the strain using infinitesimal strain theory by
		// symmetrizing the matrix.
		gradVel.Symmetrize();
		// Map the matrix into the point
		for (int j = 0; j < dim; j++) {
			for (int k = 0; k < dim; k++) {
				point.strain[j][k] = gradVel(j,k);
			}
		}
	}

}

void MFEMOlevskyLVCR::updateStress(const Kelvin::Grid & grid,
		std::vector<Kelvin::MaterialPoint> & matPoints) {

	int dim = grid.dimension();
	double twoShearMod = 2.0*shearModulus;
	double deviatoricStrainRate[dim][dim];
	double sinteringStress[dim];

	// Compute the sintering stress, which is constant for now
	for (int i = 0; i < dim; i++) {
		sinteringStress[i] = phi*(2.0*(1.0-porosity)-(1-porosity))/porosity;
	}

	int numPoints = matPoints.size();
	for (int i = 0; i < numPoints; i++) {
		double traceE = 0.0;
		auto & matPoint = matPoints[i];
		// Copy the initial strain rate into the deviatoric matrix
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				deviatoricStrainRate[i][j] = matPoint.strain[i][j];
			}
		}
		// Compute the trace of the strain rate
		for (int i = 0; i < dim; i++) {
			traceE += matPoint.strain[i][i];
		}
		// Compute the hydrostatic strain rate
		double hydrostaticStrainRate = traceE/dim;
		// Compute the deviatoric strain rate
		for (int i = 0; i < dim; i++) {
			deviatoricStrainRate[i][i] -= hydrostaticStrainRate;
		}
		// Compute the stress matrix
		for (int i = 0; i < dim; i++) {
			// "Initialize" the stress matrix with a scaled deviatoric matrix
			for (int j = 0; j < dim; j++) {
				matPoint.stress[i][j] = twoShearMod*phi*
						deviatoricStrainRate[i][j]/density;
			}
			// Add in the diagonal components
			matPoint.stress[i][i] = (twoShearMod*psi*traceE
					+ sinteringStress[i])/density;
		}
	}

	// One loop

	// Setup the sintering stress
	// Two loops

}

} /* namespace Kelvin */
