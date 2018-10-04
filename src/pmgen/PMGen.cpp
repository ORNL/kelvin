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
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <mfem.hpp>
#include <MeshContainer.h>
#include <H1FESpaceFactory.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>
#include <vector>
#include <Point.h>

using namespace Kelvin;
using namespace mfem;
using namespace std;

/**
 * This operation retrieves the points from the mesh container.
 * @param meshContainer the mesh container for the input sample mesh
 * @return a vector of particle positions
 */
vector<Kelvin::Point> getPoints(MeshContainer & meshContainer) {

	// Get the quadrature points to form the particle mesh.
	auto points = meshContainer.getQuadraturePoints();

	// Seems simple for a particle picker, but using function was desirable
	// for adding more things once the points are picked.

	return points;
}

/**
 * This operation creates a 2D quadrilateral mesh.
 * @param numBoxes the number of boxes along each side
 * @param sideLengths the side lengths of the boxes by dimension
 * @return the mesh
 */
inline Mesh * createQuadMesh(const double * numBoxes,
		const double * sideLengths) {
	return new Mesh(numBoxes[0], numBoxes[1],
			Element::Type::QUADRILATERAL, 1, sideLengths[0],
			sideLengths[1]);
}

/**
 * This operation creates a 3D hexahedral mesh.
 * @param numBoxes the number of boxes along each side
 * @param sideLengths the side lengths of the boxes by dimension
 * @return the mesh
 */
inline Mesh * createHexMesh(const double * numBoxes,
		const double * sideLengths) {
	return new Mesh(numBoxes[0], numBoxes[1], numBoxes[2],
			Element::Type::HEXAHEDRON, 1, sideLengths[0],
			sideLengths[1],	sideLengths[2]);
}

/**
 * This operation checks the bounding box and modifies it to address various
 * constraints, such as squaring a rectangular bounding box to optimize
 * particle placement with respect to grid lines.
 * @param min the minimum boundaries
 * @param max the maximum boundaries
 */
void checkBoundingBox(const int & dim, Vector & min, Vector & max) {

	int maxIndex = 0;
	double length[dim], maxLength = 0.0;

	// Compute the initial side lengths. Also find the maximum side length and
	// its index. Do a linear search since there are at most 3 entries.
	for (int i = 0; i < dim; i++) {
		length[i] = max(i) - min(i);
		if (length[i] > maxLength) {
			maxLength = length[i];
			maxIndex = i;
		}
	}

	// Square the box boundaries by shifting each maximum value by its
	// difference from the absolute max. Now the box will have equal side
	// lengths.
	for (int i = 0; i < dim; i++) {
		if (i != maxIndex) {
	       max(i) += length[maxIndex] - length[i];
		}
	}

	// FIXME! - This will square the bounding box, but it will have unequal
	// padding on one side.

	return;
}

/**
 * This function creates a reference mesh and writes it to a VTK file. This
 * function also shifts the points, if needed, to suit the position of the
 * new mesh.
 *
 * The reference mesh is either hexahedral or quadrilateral for 3D and 2D
 * meshes respectively. Elements are created using the average volume (or area
 * in 2D) of elements from the original mesh. The original mesh is padded by a
 * volume (or area) 2x greater than the original mesh. The new element volume
 * can be scaled using a multiplicative factor, (coarsenFactor).
 *
 * @param meshContainer the mesh container that contains the sample input mesh
 * @param points the points/particles picked for the sample input mesh
 * @param filename the name of the output file to which the reference mesh
 * should be written
 * @param coarsenFactor an optional multiplicative factor
 */
void createReferenceMesh(MeshContainer & meshContainer,
		vector<Kelvin::Point> & points, const char * filename,
		double coarsenFactor = 1.0) {

	// Try to create a hexahedral reference mesh.
	auto & mesh = meshContainer.getMesh();

	// Compute the average element volume of the source mesh
	int dim = meshContainer.dimension();
	int numElem = mesh.GetNE();
	double vol = 0.0;
	for (int i = 0; i < numElem; i++) {
		vol += mesh.GetElementVolume(i);
	}
	// Apply the coarsening factor when computing the average element volume.
	vol = coarsenFactor*(vol/numElem);
	// Compute the side length for a box with this volume
	double exponent = 1.0 / ((double) dim);
	double lElemSide = pow(vol, exponent);

	// Retrieve the bounds of the source mesh
	Vector min, max;
	mesh.GetBoundingBox(min, max);
	// Check the bounding box and do things like make it square/cubical, etc.
	checkBoundingBox(dim, min, max);
	double refMeshSideLengths[dim], refMeshCenter[dim], numBoxes[dim];
	// Compute the side lengths of the source bounding box. Store these
	// coordinates as the new center of the reference mesh. Double the original
	// side lengths so that the new bounding box is twice as large as that of
	// the source mesh to allow for particle movement. Also compute the number
	// of boxes on each side.
	//
	// This works by just shifting the center from the original to (lx,ly,lz).
	for (int i = 0; i < dim; i++) {
		refMeshCenter[i] = max(i) - min(i);
		refMeshSideLengths[i] = 2.0 * refMeshCenter[i];
		numBoxes[i] = refMeshSideLengths[i] / lElemSide;
	}

	// Create the new reference mesh. Using a pointer as a reference. Sometimes
	// the things I have to do for MFEM drive me nuts!
	Mesh * refMesh;
	if (dim == 2) {
		refMesh = createQuadMesh(numBoxes,refMeshSideLengths);
	} else if (dim == 3) {
		refMesh = createHexMesh(numBoxes,refMeshSideLengths);
	}
	// Write the mesh
	ofstream backgroundMeshFile;
	backgroundMeshFile.open(filename);
	refMesh->PrintVTK(backgroundMeshFile);
	// Close the mesh file
	backgroundMeshFile.close();
	// Dump the mesh
	delete refMesh;

	// Shift the particles. First, determine if the original mesh was centered
	// around the origin or another point and adjust by halving the center
	// point. This may not work for all meshes, but if covers meshes centered
	// at the origin and (lx/2,ly/2,lz/2).
	for (int i = 0; i < dim; i++) {
		if (min(i) >= 0.0) refMeshCenter[i] /= 2.0;
	}

	// For a source mesh with a bounding box centered on
	// the origin, this shift will move the center to the new center of the
	// reference mesh.
	int numPoints = points.size();
	for (int i = 0; i < numPoints; i++) {
		for (int j = 0; j < dim; j++) {
			points[i].pos[j] += refMeshCenter[j];
		}
	}

	return;
}

/**
 * This function extracts particle data from a vector and writes it to the
 * particle output file. The mesh file name is used as a reference in the
 * header (metadata) of the particle file.
 *
 * The particles are written in comma separated variables (CSV) format with a
 * short header at the top containing useful metadata.
 *
 * @param points a vector of the particles positions
 * @param particleFilename the name of the output file container particle info.
 * @param meshFilename the name of the mesh file loaded into the mesh
 * container.
 */
void writeParticles(vector<Kelvin::Point> points, const char * particleFilename,
		const char * meshFilename, const int & matId) {

	// Open the particle output file
	ofstream particleFile;
	particleFile.open(particleFilename);

	// Get the date and time in UTC
	time_t currentTime = time(0);
	tm *utc = gmtime(&currentTime);
	char * dt = asctime(utc);

	// Write the header
	particleFile << "# " << particleFilename << " generated by PMGen" << endl;
	particleFile << "# Number of particles = " << points.size() << endl;
	particleFile << "# Source mesh: " << meshFilename << endl;
	// Note that dt has a line break in it.
	particleFile << "# Created on: UTC " << dt;
	particleFile << "# x, y, z" << endl;

	// Write the particle coordinates
	auto size = points.size();
	// Assume the points are all the same dimension
	int dim = points[0].dimension();
	for (int i = 0; i < size; i++) {
		// No need for a loop over coords[j] since dim = 2 or dim = 3 and
		// branch prediction should handle any performance drop.
		if (dim == 2) {
		   particleFile << points[i].pos[0] << ", " << points[i].pos[1];
		} else if (dim == 3){
		   particleFile << points[i].pos[0] << ", " << points[i].pos[1]
					    << ", " << points[i].pos[2];
		}
		particleFile << ", " << matId << endl;
	}

	// Close the particle file
	particleFile.close();

	return;
}

/**
 * Main program
 * @param argc the number of input arguments
 * @param argv the input arguments array of argc elements
 * @return EXIT_SUCCESS if successful, otherwise another value.
 */
int main(int argc, char * argv[]) {

	const char *meshFilename = "";
	const char *backgroundMeshFilename = "background.vtk";
	const char *particleFilename = "particles.csv";
	int matId = 1;
	int order = 1;
	double coarsenFactor = 1.0;

	// Create the default command line arguments
	OptionsParser args(argc, argv);
	args.AddOption(&meshFilename, "-m", "--mesh",
			"Mesh file for sampling particle positions.");
	args.AddOption(&order, "-o", "--order",
			"Finite element order (polynomial degree) or -1 for"
					" isoparametric space.");
	args.AddOption(&coarsenFactor, "-c", "--coarsening-factor",
			"Factor for coarsening the average reference mesh element volume.");
	args.AddOption(&backgroundMeshFilename, "-b", "--background-mesh",
			"Name for the background mesh output file generated by this program.");
	args.AddOption(&particleFilename, "-p", "--particle-set",
			"Name of the particle set output file generated by this program.");
	args.AddOption(&matId, "-i", "--materialID",
			"Value of the material id that should be set by this program.");

	// Parse the arguments and do a cursory check.
	args.Parse();
	if (!args.Good()) {
		args.PrintUsage(cout);
		return EXIT_FAILURE;
	}

	// Try to load the mesh
	H1FESpaceFactory spaceFactory;
	MeshContainer meshContainer(meshFilename, order, spaceFactory);
	int dim = meshContainer.dimension();
	// Drop out if the dimensionality is unacceptable.
	if (dim != 2 && dim != 3) {
		// Only two and three dimensional meshes are supported.
		cout << "Only 2D and 3D meshes are supported. This mesh has dim = "
				<< dim << endl;
		EXIT_FAILURE;
	}

	// Get the particles
	auto points = getPoints(meshContainer);

	// Create the reference mesh. The points vector is sent along in case the
	// points need to be shifted when the new mesh is created.
	createReferenceMesh(meshContainer, points, backgroundMeshFilename,
			coarsenFactor);

	// Write particle list
	writeParticles(points, particleFilename, meshFilename, matId);

	return EXIT_SUCCESS;
}
