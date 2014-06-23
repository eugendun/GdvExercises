// ========================================================================= //
// Authors: Roman Getto, Matthias Bein	                                     //
// mailto:roman.getto@gris.informatik.tu-darmstadt.de                        //
//                                                                           //
// GRIS - Graphisch Interaktive Systeme                                      //
// Technische Universität Darmstadt                                          //
// Fraunhoferstrasse 5                                                       //
// D-64283 Darmstadt, Germany                                                //
//                                                                           //
// Creation Date: 06.06.2014                                                 //
// Content: Simple class for triangle meshes							     //
// ==========================================================================//

#ifndef VOLUMEVISUALIZATION_INCLUDED
#define VOLUMEVISUALIZATION_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include "Vec3.h"
#include "TriangleMesh.h"

class VolumeVisualization {

	typedef struct {
		Vec3f p[3];
		Vec3f g[3];
	} MC_TRIANGLE;

	typedef struct {
		Vec3f p[8];
		Vec3f n[8];
		float val[8];
	} GRIDCELL;

	typedef std::vector<float> VolumeData;
	typedef std::vector<Vec3f> VolumeNormals;

	// data
	VolumeData volumedata;
	VolumeNormals volumenormals;
	Vec3i dimension;
	Vec3f spacing;

	//mesh
	TriangleMesh mesh;

	// lookup tables for all 256 mc cases
	static const int edgeTable[256];
	static const int triTable[256][16];

public:

	//getter (returning Pointer):

	Vec3i* getDimension(){
		return &dimension;
	}

	Vec3f* getSpacing(){
		return &spacing;
	}

	std::vector<float>* getVolumeData(){
		return &volumedata;
	}

	TriangleMesh* getMesh(){
		return &mesh;
	}

	//

	void loadRAW(std::istream& in, int dimX, int dimY, int dimZ, float dx=1, float dy=1, float dz=1);
	void loadTrivariateFunction(int dimX, int dimY, int dimZ, float dx = 1, float dy = 1, float dz = 1);

	void computeMesh(float isovalue);
	void computeMeshDMC(float isovalue);


	// Given a grid cell and an isolevel, calculate the triangular
	//  facets required to represent the isosurface through the cell.
	// Return the number of triangular facets, the array "triangles"
	//  will be loaded up with the vertices at most 5 triangular facets.
	//  0 will be returned if the grid cell is either totally above
	//  of totally below the isolevel.
	unsigned int Polygonise(GRIDCELL grid,float isolevel,MC_TRIANGLE *triangles);

	/**
	 * @param dualPoints 12 vectors, i.e. one dual point (or null) for each edge
	 */
	void PolygoniseDMC(GRIDCELL grid, float isolevel, Vec3f** dualPoints);
	Vec3f generateDualPoint(MC_TRIANGLE* triangles, int indices);

	// Linearly interpolate the position where an isosurface cuts an
	//  edge between two vertices, each with their own scalar value
	Vec3f VertexInterp(float isolevel, Vec3f p1, Vec3f p2, float valp1, float valp2, float snapEpsilon = 0.0);

	float indexForCoordinates(float x, float y, float z);

private:

	float evaluateTrivariateFunction(float x, float y, float z);
	Vec3f evaluateTrivariateFunctionNormal(float x, float y, float z);
};

#endif