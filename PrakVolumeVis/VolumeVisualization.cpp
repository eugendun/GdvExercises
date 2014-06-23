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

#include "VolumeVisualization.h"
#include "mcLookUp.h"
#include "TriangleMesh.h"
#include "dualpointslist.h"
#include <map>

#define SNAP_EPSILON 0.1

void VolumeVisualization::loadRAW(std::istream& in, int dimX, int dimY, int dimZ, float dx, float dy, float dz) {
	spacing.x = dx; spacing.y = dy; spacing.z = dz; 
	dimension.x = dimX; dimension.y = dimY; dimension.z = dimZ; 
	volumedata.clear();
	volumedata.resize(dimX*dimY*dimZ);
	// we know the bit depth (8 bits by now), so all that's left is to cast each char to a float between [0,1) 
	for (int i = 0; i < dimX*dimY*dimZ; ++i) {	// TODO: include spacing: "/spacing"
		volumedata[i] = ((unsigned char)in.get())/256.0f;
	}

	// calculate a normal for each density value via central differences
	volumenormals.clear();
	volumenormals.resize(dimX*dimY*dimZ);
	for (int z = 0; z < dimension.z; z += spacing.z)	{
		for (int y = 0; y < dimension.y; y += spacing.y)	{
			for (int x = 0; x < dimension.x; x += spacing.x)	{
				float x0 = (x > 0) ? volumedata.at(z*dimension.x*dimension.y + y*dimension.x + x - 1) : 0;
				float x1 = (x < dimension.x - 1)? volumedata.at(z*dimension.x*dimension.y + y*dimension.x + x + 1) : 0;

				float y0 = (y > 0) ? volumedata.at(z*dimension.x*dimension.y + (y-1)*dimension.x + x) : 0;
				float y1 = (y < dimension.y - 1) ? volumedata.at(z*dimension.x*dimension.y + (y+1)*dimension.x + x) : 0;

				float z0 = (z > 0) ? volumedata.at((z-1)*dimension.x*dimension.y + y*dimension.x + x) : 0;
				float z1 = (z < dimension.z - 1) ? volumedata.at((z+1)*dimension.x*dimension.y + y*dimension.x + x) : 0;

				Vec3f normal = Vec3f(x1 - x0, y1 - y0, z1 - z0);
				normal *= -1;
				volumenormals[z*dimension.x*dimension.y + y*dimension.x + x] = normal.normalized();
			}
		}
	}
}

void VolumeVisualization::loadTrivariateFunction(int dimX, int dimY, int dimZ, float dx, float dy, float dz) {
	//spacing.x = dx; spacing.y = dy; spacing.z = dz;
	spacing.x = 1; spacing.y = 1; spacing.z = 1;
	dimension.x = dimX / dx; dimension.y = dimY / dy; dimension.z = dimZ / dz;
	volumedata.clear();
	volumedata.resize(dimension.x*dimension.y*dimension.z / (spacing.x*spacing.y*spacing.z));

	// calculate trivariate function
	for (int z = 0; z < dimension.z; z += spacing.z)	{
		for (int y = 0; y < dimension.y; y += spacing.y)	{
			for (int x = 0; x < dimension.x; x += spacing.x)	{
				volumedata[indexForCoordinates(x, y, z)] = evaluateTrivariateFunction(x * dx, y * dy, z * dz);
			}
		}
	}

	// calculate a normal for each density value via central differences
	volumenormals.clear();
	/*
	volumenormals.resize(dimension.x*dimension.y*dimension.z / (spacing.x*spacing.y*spacing.z));
	for (int z = 0; z < dimension.z; z += spacing.z)	{
		for (int y = 0; y < dimension.y; y += spacing.y)	{
			for (int x = 0; x < dimension.x; x += spacing.x)	{
				float x0 = (x > 0) ? volumedata.at(z*dimension.x*dimension.y + y*dimension.x + x - 1) : 0;
				float x1 = (x < dimension.x - 1) ? volumedata.at(z*dimension.x*dimension.y + y*dimension.x + x + 1) : 0;

				float y0 = (y > 0) ? volumedata.at(z*dimension.x*dimension.y + (y - 1)*dimension.x + x) : 0;
				float y1 = (y < dimension.y - 1) ? volumedata.at(z*dimension.x*dimension.y + (y + 1)*dimension.x + x) : 0;

				float z0 = (z > 0) ? volumedata.at((z - 1)*dimension.x*dimension.y + y*dimension.x + x) : 0;
				float z1 = (z < dimension.z - 1) ? volumedata.at((z + 1)*dimension.x*dimension.y + y*dimension.x + x) : 0;

				Vec3f normal = Vec3f(x1 - x0, y1 - y0, z1 - z0);
				normal *= -1;
				volumenormals[z*dimension.x*dimension.y + y*dimension.x + x] = normal.normalized();
			}
		}
	}*/
}

float VolumeVisualization::evaluateTrivariateFunction(float x, float y, float z) {
	float r = 10;
	float xm = fmodf(x + 10000, 2 * r) - r;
	float ym = fmodf(y + 10000, 2 * r) - r;
	float zm = fmodf(z + 10000, 2 * r) - r;
	float result = xm*xm + ym*ym + zm*zm - r*r;
	return -result;
}

Vec3f VolumeVisualization::evaluateTrivariateFunctionNormal(float x, float y, float z) {
	float r = 10;
	float xm = fmodf(x + 10000, 2 * r) - r;
	float ym = fmodf(y + 10000, 2 * r) - r;
	float zm = fmodf(z + 10000, 2 * r) - r;
	
	return Vec3f(xm, ym, zm).normalized();
}

float VolumeVisualization::indexForCoordinates(float x, float y, float z) {
	return (z/spacing.z) * (dimension.x/spacing.x) * (dimension.y/spacing.y) + (y/spacing.y)*(dimension.x/spacing.x) + (x/spacing.x);
}

void VolumeVisualization::computeMesh(float isovalue)	{
	MC_TRIANGLE triangles[8];
	for (int z = 0; z < dimension.z - spacing.z; z += spacing.z)	{
		for (int y = 0; y < dimension.y - spacing.y; y += spacing.y)	{
			for (int x = 0; x < dimension.x - spacing.x; x += spacing.x)	{
				volumedata.at(indexForCoordinates(x, y, z));
				GRIDCELL cell;

				// set position for each corner
				cell.p[0] = Vec3f(x,             y, z);
				cell.p[1] = Vec3f(x + spacing.x, y, z);
				cell.p[2] = Vec3f(x + spacing.x, y, z + spacing.z);
				cell.p[3] = Vec3f(x,             y, z + spacing.z);

				cell.p[4] = Vec3f(x,             y + spacing.y, z);
				cell.p[5] = Vec3f(x + spacing.x, y + spacing.y, z);
				cell.p[6] = Vec3f(x + spacing.x, y + spacing.y, z + spacing.z);
				cell.p[7] = Vec3f(x,             y + spacing.y, z + spacing.z);

				// fill with density value and normal for each corner
				for (int i = 0; i < 8; i++) {
					int dataIndex = indexForCoordinates(cell.p[i].x/spacing.x, cell.p[i].y/spacing.y, cell.p[i].z/spacing.z);
					cell.val[i] = volumedata.at(dataIndex);
					if (volumenormals.size() > 0) {
						cell.n[i] = volumenormals.at(dataIndex);
					}
				}

				int numTriangles = Polygonise(cell, isovalue, triangles);
				for (int i = 0; i < numTriangles; i++) {
					mesh.addTriangle(triangles[i].p, triangles[i].g);
				}
			}
		}

		std::cout << z << " / " << dimension.z << std::endl;
	}
}

void VolumeVisualization::computeMeshDMC(float isovalue)	{
	//std::map<Vec3f, DMC_DUAL_POINTS> dualPointEdgeMap;
	Vec3f* dualPoints[12];

	for (int z = 0; z < dimension.z - spacing.z; z += spacing.z)	{
		for (int y = 0; y < dimension.y - spacing.y; y += spacing.y)	{
			for (int x = 0; x < dimension.x - spacing.x; x += spacing.x)	{
				volumedata.at(indexForCoordinates(x, y, z));
				GRIDCELL cell;

				// set position for each corner
				cell.p[0] = Vec3f(x, y, z);
				cell.p[1] = Vec3f(x + spacing.x, y, z);
				cell.p[2] = Vec3f(x + spacing.x, y, z + spacing.z);
				cell.p[3] = Vec3f(x, y, z + spacing.z);

				cell.p[4] = Vec3f(x, y + spacing.y, z);
				cell.p[5] = Vec3f(x + spacing.x, y + spacing.y, z);
				cell.p[6] = Vec3f(x + spacing.x, y + spacing.y, z + spacing.z);
				cell.p[7] = Vec3f(x, y + spacing.y, z + spacing.z);

				// fill with density value for each corner
				for (int i = 0; i < 8; i++) {
					int dataIndex = indexForCoordinates(cell.p[i].x / spacing.x, cell.p[i].y / spacing.y, cell.p[i].z / spacing.z);
					cell.val[i] = volumedata.at(dataIndex);
				}

				// clear dual points for this cell
				for (int i = 0; i < 12; i++) {
					dualPoints[i] = NULL;
				}

				PolygoniseDMC(cell, isovalue, dualPoints);

				//----------------------
				// for i in dualPoints
				//   if dualPoints[i] != null
				//     v0 = grid.p[nodesForEdge(i)[0]]
				//     v1 = grid.p[nodesForEdge(i)[1]]
				//     edgeMap.put((v0 - v1)/2, dualPoints[i]);
				//----------------------
				for (int i = 0; i < 12; i++) {
					if (dualPoints[i] != NULL) {
						Vec3f v0 = cell.p[nodesForEdge[i][0]];
						Vec3f v1 = cell.p[nodesForEdge[i][1]];
						Vec3f key = (v0 - v1) / 2;
						Vec3f dualPoint = *dualPoints[i];
						//dualPointEdgeMap[key].push_back(dualPoint);
					}
				}
			}
		}

		std::cout << z << " / " << dimension.z << std::endl;
	}

	/*
	std::map<Vec3f, DMC_DUAL_POINTS>::iterator dual;
	for (dual = dualPointEdgeMap.begin(); dual != dualPointEdgeMap.end(); ++dual) {
		for (int i = 0; i < dual->second.size(); i++) {
			// create two triangles out of these four points
			dual->second[i];
		}
	}
	*/
}

unsigned int VolumeVisualization::Polygonise(GRIDCELL grid,float isolevel,MC_TRIANGLE *triangles) {
  int i, ntriang, cubeindex;
  Vec3f vertlist[12];
  Vec3f normallist[12];
  // Determine the index into the edge table which tells us which vertices are inside of the surface
  cubeindex = 0;
  if (grid.val[0] < isolevel) cubeindex |= 1;
  if (grid.val[1] < isolevel) cubeindex |= 2;
  if (grid.val[2] < isolevel) cubeindex |= 4;
  if (grid.val[3] < isolevel) cubeindex |= 8;
  if (grid.val[4] < isolevel) cubeindex |= 16;
  if (grid.val[5] < isolevel) cubeindex |= 32;
  if (grid.val[6] < isolevel) cubeindex |= 64;
  if (grid.val[7] < isolevel) cubeindex |= 128;
  // Cube is entirely in/out of the surface
  if (edgeTable[cubeindex] == 0)
    return(0);
  // Find the edges where the surface intersects the cube
  if (edgeTable[cubeindex] & 1) {
	  vertlist[0]   = VertexInterp(isolevel, grid.p[0], grid.p[1], grid.val[0], grid.val[1], SNAP_EPSILON);
	  normallist[0] = VertexInterp(isolevel, grid.n[0], grid.n[1], grid.val[0], grid.val[1]);
  }
  if (edgeTable[cubeindex] & 2) {
	  vertlist[1] = VertexInterp(isolevel, grid.p[1], grid.p[2], grid.val[1], grid.val[2], SNAP_EPSILON);
	  normallist[1] = VertexInterp(isolevel, grid.n[1], grid.n[2], grid.val[1], grid.val[2]);
  }
  if (edgeTable[cubeindex] & 4) {
	  vertlist[2] = VertexInterp(isolevel, grid.p[2], grid.p[3], grid.val[2], grid.val[3], SNAP_EPSILON);
	  normallist[2] = VertexInterp(isolevel, grid.n[2], grid.n[3], grid.val[2], grid.val[3]);
  }
  if (edgeTable[cubeindex] & 8) {
	  vertlist[3] = VertexInterp(isolevel, grid.p[3], grid.p[0], grid.val[3], grid.val[0], SNAP_EPSILON);
	  normallist[3] = VertexInterp(isolevel, grid.n[3], grid.n[0], grid.val[3], grid.val[0]);
  }
  if (edgeTable[cubeindex] & 16) {
	  vertlist[4] = VertexInterp(isolevel, grid.p[4], grid.p[5], grid.val[4], grid.val[5], SNAP_EPSILON);
	  normallist[4] = VertexInterp(isolevel, grid.n[4], grid.n[5], grid.val[4], grid.val[5]);
  }
  if (edgeTable[cubeindex] & 32) {
	  vertlist[5] = VertexInterp(isolevel, grid.p[5], grid.p[6], grid.val[5], grid.val[6], SNAP_EPSILON);
	  normallist[5] = VertexInterp(isolevel, grid.n[5], grid.n[6], grid.val[5], grid.val[6]);
  }
  if (edgeTable[cubeindex] & 64) {
	  vertlist[6] = VertexInterp(isolevel, grid.p[6], grid.p[7], grid.val[6], grid.val[7], SNAP_EPSILON);
	  normallist[6] = VertexInterp(isolevel, grid.n[6], grid.n[7], grid.val[6], grid.val[7]);
  }
  if (edgeTable[cubeindex] & 128) {
	  vertlist[7] = VertexInterp(isolevel, grid.p[7], grid.p[4], grid.val[7], grid.val[4], SNAP_EPSILON);
	  normallist[7] = VertexInterp(isolevel, grid.n[7], grid.n[4], grid.val[7], grid.val[4]);
  }
  if (edgeTable[cubeindex] & 256) {
	  vertlist[8] = VertexInterp(isolevel, grid.p[0], grid.p[4], grid.val[0], grid.val[4], SNAP_EPSILON);
	  normallist[8] = VertexInterp(isolevel, grid.n[0], grid.n[4], grid.val[0], grid.val[4]);
  }
  if (edgeTable[cubeindex] & 512) {
	  vertlist[9] = VertexInterp(isolevel, grid.p[1], grid.p[5], grid.val[1], grid.val[5], SNAP_EPSILON);
	  normallist[9] = VertexInterp(isolevel, grid.n[1], grid.n[5], grid.val[1], grid.val[5]);
  }
  if (edgeTable[cubeindex] & 1024) {
	  vertlist[10] = VertexInterp(isolevel, grid.p[2], grid.p[6], grid.val[2], grid.val[6], SNAP_EPSILON);
	  normallist[10] = VertexInterp(isolevel, grid.n[2], grid.n[6], grid.val[2], grid.val[6]);
  }
  if (edgeTable[cubeindex] & 2048) {
	  vertlist[11] = VertexInterp(isolevel, grid.p[3], grid.p[7], grid.val[3], grid.val[7], SNAP_EPSILON);
	  normallist[11] = VertexInterp(isolevel, grid.n[3], grid.n[7], grid.val[3], grid.val[7]);
  }
  // Create the triangles
  ntriang = 0;
  for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
    triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i  ]];
    triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i+1]];
	triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i+2]];

	if (volumenormals.size() > 0) {
		triangles[ntriang].g[0] = normallist[triTable[cubeindex][i]];
		triangles[ntriang].g[1] = normallist[triTable[cubeindex][i + 1]];
		triangles[ntriang].g[2] = normallist[triTable[cubeindex][i + 2]];
	} else {
		triangles[ntriang].g[0] = evaluateTrivariateFunctionNormal(triangles[ntriang].p[0].x, triangles[ntriang].p[0].y, triangles[ntriang].p[0].z);
		triangles[ntriang].g[1] = evaluateTrivariateFunctionNormal(triangles[ntriang].p[1].x, triangles[ntriang].p[1].y, triangles[ntriang].p[1].z);
		triangles[ntriang].g[2] = evaluateTrivariateFunctionNormal(triangles[ntriang].p[2].x, triangles[ntriang].p[2].y, triangles[ntriang].p[2].z);
	}

	ntriang++;
  }
  return(ntriang);
}


void VolumeVisualization::PolygoniseDMC(GRIDCELL grid, float isolevel, Vec3f** dualPoints) {


	int i, ntriang, cubeindex;
	Vec3f vertlist[12];
	Vec3f normallist[12];
	// Determine the index into the edge table which tells us which vertices are inside of the surface
	cubeindex = 0;
	if (grid.val[0] < isolevel) cubeindex |= 1;
	if (grid.val[1] < isolevel) cubeindex |= 2;
	if (grid.val[2] < isolevel) cubeindex |= 4;
	if (grid.val[3] < isolevel) cubeindex |= 8;
	if (grid.val[4] < isolevel) cubeindex |= 16;
	if (grid.val[5] < isolevel) cubeindex |= 32;
	if (grid.val[6] < isolevel) cubeindex |= 64;
	if (grid.val[7] < isolevel) cubeindex |= 128;
	// Cube is entirely in/out of the surface
	if (edgeTable[cubeindex] == 0)
		return;
	// Find the edges where the surface intersects the cube
	if (edgeTable[cubeindex] & 1) {
		vertlist[0] = VertexInterp(isolevel, grid.p[0], grid.p[1], grid.val[0], grid.val[1], SNAP_EPSILON);
		normallist[0] = VertexInterp(isolevel, grid.n[0], grid.n[1], grid.val[0], grid.val[1]);
	}
	if (edgeTable[cubeindex] & 2) {
		vertlist[1] = VertexInterp(isolevel, grid.p[1], grid.p[2], grid.val[1], grid.val[2], SNAP_EPSILON);
		normallist[1] = VertexInterp(isolevel, grid.n[1], grid.n[2], grid.val[1], grid.val[2]);
	}
	if (edgeTable[cubeindex] & 4) {
		vertlist[2] = VertexInterp(isolevel, grid.p[2], grid.p[3], grid.val[2], grid.val[3], SNAP_EPSILON);
		normallist[2] = VertexInterp(isolevel, grid.n[2], grid.n[3], grid.val[2], grid.val[3]);
	}
	if (edgeTable[cubeindex] & 8) {
		vertlist[3] = VertexInterp(isolevel, grid.p[3], grid.p[0], grid.val[3], grid.val[0], SNAP_EPSILON);
		normallist[3] = VertexInterp(isolevel, grid.n[3], grid.n[0], grid.val[3], grid.val[0]);
	}
	if (edgeTable[cubeindex] & 16) {
		vertlist[4] = VertexInterp(isolevel, grid.p[4], grid.p[5], grid.val[4], grid.val[5], SNAP_EPSILON);
		normallist[4] = VertexInterp(isolevel, grid.n[4], grid.n[5], grid.val[4], grid.val[5]);
	}
	if (edgeTable[cubeindex] & 32) {
		vertlist[5] = VertexInterp(isolevel, grid.p[5], grid.p[6], grid.val[5], grid.val[6], SNAP_EPSILON);
		normallist[5] = VertexInterp(isolevel, grid.n[5], grid.n[6], grid.val[5], grid.val[6]);
	}
	if (edgeTable[cubeindex] & 64) {
		vertlist[6] = VertexInterp(isolevel, grid.p[6], grid.p[7], grid.val[6], grid.val[7], SNAP_EPSILON);
		normallist[6] = VertexInterp(isolevel, grid.n[6], grid.n[7], grid.val[6], grid.val[7]);
	}
	if (edgeTable[cubeindex] & 128) {
		vertlist[7] = VertexInterp(isolevel, grid.p[7], grid.p[4], grid.val[7], grid.val[4], SNAP_EPSILON);
		normallist[7] = VertexInterp(isolevel, grid.n[7], grid.n[4], grid.val[7], grid.val[4]);
	}
	if (edgeTable[cubeindex] & 256) {
		vertlist[8] = VertexInterp(isolevel, grid.p[0], grid.p[4], grid.val[0], grid.val[4], SNAP_EPSILON);
		normallist[8] = VertexInterp(isolevel, grid.n[0], grid.n[4], grid.val[0], grid.val[4]);
	}
	if (edgeTable[cubeindex] & 512) {
		vertlist[9] = VertexInterp(isolevel, grid.p[1], grid.p[5], grid.val[1], grid.val[5], SNAP_EPSILON);
		normallist[9] = VertexInterp(isolevel, grid.n[1], grid.n[5], grid.val[1], grid.val[5]);
	}
	if (edgeTable[cubeindex] & 1024) {
		vertlist[10] = VertexInterp(isolevel, grid.p[2], grid.p[6], grid.val[2], grid.val[6], SNAP_EPSILON);
		normallist[10] = VertexInterp(isolevel, grid.n[2], grid.n[6], grid.val[2], grid.val[6]);
	}
	if (edgeTable[cubeindex] & 2048) {
		vertlist[11] = VertexInterp(isolevel, grid.p[3], grid.p[7], grid.val[3], grid.val[7], SNAP_EPSILON);
		normallist[11] = VertexInterp(isolevel, grid.n[3], grid.n[7], grid.val[3], grid.val[7]);
	}
	// Create the triangles
	ntriang = 0;
	MC_TRIANGLE triangles[8];
	for (i = 0; triTable[cubeindex][i] != -1; i += 3) {
		triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i]];
		triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i + 1]];
		triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i + 2]];

		ntriang++;
	}

	// create a dual point for each patch
	for (int i = 0; dualPointsList[cubeindex][i] != 0; i++) {
		//dualPointsList[cubeindex][i]
	}
}

/*
Algorithmus fuer DMC:
---------------------
Laufe wie bei MC ueber alle Zellen;
	Berechne Dreiecke fuer Zelle[i];
	Unterteile Dreiecke in Patches (mittels vorkommender Kanten in dualPointsList)
	Berechne DualPoint pro Patch;
	Speichere den DualPoint fuer jede Kante in dualPointsList;
Jetzt hat man eine Liste (Map), die vier DualPoints fuer jede relevante Kante enthaelt.
Erstelle fuer jede Kante aus den vier Punkten ein Viereck (Punkte mit groesstem Abstand werden mittlere Kante, danach flip)
*/

Vec3f VolumeVisualization::generateDualPoint(MC_TRIANGLE* triangles, int indices) {
	return Vec3f();
}

Vec3f VolumeVisualization::VertexInterp(float isolevel, Vec3f p1, Vec3f p2, float valp1, float valp2, float snapEpsilon) {
  if (abs(isolevel-valp1) < 0.00001)
    return(p1);
  if (abs(isolevel-valp2) < 0.00001)
    return(p2);
  if (abs(valp1-valp2) < 0.00001)
    return(p1);
  float mu;
  Vec3f p;
  mu = (isolevel - valp1) / (valp2 - valp1);
  p.x = p1.x + mu * (p2.x - p1.x);
  p.y = p1.y + mu * (p2.y - p1.y);
  p.z = p1.z + mu * (p2.z - p1.z);

  // REMARKS: BONUS EXERCISE (b)
  float snapEpsilonSqr = snapEpsilon * snapEpsilon;
  if (((p - p1).sqlength() < snapEpsilonSqr) || ((p - p2).sqlength() < snapEpsilonSqr)) {
	  if (mu < 0.5) {
		  return p1;
	  } else {
		  return p2;
	  }
  }
  return(p);
}