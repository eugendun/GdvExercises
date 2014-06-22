// ========================================================================= //
// Authors: Roman Getto, Matthias Bein                                       //
// mailto:roman.getto@gris.informatik.tu-darmstadt.de                        //
//                                                                           //
// GRIS - Graphisch Interaktive Systeme                                      //
// Technische Universität Darmstadt                                          //
// Fraunhoferstrasse 5                                                       //
// D-64283 Darmstadt, Germany                                                //
//                                                                           //
// Changed: 06.06.2014														 //
// Creation Date: 19.05.2012                                                 //
// Content: Simple class for triangle meshes							     //
// ==========================================================================//

#include "TriangleMesh.h"
#include <iostream>
#include <fstream>


// ===============================
// === CONSTRUCTOR, DESTRUCTOR ===
// ===============================

TriangleMesh::TriangleMesh() {
  clear();
}

TriangleMesh::~TriangleMesh() {
  clear();
}

void TriangleMesh::clear() {
  // clear mesh data
  vertices.clear();
  triangles.clear();
  normals.clear();
}

void TriangleMesh::addTriangle(Vec3f* p) {
	addTriangle(p, new Vec3f());
}

void TriangleMesh::addTriangle(Vec3f* p, Vec3f* n) {

	// REMARKS: BONUS EXERCISE (b)
	// return, if degenerated triangle
	if (((p[0] - p[1]).sqlength() < 0.00001) || ((p[0] - p[2]).sqlength() < 0.00001) || ((p[1] - p[2]).sqlength() < 0.00001)) {
		return;
	}


	Vec3i indices = Vec3i(-1, -1, -1);

	for (int j = 0; j < 3; j++) {

		// look, if this vertex has already been saved
		// REMARKS: BONUS EXERCISE (a) - deactivated, because too slow for debugging
		/*
		for (int i = 0; i < vertices.size(); i++) {
			if ((abs(vertices[i].x - p[j].x) < 0.0001) && (abs(vertices[i].y - p[j].y) < 0.0001) && (abs(vertices[i].z - p[j].z) < 0.0001)) {
				indices[j] = i;
				break;
			}
		}
		*/

		// save, if not found
		if (indices[j] < 0) {
			vertices.push_back(p[j]);

			if (n->sqlength() > 0.0001) {
				normals.push_back(n[j]);
			}

			indices[j] = vertices.size() - 1;
		}
	}
	
	triangles.push_back(indices);
}

// ================
// === RAW DATA ===
// ================

TriangleMesh::Vertices& TriangleMesh::getVertices() {  
	return vertices;
}
TriangleMesh::Triangles&  TriangleMesh::getTriangles() {
	return triangles;
}

TriangleMesh::Normals&  TriangleMesh::getNormals() {
  return normals;
}

  void TriangleMesh::saveAsPly(std::string outputfilename){
	  if(triangles.size() == 0){
		  return;
	  }
	std::ofstream file;
	outputfilename += ".ply";
	file.open (outputfilename);
	file << "ply \n" << "format ascii 1.0 \n" <<
		"element vertex " << vertices.size() << "\n"
		<< "property float x \n" 
		<< "property float y \n"
		<< "property float z \n";
	if(normals.size() != 0){
	file << "property float nx \n"
		<< "property float ny \n"
		<< "property float nz \n";
	}
	file << "element face " << triangles.size() << " \n";
	file << "property list uchar int vertex_indices \n" << "end_header \n";
	if(normals.size() != 0){
		for (int i = 0; i < vertices.size(); i++)
		{
			file << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z  << " "
				<< normals[i].x << " " << normals[i].y << " " << normals[i].z << " \n";
		}
	}
	else{
		for (int i = 0; i < vertices.size(); i++)
		{
			file << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z << " \n";
			if (i % 100 == 0) {
				std::cout << "vertex " << i << " / " << vertices.size() << std::endl;
			}
		}
	}
	for (int i = 0; i < triangles.size(); i++)
		{
			file << "3 " << triangles[i].x << " " << triangles[i].y << " " << triangles[i].z << " \n";
			if (i % 100 == 0) {
				std::cout << "tri " << i << " / " << triangles.size() << std::endl;
			}
		}
	
	file.close();
	}
