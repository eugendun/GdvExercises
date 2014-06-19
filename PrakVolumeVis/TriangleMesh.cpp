// ========================================================================= //
// Authors: Roman Getto, Matthias Bein                                       //
// mailto:roman.getto@gris.informatik.tu-darmstadt.de                        //
//                                                                           //
// GRIS - Graphisch Interaktive Systeme                                      //
// Technische Universit�t Darmstadt                                          //
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
		}
	}
	for (int i = 0; i < triangles.size(); i++)
		{
			file << "3 " << triangles[i].x << " " << triangles[i].y << " " << triangles[i].z << " \n";
		}
	
	file.close();
	}
