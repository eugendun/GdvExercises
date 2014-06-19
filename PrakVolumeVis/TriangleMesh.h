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
#ifndef TRIANGLEMESH_H
#define TRIANGLEMESH_H

#include <vector>
#include "Vec3.h"
#include <string>

using namespace std;

class TriangleMesh  {

public:

 // typedefs for data
  typedef Vec3i Triangle;
  typedef Vec3f Normal;
  typedef Vec3f Vertex;
  typedef vector<Triangle> Triangles;
  typedef vector<Normal> Normals;
  typedef vector<Vertex> Vertices;  

private:

  // data of TriangleMesh
  Vertices vertices;
  Normals normals;
  Triangles triangles;

public:

  // ===============================
  // === CONSTRUCTOR, DESTRUCTOR ===
  // ===============================

  TriangleMesh();
  ~TriangleMesh();

  // clears all data, sets defaults
  void clear();

  void addTriangle(Vec3f* t);

  // ================
  // === RAW DATA ===
  // ================

  // get raw data references
  Vertices& getVertices();
  Triangles& getTriangles();
  Normals& getNormals();

  void saveAsPly(std::string outputfilename);

};


#endif

