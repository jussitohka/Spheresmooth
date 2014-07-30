// Sphere prior shape model based smoothing of polygon meshes
// Copyright (C) Jussi Tohka, Institute of Signal Processing, Tampere University of
// Technology, 2004 - 2007.
// P.O. Box 553, FIN-33101, Finland
//  E-mail: jussi.tohka@tut.fi
// ----------------------------------------------------------------------------
// Permission to use, copy, modify, and distribute this software 
// for any purpose and without fee is hereby
// granted, provided that the above copyright notice appear in all
// copies.  The author and Tampere University of Technology make no representations
// about the suitability of this software for any purpose.  It is
// provided "as is" without express or implied warranty.
//  ----------------------------------------------------------------------------

// -----------------------------------------------------------------
// The method is described in 
// J. Tohka. 
// Surface Smoothing Based on a Sphere Shape Model. 
// In Proc. of 6th Nordic Signal Processing Symposium, NORSIG2004, 
// pp. 17 - 20, Espoo, Finland, June, 2004.
// ------------------------------------------------------------


#ifndef SSMOOTH_H
#define SSMOOTH_H

typedef double points[3];
typedef int intpoints[3];

#define MAX_FACE_ELEMENTS 40
#define PI 3.141592653589793
// This macro computes the number of common nodes in two triangles
#define COMMON_NODES(n1, n2, n3, m1,  m2, m3) ((((n1) == (m1)) || ((n1) == (m2)) || ((n1) == (m3))) + (((n2) == (m1)) || ((n2) == (m2)) || ((n2) == (m3))) +  (((n3) == (m1)) || ((n3) == (m2)) || ((n3) == (m3))))
#define MEAN3(n1,n2,n3) (0.3333333333*(n1 + n2 + n3))


typedef struct
{
  char no_elem; // number of nodes in the face
  int elem[MAX_FACE_ELEMENTS];     // pointer to elements  
} face;


typedef struct 
{  
  int vertices;         // number of vertices
  int triangles;        // number of triangles
  double props[5];      // properties of the surface
  points* x;            // coordinates of vertices
  intpoints* tri;     // triangles;
  points* normals;    // surface normals
  int props2[5];       // Color of polygons. Only constant color is supported and 
  int* triprops;           // I don't know what for these are. (The name triprops does not refer to anything...) 
} trimesh;     

typedef struct
{  
  int vertices;  // number of vertices
  int no_faces;     // number of faces
  points* x;     // coordinates of vertices
  intpoints* nerel;   // neighbourhood relations
  face* faces;   // faces of the mesh
} simplexmesh;  


  

// Main functions

// Reads a triangular mesh file to memory. Works with MNI style surface files. 
// Returns 0 if everything is ok.
int ReadTriMesh(char* filename, trimesh* input_mesh);
// Writes triangular mesh files.
// Returns 0 if ok.
int WriteTriMesh(char* filename, trimesh* output_mesh);
// Converts a triangular mesh into simplex. Quadratic time method. 
// Useful if meshes processed are very small.
int Tri2Simplex(trimesh* input_mesh,simplexmesh* out_mesh); 
// Converts a triangular mesh into simplex mesh. Linear time method.
int Tri2SimplexLinear(trimesh* input_mesh,simplexmesh* out_mesh);
// Converts a simplex mesh to a triangular mesh. 
// No implemented because this is not needed for surface smoothing.
trimesh* Simplex2Tri(simplexmesh* inputmesh);
// the sphere smoothing algorithm
int SphereSmooth(simplexmesh* input_mesh,double lambda, int max_iter, double multiplier);
// Changes the geometry of the triangular mesh according to a dual simplex mesh.
int Simplex2TriGeom(simplexmesh* input_smesh,trimesh* input_trimesh);

// Helper functions

// Computes the stopping threshold for the steepest descent surface smoothing.
double ComputeThreshold(simplexmesh* input_mesh,double multiplier);

// Does the same as the macro earlier on. 
// The macro is obviously faster, so this is mainly for debugging purposes 
int CommonNodes(int n1,int n2,int n3,int m1,int m2,int m3);
// This function is not required for surface smoothing.
void RemoveDoubleNeighbours(face* pface);


#endif
