
#ifndef __meshgen_h__
#define __meshgen_h__

#include "c3dlas.h"


typedef struct {
	Vector3* vertices;
	Vector3* normals;
	short* indices;
	int nElements;
	int nVertices;
	int nNormals;
	int nIndices;
	short width; // should be plenty. 8 billion triangles it too much to render anyway.
	short height;
} IndexedPatch;






typedef struct TriangleMesh {
	Vector3* vertices;
	Vector3* normals;
	unsigned int* indices;
	
	// used quantity
	int vertexCnt; // also used for normals.
	int indexCnt;
	// always GL_TRIANGLES
	
	// allocated sizes, in elements
	int szVertices;
	int szIndices;
	
	char windingDir; // c = clockwise, w = widdershins (counter-clockwise) 
} TriangleMesh;


typedef struct MeshSlice {
	
//	MeshVertex* vertices; 
	short* indices; 
	// always line strips; -1 resets
	
	int vertexCnt;
	int indexCnt;
	
	// allocated sizes, in elements
	int szVertices;
	int szIndices;
	
	Vector2 texDxDy;
} MeshSlice;



// TODO: winding and strips. currently it outputs a triangle list and the winding is probably fucked.
// in the x-y plane, 0,0 at the lower left
void mgGenFlatPatch(short width, short height, IndexedPatch* out);



float* genNoisef(short width, short height, float min, float max);


TriangleMesh* makeCube(Matrix* mat, int flat);


TriangleMesh* makeQuad(float width, float height);
TriangleMesh* makeCylinder(int sides, float height, float baseRadius, float topRadius);
TriangleMesh* makeIcosahedron(float radius);
TriangleMesh* makeSubbedIcosahedron(float radius, int subdivisions);


// uncapped cylinders
// discs
// free triangles
// tetrahedrons / cones
// shapes with planar but subdivided faces
// mesh smoothing fn
// tex coord projection

TriangleMesh* makeParallelpiped();





#endif // __meshgen_h__
