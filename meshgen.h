
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





typedef struct MeshVertex {
	Vector3 v, n;
	struct {
		unsigned short u, v;
	} t;
} MeshVertex;

typedef struct Mesh {
	MeshVertex* vertices;
	unsigned short* indices;
	
	// used quantity
	int vertexCnt;
	int indexCnt;
	// always GL_TRIANGLES
	
	// allocated sizes, in elements
	int szVertices;
	int szIndices;
	
	
} Mesh;


typedef struct MeshSlice {
	
	MeshVertex* vertices; 
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


Mesh* makeCube(Matrix* mat, int flat);








#endif // __meshgen_h__