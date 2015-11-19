
#include <stdlib.h>

#include "meshgen.h"




// in the x-y plane, 0,0 at the lower left
void mgGenFlatPatch(short width, short height, IndexedPatch* out) {
	
	int cnt; // programming is more fun if you mentally read this as "cunt"
	int x, y;
	int ii; // index index

	
	cnt = (width + 1) * (height + 1); // plus one for the edge
	out->nElements = width * height * 2;
	out->nVertices = out->nNormals = cnt; 
	out->nIndices = out->nElements * 3;
	
	out->width = width;
	out->height = height;
	
	out->vertices = (Vector*)malloc(cnt * sizeof(Vector));
	out->normals = (Vector*)malloc(cnt * sizeof(Vector));
	out->indices = (short*)malloc(cnt * sizeof(short) * 6); // too big but we'll fix that later
	
	// fill in the vertices and normals
	for(y = 0; y < height + 1; y++) {
		for(x = 0; x < width + 1; x++) {
			Vector* v, *n;
			v = &out->vertices[x+(y*(width+1))];
			n = &out->normals[x+(y*(width+1))];
			
			v->x = x; // generate it flat
			v->y = y;
			v->z = 0.0;
			
			n->x = 0.0;
			n->y = 0.0;
			n->z = 1.0;
		}
	}
	
	// fill in the indices
	ii = 0;
	for(y = 0; y < height; y++) {
		for(x = 0; x < width; x++) {
			out->indices[ii++] =  x + (y * (width + 1));
			out->indices[ii++] =  (x + 1) + (y * (width + 1));
			out->indices[ii++] =  (x + 1) + ((y + 1) * (width + 1));
	
			out->indices[ii++] =  x + (y * (width + 1));
			out->indices[ii++] =  (x + 1) + ((y + 1) * (width + 1));
			out->indices[ii++] =  x + ((y + 1) * (width + 1));
		}
	}
	
}





float* genNoisef(short width, short height, float min, float max) {
	
	int len, i;
	float* out, range;
	
	len = width * height;
	range = max - min;
	
	out = (float*)malloc(len * sizeof(float));
	
	// no need for super quality math stuff here. it's just randomish junk for textures.
	for(i = 0; i < len; i++) 
		out[i] = (((float)rand() / RAND_MAX) * range) + min;
	
	return out;
}






Mesh* extrudeAlongVector(Vector* lineStrip, int lineCount, Vector* v) {
	
	int i;
	Mesh* m;
	
	m = malloc(sizeof(Mesh));
	
	// step 1: allocate enough memory for the mesh
	m->vertexCnt = (lineCount + 1) * 2;
	m->indexCnt = 3 * vertexCnt;
	
	m->vertices = malloc(sizeof(Vector) * vertexCnt);
	m->indices = malloc(sizeof(short) * indexCnt);
	
	// fill in vertices
	for(i = 0; i <= lineCount; i++) {
		vCopy(&lineStrip[i], &m->vertices[i*2]);
		vAdd(&lineStrip[i], v, &m->vertices[(i*2)+1]);
	}
	
	// fill in indices
	for(i = 0; i <= lineCount; i++) {
		m->indices[(i*6)+0] = i;
		m->indices[(i*6)+1] = i + 1
	}
	
	
	
	
	
}







