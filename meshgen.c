
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "c3dlas.h"
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
	m->indexCnt = 3 * m->vertexCnt;
	
	m->vertices = malloc(sizeof(Vector) * m->vertexCnt);
	m->indices = malloc(sizeof(short) * m->indexCnt);
	
	// fill in vertices
	for(i = 0; i <= lineCount; i++) {
		vCopy(&lineStrip[i], &m->vertices[i*2]);
		vAdd(&lineStrip[i], v, &m->vertices[(i*2)+1]);
	}
	
	// fill in indices
	for(i = 0; i <= lineCount; i++) {
		m->indices[(i*6)+0] = i;
		m->indices[(i*6)+1] = i + 1;
	}
	
	
	
	
	
}


Mesh* allocMesh(int triangles) {
	Mesh* m;
	
	m = calloc(sizeof(Mesh), 1);
	
	if(triangles <= 0) return m;
	
	m->vertices = malloc(sizeof(MeshVertex) * triangles);
	m->indices = malloc(sizeof(unsigned short) * 3 * triangles);
	
	m->szVertices = 3 * triangles;
	m->szIndices = 3 * triangles;
	
	m->vertexCnt = 0;
	m->indexCnt = 0;
	
	return m;
}






// shitty version
Mesh* makeCube(Matrix* mat, int flat) {
	Mesh* m;
	int i;
	
	static Vector vertices[] = {
		// front face
		{-.5, -.5,  .5},
		{-.5,  .5,  .5},
		{ .5,  .5,  .5},
		{ .5, -.5,  .5},
		// back face
		{-.5, -.5, -.5},
		{-.5,  .5, -.5},
		{ .5,  .5, -.5},
		{ .5, -.5, -.5}
	};
	
	static unsigned short indices[] = {
		0,1,2, 2,3,0, 
		4,5,6, 6,7,4,
		0,1,4, 1,4,5,
		3,4,7, 3,0,4,
		2,5,6, 0,2,5,
		2,3,6, 3,6,7,
	};
	
	m = allocMesh(6 * 2);
	
	// transform the vertices
	for(i = 0; i < 8; i++)
		vMatrixMul(&vertices[i], mat, &m->vertices[i].v);
	
	// fill the index buffer
	memcpy(m->indices, indices, sizeof(indices));
	
	
	return m;
}


Mesh* makeCuboid(Vector* p1, Vector* p2) {
	
	Mesh* m;
	Vector min, max;
	int i, n;
	
	vMin(p1, p2, &min);
	vMax(p1, p2, &max);
	
	m = allocMesh(6 * 2);
	
	i = 0;
	n = 0;
	// x+ face
	vSet(max.x, min.y, min.z, &m->vertices[i++].v);
	vSet(max.x, min.y, max.z, &m->vertices[i++].v);
	vSet(max.x, max.y, min.z, &m->vertices[i++].v);
	vSet(max.x, max.y, max.z, &m->vertices[i++].v);

	// x- face
	vSet(min.x, min.y, min.z, &m->vertices[i++].v);
	vSet(min.x, min.y, max.z, &m->vertices[i++].v);
	vSet(min.x, max.y, min.z, &m->vertices[i++].v);
	vSet(min.x, max.y, max.z, &m->vertices[i++].v);

	// y+ face
	vSet(min.x, max.y, min.z, &m->vertices[i++].v);
	vSet(min.x, max.y, max.z, &m->vertices[i++].v);
	vSet(max.x, max.y, min.z, &m->vertices[i++].v);
	vSet(max.x, max.y, max.z, &m->vertices[i++].v);

	// y- face
	vSet(min.x, min.y, min.z, &m->vertices[i++].v);
	vSet(min.x, min.y, max.z, &m->vertices[i++].v);
	vSet(max.x, min.y, min.z, &m->vertices[i++].v);
	vSet(max.x, min.y, max.z, &m->vertices[i++].v);

	// z+ face
	vSet(min.x, min.y, max.z, &m->vertices[i++].v);
	vSet(max.x, min.y, max.z, &m->vertices[i++].v);
	vSet(min.x, max.y, max.z, &m->vertices[i++].v);
	vSet(max.x, max.y, max.z, &m->vertices[i++].v);

	// z- face
	vSet(min.x, min.y, min.z, &m->vertices[i++]);
	vSet(max.x, min.y, min.z, &m->vertices[i++]);
	vSet(min.x, max.y, min.z, &m->vertices[i++]);
	vSet(max.x, max.y, min.z, &m->vertices[i++]);
	
	m->vertexCnt = i;
}






// assumes vertices are not shared 
void calcFlatNormals(Mesh*) {
	int j, i1, i2, i3;
	Vector n;
	
	
	for(j = 0; j < m->indexCnt; j += 3) {
		
		i1 = m->indices[i];
		i2 = m->indices[i+1];
		i3 = m->indices[i+2];
		
		vTriFaceNormal(&m->vertices[i1].v, &m->vertices[i2].v, &m->vertices[i3].v, &n);
		vCopy(&n, &m->vertices[i1].n);
		vCopy(&n, &m->vertices[i2].n);
		vCopy(&n, &m->vertices[i3].n);
	}
}



// rewrites a mesh with no vertex reuse. required for flat shading.
// super mega naive version; git-r-done. probably more cache-friendly than fancy ones anyway.
//    a few bits of data can probably be collected in previous passes, like newVertexCnt.
void duplicateVertices(Mesh* m) {
	
	char* vusage;
	int i;
	int good_already = 1; // this variable looks silly in camelCase
	int newVertexCnt = 0;
	int newVIndex;
	
	
	vusage = calloc(1, sizeof(char) * m->vertexCnt);
	
	// see how many vertices are reused;
	for(i = 0; i < m->indexCnt; i++) {
		vusage[m->indices[i]]++;
		if(vusage[m->indices[i]] > 1) good_already = 0;
	}
	
	// return early if no vertices are used twice
	if(good_already) {
		free(vusage);
		return;
	}
	
	// calculate new vertex buffer size and realloc if necessary
	for(i = 0; i < m->vertexCnt; i++) {
		newVertexCnt += vusage[i];
	}
	
	if(newVertexCnt > szVertices) {
		realloc(m->vertices, sizeof(MeshVertex) * newVertexCnt);
		m->szVertices = newVertexCnt;
	}
	
	// "new" vertices are appended after the last "old" vertex
	newVindex = m->vertexCnt;
	for(i = 0; i < m->indexCnt; i++) {
		int vi;
		
		vi = m->indices[i]; 
		if(vusage[vi] > 1) {
			// copy vertex
			memcpy(&m->vertices[vi], &m->vertices[newVindex], sizeof(MeshVertex));
			
			// update index buffer with new vertex
			m->indices[i] = newVIndex;
			
			// mark off one of the "uses"
			vusage[vi]--;
			
			newVIndex++;
		}
	}
	
	
	m->vertexCnt = newVertexCnt;
}


