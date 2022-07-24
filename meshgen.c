
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "c3dlas.h"
#include "meshgen.h"



#if 0


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
	
	out->vertices = (Vector3*)malloc(cnt * sizeof(Vector3));
	out->normals = (Vector3*)malloc(cnt * sizeof(Vector3));
	out->indices = (short*)malloc(cnt * sizeof(short) * 6); // too big but we'll fix that later
	
	// fill in the vertices and normals
	for(y = 0; y < height + 1; y++) {
		for(x = 0; x < width + 1; x++) {
			Vector3* v, *n;
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



void genIcosahedronPointsf(float radius, Vector3* out, int* count) {
	int i;
	
	// scale factor for inscribed sphere
	float a = 0.7557613140761707304801337020250013926384447888;
	float r = radius / (a * 2); 
	
	// this one has edge length 2
	Vector3 points[] = {
		{ 1.0f, 0.0f, F_GOLDEN}, // top edge along x axis
		{ -1.0f, 0.0f, F_GOLDEN},
		
		{ 0.0f, F_GOLDEN, 1.0f}, // triangle tips off that top edge
		{ 0.0f, -F_GOLDEN, 1.0f},
		
		{ F_GOLDEN, 1.0f, 0.0f}, // waistline
		{ F_GOLDEN, -1.0f, 0.0f},
		{ -F_GOLDEN, 1.0f, 0.0f},
		{ -F_GOLDEN, -1.0f, 0.0f},
		
		{ 0.0f, F_GOLDEN, -1.0f}, // triangle tips off the bottom edge
		{ 0.0f, -F_GOLDEN, -1.0f},
		
		{ 1.0f, 0.0f, -F_GOLDEN}, // bottom edge along x axis
		{ -1.0f, 0.0f, -F_GOLDEN}
	};
	
	if(count) *count = 12;
	for(i = 0; i < 12; i++) {
		vScale3p(&points[i], r, &out[i]);
	}
}

// top half
void genHalfIcosahedronPointsf(float radius, Vector3* out, int* count) {
	int i;
	
	// scale factor for inscribed sphere
	float a = 0.7557613140761707304801337020250013926384447888;
	float r = radius / (a * 2); 
	
	// this one has edge length 2
	Vector3 points[] = {
		{ 1.0f, 0.0f, F_GOLDEN}, // top edge along x axis
		{ -1.0f, 0.0f, F_GOLDEN},
		
		{ 0.0f, F_GOLDEN, 1.0f}, // triangle tips off that top edge
		{ 0.0f, -F_GOLDEN, 1.0f},
		
		{ F_GOLDEN, 1.0f, 0.0f}, // icosahedral waistline
		{ F_GOLDEN, -1.0f, 0.0f},
		{ -F_GOLDEN, 1.0f, 0.0f},
		{ -F_GOLDEN, -1.0f, 0.0f},
		
		{ 0.0f, F_GOLDEN, 0.0f}, // half way down the vertical edges
		{ 0.0f, -F_GOLDEN, 0.0f}
	};
	
	if(count) *count = 10;
	for(i = 0; i < 10; i++) {
		vScale3p(&points[i], r, &out[i]);
	}
}

TriangleMesh* genIcosahedronMesh(float radius) {
	
	TriangleMesh* m;
	
	m = malloc(sizeof(*m));
	m->vertices = malloc(sizeof(*m->vertices) * 12);
	m->indices = malloc(sizeof(*m->indices) * 20 * 3);
	
	m->szVertices = m->vertexCnt = 12;
	m->szIndices = m->indexCnt = 20 * 3;
	
	genIcosahedronPointsf(radius, (Vector3*)m->vertices, NULL);
	
	short indices[] = {
		// top two triangles
		0, 2, 1, 
		1, 0, 3,
		
		// neighbors of those top two triangles
		// signs are shared in x and y
		2, 0, 4,
		2, 1, 6,
		3, 0, 5,
		3, 1, 7,
		
		// standing triangle end caps
		0, 5, 6,
		1, 7, 8,
		
		// sideways triangles
		// all y values share sign
		3, 9, 5,
		3, 9, 7,
		2, 8, 4,
		2, 8, 6,
		
		// inverted triangle end caps
		10, 5, 6,
		11, 7, 8,
		
		// neighbors of the two bottom triangles
		// signs are shared in x and y
		8, 10, 4,
		8, 11, 6,
		9, 10, 5,
		9, 11, 7,
		
		// bottom two triangles
		10, 8, 11, 
		11, 10, 9
	};
	
	memcpy(m->indices, indices, sizeof(short) * 20 * 3); 
	
	return m;
}


// centered, along z axis
TriangleMesh* genHollowCylinder(float innerRadius, float outerRadius, int segments, float length) {
	TriangleMesh* m;
	int i, j;
	float hlen = length / 2;
	
	m = malloc(sizeof(Mesh));

	// endcaps with two rows
	m->szVertices = m->vertexCnt = segments * 2 * 2;
	
	// two end caps, and two triangles per wall quad, two walls
	m->szIndices = m->indexCnt = 3 * (segments * (2 + 2 + 2 + 2));

	m->vertices = malloc(sizeof(Vector3) * m->szVertices);
	m->indices = malloc(sizeof(short) * m->szIndices);
	
	j = 0;
	// outside
	for(i = 0; i < segments; i++) {
		m->vertices[j].v.x = outerRadius * sin(D_2PI / i);
		m->vertices[j].v.y = outerRadius * cos(D_2PI / i);
		m->vertices[j++].v.z = -hlen;
	}
	for(i = 0; i < segments; i++) {
		m->vertices[j].v.x = outerRadius * sin(D_2PI / i);
		m->vertices[j].v.y = outerRadius * cos(D_2PI / i);
		m->vertices[j++].v.z = hlen;
	}
	// inside
	for(i = 0; i < segments; i++) {
		m->vertices[j].v.x = innerRadius * sin(D_2PI / i);
		m->vertices[j].v.y = innerRadius * cos(D_2PI / i);
		m->vertices[j++].v.z = -hlen;
	}
	for(i = 0; i < segments; i++) {
		m->vertices[j].v.x = innerRadius * sin(D_2PI / i);
		m->vertices[j].v.y = innerRadius * cos(D_2PI / i);
		m->vertices[j++].v.z = hlen;
	}
	
	// indices
	j = 0;
	
	// TODO
	return m;
}


// centered, along z axis
TriangleMesh* genCylinder(float radius, int segments, float length) {
	
	TriangleMesh* m;
	int i, j;
	float hlen = length / 2;
	
	m = malloc(sizeof(Mesh));

	// endcaps plus center points
	m->szVertices = m->vertexCnt = segments * 2 + 2;
	
	// two end caps, and two triangles per wall quad
	m->szIndices = m->indexCnt = 3 * (segments * (2 + 2));

	m->vertices = malloc(sizeof(Vector3) * m->szVertices);
	m->indices = malloc(sizeof(short) * m->szIndices);

	// vertices
	for(i = 0; i < segments; i++) {
		m->vertices[i].v.x = radius * sin(D_2PI / i);
		m->vertices[i].v.y = radius * cos(D_2PI / i);
		m->vertices[i].v.z = -hlen;
	}
	m->vertices[segments].v.x = 0.0f;
	m->vertices[segments].v.y = 0.0f;
	m->vertices[segments].v.z = -hlen;

	for(i = 0; i < segments; i++) {
		m->vertices[segments + 1 + i].v.x = radius * sin(D_2PI / i);
		m->vertices[segments + 1 + i].v.y = radius * cos(D_2PI / i);
		m->vertices[segments + 1 + i].v.z = hlen;
	}
	m->vertices[segments * 2 + 1].v.x = 0.0f;
	m->vertices[segments * 2 + 1].v.y = 0.0f;
	m->vertices[segments * 2 + 1].v.z = hlen;

	// indices
	j = 0;
	// first the starting disk
	for(i = 0; i < segments; i++) {
		m->indices[j++] = i;
		m->indices[j++] = (i + 1) % segments;
		m->indices[j++] = segments;
	}
	// then the second disk
	for(i = 0; i < segments; i++) {
		m->indices[j++] = segments + 1 + i;
		m->indices[j++] = segments + 1 + ((i + 1) % segments);
		m->indices[j++] = segments + 1 + segments;
	}
	//finally the walls
	for(i = 0; i < segments; i++) {
		m->indices[j++] = i;
		m->indices[j++] = (i + 1) % segments;
		m->indices[j++] = segments + 1 + ((i + 1) % segments);
		
		m->indices[j++] = i;
		m->indices[j++] = segments + 1 + ((i + 1) % segments);
		m->indices[j++] = segments + 1 + segments;
	}
	
	return m;
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






TriangleMesh* extrudeAlongVector(Vector3* lineStrip, int lineCount, Vector3* v) {
	
	int i;
	Mesh* m;
	
	m = malloc(sizeof(Mesh));
	
	// step 1: allocate enough memory for the mesh
	m->vertexCnt = (lineCount + 1) * 2;
	m->indexCnt = 3 * m->vertexCnt;
	
	m->vertices = malloc(sizeof(Vector3) * m->vertexCnt);
	m->indices = malloc(sizeof(short) * m->indexCnt);
	
	// fill in vertices
	for(i = 0; i <= lineCount; i++) {
//		vCopy3p(&lineStrip[i], &m->vertices[i*2].v);
		vAdd3p(&lineStrip[i], v, &m->vertices[(i*2)+1].v);
	}
	
	// fill in indices
	for(i = 0; i <= lineCount; i++) {
		m->indices[(i*6)+0] = i;
		m->indices[(i*6)+1] = i + 1;
	}
	
	
	
	
	return m;
}


TriangleMesh* allocMesh(int triangles) {
	TriangleMesh* m;
	
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



void growMeshVertices(TriangleMesh* m, int count) {
	if(count < m->szVertices) return;
	
	m->vertices = realloc(m->vertices, sizeof(MeshVertex) * count);
	m->szVertices = count;
}

void growMeshIndices(TriangleMesh* m, int count) {
	if(count < m->szIndices) return;
	
	m->indices = realloc(m->indices, sizeof(MeshVertex) * count);
	m->szIndices = count;
}

void checkGrowMesh(TriangleMesh* m, int newVertices, int newIndices) {
	if(m->szVertices - m->vertexCnt < newVertices) {
		if(m->szVertices < 16) m->szVertices = 16;
		// grow by the next multiple of szVertices large enough to fit the requested quantity
		growMeshVertices(m, (((m->szVertices + newVertices) / m->szVertices) + 1) * m->szVertices);
	}
	
	if(m->szIndices - m->indexCnt < newIndices) {
		if(m->szIndices < 16) m->szIndices = 16;
		// same as above
		growMeshIndices(m, (((m->szIndices + newIndices) / m->szIndices) + 1) * m->szIndices);
	}
	
}



// doesn't check mesh allocation size. Real Programmers(TM) already know how much memory there is.
void appendTriangleFast(TriangleMesh* m, Vector3* v0, Vector3* v1, Vector3* v2) {
	int i = m->vertexCnt;
	int j = m->indexCnt;
	
	m->indices[j++] = i;
	m->indices[j++] = i+1;
	m->indices[j++] = i+2;
	
//	vCopy3p(v0, &m->vertices[i++].v);
//	vCopy3p(v1, &m->vertices[i++].v);
//	vCopy3p(v2, &m->vertices[i++].v);
	
	m->vertexCnt = i;
	m->indexCnt = j;
}

void appendTriangle(TriangleMesh* m, Vector3* v0, Vector3* v1, Vector3* v2) {
	
	checkGrowMesh(m, 3, 3);
	
	appendTriangleFast(m, v0, v1, v2);
}


// appends a face of vnct vertices to the mesh
// CW winding, convex polygons only. will be tesselated as a fan pivoting around v[0].
void appendFace(TriangleMesh* m, Vector3* v, int vcnt) {
	int i;
	
	checkGrowMesh(m, vcnt, 3 * (vcnt - 2));
	
	for(i = 0; i < vcnt - 2; i++) {
		appendTriangleFast(m, &v[0], &v[i+1], &v[i+2]);
	}
}


// shitty version
TriangleMesh* makeCube(Matrix* mat, int flat) {
	Mesh* m;
	int i;
	
	static Vector3 vertices[] = {
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
		vMatrixMul3p(&vertices[i], mat, &m->vertices[i].v);
	
	// fill the index buffer
	memcpy(m->indices, indices, sizeof(indices));
	
	
	return m;
}


TriangleMesh* makeCuboid(Vector3* p1, Vector3* p2) {
/*	
	Mesh* m;
	Vector3 min, max;
	int i, n;
	
	vMin3p(p1, p2, &min);
	vMax3p(p1, p2, &max);
	
	m = allocMesh(6 * 2);
	
	i = 0;
	n = 0;
	// x+ face
	vSet3p(max.x, min.y, min.z, &m->vertices[i++].v);
	vSet3p(max.x, min.y, max.z, &m->vertices[i++].v);
	vSet3p(max.x, max.y, min.z, &m->vertices[i++].v);
	vSet3p(max.x, max.y, max.z, &m->vertices[i++].v);

	// x- face
	vSet3p(min.x, min.y, min.z, &m->vertices[i++].v);
	vSet3p(min.x, min.y, max.z, &m->vertices[i++].v);
	vSet3p(min.x, max.y, min.z, &m->vertices[i++].v);
	vSet3p(min.x, max.y, max.z, &m->vertices[i++].v);

	// y+ face
	vSet3p(min.x, max.y, min.z, &m->vertices[i++].v);
	vSet3p(min.x, max.y, max.z, &m->vertices[i++].v);
	vSet3p(max.x, max.y, min.z, &m->vertices[i++].v);
	vSet3p(max.x, max.y, max.z, &m->vertices[i++].v);

	// y- face
	vSet3p(min.x, min.y, min.z, &m->vertices[i++].v);
	vSet3p(min.x, min.y, max.z, &m->vertices[i++].v);
	vSet3p(max.x, min.y, min.z, &m->vertices[i++].v);
	vSet3p(max.x, min.y, max.z, &m->vertices[i++].v);

	// z+ face
	vSet3p(min.x, min.y, max.z, &m->vertices[i++].v);
	vSet3p(max.x, min.y, max.z, &m->vertices[i++].v);
	vSet3p(min.x, max.y, max.z, &m->vertices[i++].v);
	vSet3p(max.x, max.y, max.z, &m->vertices[i++].v);

	// z- face
	vSet3p(min.x, min.y, min.z, &m->vertices[i++].v);
	vSet3p(max.x, min.y, min.z, &m->vertices[i++].v);
	vSet3p(min.x, max.y, min.z, &m->vertices[i++].v);
	vSet3p(max.x, max.y, min.z, &m->vertices[i++].v);
	
	m->vertexCnt = i;
	
	return m;
*/	
	return NULL;
}






// assumes vertices are not shared
void calcFlatNormals(TriangleMesh* m) {
	int j, i1, i2, i3;
	Vector3 n;
	
	
	for(j = 0; j < m->indexCnt; j += 3) {
		
		i1 = m->indices[j];
		i2 = m->indices[j+1];
		i3 = m->indices[j+2];
		
		vTriFaceNormal3p(&m->vertices[i1].v, &m->vertices[i2].v, &m->vertices[i3].v, &n);
//		vCopy3p(&n, &m->vertices[i1].n);
//		vCopy3p(&n, &m->vertices[i2].n);
//		vCopy3p(&n, &m->vertices[i3].n);
	}
}


// only works if vertices are welded first.
void calcSmoothNormals(TriangleMesh* m) {
	
	Vector3* sums;
	int i, vi0, vi1, vi2;
	Vector3 n;
	
	sums = calloc(1, m->vertexCnt * sizeof(Vector3));
	
	for(i = 0; i < m->indexCnt; i += 3) {
		vi0 = m->indices[i];
		vi1 = m->indices[i];
		vi2 = m->indices[i];
		
		vTriFaceNormal3p(&m->vertices[vi0].v, &m->vertices[vi1].v, &m->vertices[vi2].v, &n);
		
		vAdd3p(&n, &sums[vi0], &sums[vi0]);
		vAdd3p(&n, &sums[vi1], &sums[vi1]);
		vAdd3p(&n, &sums[vi2], &sums[vi2]);
	};
	
	for(i = 0; i < m->vertexCnt; i++) {
		vNorm3p(&sums[i], &m->vertices[i].n);
	}
	
	free(sums);
}



void weldVertices(TriangleMesh* m, float epsilon) {
	
	
	
}


// rewrites a mesh with no vertex reuse. required for flat shading.
// super mega naive version; git-r-done. probably more cache-friendly than fancy ones anyway.
//    a few bits of data can probably be collected in previous passes, like newVertexCnt.
void unweldVertices(TriangleMesh* m) {
	
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
	
	if(newVertexCnt > m->szVertices) {
		realloc(m->vertices, sizeof(MeshVertex) * newVertexCnt);
		m->szVertices = newVertexCnt;
	}
	
	// "new" vertices are appended after the last "old" vertex
	newVIndex = m->vertexCnt;
	for(i = 0; i < m->indexCnt; i++) {
		int vi;
		
		vi = m->indices[i];
		if(vusage[vi] > 1) {
			// copy vertex
			memcpy(&m->vertices[vi], &m->vertices[newVIndex], sizeof(MeshVertex));
			
			// update index buffer with new vertex
			m->indices[i] = newVIndex;
			
			// mark off one of the "uses"
			vusage[vi]--;
			
			newVIndex++;
		}
	}
	
	
	m->vertexCnt = newVertexCnt;

	free(vusage);
}



MeshSlice* allocMeshSlice(int vcnt, int icnt) {
	
	MeshSlice* ms;
	
	ms = calloc(1, sizeof(MeshSlice));
	
	ms->vertices = calloc(1, vcnt * sizeof(MeshVertex));
	ms->indices = calloc(1, icnt * sizeof(short));
	
	ms->szVertices = vcnt;
	ms->szIndices = icnt;
	
	return ms;
}



MeshSlice* makeCircle(float radius, int divisions) {
	
	MeshSlice* ms;
	int i;
	float divf = (float)divisions;
	
	ms = allocMeshSlice(divisions, divisions + 1);
	
	for(i = 0; i < divisions; i++) {
		ms->vertices[i].v.x = radius * sin((i / divf) * F_2PI);
		ms->vertices[i].v.y = radius * cos((i / divf) * F_2PI);
		ms->vertices[i].v.z = 0;
//		vCopy3p(&ms->vertices[i].v, &ms->vertices[i].n);
		ms->vertices[i].t.u = i / divf;
		ms->vertices[i].t.v = 0;
		
		ms->indices[i] = i;
	}
	
	ms->indices[divisions] = 0;
	
	ms->vertexCnt = divisions;
	ms->indexCnt = divisions + 1;
	
	ms->texDxDy.x = 0;
	ms->texDxDy.y = 1;
	
	return ms;
}


void appendBezierSpline(MeshSlice* ms, BezierSpline3* bs, float tolerance, int connectToEnd) {
	
	
	
}

#endif
