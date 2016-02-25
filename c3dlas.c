
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "c3dlas.h"




// vector operations

int vEq(Vector* a, Vector* b) {
	return vEqEp(a, b, FLT_EPSILON);
}

int vEqEp(Vector* a, Vector* b, float epsilon) {
	float x, y, z, n;
	
	x = a->x - b->x;
	x = a->y - b->y;
	x = a->z - b->z;
	
	n = fabs(x * x + y * y + z * z);
	
	return n <= epsilon * epsilon; 
}


void vCopy(const Vector* src, Vector* dst) {
	dst->x = src->x;
	dst->y = src->y;
	dst->z = src->z;
}

void vSwap(Vector* a, Vector* b) { // swap two vectors
	float x, y, z;
	x = a->x;
	y = a->y;
	z = a->z;
	a->x = b->x;
	a->y = b->y;
	a->z = b->z;
	b->x = x;
	b->y = y;
	b->z = z;
}

void vAdd(Vector* a, Vector* b, Vector* out) {
	out->x = a->x + b->x;
	out->y = a->y + b->y;
	out->z = a->z + b->z;
}

void vSub(Vector* from, Vector* what, Vector* diff) { // diff = from - what
	diff->x = from->x - what->x;
	diff->y = from->y - what->y;
	diff->z = from->z - what->z;
}

void vScale(Vector* v, float scalar, Vector* out) {
	out->x = v->x * scalar;
	out->y = v->y * scalar;
	out->z = v->z * scalar;
}

void vInverse(Vector* v, Vector* out) {
	// yeah yeah yeah shut up. in games, pesky details are just annoying. this function does what you mean rather than sucking Gauss and Euclid.
	out->x = v->x == 0.0f ? FLT_MAX : 1.0f / v->x;
	out->y = v->y == 0.0f ? FLT_MAX : 1.0f / v->y;
	out->z = v->z == 0.0f ? FLT_MAX : 1.0f / v->z;
}

float vMag(Vector* v) {
	return sqrt((float)((v->x * v->x) + (v->y * v->y) + (v->z * v->z)));
}

float vDot(Vector* a, Vector* b) {
	return (float)((a->x * b->x) + (a->y * b->y) + (a->z * b->z));
}


void vNorm(Vector* v, Vector* out) {
	vUnit(v, out);
}

void vUnit(Vector* v, Vector* out) {
	float n;
	n = (v->x * v->x) + (v->y * v->y) + (v->z * v->z);
	
	if(n >= 1.0f - FLT_EPSILON || n >= 1.0f + FLT_EPSILON) return; // very exact here
	
	n = 1.0f / sqrtf(n);
	out->x = v->x * n;
	out->y = v->y * n;
	out->z = v->z * n;
}

void vCross(Vector* a, Vector* b, Vector* out) { // c = a x b
	out->x = (a->y * b->z) - (a->z * b->y);
	out->y = (a->z * b->x) - (a->x * b->z);
	out->z = (a->x * b->y) - (a->y * b->x);
}

float vScalarTriple(Vector* a, Vector* b, Vector* c) { // a . (b x c)
	return (float)((a->x * b->y * c->z) - (a->x * b->z * c->y) - (a->y * b->x * c->z)
				 + (a->z * b->x * c->y) + (a->y * b->z * c->x) - (a->z * b->y * c->x));
}


// feeding a zero vector into this will cause div/0 and you will deserve it
void  vProject(Vector* what, Vector* onto, Vector* out) { // slower; onto may not be normalized
	float wdo = vDot(what, onto);
	float odo = vDot(onto, onto);
	vScale(onto, wdo / odo, out);
}

void  vProjectNorm(Vector* what, Vector* onto, Vector* out) { // faster; onto must be normalized
	float wdo = vDot(what, onto);
	vScale(onto, wdo, out);
}


// returns the minimum values of each component
void  vMin(Vector* a, Vector* b, Vector* out) {
	out->x = fmin(a->x, b->x);
	out->y = fmin(a->y, b->y);
	out->z = fmin(a->z, b->z);
}

// returns the maximum values of each component
void  vMax(Vector* a, Vector* b, Vector* out) {
	out->x = fmax(a->x, b->x);
	out->y = fmax(a->y, b->y);
	out->z = fmax(a->z, b->z);
}

void inline vSet(float x, float y, float z, Vector* out) {
	out->x = x;
	out->y = y;
	out->z = x;
}


// reflects the distance from v to pivot across pivot. 
// out, pivot, and v will form a straight line with pivot exactly in the middle.
void vReflectAcross(Vector* v, Vector* pivot, Vector* out) {
	Vector diff;
	
	vSub(pivot, v, &diff);
	vAdd(pivot, &diff, out);
}

// calculate a unit vector normal to a triangle's face.
void  vTriFaceNormal(Vector* a, Vector* b, Vector* c, Vector* out) {
	Vector a_b, a_c;
	
	vSub(a, b, &a_b);
	vSub(a, c, &a_c);
	vCross(&a_b, &a_b, out);
	vNorm(out, out);
}




// 2d vector stuff

void vSwap2(Vector2* a, Vector2* b) { // swap two vectors
	float x, y;
	x = a->x;
	y = a->y;
	a->x = b->x;
	a->y = b->y;
	b->x = x;
	b->y = y;
}

void vAdd2(Vector2* a, Vector2* b, Vector2* out) {
	out->x = a->x + b->x;
	out->y = a->y + b->y;
}

void vSub2(Vector2* from, Vector2* what, Vector2* diff) { // diff = from - what
	diff->x = from->x - what->x;
	diff->y = from->y - what->y;
}

void vScale2(Vector2* v, float scalar, Vector2* out) {
	out->x = v->x * scalar;
	out->y = v->y * scalar;
}

void vNorm2(Vector2* v, Vector2* out) {
	vUnit2(v, out);
}

void vUnit2(Vector2* v, Vector2* out) {
	float n;
	n = (v->x * v->x) + (v->y * v->y);
	
	if(n >= 1.0f - FLT_EPSILON || n >= 1.0f + FLT_EPSILON) return; // very exact here
	
	n = 1.0f / sqrtf(n);
	out->x = v->x * n;
	out->y = v->y * n;
}


// reflects the distance from v to pivot across pivot. 
// out, pivot, and v will form a straight line with pivot exactly in the middle.
void vReflectAcross2(Vector2* v, Vector2* pivot, Vector2* out) {
	Vector2 diff;
	
	vSub2(pivot, v, &diff);
	vAdd2(pivot, &diff, out);
}





void vSwap2i(Vector2i* a, Vector2i* b) { // swap two vectors
	int x, y;
	x = a->x;
	y = a->y;
	a->x = b->x;
	a->y = b->y;
	b->x = x;
	b->y = y;
}


// degenerate cases may not give desired results. GIGO.
void vRoundAway2(const Vector2* in, const Vector2* center, Vector2i* out) {
	
	if(in->x > center->x) out->x = ceilf(in->x);
	else out->x = floorf(in->x);
	
	if(in->y > center->y) out->y = ceilf(in->y);
	else out->y = floorf(in->y);
}

// degenerate cases may not give desired results. GIGO.
void vRoundToward2(const Vector2* in, const Vector2* center, Vector2i* out) {
	
	if(in->x > center->x) out->x = floorf(in->x);
	else out->x = ceilf(in->x);
	
	if(in->y > center->y) out->y = floorf(in->y);
	else out->y = ceilf(in->y);
}

// plane-vector operations

// distance from point to plane
float pvDist(Plane* p, Vector* v) {
	return vDot(v, &p->n) + p->d;
}



// matrix-vector operations


// multiply a vector by a matrix
void vMatrixMul(Vector* in, Matrix* m, Vector* out) {
	vMatrixMulf(in->x, in->y, in->z, m, out);
}

void vMatrixMulf(float x, float y, float z, Matrix* m, Vector* out) { 
	Vector4 v;

	v.x = x * m->m[0+0] + y * m->m[0+1] + z * m->m[0+2] + 1 * m->m[0+3];
	v.y = x * m->m[4+0] + y * m->m[4+1] + z * m->m[4+2] + 1 * m->m[4+3];
	v.z = x * m->m[8+0] + y * m->m[8+1] + z * m->m[8+2] + 1 * m->m[8+3];
	v.w = x * m->m[12+0] + y * m->m[12+1] + z * m->m[12+2] + 1 * m->m[12+3];
	
	out->x = v.x / v.w;
	out->y = v.y / v.w;
	out->z = v.z / v.w;
}



// matrix operations






const Matrix IDENT_MATRIX = { { 1, 0, 0, 0, 
                                0, 1, 0, 0, 
                                0, 0, 1, 0, 
                                0, 0, 0, 1 } };


void mIdent(Matrix* m) {
	*m = IDENT_MATRIX;
}

void mCopy(Matrix* in, Matrix* out) { 
	memcpy(in->m, out->m, sizeof(out->m));
}


void mFastMul(Matrix* a, Matrix* b, Matrix* out) {
	int r, c;
	
	for(r = 0; r < 4; r++) {
		for(c = 0; c < 4; c++) {
			out->m[c + r * 4] = 
				(a->m[r * 4 + 0] * b->m[c + 0]) + 
				(a->m[r * 4 + 1] * b->m[c + 4]) + 
				(a->m[r * 4 + 2] * b->m[c + 8]) + 
				(a->m[r * 4 + 3] * b->m[c + 12]);
		}
	}
}

void mMul(Matrix* a, Matrix* out) {
	Matrix b;
	
	mCopy(&b, out);
	
	mFastMul(a, &b, out);
}



void mTransv(Vector* v, Matrix* out) {
	mTrans3f(v->x, v->y, v->z, out);
}

void mTrans3f(float x, float y, float z, Matrix* out) {
	Matrix t;
	
	t = IDENT_MATRIX;
	t.m[12] = x;
	t.m[13] = y;
	t.m[14] = z;
	
	mMul(&t, out);
}



void mScalev(Vector* v, Matrix* out) {
	mTrans3f(v->x, v->y, v->z, out);
}

void mScale3f(float x, float y, float z, Matrix* out) {
	Matrix t;
	
	t = IDENT_MATRIX;
	t.m[0] = x;
	t.m[5] = y;
	t.m[10] = z;
	
	mMul(&t, out);
}



void mRotv(Vector* v, float theta, Matrix* out) {
	mRot3f(v->x, v->y, v->z, theta, out);
}

void mRot3f(float x, float y, float z, float theta, Matrix* out) {
	
	float costh, omcosth;
	float sinth;
	Matrix r;
	
	costh = cos(theta);
	omcosth = 1 - costh;
	sinth = sin(theta);
	
	r.m[0] = costh + (x * x * omcosth);
	r.m[1] = (z * sinth) + (x * y * omcosth);
	r.m[2] = (-y * sinth) + (x * z * omcosth);
	r.m[3] = 0;
	
	r.m[4] = (x * y * omcosth) - (z * sinth);
	r.m[5] = costh + (y * y * omcosth);
	r.m[6] = (x * sinth) + (y * z * omcosth);
	r.m[7] = 0;
	
	r.m[8]  = (y * sinth) + (x * z * omcosth);
	r.m[9]  = (-x * sinth) + (y * z * omcosth);
	r.m[10] = costh + (z * z * omcosth);
	r.m[11] = 0;
	
	r.m[12] = 0;
	r.m[13] = 0;
	r.m[14] = 0;
	r.m[15] = 1;
	
	mMul(&r, out);
}


void mRotX(float theta, Matrix* out) {
	mRot3f(1,0,0, theta, out);
}

void mRotY(float theta, Matrix* out) {
	mRot3f(0,1,0, theta, out);
}

void mRotZ(float theta, Matrix* out) {
	mRot3f(0,0,1, theta, out);
}


void mTranspose(Matrix* in, Matrix* out) {
	Matrix t;
	int i;

	for(i = 0; i < 4; i++) {
		t.m[i]      = in->m[i * 4];
		t.m[i + 4]  = in->m[(i * 4) + 1];
		t.m[i + 8]  = in->m[(i * 4) + 2];
		t.m[i + 12] = in->m[(i * 4) + 3];
	}
	
	mCopy(&t, out);
}

void mTransposeFast(Matrix* in, Matrix* out) {
	int i;
	for(i = 0; i < 4; i++) {
		out->m[i]      = in->m[i * 4];
		out->m[i + 4]  = in->m[(i * 4) + 1];
		out->m[i + 8]  = in->m[(i * 4) + 2];
		out->m[i + 12] = in->m[(i * 4) + 3];
	}
}


float mDeterminate(Matrix* m) {
	return
		m->m[3] * m->m[6] * m->m[9]  * m->m[12] - m->m[2] * m->m[7] * m->m[9]  * m->m[12] -
		m->m[3] * m->m[5] * m->m[10] * m->m[12] + m->m[1] * m->m[7] * m->m[10] * m->m[12] +
		m->m[2] * m->m[5] * m->m[11] * m->m[12] - m->m[1] * m->m[6] * m->m[11] * m->m[12] -
		m->m[3] * m->m[6] * m->m[8]  * m->m[13] + m->m[2] * m->m[7] * m->m[8]  * m->m[13] +
		m->m[3] * m->m[4] * m->m[10] * m->m[13] - m->m[0] * m->m[7] * m->m[10] * m->m[13] -
		m->m[2] * m->m[4] * m->m[11] * m->m[13] + m->m[0] * m->m[6] * m->m[11] * m->m[13] +
		m->m[3] * m->m[5] * m->m[8]  * m->m[14] - m->m[1] * m->m[7] * m->m[8]  * m->m[14] -
		m->m[3] * m->m[4] * m->m[9]  * m->m[14] + m->m[0] * m->m[7] * m->m[9]  * m->m[14] +
		m->m[1] * m->m[4] * m->m[11] * m->m[14] - m->m[0] * m->m[5] * m->m[11] * m->m[14] -
		m->m[2] * m->m[5] * m->m[8]  * m->m[15] + m->m[1] * m->m[6] * m->m[8]  * m->m[15] +
		m->m[2] * m->m[4] * m->m[9]  * m->m[15] - m->m[0] * m->m[6] * m->m[9]  * m->m[15] -
		m->m[1] * m->m[4] * m->m[10] * m->m[15] + m->m[0] * m->m[5] * m->m[10] * m->m[15];
}


// shamelessly lifted from this SO post and modified into C, added error checking, and transposed for column major ordering:
// http://stackoverflow.com/a/7596981
// matrix inversions suck. maybe one day i'll lift intel's super fast SSE one instead. 
// functions returns 0 if sucessful, 1 if there is no inverse
int mInverse(Matrix* in, Matrix* out) {
	
	float s0, s1, s2, s3, s4, s5;
	float c0, c1, c2, c3, c4, c5;
	float invdet;

	s0 = in->m[0] * in->m[5]  - in->m[1] * in->m[4];
	s1 = in->m[0] * in->m[9]  - in->m[1] * in->m[8];
	s2 = in->m[0] * in->m[13] - in->m[1] * in->m[12];
	s3 = in->m[4] * in->m[9]  - in->m[5] * in->m[8];
	s4 = in->m[4] * in->m[13] - in->m[5] * in->m[12];
	s5 = in->m[8] * in->m[13] - in->m[9] * in->m[12];

	c5 = in->m[10] * in->m[15] - in->m[11] * in->m[14];
	c4 = in->m[6]  * in->m[15]  - in->m[7]  * in->m[14];
	c3 = in->m[6]  * in->m[11]  - in->m[7]  * in->m[10];
	c2 = in->m[2]  * in->m[15]  - in->m[3]  * in->m[14];
	c1 = in->m[2]  * in->m[11]  - in->m[3]  * in->m[10];
	c0 = in->m[2]  * in->m[7]   - in->m[3]  * in->m[6];

	// Should check for 0 determinant

	invdet = (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0);
	if(invdet == 0.0) {
		fprintf(stderr, "ERROR: Matrix has no inverse!!!\n");
		return 1;
	}
	invdet = 1.0 / invdet;
	
	out->m[0]  = ( in->m[5] * c5 - in->m[9]  * c4 + in->m[13] * c3) * invdet;
	out->m[4]  = (-in->m[4] * c5 + in->m[8]  * c4 - in->m[12] * c3) * invdet;
	out->m[8]  = ( in->m[7] * s5 - in->m[11] * s4 + in->m[15] * s3) * invdet;
	out->m[12] = (-in->m[6] * s5 + in->m[10] * s4 - in->m[14] * s3) * invdet;

	out->m[1]  = (-in->m[1] * c5 + in->m[9]  * c2 - in->m[13] * c1) * invdet;
	out->m[5]  = ( in->m[0] * c5 - in->m[8]  * c2 + in->m[12] * c1) * invdet;
	out->m[9]  = (-in->m[3] * s5 + in->m[11] * s2 - in->m[15] * s1) * invdet;
	out->m[13] = ( in->m[2] * s5 - in->m[10] * s2 + in->m[14] * s1) * invdet;

	out->m[2]  = ( in->m[1] * c4 - in->m[5] * c2 + in->m[13] * c0) * invdet;
	out->m[6]  = (-in->m[0] * c4 + in->m[4] * c2 - in->m[12] * c0) * invdet;
	out->m[10] = ( in->m[3] * s4 - in->m[7] * s2 + in->m[15] * s0) * invdet;
	out->m[14] = (-in->m[2] * s4 + in->m[6] * s2 - in->m[14] * s0) * invdet;

	out->m[3]  = (-in->m[1] * c3 + in->m[5] * c1 - in->m[9]  * c0) * invdet;
	out->m[7]  = ( in->m[0] * c3 - in->m[4] * c1 + in->m[8]  * c0) * invdet;
	out->m[11] = (-in->m[3] * s3 + in->m[7] * s1 - in->m[11] * s0) * invdet;
	out->m[15] = ( in->m[2] * s3 - in->m[6] * s1 + in->m[10] * s0) * invdet;

	return 0;
}



// analogous to glFrustum
// no div/0 checking here for right == left etc. just don't be an idiot.
void mFrustum(float left, float right, float top, float bottom, float near, float far, Matrix* out) {
	
	Matrix m;
	
	m = IDENT_MATRIX;
	
	m.m[0] = (2 * near) / (right - left);
	m.m[5] = (2 * near) / (top - bottom);
	m.m[8] = (right + left) / (right - left);
	m.m[9] = (top + bottom) / (top - bottom);
	m.m[10] = -(far + near) / (far - near);
	m.m[11] = -1;
	m.m[14] = (-2 * far * near) / (far - near);
	m.m[15] = 0;
	
	mMul(&m, out);
}


// analogous to gluPerspective
// same div/0 warnings apply. if you get an FP exception you deserve it. 
// use a double for fov; the precision matters often.
// https://www.opengl.org/archives/resources/faq/technical/transformations.htm
// https://www.opengl.org/sdk/docs/man2/xhtml/gluPerspective.xml
void mPerspective(double fov, float aspect, float near, float far, Matrix* out) {
	
	Matrix m;
	double f;
	
	m = IDENT_MATRIX;
	f = 1.0 / tan(fov * DEG2RAD / 2);
	
	m.m[0] = f / aspect;
	m.m[5] = f;
	m.m[10] = (far + near) / (near - far);
	m.m[11] = -1;
	m.m[14] = (2 * far * near) / (near - far);
	m.m[15] = 1;
	
	mMul(&m, out);
}




// orthographic projection. use this for a "2D" look.
void mOrtho(float left, float right, float top, float bottom, float near, float far, Matrix* out) {
	
	Matrix m;
	
	m = IDENT_MATRIX;
	
	m.m[0] = 2 / (right - left);
	m.m[5] = 2 / (top - bottom);
	m.m[10] = -2 / (far - near);
	m.m[12] = -(right + left) / (right - left);
	m.m[13] = -(top + bottom) / (top - bottom);		
	m.m[14] = -(far + near) / (far - near);
	m.m[15] = 1;
	
	mMul(&m, out);
}

// analgous to gluLookAt
// https://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
void mLookAt(Vector* eye, Vector* center, Vector* up, Matrix* out) {
	
	Vector f, upn, s, u, sn;
	Matrix m, m2;
	
	vSub(center, eye, &f);
	vNorm(&f, &f);
	
	vNorm(up, &upn);
	
	vCross(&f, &upn, &s);
	vNorm(&s, &sn);
	
	vCross(&sn, &f, &u);
	
	m.m[0] = s.x;
	m.m[1] = u.x;
	m.m[2] = -f.x;
	m.m[3] = 0;

	m.m[4] = s.y;
	m.m[5] = u.y;
	m.m[6] = -f.y;
	m.m[7] = 0;
	
	m.m[8] = s.z;
	m.m[9] = u.z;
	m.m[10] = -f.z;
	m.m[11] = 0;
	
	m.m[12] = 0;
	m.m[13] = 0;
	m.m[14] = 0;
	m.m[15] = 1;

	mTrans3f(-eye->x, -eye->y, -eye->z, &m2);
	mFastMul(&m, &m2, out);
}


void mPrint(Matrix* m, FILE* f) {
	int r, c;
	
	if(!f) f = stdout;
	
	for(r = 0; r < 4; r++) {
		fprintf(f, "% .3e % .3e % .3e % .3e\n", m->m[r*4], m->m[r*4+1], m->m[r*4+2], m->m[r*4+3]);
	}
	
	fprintf(f, "\n");
}




// make sure you allocate enough. when it's out, it's out. no surprise mallocs later on. (yet)
void msAlloc(int size, MatrixStack* ms) {
	
	ms->stack = (Matrix*)malloc(size * sizeof(Matrix));
	
	ms->size = size;
	ms->top = 0;
};
	
void msFree(MatrixStack* ms) {
	free(ms->stack);
	ms->stack = NULL;
}

// push a new matrix on the stack. if m is null, push an identity matrix
int msPush(MatrixStack* ms) {
	if(ms->top == ms->size - 1) {
		fprintf(stderr, "Matrix Stack overflowed.\n");
		return 1;
	}
	
	ms->top++;
	
	mCopy(&ms->stack[ms->top], &ms->stack[ms->top - 1]);
	
	return 0;
}

void msPop(MatrixStack* ms) {
	if(ms->top == 0) return;
	ms->top--;
}


Matrix* msGetTop(MatrixStack* ms) {
	return &ms->stack[ms->top];
}

void msPrintAll(MatrixStack* ms, FILE* f) {
	int i;
	
	if(!f) f = stdout;
	
	for(i = 0; i <= ms->top; i++) {
		fprintf(f, "--%3d--------------------------\n", i);
		mPrint(&ms->stack[i], f);
	}
	
	fprintf(f, "-------------------------------\n");
}


void msIdent(MatrixStack* ms) {
	mIdent(msGetTop(ms));
}

void msCopy(Matrix* in, MatrixStack* ms) {
	mCopy(in, msGetTop(ms));
}

void msMul(Matrix* a, MatrixStack* ms) { // makes a copy of out before multiplying over it
	mMul(a, msGetTop(ms));
}

void msTransv(Vector* v, MatrixStack* ms) { // translation
	mTransv(v, msGetTop(ms));
}

void msTrans3f(float x, float y, float z, MatrixStack* ms) { // translation
	mTrans3f(x, y, z, msGetTop(ms));
}

void msScalev(Vector* v, MatrixStack* ms) {
	mScalev(v, msGetTop(ms));
}

void msScale3f(float x, float y, float z, MatrixStack* ms) {
	mScale3f(x, y, z, msGetTop(ms));
}

void msRotv(Vector* v, float theta, MatrixStack* ms) { // rotate about a vector
	mRotv(v, theta, msGetTop(ms));
}

void msRot3f(float x, float y, float z, float theta, MatrixStack* ms) { // rotate about a vector
	mRot3f(x, y, z, theta, msGetTop(ms));
}

void msFrustum(float left, float right, float top, float bottom, float near, float far, MatrixStack* ms) { 
	mFrustum(left, right, top, bottom, near, far, msGetTop(ms));
}

void msPerspective(double fov, float aspect, float near, float far, MatrixStack* ms) {
	mPerspective(fov, aspect, near, far, msGetTop(ms));
}

void msOrtho(float left, float right, float top, float bottom, float near, float far, MatrixStack* ms) {
	mOrtho(left, right, top, bottom, near, far, msGetTop(ms));
}

void msLookAt(Vector* eye, Vector* center, Vector* up, MatrixStack* ms) {
	mLookAt(eye, center, up, msGetTop(ms));
}



void evalBezier(Vector* e1, Vector* e2, Vector* c1, Vector* c2, float t, Vector* out) {
	out->x = evalBezier1D(e1->x, e2->x, c1->x, c2->x, t);
	out->y = evalBezier1D(e1->y, e2->y, c1->y, c2->y, t);
	out->z = evalBezier1D(e1->z, e2->z, c1->z, c2->z, t);
}

void evalBezier2D(Vector2* e1, Vector2* e2, Vector2* c1, Vector2* c2, float t, Vector2* out) {
	out->x = evalBezier1D(e1->x, e2->x, c1->x, c2->x, t);
	out->y = evalBezier1D(e1->y, e2->y, c1->y, c2->y, t);
}

float evalBezier1D(float e1, float e2, float c1, float c2, float t) {
	float mt, mt2, t2;
	mt = 1 - t;
	mt2 = mt * mt;
	t2 = t * t;
	return (mt2 * mt * e1) + (3 * mt2 * t * c1) + (3 * t2 * mt * c2) + (t2 * t * e2);
}



void evalBezierTangent(Vector* e1, Vector* e2, Vector* c1, Vector* c2, float t, Vector* out) {
	out->x = evalBezier1D_dt(e1->x, e2->x, c1->x, c2->x, t);
	out->y = evalBezier1D_dt(e1->y, e2->y, c1->y, c2->y, t);
	out->z = evalBezier1D_dt(e1->z, e2->z, c1->z, c2->z, t);
}

float evalBezier1D_dt(float e1, float e2, float c1, float c2, float t) {
	float mt, mt2, t2;
	mt = 1 - t;
	mt2 = mt * mt;
	t2 = t * t;
	
	return (3 * mt2 * (c1 - e1)) + (6 * mt * t * (c2 - c1)) + (3 * t2 * (e2 - c2));
}


	
float evalBezier1D_ddt(float e1, float e2, float c1, float c2, float t) {
	return (6 * (1 - t) * (c2 - c1 - c1 + e1)) + (6 * t * (e2 - c2 - c2 - c1));
}

void evalBezierNorm(Vector* e1, Vector* e2, Vector* c1, Vector* c2, float t, Vector* out) {
	out->x = evalBezier1D_ddt(e1->x, e2->x, c1->x, c2->x, t);
	out->y = evalBezier1D_ddt(e1->y, e2->y, c1->y, c2->y, t);
	out->z = evalBezier1D_ddt(e1->z, e2->z, c1->z, c2->z, t);
}




///// bounding box functions


// 3D versions

int boxDisjoint(const AABB* a, const AABB* b) {

	return a->max.x < b->min.x || b->max.x < a->min.x
		|| a->max.y < b->min.y || b->max.y < a->min.y
		|| a->max.z < b->min.z || b->max.z < a->min.z;
}

int boxOverlaps(const AABB* a, const AABB* b) {
	return !boxDisjoint(a, b);
}



int boxContainsPoint(const AABB* b, const Vector* p) {
	return b->min.x <= p->x && b->max.x >= p->x
		&& b->min.y <= p->y && b->max.y >= p->y
		&& b->min.z <= p->z && b->max.z >= p->z;
}


void boxCenter(const AABB* b, Vector* out) {
	out->x = (b->max.x + b->min.x) / 2;
	out->y = (b->max.y + b->min.y) / 2;
	out->z = (b->max.z + b->min.z) / 2;
}

void boxSize(const AABB* b, Vector* out) {
	out->x = b->max.x - b->min.x;
	out->y = b->max.y - b->min.y;
	out->z = b->max.z - b->min.z;
}



// 2D versions

int boxDisjoint2(const AABB2* a, const AABB2* b) {

	return a->max.x < b->min.x || b->max.x < a->min.x
		|| a->max.y < b->min.y || b->max.y < a->min.y;
}

int boxOverlaps2(const AABB2* a, const AABB2* b) {
	return !boxDisjoint2(a, b);
}



int boxContainsPoint2(const AABB2* b, const Vector2* p) {
	return b->min.x <= p->x && b->max.x >= p->x
		&& b->min.y <= p->y && b->max.y >= p->y;
}


void boxCenter2(const AABB2* b, Vector2* out) {
	out->x = (b->max.x + b->min.x) / 2;
	out->y = (b->max.y + b->min.y) / 2;
}

void boxSize2(const AABB2* b, Vector2* out) {
	out->x = b->max.x - b->min.x;
	out->y = b->max.y - b->min.y;
}

void boxQuadrant2(const AABB2* in, char ix, char iy, AABB2* out) {
	Vector2 sz, c;
	
	boxCenter2(in, &c);
	boxSize2(in, &sz);
	sz.x *= .5;
	sz.y *= .5;
	
	out->min.x = c.x - (ix ? 0.0f : sz.x);
	out->min.y = c.y - (iy ? 0.0f : sz.y);
	out->max.x = c.x + (ix ? sz.x : 0.0f);
	out->max.y = c.y + (iy ? sz.y : 0.0f);
}


// 2D integer versions

int boxDisjoint2i(const AABB2i* a, const AABB2i* b) {

	return a->max.x < b->min.x || b->max.x < a->min.x
		|| a->max.y < b->min.y || b->max.y < a->min.y;
}

int boxOverlaps2i(const AABB2i* a, const AABB2i* b) {
	return !boxDisjoint2i(a, b);
}



int boxContainsPoint2i(const AABB2i* b, const Vector2i* p) {
	return b->min.x <= p->x && b->max.x >= p->x
		&& b->min.y <= p->y && b->max.y >= p->y;
}


void boxCenter2i(const AABB2i* b, Vector2* out) {
	out->x = (b->max.x + b->min.x) / 2.0f;
	out->y = (b->max.y + b->min.y) / 2.0f;
}

void boxSize2i(const AABB2i* b, Vector2i* out) {
	out->x = b->max.x - b->min.x;
	out->y = b->max.y - b->min.y;
}

// BUG: needs some fancy math work to keep everything tight. integers don't split nicely
void boxQuadrant2i(const AABB2i* in, char ix, char iy, AABB2i* out) {
	Vector2 sz, c;
	
	printf("fix me: %s:%d", __FILE__, __LINE__);
	exit(666);
	
	boxCenter2(in, &c);
	boxSize2(in, &sz);
	sz.x *= .5;
	sz.y *= .5;
	
	out->min.x = c.x - (ix ? 0.0f : sz.x);
	out->min.y = c.y - (iy ? 0.0f : sz.y);
	out->max.x = c.x + (ix ? sz.x : 0.0f);
	out->max.y = c.y + (iy ? sz.y : 0.0f);
}



// find the center of a quad
void quadCenter2(const Quad2* q, Vector2* out) {
	Vector2 c;
	int i;
	
	for(i = 0; i < 4; i++) {
		c.x += q->v[i].x;
		c.y += q->v[i].y;
	}
	
	out->x = c.x / 4;
	out->y = c.y / 4;
}


void quadRoundOutward2(const Quad2* in, Quad2i* out) {
	Vector2 c;
	int i;
	
	quadCenter2(in, &c);
	
	for(i = 0; i < 4; i++)
		vRoundAway2(&in->v[i], &c, &out->v[i]);
}

void quadRoundInward2(const Quad2* in, Quad2i* out) {
	Vector2 c;
	int i;
	
	quadCenter2(in, &c);
	
	for(i = 0; i < 4; i++)
		vRoundToward2(&in->v[i], &c, &out->v[i]);
}


int quadIsPoint2i(const Quad2i* q) {
	return (
		q->v[0].x == q->v[1].x == q->v[2].x == q->v[3].x && 
		q->v[0].y == q->v[1].y == q->v[2].y == q->v[3].y);
}

int quadIsAARect2i(const Quad2i* q) {
	return (
		q->v[0].x == q->v[3].x && q->v[1].x == q->v[2].x && 
		q->v[0].y == q->v[1].y && q->v[2].y == q->v[3].y);
}


// ray stuff
void makeRay(Vector* origin, Vector* direction, Ray* out) {
	
	out->o.x = origin->x;
	out->o.y = origin->y;
	out->o.z = origin->z;
	
	vNorm(direction, &out->d);
	vInverse(&out->d, &out->id);
}

// this version has no branching, but only answers yes or no. 
// algorithm explanation here. hopefully my extrapolation into 3 dimensions is correct.
// http://tavianator.com/fast-branchless-raybounding-box-intersections/
int boxRayIntersectFast(const AABB* b, const Ray* r) {
	Vector t1, t2;
	float tmin, tmax;
	
	t1.x = (b->min.x - r->o.x) * r->id.x;
	t2.x = (b->max.x - r->o.x) * r->id.x;
	tmin = fmin(t1.x, t2.x);
	tmax = fmax(t1.x, t2.x);
	
	t1.y = (b->min.y - r->o.y) * r->id.y;
	t2.y = (b->max.y - r->o.y) * r->id.y;
	tmin = fmax(tmin, fmin(t1.y, t2.y));
	tmax = fmin(tmax, fmax(t1.y, t2.y));
	
	t1.z = (b->min.z - r->o.z) * r->id.z;
	t2.z = (b->max.z - r->o.z) * r->id.z;
	tmin = fmax(tmin, fmin(t1.z, t2.z));
	tmax = fmin(tmax, fmax(t1.z, t2.z));

	return tmax >= tmin;
}


// this version gives the point of intersection as well as distance
// algorithm explanation here. hopefully my extrapolation into 3 dimensions is correct.
// http://tavianator.com/fast-branchless-raybounding-box-intersections/
int boxRayIntersect(const AABB* b, const Ray* r, Vector* ipoint, float* idist) {
	Vector t1, t2;
	float tmin, tmax;
	
	t1.x = (b->min.x - r->o.x) * r->id.x;
	t2.x = (b->max.x - r->o.x) * r->id.x;
	tmin = fmin(t1.x, t2.x);
	tmax = fmax(t1.x, t2.x);
	
	t1.y = (b->min.y - r->o.y) * r->id.y;
	t2.y = (b->max.y - r->o.y) * r->id.y;
	tmin = fmax(tmin, fmin(t1.y, t2.y));
	tmax = fmin(tmax, fmax(t1.y, t2.y));
	
	t1.z = (b->min.z - r->o.z) * r->id.z;
	t2.z = (b->max.z - r->o.z) * r->id.z;
	tmin = fmax(tmin, fmin(t1.z, t2.z));
	tmax = fmin(tmax, fmax(t1.z, t2.z));
	
	if(tmax < tmin) return 0;
	
	if(idist) *idist = tmin;
	
	if(ipoint) {
		ipoint->x = r->o.x + (r->d.x * tmin);
		ipoint->y = r->o.y + (r->d.y * tmin);
		ipoint->z = r->o.z + (r->d.z * tmin);
	}
	
	return 1;
}


// returns a local t value for the segment the normalized value falls into
static float bsNormalToLocalT2(BezierSpline2* bs, float normalT, int* segNum) {
	
	float segT = 1.0 / (bs->length - (!bs->isLoop));
	if(segNum) *segNum = (int)(normalT / segT);
	return fmod(normalT, segT) / segT;
}

// returns which segment a normalized t falls in
static int bsSegNum2(BezierSpline2* bs, float normalT) {
	
	float segT = 1.0 / (bs->length - (!bs->isLoop));
	return (int)(normalT / segT);
}


// returns a full set of control points for the segment t falls into
// out's order is e1, c1, c2, e2
static void bsSegmentForT2(BezierSpline2* bs, float normalT, Vector2* out) {
	BezierSplineSegment2* p, *n;
	int segN, i;
	
	segN = i = bsSegNum2(bs, normalT);
	
	p = bs->segments;
	while(i--) p = p->next;
	
	// a loop wraps around at the end
	n = (bs->isLoop && (segN == (bs->length - 1))) ? bs->segments : p->next;
	
	// end 1
	out[0].x = p->e.x;
	out[0].y = p->e.y;
	
	// control 1
	out[1].x = p->c.x;
	out[1].y = p->c.y;
	
	// control 2 - this one is reflected across e2
	vReflectAcross(&n->c, &n->e, &out[2]); 
	
	// end 2
	out[3].x = n->e.x;
	out[3].y = n->e.y;
}



	


// BUG: this function is probably horribly broken
// this function evaluates a spline with a [0,1] normalized value of t across the whole spline
void bsEvalPoint2(BezierSpline2* bs, float normalT, Vector2* out) {
	
	int segN;
	float localT;
	
	// find which spline segment the t is in
	localT = bsNormalToLocalT2(bs, normalT, &segN);  
	
	// find the local value of t
	Vector2 cp[4];
	bsSegmentForT2(bs, normalT, cp);
	
	evalBezier2D(&cp[0], &cp[3], &cp[1], &cp[2], localT, out);
}


// subdivide a spline into a set of line segments. linePoints must be allocated externally.
// this function is faster than more accurate ones.
void bsEvenLines(BezierSpline2* bs, int lineCount, Vector2* linePoints) {
	
	float tIncrement;
	
	
	tIncrement = (float)(bs->length - (!bs->isLoop)) / (float)lineCount;
	
	
	
	
	
}


