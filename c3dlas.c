
#include "c3dlas.h"
#include "string.h"
#include "math.h"




// vector operations


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
	n = (float)1.0 / sqrt((float)((v->x * v->x) + (v->y * v->y) + (v->z * v->z)));
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



// matrix operations



const Matrix IDENT_MATRIX = { { 1, 0, 0, 0, 
                                0, 1, 0, 0, 
                                0, 0, 1, 0, 
                                0, 0, 0, 1 } };


void mIdent(Matrix* m) {
	*m = IDENT_MATRIX;
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
	
	memcpy(b.m, out->m, sizeof(out->m));
	
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
	f = 1.0 / tan(fov / 2);
	
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




