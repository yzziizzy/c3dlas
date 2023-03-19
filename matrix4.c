


const Matrix IDENT_MATRIX = { { 1, 0, 0, 0,
                                0, 1, 0, 0,
                                0, 0, 1, 0,
                                0, 0, 0, 1 } };


void mIdent(Matrix* m) {
	*m = IDENT_MATRIX;
}

void mCopy(Matrix* in, Matrix* out) {
	*out = *in;
}


// out cannot overlap with a or b
// with restrict and -O2, this vectorizes nicely.
void mFastMul(Matrix* restrict a, Matrix* restrict b, Matrix* restrict out) {
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

// out cannot overlap with a. make a copy first if you want to do weird stuff.
void mMul(Matrix* restrict a, Matrix* restrict out) {
	Matrix b;
	
	mCopy(out, &b);
	
	mFastMul(a, &b, out);
}



void mTransv(Vector3* v, Matrix* out) {
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



void mScalev(Vector3* v, Matrix* out) {
	mScale3f(v->x, v->y, v->z, out);
}

void mScale3f(float x, float y, float z, Matrix* out) {
	Matrix t;
	
	t = IDENT_MATRIX;
	t.m[0] = x;
	t.m[5] = y;
	t.m[10] = z;

	mMul(&t, out);
}



void mRotv(Vector3* v, float theta, Matrix* out) {
	mRot3f(v->x, v->y, v->z, theta, out);
}

void mRot3f(float x, float y, float z, float theta, Matrix* out) {
	
	float costh, omcosth;
	float sinth;
	Matrix r;
	
	sincosf(theta, &sinth, &costh);
	omcosth = 1 - costh;
	
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


void mTransposeFast(Matrix* in, Matrix* out) {
	int i;
	
	for(i = 0; i < 4; i++) {
		out->m[i]      = in->m[i * 4];
		out->m[i + 4]  = in->m[(i * 4) + 1];
		out->m[i + 8]  = in->m[(i * 4) + 2];
		out->m[i + 12] = in->m[(i * 4) + 3];
	}
}

void mTranspose(Matrix* in, Matrix* out) {
	Matrix t;
	
	mTransposeFast(in, &t);
	
	*out = t;
}


void mAdd(Matrix* a, Matrix* b, Matrix* out) {
	for(int i = 0; i < 16; i++) {
		out->m[i] = a->m[i] + b->m[i];
	}
}

void mScalarMul(Matrix* a, float f, Matrix* out) {
	for(int i = 0; i < 16; i++) {
		out->m[i] = a->m[i] * f;
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
int mInverse(Matrix* restrict in, Matrix* restrict out) {
	
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
	f = 1.0 / tan(fov * DEG2RAD / 2.0);
	
	m.m[0] = f / aspect;
	m.m[5] = f;
	m.m[10] = (far + near) / (near - far);
	m.m[11] = -1.0;
	m.m[14] = (2.0 * far * near) / (near - far);
	m.m[15] = 0.0;

	mMul(&m, out);
}



// extract the near and far planes from a prespective matrix
void mPerspExtractNF(Matrix* m, double* near, double* far) {
	
	// perspective computations can be sensitive to precision.
	double a = m->m[10];
	double b = m->m[14];
	
	*far = b / (1.0 + a);
	*near = (-a * b) / (a - 1.0);
// 	printf("a: %f, b: %f, f: %f, n: %f\n", a, b, f, n);
}


// set the near and far planes for an existing prespective matrix
void mPerspSetNF_ZUp(Matrix* m, float near, float far) { // Z-up
	m->m[6] = (far + near) / (near - far);
	m->m[14] = (2.0 * far * near) / (near - far);
}

void mPerspSetNF(Matrix* m, float near, float far) { // Y-up	
	m->m[10] = (far + near) / (far - near);
	m->m[14] = (2.0 * far * near) / (far - near);
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



void mOrthoFromRadius(float r, Matrix* out) {
	Matrix m;
	
	float right = -r ;
	float left = r;
	float top = r;
	float bottom = -r;
	float near = -r;
	float far = r;
	
	*out = IDENT_MATRIX;
	out->m[0 + (0*4)] = 2. / (right - left);
	out->m[1 + (1*4)] = 2. / (top - bottom);
	out->m[2 + (2*4)] = -2. / (far - near);
	
	out->m[0 + (3*4)] = -(top + bottom) / (top - bottom);
	out->m[1 + (3*4)] = -(right + left) / (right - left);
	out->m[2 + (3*4)] = -(far + near) / (far - near);
	out->m[3 + (3*4)] = 1;
}	

//void mOrthoFromSphere(Sphere* s, Vector3* eyePos, Matrix* out) {
void mOrthoFromSphere(Sphere s, Vector3 eyePos, Vector3 up, Matrix* out) {
	Matrix m;

	mOrthoFromRadius(s.r, &m);

	Vector3 d = vNorm3(vSub3(s.center, eyePos));
	
	Matrix m2;
	mLookAt(eyePos, s.center, up, &m2);
	
	mFastMul(&m2, &m, out);
}

// extract the planes from an orthographic projection matrix.
void mOrthoExtractPlanes(Matrix* m, float* left, float* right, float* top, float* bottom, float* near, float* far) {
	
	*left = (m->m[12] + 1) / -m->m[0]; 
	*right = (-m->m[12] + 1) / m->m[0]; 
	*bottom = (m->m[13] + 1) / m->m[5]; 
	*top = (-m->m[13] + 1) / m->m[5]; 
	*near = (m->m[14] + 1) / m->m[10];
	*far = (m->m[14] - 1) / m->m[10];
}


void mOrthoSetNF(Matrix* m, float near, float far) {
	m->m[2 + (2*4)] = -2. / (far - near);
	m->m[2 + (3*4)] = -(far + near) / (far - near);
}

// analgous to gluLookAt
// up is not required to be orthogonal to anything, so long as it's not parallel to anything
// http://www.songho.ca/opengl/gl_camera.html#lookat
void mLookAt(Vector3 eye, Vector3 center, Vector3 up, Matrix* out) {
	
	Vector3 forward, left, upn;
	Matrix m;
	
	
	forward = vNorm3(vSub3(eye, center));
	left = vNorm3(vCross3(forward, up));
	upn = vCross3(left, forward);
	
	m.m[0+(0*4)] = left.x;
	m.m[1+(0*4)] = upn.x;
	m.m[2+(0*4)] = forward.x;
	m.m[3+(0*4)] = 0;

	m.m[0+(1*4)] = left.y;
	m.m[1+(1*4)] = upn.y;
	m.m[2+(1*4)] = forward.y;
	m.m[3+(1*4)] = 0;
	
	m.m[0+(2*4)] = left.z;
	m.m[1+(2*4)] = upn.z;
	m.m[2+(2*4)] = forward.z;
	m.m[3+(2*4)] = 0;
	
	m.m[0+(3*4)] = -left.x * center.x    - left.y * center.y    - left.z * center.z;
	m.m[1+(3*4)] = -upn.x * center.x     - upn.y * center.y     - upn.z * center.z;
	m.m[2+(3*4)] = -forward.x * center.x - forward.y * center.y - forward.z * center.z;
	m.m[3+(3*4)] = 1;
	
	*out = m;
}

void mLookDir(Vector3 eye, Vector3 center, Vector3 dir, Vector3 up, Matrix* out) {
	
	Vector3 forward, left, upn;
	Matrix m;
	
	
	forward = dir; //vNorm3(vSub3(eye, center));
	left = vNorm3(vCross3(forward, up));
	upn = vCross3(left, forward);
	
	m.m[0+(0*4)] = left.x;
	m.m[1+(0*4)] = upn.x;
	m.m[2+(0*4)] = forward.x;
	m.m[3+(0*4)] = 0;

	m.m[0+(1*4)] = left.y;
	m.m[1+(1*4)] = upn.y;
	m.m[2+(1*4)] = forward.y;
	m.m[3+(1*4)] = 0;
	
	m.m[0+(2*4)] = left.z;
	m.m[1+(2*4)] = upn.z;
	m.m[2+(2*4)] = forward.z;
	m.m[3+(2*4)] = 0;
	
	m.m[0+(3*4)] = -left.x * center.x    - left.y * center.y    - left.z * center.z;
	m.m[1+(3*4)] = -upn.x * center.x     - upn.y * center.y     - upn.z * center.z;
	m.m[2+(3*4)] = -forward.x * center.x - forward.y * center.y - forward.z * center.z;
	m.m[3+(3*4)] = 1;
	
	*out = m;
}


void mPrint(Matrix* m, FILE* f) {
	int r;
	
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
	
	mCopy(&ms->stack[ms->top - 1], &ms->stack[ms->top]);
	
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

void msTransv(Vector3* v, MatrixStack* ms) { // translation
	mTransv(v, msGetTop(ms));
}

void msTrans3f(float x, float y, float z, MatrixStack* ms) { // translation
	mTrans3f(x, y, z, msGetTop(ms));
}

void msScalev(Vector3* v, MatrixStack* ms) {
	mScalev(v, msGetTop(ms));
}

void msScale3f(float x, float y, float z, MatrixStack* ms) {
	mScale3f(x, y, z, msGetTop(ms));
}

void msRotv(Vector3* v, float theta, MatrixStack* ms) { // rotate about a vector
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

void msLookAt(Vector3* eye, Vector3* center, Vector3* up, MatrixStack* ms) {
	mLookAt(*eye, *center, *up, msGetTop(ms));
}



void msymPlaneFromTrip(Vector3* a, Vector3* b, Vector3* c, MatrixSym* m) {

	Vector3 n = vTriFaceNormal3(*a, *b, *c);
	float d = vDot3p(&n, a);
	
	m->m[0] = n.x * n.x;
	m->m[1] = n.x * n.y;
	m->m[2] = n.y * n.y;
	m->m[3] = n.x * n.z;
	m->m[4] = n.y * n.z;
	m->m[5] = n.z * n.z;
	m->m[6] = n.x * d;
	m->m[7] = n.y * d;
	m->m[8] = n.z * d;
	m->m[9] = d * d;
}


void msymAddp(MatrixSym* a, MatrixSym* b, MatrixSym* c) {
	for(int i = 0; i < 10; i++) c->m[i] = a->m[i] + b->m[i];
} 



/* 
 _            _
|  0  1  3  6  |
|  1  2  4  7  |
|  3  4  5  8  |
|_ 6  7  8  9 _|

*/

// v * m * v^t
float msymConjugateMulV3p(MatrixSym* m, Vector3* v) {
	float f0 = v->x * m->m[0] + v->x * m->m[1] + v->x * m->m[3] + v->x * m->m[6];
	float f1 = v->y * m->m[1] + v->y * m->m[2] + v->y * m->m[4] + v->y * m->m[7];
	float f2 = v->z * m->m[3] + v->z * m->m[4] + v->z * m->m[5] + v->z * m->m[8];
	float f3 =  1.0 * m->m[6] +  1.0 * m->m[7] +  1.0 * m->m[8] +  1.0 * m->m[9];
	
	return f0 * v->x + f1 * v->y + f2 * v->z + f3 * 1.0;
} 

