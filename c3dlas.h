
#ifndef __c3dlas_h__
#define __c3dlas_h__


#include <stdlib.h> // rand() et al.
#include <stdint.h> 
#include <math.h> // fmin/fmax

//#define C3DLAS_USE_SIMD 1

#define _0000b 0x00
#define _0001b 0x01
#define _0010b 0x02
#define _0011b 0x03
#define _0100b 0x04
#define _0101b 0x05
#define _0110b 0x06
#define _0111b 0x07
#define _1000b 0x08
#define _1001b 0x09
#define _1010b 0x0a
#define _1011b 0x0b
#define _1100b 0x0c
#define _1101b 0x0d
#define _1110b 0x0e
#define _1111b 0x0f

#define F_PI ((float)3.1415926535897932384626433832795028841971693993751)
#define D_PI ((double)3.1415926535897932384626433832795028841971693993751)
#define F_2PI ((float)6.2831853071795864769252867665590057683943387987502)
#define D_2PI ((double)6.2831853071795864769252867665590057683943387987502)
#define F_1_PI ((float)0.3183098861837906715377675267450287240689192914809)
#define D_1_PI ((double)0.3183098861837906715377675267450287240689192914809)
#define F_PI_2 ((float)1.5707963267948966192313216916397514420985846996875)
#define D_PI_2 ((double)1.5707963267948966192313216916397514420985846996875)
#define F_3PI_2 ((float)4.7123889803846898576939650749192543262957540990626)
#define D_3PI_2 ((double)4.7123889803846898576939650749192543262957540990626)

#define F_GOLDEN ((float)1.61803398874989484820458683436563811772030917980576f)
#define D_GOLDEN ((double)1.61803398874989484820458683436563811772030917980576)

#define RAD2DEG (57.29577951308232087679815481410517033240547246656432154916024386)
#define DEG2RAD (0.0174532925199432957692369076848861271344287188854172545609719144)

#define FLT_CMP_EPSILON 0.000001
#define FLT_CMP_EPSILON_SQ (FLT_CMP_EPSILON * FLT_CMP_EPSILON)

#define C3DLAS_COPLANAR  (0)
#define C3DLAS_FRONT     (1)
#define C3DLAS_BACK      (2)
#define C3DLAS_INTERSECT (3)
#define C3DLAS_DISJOINT  (4)
#define C3DLAS_PARALLEL  (5)

static const char* c3dlas_EnumString(int e) {
	switch(e) {
		case 0: return "C3DLAS_COPLANAR";
		case 1: return "C3DLAS_FRONT";
		case 2: return "C3DLAS_BACK";
		case 3: return "C3DLAS_INTERSECT";
		case 4: return "C3DLAS_DISJOINT";
		case 5: return "C3DLAS_PARALLEL";
		default: return "Unknown Code";
	}
}


#define MAX(a,b) ({ \
	__typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b); \
	_a > _b ? _a : _b; \
})
#define MIN(a,b) ({ \
	__typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b); \
	_a < _b ? _a : _b; \
})
#define MAXE(a,b) ({ \
	__typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b); \
	_a >= _b ? _a : _b; \
})
#define MINE(a,b) ({ \
	__typeof__ (a) _a = (a); \
	__typeof__ (b) _b = (b); \
	_a <= _b ? _a : _b; \
})


typedef struct {
	float x,y;
} Vector2;

typedef struct {
	float x,y,z;
} Vector;

typedef struct {
	float x,y,z,w;
} Vector4;

typedef struct Vector2i {
	int x,y;
} Vector2i;

typedef struct {
	Vector o; // origin
	Vector d; // normalized direction
} Ray;

typedef struct {
	Vector2 o; // origin
	Vector2 d; // normalized direction
} Ray2;

typedef struct {
	Vector start, end;
} LineSegment;


typedef struct BezierSplineSegment {
	Vector e, c; // end and control
	struct BezierSplineSegment* next;
} BezierSplineSegment;

typedef struct {
	int length;
	unsigned char isLoop;
	BezierSplineSegment* segments;
} BezierSpline;

typedef struct BezierSplineSegment2 {
	Vector2 e, c; // end and control
	struct BezierSplineSegment2* next;
} BezierSplineSegment2;

typedef struct {
	int length;
	unsigned char isLoop;
	BezierSplineSegment2* segments;
} BezierSpline2;

typedef struct {
	Vector n; // normal
	float d; // distance along normal to the origin
} Plane;

typedef struct {
	Vector center;
	float r;
} Sphere;

typedef struct {
	Plane planes[6]; // near, far, sides[4]
	Vector points[8]; // near then far
} Frustum;

typedef struct { // does not have to be coplanar
	Vector v[4];
} Quad;

typedef struct {
	Vector2 v[4];
} Quad2;

typedef struct {
	Vector2i v[4];
} Quad2i;

typedef struct {
	float m[16];
} Matrix;


typedef struct MatrixStack {
	short size;
	short top;
	Matrix* stack;
} MatrixStack;


// axis-aligned bounding box
typedef struct AABB {
	Vector min;
	Vector max;
} AABB;

typedef struct AABB2 {
	Vector2 min;
	Vector2 max;
} AABB2;


typedef struct AABB2i {
	Vector2i min;
	Vector2i max;
} AABB2i;


extern const Matrix IDENT_MATRIX;

// utilities

uint32_t bitReverse32(uint32_t x);
uint32_t reverseBits(uint32_t n, int len);

float pcg_f(uint64_t* state, uint64_t stream);

static inline float frand(float low, float high) {
	return low + ((high - low) * ((float)rand() / (float)RAND_MAX));
}

static inline float frandNorm(void) {
	return ((float)rand() / (float)RAND_MAX);
}

static inline double drand(double low, double high) {
	return low + ((high - low) * ((double)rand() / (double)RAND_MAX));
}

static inline double drandNorm(void) {
	return ((double)rand() / (double)RAND_MAX);
}

static inline float fclamp(float val, float min, float max) {
	return fmin(max, fmax(min, val));
}

static inline float fclampNorm(float val) {
	return fclamp(val, 0.0f, 1.0f);
}

static inline int iclamp(int val, int min, int max) {
	return MIN(max, MAX(min, val));
}

static inline float flerp(float a, float b, float t) {
	return a  + ((b - a) * t);
}

static inline float flerp2D(float xx, float xy, float yx, float yy, float xt, float yt) {
	float a = xx  + ((xy - xx) * yt);
	float b = yx  + ((yy - yx) * yt);
	return a  + ((b - a) * xt);
}


// vectors
int   vEq(Vector* a, Vector* b); // safe equivalence, to FLT_CMP_EPSILON
int   vEqEp(Vector* a, Vector* b, float epsilon); // safe equivalence, to arbitrary epsilon
void  vCopy(const Vector* src, Vector* dst); // copy vector values
void  vSwap(Vector* a, Vector* b); // swap two vectors
void  vAdd(Vector* a, Vector* b, Vector* out); // add two vectors
void  vSub(Vector* from, Vector* what, Vector* diff); // diff = from - what
void  vScale(Vector* v, float scalar, Vector* out); // scalar muliplication
void  vLerp(Vector* a, Vector* b, float t, Vector* out); // Linear interpolation between two vectors
void  vInverse(const Vector* v, Vector* out); // inverse
float vMag(const Vector* v); // return the magnitude
float vDot(const Vector* a, const Vector* b); // dot product
float vDist(Vector* from, Vector* to); // distance from one point to another
float vDistSq(Vector* from, Vector* to); // squared distance from one point to another
void  vNorm(Vector* v, Vector* out); // normalize the vector
void  vUnit(Vector* v, Vector* out); // normalise the vector, alternate name
void  vCross(Vector* a, Vector* b, Vector* out); // cross product: out = a x b
float vScalarTriple(Vector* a, Vector* b, Vector* c); // scalar triple product: a . (b x c)
void  vProject(Vector* what, Vector* onto, Vector* out); // slower; onto may not be normalized
void  vProjectNorm(Vector* what, Vector* onto, Vector* out); // faster; onto must be normalized
void  vMin(Vector* a, Vector* b, Vector* out); // returns the minimum values of each component
void  vMax(Vector* a, Vector* b, Vector* out); // returns the maximum values of each component
void  vSet(float x, float y, float z, Vector* out);
void  vPointAvg(Vector* a, Vector* b, Vector* out);

void vRandom(Vector* end1, Vector* end2, Vector* out);
void vRandomNorm(Vector* out);

void  vLerp4(Vector4* a, Vector4* b, float t, Vector4* out); // Linear interpolation between two vectors


// http://geomalgorithms.com/a07-_distance.html
// _PARALLEL with no output on parallel lines
// _INTERSECT with one point of output on intersection
// _DISJOINT with two outputs otherwise
int shortestLineFromRayToRay(Ray* r1, Ray* r2, Vector* pOut);

// reflects the distance from v to pivot across pivot.
// out, pivot, and v will form a straight line with pivot exactly in the middle.
void  vReflectAcross(Vector* v, Vector* pivot, Vector* out);
void  vTriFaceNormal(Vector* a, Vector* b, Vector* c, Vector* out); // returns a normalized face normal for the given triangle
void  vpTriFaceNormal(Vector* tri, Vector* out); // returns a normalized face normal for the given triangle

void  vProjectOntoPlane(Vector* v, Plane* p, Vector* out);
void  vProjectOntoPlaneNormalized(Vector* v, Plane* p, Vector* out);

void  planeFromTriangle(Vector* v1, Vector* v2, Vector* v3, Plane* out); // calculates a plane form a triangle
void  planeCopy(Plane* in, Plane* out); // copy a plane
void  planeInverse(Plane* in, Plane* out); // flips the plane's direction
int   planeClassifyPoint(Plane* p, Vector* pt); // classifies a point by which side of the plane it's on, default espilon
int   planeClassifyPointEps(Plane* p, Vector* pt, float epsilon); // classifies a point by which side of the plane it's on, custom espilon
// closest distance from an arbitrary point to the plane 
float planePointDist(Plane* pl, Vector* p);
// signed closest distance from an arbitrary point to the plane 
float planePointDistSigned(Plane* pl, Vector* p);

// C3DLAS_INTERSECT, _COPLANAR or _DISJOINT
int planeLineFindIntersect(Plane* pl, Vector* la, Vector* lb, Vector* out);

// Assumes full proper intersection.
// C3DLAS_INTERSECT
int planeLineFindIntersectFast(Plane* pl, Vector* la, Vector* lb, Vector* out);

// C3DLAS_COPLANAR, _INTERSECT, or _DISJOINT
int triPlaneTestIntersect(Vector* pTri, Plane* pl);

// C3DLAS_COPLANAR, _INTERSECT, or _DISJOINT
int triPlaneClip(
	Vector* pTri, 
	Plane* pl, 
	Vector* aboveOut, 
	Vector* belowOut, 
	int* aboveCnt,
	int* belowCnt
);

// C3DLAS_COPLANAR, _PARALLEL, _INTERSECT, or _DISJOINT
// aboveCnt and belowCnt are always set.
int linePlaneClip(
	Vector* la, 
	Vector* lb, 
	Plane* pl, 
	Vector* aboveOut, 
	Vector* belowOut,
	int* aboveCnt,
	int* belowCnt
);

void frustumCenter(Frustum* f, Vector* out);
void frustumBoundingSphere(Frustum* f, Sphere* out);

void quadCenterp(Vector* a, Vector* b, Vector* c, Vector* d, Vector* out);

void frustumFromMatrix(Matrix* m, Frustum* out);
void frustumInnerBoundingSphere(Frustum* f, Sphere* out);
void frustumOuterBoundingSphere(Frustum* f, Sphere* out);

// 2d vector stuff, same as 3d except one less d
int   vEq2(Vector2* a, Vector2* b); // safe equivalence, to FLT_CMP_EPSILON
int   vEqEp2(Vector2* a, Vector2* b, float epsilon); // safe equivalence, to arbitrary epsilon
void  vCopy2(const Vector2* src, Vector2* dst); // copy vector values
void  vSwap2(Vector2* a, Vector2* b); // swap two vectors
void  vAdd2(Vector2* a, Vector2* b, Vector2* out); // add two vectors
void  vSub2(Vector2* from, Vector2* what, Vector2* diff); // diff = from - what
void  vScale2(Vector2* v, float scalar, Vector2* out); // scalar muliplication
float  vDist2(Vector2* a, Vector2* b); // distance between points
void  vLerp2(Vector2* a, Vector2* b, float t, Vector2* out); // linear interpolation
void  vInverse2(const Vector2* v, Vector2* out); // inverse
float vMag2(Vector2* v); // return the magnitude
float vDot2(Vector2* a, Vector2* b); // dot product
void  vNorm2(Vector2* v, Vector2* out); // normalize the vector
void  vUnit2(Vector2* v, Vector2* out); // normalise the vector, alternate name
void  vMin2(Vector2* a, Vector2* b, Vector2* out); // returns the minimum values of each component
void  vMax2(Vector2* a, Vector2* b, Vector2* out); // returns the maximum values of each component
void  vSet2(float x, float y, Vector2* out);

// reflects the distance from v to pivot across pivot.
// out, pivot, and v will form a straight line with pivot exactly in the middle.
void  vReflectAcross2(Vector2* v, Vector2* pivot, Vector2* out);

// degenerate cases may not give desired results. GIGO.
void  vRoundAway2(const Vector2* in, const Vector2* center, Vector2i* out);
void  vRoundToward2(const Vector2* in, const Vector2* center, Vector2i* out);

// returns the *signed* area of a triangle. useful for determining winding
// positive values mean a clockwise triangle
float triArea2(Vector2* a, Vector2* b, Vector2* c);

// determines if a point is inside a triangle
int triPointInside2(Vector2* p, Vector2* a, Vector2* b, Vector2* c);

// 2d integer vector stuff
int   vEq2i(Vector2i* a, Vector2i* b);
void  vCopy2i(const Vector2i* src, Vector2i* dst); // copy vector values
void  vSwap2i(Vector2i* a, Vector2i* b); // swap two vectors
void  vAdd2i(Vector2i* a, Vector2i* b, Vector2i* out); // add two vectors
void  vSub2i(Vector2i* from, Vector2i* what, Vector2i* diff); // diff = from - what
void  vScale2i(Vector2i* v, int scalar, Vector2i* out); // scalar muliplication
int   vDot2i(Vector2i* a, Vector2i* b); // dot product
void  vMin2i(Vector2i* a, Vector2i* b, Vector2i* out); // returns the minimum values of each component
void  vMax2i(Vector2i* a, Vector2i* b, Vector2i* out); // returns the maximum values of each component
void  vSet2i(int x, int y, Vector2i* out);
float vDist2i(Vector2i* a, Vector2i* b); // returns the absolute distance between two vectors





float pvDist(Plane* p, Vector* v);

void vMatrixMul(Vector* in, Matrix* m, Vector* out); // multiply a vector by a matrix
void vMatrixMulf(float x, float y, float z, Matrix* m, Vector* out); // multiply a vector by a matrix


void mIdent(Matrix* m); // set m to the identity matrix
void mCopy(Matrix* in, Matrix* out);
void mFastMul(Matrix* a, Matrix* b, Matrix* out); // a and b cannot also be out. mostly internal use.
void mMul(Matrix* a, Matrix* out); // makes a copy of out before multiplying over it
void mTransv(Vector* v, Matrix* out); // translation
void mTrans3f(float x, float y, float z, Matrix* out); // translation
void mScalev(Vector* v, Matrix* out);
void mScale3f(float x, float y, float z, Matrix* out);
void mRotv(Vector* v, float theta, Matrix* out); // rotate about a vector
void mRot3f(float x, float y, float z, float theta, Matrix* out); // rotate about a vector
void mRotX(float theta, Matrix* out); //
void mRotY(float theta, Matrix* out); // rotate about axes
void mRotZ(float theta, Matrix* out); //
void mTranspose(Matrix* in, Matrix* out);
void mTransposeFast(Matrix* in, Matrix* out); // in cannot be out
float mDeterminate(Matrix* m);
int mInverse(Matrix* in, Matrix* out); // returns 0 on success, 1 if there is no inverse; out remains unchanged


// analogous to glFrustum
// no div/0 checking here for right == left etc. just don't be an idiot.
void mFrustum(float left, float right, float top, float bottom, float near, float far, Matrix* out);

// analogous to gluPerspective
// same div/0 warnings apply. if you get an FP exception you deserve it.
// use a double for fov; the precision matters often.
// https://www.opengl.org/archives/resources/faq/technical/transformations.htm
// https://www.opengl.org/sdk/docs/man2/xhtml/gluPerspective.xml
void mPerspective(double fov, float aspect, float near, float far, Matrix* out);

// extract the near and far planes from a prespective matrix
void mPerspExtractNF(Matrix* m, double* near, double* far);

// set the near and far planes for an existing prespective matrix
void mPerspSetNF(Matrix* m, float near, float far);

// orthographic projection. use this for a "2D" look.
// same div/0 warnings.
void mOrtho(float left, float right, float top, float bottom, float near, float far, Matrix* out);

// calculates an orthographic matrix that encloses the sphere, looking from eyePos
void mOrthoFromSphere(Sphere* s, Vector* eyePos, Matrix* out);

// extract the planes from an orthographic projection matrix.
void mOrthoExtractPlanes(Matrix* m, float* left, float* right, float* top, float* bottom, float* near, float* far);


// analgous to gluLookAt
// https://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
void mLookAt(Vector* eye, Vector* center, Vector* up, Matrix* out);

void mPrint(Matrix* m, FILE* f);



// matrix stack functions

// make sure you allocate enough. when it's out, it's out. no surprise mallocs later on. (yet)
void msAlloc(int size, MatrixStack* ms);
void msFree(MatrixStack* ms);

int msPush(MatrixStack* ms);
void msPop(MatrixStack* ms);
Matrix* msGetTop(MatrixStack* ms);

void msPrintAll(MatrixStack* ms, FILE* f);

// these are all wrappers around the functions listed above
void msIdent(MatrixStack* ms); // set to the identity matrix
void msCopy(Matrix* in, MatrixStack* ms);
void msMul(Matrix* a, MatrixStack* ms); // makes a copy of out before multiplying over it
void msTransv(Vector* v, MatrixStack* ms); // translation
void msTrans3f(float x, float y, float z, MatrixStack* ms); // translation
void msScalev(Vector* v, MatrixStack* ms);
void msScale3f(float x, float y, float z, MatrixStack* ms);
void msRotv(Vector* v, float theta, MatrixStack* ms); // rotate about a vector
void msRot3f(float x, float y, float z, float theta, MatrixStack* ms); // rotate about a vector
void msFrustum(float left, float right, float top, float bottom, float near, float far, MatrixStack* ms);
void msPerspective(double fov, float aspect, float near, float far, MatrixStack* ms);
void msOrtho(float left, float right, float top, float bottom, float near, float far, MatrixStack* ms);
void msLookAt(Vector* eye, Vector* center, Vector* up, MatrixStack* ms);


// cubic Bezier curves
void evalBezier(Vector* e1, Vector* e2, Vector* c1, Vector* c2, float t, Vector* out);
void evalBezierTangent(Vector* e1, Vector* e2, Vector* c1, Vector* c2, float t, Vector* out); // tangent vector; not normalized
void evalBezierNorm(Vector* e1, Vector* e2, Vector* c1, Vector* c2, float t, Vector* out); // normal vector; not normalized
float evalBezier1D(float e1, float e2, float c1, float c2, float t);
float evalBezier1D_dt(float e1, float e2, float c1, float c2, float t); // first derivative with respect to t
float evalBezier1D_ddt(float e1, float e2, float c1, float c2, float t); // second derivative with respect to t

// quadratic Bezier curves
float evalQBezier1D(float e1, float e2, float c1, float t);
void evalQBezier2D(Vector2* e1, Vector2* e2, Vector2* c1, float t, Vector2* out);
void evalQBezier(Vector* e1, Vector* e2, Vector* c1, float t, Vector* out);

///// bounding box functions

// 3D versions
int boxDisjoint(const AABB* a, const AABB* b);
int boxOverlaps(const AABB* a, const AABB* b);
int boxContainsPoint(const AABB* b, const Vector* p);

void boxCenter(const AABB* b, Vector* out); // calculates the center of the box
void boxSize(const AABB* b, Vector* out); // calculates the size of the box


void makeRay(Vector* origin, Vector* direction, Ray* out);
int boxRayIntersectFast(const AABB* b, const Ray* r);
int boxRayIntersect(const AABB* b, const Ray* r, Vector* ipoint, float* idist);


// 2D versions
int boxDisjoint2(const AABB2* a, const AABB2* b);
int boxOverlaps2(const AABB2* a, const AABB2* b);
int boxContainsPoint2(const AABB2* b, const Vector2* p);

void boxCenter2(const AABB2* b, Vector2* out); // calcuates the center of the box
void boxSize2(const AABB2* b, Vector2* out); // calculates the size of the box
void boxQuadrant2(const AABB2* in, char ix, char iy, AABB2* out);

// 2D integer versions
int boxDisjoint2i(const AABB2i* a, const AABB2i* b);
int boxOverlaps2i(const AABB2i* a, const AABB2i* b);
int boxContainsPoint2i(const AABB2i* b, const Vector2i* p);

void boxCenter2i(const AABB2i* b, Vector2* out); // calcuates the center of the box
void boxSize2i(const AABB2i* b, Vector2* out); // calculates the size of the box
void boxQuadrant2i(const AABB2i* in, char ix, char iy, AABB2i* out);

// find the center of a quad
void quadCenter2(const Quad2* in, Vector2* out);
void quadRoundOutward2(const Quad2* in, Quad2i* out);
void quadRoundInward2(const Quad2* in, Quad2i* out);


float   evalCubicHermite1D(float t, float p0, float p1, float m0, float m1);
Vector2 evalCubicHermite2D(float t, Vector2 p0, Vector2 p1, Vector2 m0, Vector2 m1); 
Vector  evalCubicHermite3D(float t, Vector p0, Vector p1, Vector m0, Vector m1);

#endif // __c3dlas_h__

