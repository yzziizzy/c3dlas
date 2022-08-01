
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


typedef struct Vector2i {
	int x,y;
} Vector2i;

typedef struct {
	//    cartesian  polar
	union { float x, rho; };
	union { float y, theta; };
} Vector2;


// The spherical coordinates use the "mathematical" conventions as opposed to the "physics" conventions
// This is because of the struct layout convenience for overlapping with 2D polar coordinates. 

typedef struct {
	//    cartesian  spherical  color
	union { float x, rho,       r; }; // rho is the radius
	union { float y, theta,     g; }; // rotation in the X-Y plane
	union { float z, phi,       b; }; // rotation in the plane passing through the Z axis
} Vector3;

typedef struct {
	//    cartesian  spherical  color
	union { float x, rho,       r; }; // rho is the radius
	union { float y, theta,     g; }; // rotation in the X-Y plane
	union { float z, phi,       b; }; // rotation in the plane passing through the Z axis
	union { float w,            a; };
} Vector4;


// Rays, but also infinite lines (mathematical lines) because there is no practical
//   difference besides whether you conside the ray to be one-sided or not.
typedef struct {
	Vector3 o; // origin
	Vector3 d; // normalized direction
} Ray3;

typedef struct {
	Vector2 o; // origin
	Vector2 d; // normalized direction
} Ray2;

// Line *segments*
typedef struct {
	Vector2 start, end;
} Line2;

typedef struct {
	Vector3 start, end;
} Line3;

// Polygons are 2-dimensional by definition
typedef struct Polygon {
	Vector2* points;
	long pointAlloc;
	long pointCount;
	
	Vector2 centroid;
	float maxRadiusSq; // squared distance from the centroid to the furthest point
} Polygon;

typedef struct BezierSplineSegment3 {
	Vector3 e, c; // end and control
	struct BezierSplineSegment3* next;
} BezierSplineSegment3;

typedef struct {
	int length;
	unsigned char isLoop;
	BezierSplineSegment3* segments;
} BezierSpline3;

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
	Vector3 n; // normal
	float d; // distance along normal to the origin
} Plane;

typedef struct {
	Vector3 center;
	float r;
} Sphere;

typedef struct {
	Plane planes[6]; // near, far, sides[4]
	Vector3 points[8]; // near then far
} Frustum;

typedef struct {
	Vector2i v[4];
} Quad2i;

typedef struct {
	Vector2 v[4];
} Quad2;

typedef struct { // does not have to be coplanar
	Vector3 v[4];
} Quad3;


typedef struct {
	float m[16];
} Matrix;

typedef struct MatrixStack {
	short size;
	short top;
	Matrix* stack;
} MatrixStack;


// axis-aligned bounding box
typedef struct AABB2i {
	Vector2i min;
	Vector2i max;
} AABB2i;

typedef struct AABB2 {
	Vector2 min;
	Vector2 max;
} AABB2;

typedef struct AABB3 {
	Vector3 min;
	Vector3 max;
} AABB3;





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

static inline float frandNorm2(void) {
	return ((float)rand() / (float)RAND_MAX) * 2.0 - 1.0;
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



//
// Vectors
//


// exact equivalence
int vEq2i(Vector2i a, Vector2i b); 
int vEqExact2i(Vector2i a, Vector2i b); 
int vEqExact2(Vector2 a, Vector2 b); 
int vEqExact3(Vector3 a, Vector3 b); 
int vEqExact4(Vector4 a, Vector4 b); 
int vEq2ip(Vector2i* a, Vector2i* b); 
int vEqExact2ip(Vector2i* a, Vector2i* b); 
int vEqExact2p(Vector2* a, Vector2* b); 
int vEqExact3p(Vector3* a, Vector3* b); 
int vEqExact4p(Vector4* a, Vector4* b); 

// safe equivalence, to FLT_CMP_EPSILON
int vEq2(Vector2 a, Vector2 b); 
int vEq3(Vector3 a, Vector3 b); 
int vEq4(Vector4 a, Vector4 b); 
int vEq2p(Vector2* a, Vector2* b); 
int vEq3p(Vector3* a, Vector3* b); 
int vEq4p(Vector4* a, Vector4* b); 

// safe equivalence, to arbitrary epsilon
int vEqEp2i(Vector2i a, Vector2i b, double epsilon);
int vEqEp2(Vector2 a, Vector2 b, float epsilon);
int vEqEp3(Vector3 a, Vector3 b, float epsilon);
int vEqEp4(Vector4 a, Vector4 b, float epsilon);
int vEqEp2ip(Vector2i* a, Vector2i* b, double epsilon);
int vEqEp2p(Vector2* a, Vector2* b, float epsilon);
int vEqEp3p(Vector3* a, Vector3* b, float epsilon);
int vEqEp4p(Vector4* a, Vector4* b, float epsilon);

// Swap two vectors
void vSwap2ip(Vector2i* a, Vector2i* b);
void vSwap2p(Vector2* a, Vector2* b);
void vSwap3p(Vector3* a, Vector3* b);
void vSwap4p(Vector4* a, Vector4* b);

// Vector addition
Vector2i vAdd2i(Vector2i a, Vector2i b);
Vector2  vAdd2(Vector2 a, Vector2 b);
Vector3  vAdd3(Vector3 a, Vector3 b);
Vector4  vAdd4(Vector4 a, Vector4 b);
void     vAdd2ip(Vector2i* a, Vector2i* b, Vector2i* out);
void     vAdd2p(Vector2* a, Vector2* b, Vector2* out);
void     vAdd3p(Vector3* a, Vector3* b, Vector3* out);
void     vAdd4p(Vector4* a, Vector4* b, Vector4* out);

// Vector subtraction. diff = from - what
Vector2i vSub2i(Vector2i from, Vector2i what);
Vector2  vSub2(Vector2 from, Vector2 what);
Vector3  vSub3(Vector3 from, Vector3 what);
Vector4  vSub4(Vector4 from, Vector4 what);
void     vSub2ip(Vector2i* from, Vector2i* what, Vector2i* diff);
void     vSub2p(Vector2* from, Vector2* what, Vector2* diff);
void     vSub3p(Vector3* from, Vector3* what, Vector3* diff);
void     vSub4p(Vector4* from, Vector4* what, Vector4* diff);

// Scalar muliplication
Vector2  vScale2(Vector2 v, float scalar);
Vector3  vScale3(Vector3 v, float scalar);
void     vScale2ip(Vector2i* v, float scalar, Vector2i* out);
void     vScale2p(Vector2* v, float scalar, Vector2* out);
void     vScale3p(Vector3* v, float scalar, Vector3* out);

// Dot product (inner product)
double vDot2i(const Vector2i a, const Vector2i b);
float  vDot2(const Vector2 a, const Vector2 b);
float  vDot3(const Vector3 a, const Vector3 b);
float  vDot4(const Vector4 a, const Vector4 b);
double vDot2ip(const Vector2i* a, const Vector2i* b);
float  vDot2p(const Vector2* a, const Vector2* b);
float  vDot3p(const Vector3* a, const Vector3* b);
float  vDot4p(const Vector4* a, const Vector4* b);

// Cross product: out = a x b
Vector3 vCross3(Vector3 a, Vector3 b);
void  vCross3p(Vector3* a, Vector3* b, Vector3* out); 

// Scalar triple product: a . (b x c)
float vScalarTriple3(Vector3 a, Vector3 b, Vector3 c);
float vScalarTriple3p(Vector3* a, Vector3* b, Vector3* c);

// Linear interpolation between two vectors
Vector2  vLerp2(Vector2 a, Vector2 b, float t);
Vector3  vLerp3(Vector3 a, Vector3 b, float t);
Vector4  vLerp4(Vector4 a, Vector4 b, float t); 
void     vLerp2p(Vector2* a, Vector2* b, float t, Vector2* out);
void     vLerp3p(Vector3* a, Vector3* b, float t, Vector3* out);
void     vLerp4p(Vector4* a, Vector4* b, float t, Vector4* out);

// Vector Inverse. Returns FLT_MAX on div/0
Vector2 vInverse2(const Vector2 v); 
Vector3 vInverse3(const Vector3 v); 
Vector4 vInverse4(const Vector4 v); 
void    vInverse2p(const Vector2* v, Vector2* out); 
void    vInverse3p(const Vector3* v, Vector3* out); 
void    vInverse4p(const Vector4* v, Vector4* out); 

// Vector magnitude (length)
double vMag2i(const Vector2i v);
float  vMag2(const Vector2 v);
float  vMag3(const Vector3 v);
float  vMag4(const Vector4 v);
double vMag2ip(const Vector2i* v);
float  vMag2p(const Vector2* v);
float  vMag3p(const Vector3* v);
float  vMag4p(const Vector4* v);

// Squared distance from one point to another
double vDistSq2i(Vector2i a, Vector2i b); 
float  vDistSq2(Vector2 a, Vector2 b); 
float  vDistSq3(Vector3 a, Vector3 b); 
float  vDistSq4(Vector4 a, Vector4 b); 
double vDistSq2ip(Vector2i* a, Vector2i* b); 
float  vDistSq2p(Vector2* a, Vector2* b); 
float  vDistSq3p(Vector3* a, Vector3* b); 
float  vDistSq4p(Vector4* a, Vector4* b); 

// Distance from one point to another
double vDist2i(Vector2i a, Vector2i b); 
float  vDist2(Vector2 a, Vector2 b); 
float  vDist3(Vector3 a, Vector3 b); 
float  vDist4(Vector4 a, Vector4 b); 
double vDist2ip(Vector2i* a, Vector2i* b); 
float  vDist2p(Vector2* a, Vector2* b); 
float  vDist3p(Vector3* a, Vector3* b); 
float  vDist4p(Vector4* a, Vector4* b); 

// Vector normalize (scale to unit length)
Vector2 vNorm2(Vector2 v);
Vector3 vNorm3(Vector3 v);
Vector4 vNorm4(Vector4 v);
void    vNorm2p(const Vector2* v, Vector2* out);
void    vNorm3p(const Vector3* v, Vector3* out);
void    vNorm4p(const Vector4* v, Vector4* out);
// alternate name
Vector2 vUnit2(Vector2 v);
Vector3 vUnit3(Vector3 v);
Vector4 vUnit4(Vector4 v);
void    vUnit2p(const Vector2* v, Vector2* out); 
void    vUnit3p(const Vector3* v, Vector3* out); 
void    vUnit4p(const Vector4* v, Vector4* out); 

// Returns the minimum values of each component
Vector2i vMin2i(Vector2i a, Vector2i b);
Vector2  vMin2(Vector2 a, Vector2 b);
Vector3  vMin3(Vector3 a, Vector3 b);
Vector4  vMin4(Vector4 a, Vector4 b);
void     vMin2ip(Vector2i* a, Vector2i* b, Vector2i* out);
void     vMin2p(Vector2* a, Vector2* b, Vector2* out);
void     vMin3p(Vector3* a, Vector3* b, Vector3* out);
void     vMin4p(Vector4* a, Vector4* b, Vector4* out);

// Returns the maximum values of each component
Vector2i vMax2i(Vector2i a, Vector2i b);
Vector2  vMax2(Vector2 a, Vector2 b);
Vector3  vMax3(Vector3 a, Vector3 b);
Vector4  vMax4(Vector4 a, Vector4 b);
void     vMax2ip(Vector2i* a, Vector2i* b, Vector2i* out);
void     vMax2p(Vector2* a, Vector2* b, Vector2* out);
void     vMax3p(Vector3* a, Vector3* b, Vector3* out);
void     vMax4p(Vector4* a, Vector4* b, Vector4* out);

// Cartesian to Spherical
Vector3 vC2S3(Vector3 cart);

// Spherical to Cartesian
Vector3 vS2C3(Vector3 s);


// Distance from a point to a line segment 
float vDistPointLine2(Vector2 p, Line2 ls);
float vDistPointLine3(Vector3 p, Line3 ls);

// Also returns the normalized distance along the line to the closest point
float vDistTPointLine2(Vector2 p, Line2 ls, float* T);
float vDistTPointLine3(Vector3 p, Line3 ls, float* T);


int vInsidePolygon(Vector2 p, Polygon* poly);

// Returns the distance from p to the closest point on the polygon.
// Interior distances are negative
float vDistPolygon(Vector2 p, Polygon* poly);
void polyCalcCentroid(Polygon* poly);
void polyCalcRadiusSq(Polygon* poly) ;

void  vProject3p(Vector3* what, Vector3* onto, Vector3* out); // slower; onto may not be normalized
void  vProjectNorm3p(Vector3* what, Vector3* onto, Vector3* out); // faster; onto must be normalized
void  vPointAvg3p(Vector3* a, Vector3* b, Vector3* out);

void vRandom3p(Vector3* end1, Vector3* end2, Vector3* out);
void vRandomNorm3p(Vector3* out);



// http://geomalgorithms.com/a07-_distance.html
// _PARALLEL with no output on parallel lines
// _INTERSECT with one point of output on intersection
// _DISJOINT with two outputs otherwise
int shortestLineFromRayToRay3p(Ray3* r1, Ray3* r2, Vector3* pOut);

// reflects the distance from v to pivot across pivot.
// out, pivot, and v will form a straight line with pivot exactly in the middle.
void  vReflectAcross3p(Vector3* v, Vector3* pivot, Vector3* out);
void  vTriFaceNormal3p(Vector3* a, Vector3* b, Vector3* c, Vector3* out); // returns a normalized face normal for the given triangle
void  vpTriFaceNormal3p(Vector3* tri, Vector3* out); // returns a normalized face normal for the given triangle

void  vProjectOntoPlane3p(Vector3* v, Plane* p, Vector3* out);
void  vProjectOntoPlaneNormalized3p(Vector3* v, Plane* p, Vector3* out);

void  planeFromTriangle3p(Vector3* v1, Vector3* v2, Vector3* v3, Plane* out); // calculates a plane form a triangle
void  planeCopy3p(Plane* in, Plane* out); // copy a plane
void  planeInverse3p(Plane* in, Plane* out); // flips the plane's direction
int   planeClassifyPoint3p(Plane* p, Vector3* pt); // classifies a point by which side of the plane it's on, default espilon
int   planeClassifyPointEps3p(Plane* p, Vector3* pt, float epsilon); // classifies a point by which side of the plane it's on, custom espilon
// closest distance from an arbitrary point to the plane 
float planePointDist3p(Plane* pl, Vector3* p);
// signed closest distance from an arbitrary point to the plane 
float planePointDistSigned3p(Plane* pl, Vector3* p);

// C3DLAS_INTERSECT, _COPLANAR or _DISJOINT
int planeLineFindIntersect3p(Plane* pl, Vector3* la, Vector3* lb, Vector3* out);

// Assumes full proper intersection.
// C3DLAS_INTERSECT
int planeLineFindIntersectFast3p(Plane* pl, Vector3* la, Vector3* lb, Vector3* out);

// C3DLAS_INTERSECT, or _DISJOINT
int IntersectPlaneRay3p(Plane* p, Ray3* r, Vector3* out);


// C3DLAS_COPLANAR, _INTERSECT, or _DISJOINT
int triPlaneTestIntersect3p(Vector3* pTri, Plane* pl);

// C3DLAS_COPLANAR, _INTERSECT, or _DISJOINT
int triPlaneClip3p(
	Vector3* pTri, 
	Plane* pl, 
	Vector3* aboveOut, 
	Vector3* belowOut, 
	int* aboveCnt,
	int* belowCnt
);

// C3DLAS_COPLANAR, _PARALLEL, _INTERSECT, or _DISJOINT
// aboveCnt and belowCnt are always set.
int linePlaneClip3p(
	Vector3* la, 
	Vector3* lb, 
	Plane* pl, 
	Vector3* aboveOut, 
	Vector3* belowOut,
	int* aboveCnt,
	int* belowCnt
);

void frustumCenter(Frustum* f, Vector3* out);
void frustumBoundingSphere(Frustum* f, Sphere* out);

void quadCenterp3p(Vector3* a, Vector3* b, Vector3* c, Vector3* d, Vector3* out);

void frustumFromMatrix(Matrix* m, Frustum* out);
void frustumInnerBoundingSphere(Frustum* f, Sphere* out);
void frustumOuterBoundingSphere(Frustum* f, Sphere* out);


// reflects the distance from v to pivot across pivot.
// out, pivot, and v will form a straight line with pivot exactly in the middle.
void  vReflectAcross2p(Vector2* v, Vector2* pivot, Vector2* out);

// degenerate cases may not give desired results. GIGO.
void  vRoundAway2p(const Vector2* in, const Vector2* center, Vector2i* out);
void  vRoundToward2p(const Vector2* in, const Vector2* center, Vector2i* out);

// returns the *signed* area of a triangle. useful for determining winding
// positive values mean a clockwise triangle
float triArea2p(Vector2* a, Vector2* b, Vector2* c);

// determines if a point is inside a triangle
int triPointInside2p(Vector2* p, Vector2* a, Vector2* b, Vector2* c);






float pvDist3p(Plane* p, Vector3* v);

void vMatrixMul3p(Vector3* in, Matrix* m, Vector3* out); // multiply a vector by a matrix
void vMatrixMulf3p(float x, float y, float z, Matrix* m, Vector3* out); // multiply a vector by a matrix
Vector3 vMatrixMul3(Vector3 in, Matrix* m);
Vector4 vMatrixMul4(Vector4 in, Matrix* m);
Vector3 vMatrixMulProjectedMagic3(Vector3 in, Matrix* m);

void mIdent(Matrix* m); // set m to the identity matrix
void mCopy(Matrix* in, Matrix* out);
void mFastMul(Matrix* a, Matrix* b, Matrix* out); // a and b cannot also be out. mostly internal use.
void mMul(Matrix* a, Matrix* out); // makes a copy of out before multiplying over it
void mTransv(Vector3* v, Matrix* out); // translation
void mTrans3f(float x, float y, float z, Matrix* out); // translation
void mScalev(Vector3* v, Matrix* out);
void mScale3f(float x, float y, float z, Matrix* out);
void mRotv(Vector3* v, float theta, Matrix* out); // rotate about a vector
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
void mOrthoFromSphere(Sphere* s, Vector3* eyePos, Matrix* out);

// extract the planes from an orthographic projection matrix.
void mOrthoExtractPlanes(Matrix* m, float* left, float* right, float* top, float* bottom, float* near, float* far);


// analgous to gluLookAt
// https://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
void mLookAt(Vector3* eye, Vector3* center, Vector3* up, Matrix* out);

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
void msTransv(Vector3* v, MatrixStack* ms); // translation
void msTrans3f(float x, float y, float z, MatrixStack* ms); // translation
void msScalev(Vector3* v, MatrixStack* ms);
void msScale3f(float x, float y, float z, MatrixStack* ms);
void msRotv(Vector3* v, float theta, MatrixStack* ms); // rotate about a vector
void msRot3f(float x, float y, float z, float theta, MatrixStack* ms); // rotate about a vector
void msFrustum(float left, float right, float top, float bottom, float near, float far, MatrixStack* ms);
void msPerspective(double fov, float aspect, float near, float far, MatrixStack* ms);
void msOrtho(float left, float right, float top, float bottom, float near, float far, MatrixStack* ms);
void msLookAt(Vector3* eye, Vector3* center, Vector3* up, MatrixStack* ms);


// cubic Bezier curves
void evalBezier3p(Vector3* e1, Vector3* e2, Vector3* c1, Vector3* c2, float t, Vector3* out);
void evalBezierTangent3p(Vector3* e1, Vector3* e2, Vector3* c1, Vector3* c2, float t, Vector3* out); // tangent vector; not normalized
void evalBezierNorm3p(Vector3* e1, Vector3* e2, Vector3* c1, Vector3* c2, float t, Vector3* out); // normal vector; not normalized
float evalBezier1D(float e1, float e2, float c1, float c2, float t);
float evalBezier1D_dt(float e1, float e2, float c1, float c2, float t); // first derivative with respect to t
float evalBezier1D_ddt(float e1, float e2, float c1, float c2, float t); // second derivative with respect to t

// quadratic Bezier curves
float evalQBezier1D(float e1, float e2, float c1, float t);
void evalQBezier2D3p(Vector2* e1, Vector2* e2, Vector2* c1, float t, Vector2* out);
void evalQBezier3p(Vector3* e1, Vector3* e2, Vector3* c1, float t, Vector3* out);

///// bounding box functions

// 3D versions
int boxDisjoint3p(const AABB3* a, const AABB3* b);
int boxOverlaps3p(const AABB3* a, const AABB3* b);
int boxContainsPoint3p(const AABB3* b, const Vector3* p);

Vector3 boxCenter3(const AABB3 b); // calculates the center of the box
void boxCenter3p(const AABB3* b, Vector3* out); // calculates the center of the box
Vector2 boxSize2(const AABB2 b); // calculates the size of the box
Vector3 boxSize3(const AABB3 b); // calculates the size of the box
void boxSize2p(const AABB2* b, Vector2* out); // calculates the size of the box
void boxSize3p(const AABB3* b, Vector3* out); // calculates the size of the box
void boxExpandTo3p(AABB3* b, Vector3* p);
void boxExpandTo3(AABB3* b, Vector3 p);


void makeRay3p(Vector3* origin, Vector3* direction, Ray3* out);
int boxRayIntersectFast3p(const AABB3* b, const Ray3* r);
int boxRayIntersect3p(const AABB3* b, const Ray3* r, Vector3* ipoint, float* idist);


// 2D versions
int boxDisjoint2p(const AABB2* a, const AABB2* b);
int boxOverlaps2p(const AABB2* a, const AABB2* b);
int boxContainsPoint2p(const AABB2* b, const Vector2* p);

void boxCenter2p(const AABB2* b, Vector2* out); // calcuates the center of the box
void boxSize2p(const AABB2* b, Vector2* out); // calculates the size of the box
void boxQuadrant2p(const AABB2* in, char ix, char iy, AABB2* out);

// 2D integer versions
int boxDisjoint2ip(const AABB2i* a, const AABB2i* b);
int boxOverlaps2ip(const AABB2i* a, const AABB2i* b);
int boxContainsPoint2ip(const AABB2i* b, const Vector2i* p);

void boxCenter2ip(const AABB2i* b, Vector2* out); // calcuates the center of the box
void boxSize2ip(const AABB2i* b, Vector2* out); // calculates the size of the box
void boxQuadrant2ip(const AABB2i* in, char ix, char iy, AABB2i* out);

// find the center of a quad
void quadCenter2p(const Quad2* in, Vector2* out);
void quadRoundOutward2p(const Quad2* in, Quad2i* out);
void quadRoundInward2p(const Quad2* in, Quad2i* out);


float   evalCubicHermite1D(float t, float p0, float p1, float m0, float m1);
Vector2 evalCubicHermite2D(float t, Vector2 p0, Vector2 p1, Vector2 m0, Vector2 m1); 
Vector3  evalCubicHermite3D(float t, Vector3 p0, Vector3 p1, Vector3 m0, Vector3 m1);

#endif // __c3dlas_h__

