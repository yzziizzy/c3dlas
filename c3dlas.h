#ifndef __c3dlas_h__
#define __c3dlas_h__


#include <stdlib.h> // rand() et al.
#include <stdint.h> 
#include <math.h> // fmin/fmax
#include <tgmath.h> // fmin/fmax

#undef I // because of some bullshit in <complex.h>

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

#ifndef FLT_CMP_EPSILON
	#define FLT_CMP_EPSILON 0.000001
#endif

#ifndef DBL_CMP_EPSILON
	#define DBL_CMP_EPSILON 0.00000000000001
#endif

#define FLT_CMP_EPSILON_SQ (FLT_CMP_EPSILON * FLT_CMP_EPSILON)
#define DBL_CMP_EPSILON_SQ (DBL_CMP_EPSILON_SQ * DBL_CMP_EPSILON_SQ)

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

#ifndef MAX
	#define MAX(a,b) ({ \
		__typeof__ (a) _a = (a); \
		__typeof__ (b) _b = (b); \
		_a > _b ? _a : _b; \
	})
#endif

#ifndef MIN
	#define MIN(a,b) ({ \
		__typeof__ (a) _a = (a); \
		__typeof__ (b) _b = (b); \
		_a < _b ? _a : _b; \
	})
#endif


// suffix, type, float
//               type
#define C3DLAS_VECTOR_TYPE_LIST(X, ...) \
	X(,  float,  float,  __VA_ARGS__) \
	X(d, double, double, __VA_ARGS__) \
	X(i, int,    double, __VA_ARGS__) \
	X(l, long,   double, __VA_ARGS__) \


// suffix, type, float
//               type
#define C3DLAS_VECTOR_LIST(X, ...) \
	X(2, 2,  float,  float,  2,  FLT __VA_OPT__(,) __VA_ARGS__) \
	X(3, 3,  float,  float,  3,  FLT __VA_OPT__(,) __VA_ARGS__) \
	X(4, 4,  float,  float,  4,  FLT __VA_OPT__(,) __VA_ARGS__) \
	X(2, 2d, double, double, 2d, DBL __VA_OPT__(,) __VA_ARGS__) \
	X(3, 3d, double, double, 3d, DBL __VA_OPT__(,) __VA_ARGS__) \
	X(4, 4d, double, double, 4d, DBL __VA_OPT__(,) __VA_ARGS__) \
	X(2, 2i, int,    double, 2d, DBL __VA_OPT__(,) __VA_ARGS__) \
	X(3, 3i, int,    double, 3d, DBL __VA_OPT__(,) __VA_ARGS__) \
	X(4, 4i, int,    double, 4d, DBL __VA_OPT__(,) __VA_ARGS__) \
	X(2, 2l, long,   double, 2d, DBL __VA_OPT__(,) __VA_ARGS__) \
	X(3, 3l, long,   double, 3d, DBL __VA_OPT__(,) __VA_ARGS__) \
	X(4, 4l, long,   double, 4d, DBL __VA_OPT__(,) __VA_ARGS__) \


#define C3DLAS_GEN_HELPER(sz, suf, ty, ft, sufft, pref, name, ...) Vector##suf: name##suf,




#define X(suf, t, ...) \
	typedef struct Vector2 ## suf { \
		/*    cartesian polar */ \
		union { t x, rho; }; \
		union { t y, theta; }; \
	} Vector2 ## suf;
	
	C3DLAS_VECTOR_TYPE_LIST(X)
#undef X

// The spherical coordinates use the "mathematical" conventions as opposed to the "physics" conventions
// This is because of the struct layout convenience for overlapping with 2D polar coordinates. 

#define X(suf, t, ...) \
	typedef struct Vector3 ## suf { \
		/*    cartesian  spherical  color */ \
		union { t x, rho,       r; }; /* rho is the radius */ \
		union { t y, theta,     g; }; /* rotation in the X-Y plane */ \
		union { t z, phi,       b; }; /* rotation in the plane passing through the Z axis */ \
	} Vector3 ## suf;
	
	C3DLAS_VECTOR_TYPE_LIST(X)
#undef X

#define X(suf, t, ...) \
	typedef struct Vector4 ## suf { \
		/*    cartesian  spherical  color  quaternion basis */ \
		union { t x, rho,       r,     i; }; /* rho = the radius */ \
		union { t y, theta,     g,     j; }; /* theta = rotation in the X-Y plane */ \
		union { t z, phi,       b,     k; }; /* phi = rotation in the plane passing through the Z axis */ \
		union { t w,            a,     real; }; /* real is last for memory alignment with 4d cartesian vectors */ \
	} Vector4 ## suf;
	
	C3DLAS_VECTOR_TYPE_LIST(X)
#undef X



// Best quaternion site on the internet so far: http://www.songho.ca/math/quaternion/quaternion.html
typedef struct Vector4 Quaternion;
typedef struct Vector4d Quaterniond;

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
	union { Vector3 start, a; };
	union { Vector3 end, b; };
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
	float d; // distance along normal from the origin
} Plane;

typedef struct {
	Vector3 n; // normal
	Vector3 p; // an arbitrary point on the plane
} PlaneP;

typedef struct {
	Vector3 center;
	float r;
} Sphere;

typedef struct {
	Vector3 a, b; // the two centroids
	float r; // radius
} Capsule;

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


/* Column-major, for OpenGL compatibility:
 _            _
|  0  4  8 12  |
|  1  5  9 13  |
|  2  6 10 14  |
|_ 3  7 11 15 _|

*/

typedef struct {
	float m[16];
} Matrix;

typedef struct MatrixStack {
	short size;
	short top;
	Matrix* stack;
} MatrixStack;


// Symmetric matrix
/* Column-major, for OpenGL compatibility:
 _            _
|  0  1  3  6  |
|     2  4  7  |
|        5  8  |
|_          9 _|

       ==
 _            _
|  0  1  3  6  |
|  1  2  4  7  |
|  3  4  5  8  |
|_ 6  7  8  9 _|

*/

typedef struct {
	float m[10];
} MatrixSym;


// 3x3 matrix
/* Column-major, for OpenGL compatibility:
 _         _
|  0  3  6  |
|  1  4  7  |
|_ 2  5  8 _|

*/
typedef struct {
	float m[9];
} Matrix3;

// 2x2 matrix
/* Column-major, for OpenGL compatibility:
 _      _
|  0  2  |
|_ 1  3 _|

*/
typedef struct {
	float m[4];
} Matrix2;


// axis-aligned bounding box
#define X(sz, suf, t, ...) \
	typedef struct AABB ## suf { \
		Vector ## suf min, max; \
	} AABB ## suf;
	
	C3DLAS_VECTOR_LIST(X)
#undef X


typedef struct {
	uint64_t state, stream;
} PCG;



extern const Matrix IDENT_MATRIX;
extern const Matrix3 IDENT_MATRIX3;

// utilities

uint32_t bitReverse32(uint32_t x);
uint32_t reverseBits(uint32_t n, int len);

float pcg_f(uint64_t* state, uint64_t stream);
float frandPCG(float low, float high, PCG* pcg);

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
	return fminf(max, fmaxf(min, val));
}

static inline double dclamp(double val, double min, double max) {
	return fmin(max, fmax(min, val));
}

static inline float fclampNorm(float val) {
	return fclamp(val, 0.0f, 1.0f);
}

static inline double dclampNorm(double val) {
	return dclamp(val, 0.0, 1.0);
}

static inline int iclamp(int val, int min, int max) {
	return MIN(max, MAX(min, val));
}

static inline long lclamp(long val, long min, long max) {
	return MIN(max, MAX(min, val));
}

static inline float flerp(float a, float b, float t) {
	return a  + ((b - a) * t);
}
static inline double dlerp(double a, double b, double t) {
	return a  + ((b - a) * t);
}

static inline float flerp2D(float xx, float xy, float yx, float yy, float xt, float yt) {
	float a = xx  + ((xy - xx) * yt);
	float b = yx  + ((yy - yx) * yt);
	return a  + ((b - a) * xt);
}

static inline double dlerp2D(double xx, double xy, double yx, double yy, double xt, double yt) {
	double a = xx  + ((xy - xx) * yt);
	double b = yx  + ((yy - yx) * yt);
	return a  + ((b - a) * xt);
}


// Returns an arbitrary unit vector perpendicular to the input
// The input vector does not need to be normalized
void vPerp2p(Vector2* n, Vector2* out);
Vector2 vPerp2(Vector2 n);
void vPerp3p(Vector3* n, Vector3* out);
Vector3 vPerp3(Vector3 n);

//
// Vectors
//

#define X(sz, suf, ty, ft, sufft, ...) \
	int vEq##suf(const Vector##suf a, const Vector##suf b); \
	int vEq##suf##p(const Vector##suf* a, const Vector##suf* b); \
	\
	int vEqEp##suf(const Vector##suf a, const Vector##suf b, ft epsilon); \
	int vEqEp##suf##p(const Vector##suf* a, const Vector##suf* b, ft epsilon); \
	\
	int vEqExact##suf(const Vector##suf a, const Vector##suf b); \
	int vEqExact##suf##p(const Vector##suf* a, const Vector##suf* b); \
	\
	Vector##suf vAdd##suf(const Vector##suf a, const Vector##suf b); \
	void vAdd##suf##p(const Vector##suf* a, const Vector##suf* b, Vector##suf* out); \
	\
	Vector##suf vSub##suf(const Vector##suf from, const Vector##suf what); \
	void vSub##suf##p(const Vector##suf* from, const Vector##suf* what, Vector##suf* out); \
	\
	Vector##suf vMul##suf(const Vector##suf a, const Vector##suf b); \
	void vMul##suf##p(const Vector##suf* a, const Vector##suf* b, Vector##suf* out); \
	\
	Vector##sufft vScale##suf(const Vector##suf v, ft scalar); \
	void vScale##suf##p(const Vector##suf* v, ft scalar, Vector##sufft* out); \
	\
	Vector##suf vMin##suf(const Vector##suf a, const Vector##suf b); \
	void vMin##suf##p(const Vector##suf* a, const Vector##suf* b, Vector##suf* out); \
	\
	Vector##suf vMax##suf(const Vector##suf a, const Vector##suf b); \
	void vMax##suf##p(const Vector##suf* a, const Vector##suf* b, Vector##suf* out); \
	\
	Vector##suf vClamp##suf(const Vector##suf in, const Vector##suf min, const Vector##suf max); \
	void vClamp##suf##p(const Vector##suf* in, const Vector##suf* min, const Vector##suf* max, Vector##suf* out); \
	\
	int vMinComp##suf(const Vector##suf a); \
	int vMinComp##suf##p(const Vector##suf* a); \
	\
	int vMaxComp##suf(const Vector##suf a); \
	int vMaxComp##suf##p(const Vector##suf* a); \
	\
	ft vDot##suf(const Vector##suf a, const Vector##suf b); \
	ft vDot##suf##p(const Vector##suf* a, const Vector##suf* b); \
	\
	Vector##sufft vAvg##suf(const Vector##suf a, const Vector##suf b); \
	void vAvg##suf##p(const Vector##suf* a, const Vector##suf* b, Vector##sufft* out); \
	\
	ft vDist##suf(const Vector##suf a, const Vector##suf b); \
	ft vDist##suf##p(const Vector##suf* a, const Vector##suf* b); \
	\
	ft vDistSq##suf(const Vector##suf a, const Vector##suf b); \
	ft vDistSq##suf##p(const Vector##suf* a, const Vector##suf* b); \
	\
	Vector##sufft vNorm##suf(const Vector##suf v); \
	void vNorm##suf##p(const Vector##suf* v, Vector##sufft* out); \
	\
	Vector##sufft vUnit##suf(const Vector##suf v); \
	void vUnit##suf##p(const Vector##suf* v, Vector##sufft* out); \
	\
	ft vLen##suf(const Vector##suf v); \
	ft vLen##suf##p(const Vector##suf* v); \
	ft vMag##suf(const Vector##suf v); \
	ft vMag##suf##p(const Vector##suf* v); \
	\
	ft vLenSq##suf(const Vector##suf v); \
	ft vLenSq##suf##p(const Vector##suf* v); \
	\
	ft vInvLen##suf(const Vector##suf v); \
	ft vInvLen##suf##p(const Vector##suf* v); \
	\
	Vector##suf vAbs##suf(const Vector##suf v); \
	void vAbs##suf##p(const Vector##suf* v, Vector##suf* out); \
	\
	Vector##sufft vRecip##suf(const Vector##suf v); \
	void vRecip##suf##p(const Vector##suf* v, Vector##sufft* out); \
	Vector##sufft vInv##suf(const Vector##suf v); \
	void vInv##suf##p(const Vector##suf* v, Vector##sufft* out); \
	\
	Vector##suf vNeg##suf(const Vector##suf v); \
	void vNeg##suf##p(const Vector##suf* v, Vector##suf* out); \
	\
	Vector##suf vSign##suf(const Vector##suf v); \
	void vSign##suf##p(const Vector##suf* v, Vector##suf* out); \
	\
	Vector##suf vStep##suf(const Vector##suf edge, const Vector##suf v); \
	void vStep##suf##p(const Vector##suf* edge, const Vector##suf* v, Vector##suf* out); \
	\
	Vector##sufft vLerp##suf(const Vector##suf a, const Vector##suf b, ft t); \
	void vLerp##suf##p(const Vector##suf* a, const Vector##suf* b, ft t, Vector##sufft* out); \
	\
	
		
	C3DLAS_VECTOR_LIST(X)

#undef X

	

#define vAdd(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vAdd) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vSub(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vSub) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMul(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMul) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vScale(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vScale) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vLerp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vLerp) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vRecip(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vRecip) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vInv(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vInv) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMag(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMag) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vInvLen(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vInvLen) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vLen(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vLen) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vLenSq(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vLenSq) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vDot(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vDot) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vDist(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vDist) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vDistSq(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vDistSq) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vNorm(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vNorm) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vNeg(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vNeg) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vSign(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vSign) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vStep(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vStep) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMin(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMin) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMax(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMax) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMinComp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMinComp) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMaxComp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMaxComp) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vAbs(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vAbs) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vAvg(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vAvg) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vClamp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vClamp) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vEq(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vEq) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vEqEp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vEqEp) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vEqExact(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vEqExact) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)





// Swap two vectors
void vSwap2ip(Vector2i* a, Vector2i* b);
void vSwap2p(Vector2* a, Vector2* b);
void vSwap3p(Vector3* a, Vector3* b);
void vSwap4p(Vector4* a, Vector4* b);

// Vector addition


// Vector subtraction. diff = from - what

// Scalar muliplication
//Vector2  vScale2(Vector2 v, float scalar);
//Vector3  vScale3(Vector3 v, float scalar);
//void     vScale2ip(Vector2i* v, float scalar, Vector2i* out);
//void     vScale2p(Vector2* v, float scalar, Vector2* out);
//void     vScale3p(Vector3* v, float scalar, Vector3* out);

// Component-wise vector muliplication


// Dot product (inner product)

// Cross product: out = a x b
Vector3 vCross3(Vector3 a, Vector3 b);
Vector3 vCross3(Vector3 a, Vector3 b);
void  vCross3p(Vector3* a, Vector3* b, Vector3* out); 

// Scalar triple product: a . (b x c)
float vScalarTriple3(Vector3 a, Vector3 b, Vector3 c);
float vScalarTriple3p(Vector3* a, Vector3* b, Vector3* c);

// Linear interpolation between two vectors

// Vector Inverse. Returns FLT/DBL_MAX on div/0. Integer functions return double vectors

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




Vector2 vClamp2f(Vector2 in, float min, float max);
Vector3 vClamp3f(Vector3 in, float min, float max);

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

float distLineLine3(Line3* a, Line3* b);

float distTPointRay3(Vector3 p, Ray3 r, float* T);
float dist2TPointRay3(Vector3 p, Ray3 r, float* T);

int vInsidePolygon(Vector2 p, Polygon* poly);

// Returns the distance from p to the closest point on the polygon.
// Interior distances are negative
float vDistPolygon(Vector2 p, Polygon* poly);
void polyCalcCentroid(Polygon* poly);
void polyCalcRadiusSq(Polygon* poly) ;

void  vProject3p(Vector3* what, Vector3* onto, Vector3* out); // slower; onto may not be normalized
void  vProjectNorm3p(Vector3* what, Vector3* onto, Vector3* out); // faster; onto must be normalized
void  vPointAvg3p(Vector3* a, Vector3* b, Vector3* out);



void vRandomPCG3p(Vector3* end1, Vector3* end2, PCG* pcg, Vector3* out);
Vector3 vRandomPCG3(Vector3 end1, Vector3 end2, PCG* pcg);

void vRandom3p(Vector3* end1, Vector3* end2, Vector3* out);
Vector3 vRandom3(Vector3 end1, Vector3 end2);
void vRandomNorm3p(Vector3* out);



// http://geomalgorithms.com/a07-_distance.html
// _PARALLEL with no output on parallel lines
// _INTERSECT with one point of output on intersection
// _DISJOINT with two outputs otherwise
int shortestLineFromRayToRay3p(Ray3* r1, Ray3* r2, Vector3* pOut);

// reflects the distance from v to pivot across pivot.
// out, pivot, and v will form a straight line with pivot exactly in the middle.
void    vReflectAcross3p(Vector3* v, Vector3* pivot, Vector3* out);
Vector3 vReflectAcross3(Vector3 v, Vector3 pivot);
void  vTriFaceNormal3p(Vector3* a, Vector3* b, Vector3* c, Vector3* out); // returns a normalized face normal for the given triangle
Vector3 vTriFaceNormal3(Vector3 a, Vector3 b, Vector3 c);
void  vpTriFaceNormal3p(Vector3* tri, Vector3* out); // returns a normalized face normal for the given triangle

void  vProjectOntoPlane3p(Vector3* v, Plane* p, Vector3* out);
void  vProjectOntoPlaneNormalized3p(Vector3* v, Plane* p, Vector3* out);

void  planeFromPointNormal(Vector3* p, Vector3* norm, Plane* out);
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

// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
// returns _INTERSECT or _DISJOINT
int rayTriangleIntersect(
	Vector3* a, Vector3* b, Vector3* c, // triangle
	Vector3* ray_origin, Vector3* ray_dir, // ray
	float* u, float* v, float* t // barycentric out coords, t of intersection point along ray 
);

Vector3 triangleClosestPoint(
	Vector3* a, Vector3* b, Vector3* c, // triangle
	Vector3* p, // test point
	float* out_u, float* out_v // barycentric out coords of closest point 
);

Vector3 triangleClosestPoint_Reference(
	Vector3* a, Vector3* b, Vector3* c, // triangle
	Vector3* p, // test point
	float* out_u, float* out_v // barycentric out coords of closest point 
);



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



void mIdent3(Matrix3* m);

// out cannot overlap with a or b
// with restrict and -O2, this vectorizes nicely.
void mFastMul3(Matrix3* restrict a, Matrix3* restrict b, Matrix3* restrict out);
void mMul3(Matrix3* a, Matrix3* out);
float mDeterminate3(Matrix3* m);
void mInverse3(Matrix3* m, Matrix3* out);
void mScalarMul3(Matrix3* m, float scalar, Matrix3* out);
Vector3 vMatrix3Mul(Vector3 v, Matrix3* restrict m);
void mTranspose3(Matrix3* m, Matrix3* out);
float mTrace3(Matrix3* m); // sum of the diagonal elements



float pvDist3p(Plane* p, Vector3* v);

void vMatrixMul3p(Vector3* in, Matrix* m, Vector3* out); // multiply a vector by a matrix
void vMatrixMulf3p(float x, float y, float z, Matrix* m, Vector3* out); // multiply a vector by a matrix
Vector3 vMatrixMul3(Vector3 in, Matrix* m);
Vector4 vMatrixMul4(Vector4 in, Matrix* m);
Vector3 vMatrixMulProjectedMagic3(Vector3 in, Matrix* m);

// These are 3d spatial operations
void mIdent(Matrix* m); // set m to the identity matrix
void mCopy(Matrix* in, Matrix* out);
void mFastMul(Matrix* a, Matrix* b, Matrix* out); // a and b cannot also be out. mostly internal use.
void mMul(Matrix* a, Matrix* out); // makes a copy of out before multiplying over it
void mTransv(Vector3* v, Matrix* out); // translation
void mTrans3f(float x, float y, float z, Matrix* out); // translation
void mScalev(Vector3* v, Matrix* out);
void mScale3f(float x, float y, float z, Matrix* out);
void mRotv(Vector3* v, float theta, Matrix* out); // rotate about a vector
//void mRotq(Vector3* v, float theta, Matrix* out); // rotate by a Quaternion
void mRot3f(float x, float y, float z, float theta, Matrix* out); // rotate about a vector
void mRotX(float theta, Matrix* out); //
void mRotY(float theta, Matrix* out); // rotate about axes
void mRotZ(float theta, Matrix* out); //
void mTranspose(Matrix* in, Matrix* out);
void mTransposeFast(Matrix* in, Matrix* out); // in cannot be out
float mDeterminate(Matrix* m);
int mInverse(Matrix* in, Matrix* out); // returns 0 on success, 1 if there is no inverse; out remains unchanged
float mTrace(Matrix* m); // sum of the diagonal elements

// simple component-wise mathematical operations
void mAdd(Matrix* a, Matrix* b, Matrix* out); 
void mScalarMul(Matrix* a, float f, Matrix* out);


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
void mPerspSetNF_ZUp(Matrix* m, float near, float far);

// orthographic projection. use this for a "2D" look.
// same div/0 warnings.
void mOrtho(float left, float right, float top, float bottom, float near, float far, Matrix* out);

// calculates a cubical orthographic matrix with a side length of 2*r
void mOrthoFromRadius(float r, Matrix* out);

// calculates an orthographic matrix that encloses the sphere, looking from eyePos
void mOrthoFromSphere(Sphere s, Vector3 eyePos, Vector3 up, Matrix* out);

// extract the planes from an orthographic projection matrix.
void mOrthoExtractPlanes(Matrix* m, float* left, float* right, float* top, float* bottom, float* near, float* far);

void mOrthoSetNF(Matrix* m, float near, float far);


// analgous to gluLookAt
// up is not required to be orthogonal to anything, so long as it's not parallel to anything
// http://www.songho.ca/opengl/gl_camera.html#lookat
void mLookAt(Vector3 eye, Vector3 center, Vector3 up, Matrix* out);


void mLookDir(Vector3 eye, Vector3 center, Vector3 dir, Vector3 up, Matrix* out);


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
int intersectBoxLine3p(const AABB3* b, const Line3* l, Vector3* ipoint, float* idist);
int intersectBoxLine3(AABB3 b, Line3 l, Vector3* ipoint, float* idist);

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


Quaternion qFromRTheta(Vector3 r, float theta);

Quaternion qAdd(Quaternion l, Quaternion r);
Quaternion qSub(Quaternion l, Quaternion r);
Quaternion qScale(Quaternion q, float s);
Quaternion qMul(Quaternion l, Quaternion r);
Quaternion qDiv(Quaternion n, Quaternion d);
Quaternion qRot(Quaternion a, Quaternion r);
Vector3 qRot3(Vector3 a, Quaternion r);
Quaternion qConj(Quaternion q);
Quaternion qInv(Quaternion q);
Quaternion qNorm(Quaternion q);

Quaternion qSlerp(Quaternion a, Quaternion b, float t);
Quaternion qNlerp(Quaternion a, Quaternion b, float t);

// these appear to all mean the same thing for quaternions.
float qMod(Quaternion q);
float qMag(Quaternion q);
float qLen(Quaternion q);


Quaternion qFromBasis(Vector3 bx, Vector3 by, Vector3 bz);

// Applies the full conjugate multiplication qvq*
void qNonUnitToMatrix(Quaternion q, Matrix* out);

// faster
void qUnitToMatrix3(Quaternion q, Matrix3* out);
void qUnitToMatrix(Quaternion q, Matrix* out);

#endif // __c3dlas_h__

