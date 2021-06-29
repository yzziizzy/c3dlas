

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <x86intrin.h>

#include "c3dlas.h"




// utilities

// reverses the argument
uint32_t bitReverse32(uint32_t x) {
	x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
	x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
	x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
	x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
	return ((x >> 16) | (x << 16));
}


// reverses the least significant (len) bits, zeroing the top
uint32_t reverseBits(uint32_t n, int len) {
	uint32_t rn = bitReverse32(n);
	
	return rn >> (32 - len);
}


// random numbers

// returns a random number in (-1, 1) uninclusive
float pcg_f(uint64_t* state, uint64_t stream) {
	union {
		uint32_t fu;
		float ff;
	} u;
	
	
	uint64_t last = *state;
	*state = (last * 6364136223846793005ULL) + (stream | 1);
	uint32_t xs = ((last >> 18) ^ last) >> 27;
	uint32_t rot = last >> 59;
	u.fu = (((xs >> rot) | (xs << ((-rot) & 31))) & 0x807fffff) | 0x3f000000;
	return u.ff;
}


// BUG: totally untested
// SIMD and C versions do not return the same values. 
void pcg_f8(uint64_t* state, uint64_t stream, float* out) {
	
#if defined(C3DLAS_USE_SIMD)
	__m256i s1, s2, xs1, xs2, xs, r, nra, q, f;
	
	s1 = _mm256_add_epi64(_mm256_set1_epi64x(*state), _mm256_set_epi64x(1,2,3,4));
	s2 = _mm256_add_epi64(_mm256_set1_epi64x(*state), _mm256_set_epi64x(5,6,7,8));
	
	// cycle the state
	*state = (*state * 6364136223846793005ULL) + (stream | 1);
	
	xs1 = _mm256_srli_epi64(_mm256_xor_si256(_mm256_srli_epi64(s1, 18), s1), 27);
	xs2 = _mm256_srli_epi64(_mm256_xor_si256(_mm256_srli_epi64(s2, 18), s2), 27);
	
	xs = _mm256_unpacklo_epi32(xs1, xs2);
	
	r = _mm256_srai_epi32(xs, 59);
	nra = _mm256_and_si256(_mm256_sign_epi32(r, _mm256_set1_epi32(-1)), _mm256_set1_epi32(31));
	
	q = _mm256_or_si256(_mm256_srav_epi32(xs, r), _mm256_sllv_epi32(xs, nra)); 
	
	// q is full of random 32bit integers now
	
	// convert to (-1, 1) floats by jamming in some exponent info
	f = _mm256_or_si256(_mm256_and_si256(q, _mm256_set1_epi32(0x807fffff)), _mm256_set1_epi32(0x3f000000));
	 
	_mm256_storeu_si256((__m256i*)out, f); 
	
#else
	out[0] = pcg_f(state, stream);
	out[1] = pcg_f(state, stream);
	out[2] = pcg_f(state, stream);
	out[3] = pcg_f(state, stream);
	out[4] = pcg_f(state, stream);
	out[5] = pcg_f(state, stream);
	out[6] = pcg_f(state, stream);
	out[7] = pcg_f(state, stream);
#endif
}




int vEq2i(Vector2i a, Vector2i b) { return vEq2ip(&a, &b); }
int vEqExact2i(Vector2i a, Vector2i b) { return vEqExact2ip(&a, &b); } 
int vEqExact2(Vector2 a, Vector2 b) { return vEqExact2p(&a, &b); } 
int vEqExact3(Vector3 a, Vector3 b) { return vEqExact3p(&a, &b); } 
int vEqExact4(Vector4 a, Vector4 b) { return vEqExact4p(&a, &b); } 

int vEq2ip(Vector2i* a, Vector2i* b) { return vEqExact2ip(a, b); } 
int vEqExact2ip(Vector2i* a, Vector2i* b)  {
	return a->x == b->x && a->y == b->y;
}
 
int vEqExact2p(Vector2* a, Vector2* b) {
	return a->x == b->x && a->y == b->y;
}
 
int vEqExact3p(Vector3* a, Vector3* b) {
	return a->x == b->x && a->y == b->y && a->z == b->z;
}
 
int vEqExact4p(Vector4* a, Vector4* b) {
	return a->x == b->x && a->y == b->y && a->z == b->z && a->w == b->w;
}

int vEq2(Vector2 a, Vector2 b) { return vEq2p(&a, &b); } 
int vEq3(Vector3 a, Vector3 b) { return vEq3p(&a, &b); } 
int vEq4(Vector4 a, Vector4 b) { return vEq4p(&a, &b); }

int vEq2p(Vector2* a, Vector2* b) {
	return vEqEp2p(a, b, FLT_CMP_EPSILON);
}

int vEq3p(Vector3* a, Vector3* b) {
	return vEqEp3p(a, b, FLT_CMP_EPSILON);
}

int vEq4p(Vector4* a, Vector4* b) {
	return vEqEp4p(a, b, FLT_CMP_EPSILON);
}

int vEqEp2i(Vector2i a, Vector2i b, double epsilon) { return vEqEp2ip(&a, &b, epsilon); }
int vEqEp2(Vector2 a, Vector2 b, float epsilon) { return vEqEp2p(&a, &b, epsilon); }
int vEqEp3(Vector3 a, Vector3 b, float epsilon) { return vEqEp3p(&a, &b, epsilon); }
int vEqEp4(Vector4 a, Vector4 b, float epsilon) { return vEqEp4p(&a, &b, epsilon); }

int vEqEp2ip(Vector2i* a, Vector2i* b, double epsilon) {
	double x, y, n;
	
	x = (double)a->x - (double)b->x;
	y = (double)a->y - (double)b->y;
	
	n = fabs(x * x + y * y);
	
	return n <= epsilon * epsilon;
}

int vEqEp2p(Vector2* a, Vector2* b, float epsilon) {
	float x, y, n;
	
	x = a->x - b->x;
	y = a->y - b->y;
	
	n = fabsf(x * x + y * y);
	
	return n <= epsilon * epsilon;
}

int vEqEp3p(Vector3* a, Vector3* b, float epsilon) {
	float x, y, z, n;
	
	x = a->x - b->x;
	y = a->y - b->y;
	z = a->z - b->z;
	
	n = fabsf(x * x + y * y + z * z);
	
	return n <= epsilon * epsilon;
}

int vEqEp4p(Vector4* a, Vector4* b, float epsilon) {
	float x, y, z, w, n;
	
	x = a->x - b->x;
	y = a->y - b->y;
	z = a->z - b->z;
	w = a->w - b->w;
	
	n = fabsf(x * x + y * y + z * z + w * w);
	
	return n <= epsilon * epsilon;
}


void vCopy3p(const Vector3* src, Vector3* dst) {
	*dst = *src;
}




// swap two vectors

void vSwap2ip(Vector2i* a, Vector2i* b) {
	Vector2i t;
	t = *b;
	*b = *a;
	*a = t;
}

void vSwap2p(Vector2* a, Vector2* b) {
	Vector2 t;
	t = *b;
	*b = *a;
	*a = t;
}

void vSwap3p(Vector3* a, Vector3* b) {
	Vector3 t;
	t = *b;
	*b = *a;
	*a = t;
}

void vSwap4p(Vector4* a, Vector4* b) {
	Vector4 t;
	t = *b;
	*b = *a;
	*a = t;
}




// Vector addition

Vector2i vAdd2i(Vector2i a, Vector2i b) {
	Vector2i out;
	vSub2ip(&a, &b, &out);
	return out;
}

void vAdd2ip(Vector2i* a, Vector2i* b, Vector2i* sum) {
	sum->x = a->x + b->x;
	sum->y = a->y + b->y;
}

Vector2 vAdd2(Vector2 a, Vector2 b) {
	Vector2 out;
	vAdd2p(&a, &b, &out);
	return out;
}

void vAdd2p(Vector2* a, Vector2* b, Vector2* sum) {
	sum->x = a->x + b->x;
	sum->y = a->y + b->y;
}

Vector3 vAdd3(Vector3 a, Vector3 b) {
	Vector3 out;
	vAdd3p(&a, &b, &out);
	return out;
}

void vAdd3p(Vector3* a, Vector3* b, Vector3* sum) {
#ifdef C3DLAS_USE_SIMD
	__m128 f_ = _mm_loadu_ps((float*)a);
	__m128 w_ = _mm_loadu_ps((float*)b);
	w_ = _mm_add_ps(f_, w_);
	
	_mm_maskstore_ps((float*)sum, _mm_set_epi32(0, -1, -1, -1), w_);
#else
	sum->x = a->x + b->x;
	sum->y = a->y + b->y;
	sum->z = a->z + b->z;
#endif
}

Vector4 vAdd4(Vector4 a, Vector4 b) {
	Vector4 out;
	vAdd4p(&a, &b, &out);
	return out;
}

void vAdd4p(Vector4* a, Vector4* b, Vector4* sum) {
	sum->x = a->x + b->x;
	sum->y = a->y + b->y;
	sum->z = a->z + b->z;
	sum->w = a->w + b->w;
}




// Vector subtraction. diff = from - what

Vector2i vSub2i(Vector2i from, Vector2i what) {
	Vector2i out;
	vSub2ip(&from, &what, &out);
	return out;
}

void vSub2ip(Vector2i* from, Vector2i* what, Vector2i* diff) {
	diff->x = from->x - what->x;
	diff->y = from->y - what->y;
}

Vector2 vSub2(Vector2 from, Vector2 what) {
	Vector2 out;
	vSub2p(&from, &what, &out);
	return out;
}

void vSub2p(Vector2* from, Vector2* what, Vector2* diff) {
	diff->x = from->x - what->x;
	diff->y = from->y - what->y;
}

Vector3 vSub3(Vector3 from, Vector3 what) {
	Vector3 out;
	vSub3p(&from, &what, &out);
	return out;
}

void vSub3p(Vector3* from, Vector3* what, Vector3* diff) {
#ifdef C3DLAS_USE_SIMD
	__m128 f_ = _mm_loadu_ps((float*)from);
	__m128 w_ = _mm_loadu_ps((float*)what);
	w_ = _mm_sub_ps(f_, w_);
	
	_mm_maskstore_ps((float*)diff, _mm_set_epi32(0, -1, -1, -1), w_);
#else
	diff->x = from->x - what->x;
	diff->y = from->y - what->y;
	diff->z = from->z - what->z;
#endif
}

Vector4 vSub4(Vector4 from, Vector4 what) {
	Vector4 out;
	vSub4p(&from, &what, &out);
	return out;
}

void vSub4p(Vector4* from, Vector4* what, Vector4* diff) {
	diff->x = from->x - what->x;
	diff->y = from->y - what->y;
	diff->z = from->z - what->z;
	diff->w = from->w - what->w;
}



// scalar muliplication

Vector2i vScale2i(Vector2i v, float scalar) {
	Vector2i out;
	vScale2ip(&v, scalar, &out);
	return out;
}

void vScale2ip(Vector2i* v, float scalar, Vector2i* out) {
	out->x = (float)v->x * scalar;
	out->y = (float)v->y * scalar;
}

Vector2 vScale2(Vector2 v, float scalar) {
	Vector2 out;
	vScale2p(&v, scalar, &out);
	return out;
}

void vScale2p(Vector2* v, float scalar, Vector2* out) {
	out->x = v->x * scalar;
	out->y = v->y * scalar;
}

Vector3 vScale3(Vector3 v, float scalar) {
	Vector3 out;
	vScale3p(&v, scalar, &out);
	return out;
}

void vScale3p(Vector3* v, float scalar, Vector3* out) {
// #ifdef C3DLAS_USE_SIMD
// 	__m128 v_ = _mm_loadu_ps((float*)v);
// 	       v_ = _mm_mul_ps(v_, _mm_set_ps1(scalar));
// 	
// 	_mm_maskstore_ps((float*)out, _mm_set_epi32(0, -1, -1, -1), v_);
// #else
	out->x = v->x * scalar;
	out->y = v->y * scalar;
	out->z = v->z * scalar;
// #endif
}



// Dot product (inner product)

double vDot2i(const Vector2i a, const Vector2i b) { return vDot2ip(&a, &b); }
float vDot2(const Vector2 a, const Vector2 b) { return vDot2p(&a, &b); }
float vDot3(const Vector3 a, const Vector3 b) { return vDot3p(&a, &b); }
float vDot4(const Vector4 a, const Vector4 b) { return vDot4p(&a, &b); }

double vDot2ip(const Vector2i* a, const Vector2i* b) {
	return ((double)a->x * (double)b->x) + ((double)a->y * (double)b->y);
}
float vDot2p(const Vector2* a, const Vector2* b) {
	return (a->x * b->x) + (a->y * b->y);
}
float vDot3p(const Vector3* a, const Vector3* b) {
	return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}
float vDot4p(const Vector4* a, const Vector4* b) {
	return (a->x * b->x) + (a->y * b->y) + (a->z * b->z) + (a->w * b->w);
}




// Cross product: out = a x b
// Cross products only exist in 3 and 7 dimensions
Vector3 vCross3(Vector3 a, Vector3 b) {
	Vector3 out;
	vCross3p(&a, &b, &out);
	return out;
}

void vCross3p(Vector3* a, Vector3* b, Vector3* out) {
	out->x = (a->y * b->z) - (a->z * b->y);
	out->y = (a->z * b->x) - (a->x * b->z);
	out->z = (a->x * b->y) - (a->y * b->x);
}




// Scalar triple product: a . (b x c)
float vScalarTriple3(Vector3 a, Vector3 b, Vector3 c) {
	return vScalarTriple3p(&a, &b, &c);
}

float vScalarTriple3p(Vector3* a, Vector3* b, Vector3* c) {
	return (float)((a->x * b->y * c->z) - (a->x * b->z * c->y) - (a->y * b->x * c->z)
				 + (a->z * b->x * c->y) + (a->y * b->z * c->x) - (a->z * b->y * c->x));
}




// Linear interpolation between two vectors

Vector2 vLerp2(Vector2 a, Vector2 b, float t) {
	Vector2 out;
	vLerp2p(&a, &b, t, &out);
	return out;
}

void vLerp2p(Vector2* a, Vector2* b, float t, Vector2* out) {
	out->x = a->x + ((b->x - a->x) * t);
	out->y = a->y + ((b->y - a->y) * t);
}

Vector3 vLerp3(Vector3 a, Vector3 b, float t) {
	Vector3 out;
	vLerp3p(&a, &b, t, &out);
	return out;	
}

void vLerp3p(Vector3* a, Vector3* b, float t, Vector3* out) {
#ifdef C3DLAS_USE_SIMD
	__m128 a_ = _mm_loadu_ps((float*)a);
	__m128 b_ = _mm_loadu_ps((float*)b);
	__m128 t_ = _mm_set_ps1(t);
	
	__m128 q = _mm_sub_ps(b_, a_);
	       q = _mm_mul_ps(q, t_);
	       q = _mm_add_ps(q, a_);
	
	_mm_maskstore_ps((float*)out, _mm_set_epi32(0, -1, -1, -1), q);
#else
	out->x = a->x + ((b->x - a->x) * t);
	out->y = a->y + ((b->y - a->y) * t);
	out->z = a->z + ((b->z - a->z) * t);
#endif
}


Vector4 vLerp4(Vector4 a, Vector4 b, float t) {
	Vector4 out; // the compiler seems to properly inline this way but not the other way around
	vLerp4p(&a, &b, t, &out);
	return out;
}

void vLerp4p(Vector4* a, Vector4* b, float t, Vector4* out) {

#ifdef C3DLAS_HAS_SSE1
	// seems to be equivalent to the best compiler-generated code. 
	__m128 a_ = _mm_loadu_ps((float*)a);
	__m128 b_ = _mm_loadu_ps((float*)b);
	__m128 t_ = _mm_set_ps1(t);
	
	__m128 q = _mm_sub_ps(b_, a_);
	       q = _mm_mul_ps(q, t_);
	       q = _mm_add_ps(q, a_);
	
	_mm_storeu_ps((float*)out, q);
#else
	out->x = a->x + ((b->x - a->x) * t);
	out->y = a->y + ((b->y - a->y) * t);
	out->z = a->z + ((b->z - a->z) * t);
	out->w = a->w + ((b->w - a->w) * t);
#endif
}




// Vector Inverse. Returns FLT_MAX on div/0

Vector2 vInverse2(const Vector2 v) {
	Vector2 out;
	vInverse2p(&v, &out);
	return out;
}

void vInverse2p(const Vector2* v, Vector2* out) {
	out->x = v->x == 0.0f ? FLT_MAX : 1.0f / v->x;
	out->y = v->y == 0.0f ? FLT_MAX : 1.0f / v->y;
}

Vector3 vInverse3(const Vector3 v) {
	Vector3 out;
	vInverse3p(&v, &out);
	return out;
}

void vInverse3p(const Vector3* v, Vector3* out) {
	// yeah yeah yeah shut up. In games, pedeantic details are just annoying. 
	// This function does what you mean rather than sucking off the IEEE Standards Committee
	out->x = v->x == 0.0f ? FLT_MAX : 1.0f / v->x;
	out->y = v->y == 0.0f ? FLT_MAX : 1.0f / v->y;
	out->z = v->z == 0.0f ? FLT_MAX : 1.0f / v->z;
}

Vector4 vInverse4(const Vector4 v) {
	Vector4 out;
	vInverse4p(&v, &out);
	return out;
}

void vInverse4p(const Vector4* v, Vector4* out) {
	out->x = v->x == 0.0f ? FLT_MAX : 1.0f / v->x;
	out->y = v->y == 0.0f ? FLT_MAX : 1.0f / v->y;
	out->z = v->z == 0.0f ? FLT_MAX : 1.0f / v->z;
	out->w = v->w == 0.0f ? FLT_MAX : 1.0f / v->w;
}



// Vector magnitude (length)

double vMag2i(const Vector2i v) { return vMag2ip(&v); }
float  vMag2(const Vector2 v) { return vMag2p(&v); }
float  vMag3(const Vector3 v) { return vMag3p(&v); }
float  vMag4(const Vector4 v) { return vMag4p(&v); }

double vMag2ip(const Vector2i* v) {
	return sqrt((double)v->x * (double)v->x + (double)v->y * (double)v->y);
}

float vMag2p(const Vector2* v) {
	return sqrtf((v->x * v->x) + (v->y * v->y));
}

float vMag3p(const Vector3* v) {
	return sqrtf((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
}

float vMag4p(const Vector4* v) {
	return sqrtf((v->x * v->x) + (v->y * v->y) + (v->z * v->z) + (v->w * v->w));
}





// Squared distance from one point to another

double vDistSq2i(Vector2i a, Vector2i b) { return vDistSq2ip(&a, &b); } 
float  vDistSq2(Vector2 a, Vector2 b) { return vDistSq2p(&a, &b); } 
float  vDistSq3(Vector3 a, Vector3 b) { return vDistSq3p(&a, &b); }
float  vDistSq4(Vector4 a, Vector4 b) { return vDistSq4p(&a, &b); }

double vDistSq2ip(Vector2i* a, Vector2i* b) {
	double x, y;
	
	x = a->x - b->x;
	y = a->y - b->y;
	
	return x*x + y*y;
}

float vDistSq2p(Vector2* a, Vector2* b) {
	float x, y;
	
	x = a->x - b->x;
	y = a->y - b->y;
	
	return x*x + y*y;
}

float vDistSq3p(Vector3* a, Vector3* b) {
	float x, y, z;
	
	x = a->x - b->x;
	y = a->y - b->y;
	z = a->z - b->z;
	
	return x*x + y*y + z*z;
}

float vDistSq4p(Vector4* a, Vector4* b) {
	float x, y, z, w;
	
	x = a->x - b->x;
	y = a->y - b->y;
	z = a->z - b->z;
	w = a->w - b->w;
	
	return x*x + y*y + z*z + w*w;
}



// Distance from one point to another

double vDist2i(Vector2i a, Vector2i b) { return vDist2ip(&a, &b); } 
float  vDist2(Vector2 a, Vector2 b) { return vDist2p(&a, &b); } 
float  vDist3(Vector3 a, Vector3 b) { return vDist3p(&a, &b); }
float  vDist4(Vector4 a, Vector4 b) { return vDist4p(&a, &b); }
double vDist2ip(Vector2i* a, Vector2i* b) { return sqrt(vDistSq2ip(a, b)); } 
float  vDist2p(Vector2* a, Vector2* b) { return sqrtf(vDistSq2p(a, b)); } 
float  vDist3p(Vector3* a, Vector3* b) { return sqrtf(vDistSq3p(a, b)); }
float  vDist4p(Vector4* a, Vector4* b) { return sqrtf(vDistSq4p(a, b)); }




// Vector normalize (scale to unit length)
Vector2 vNorm2(Vector2 v) { 
	Vector2 out;
	vNorm2p(&v, &out);
	return out;
}

Vector3 vNorm3(Vector3 v) { 
	Vector3 out;
	vNorm3p(&v, &out);
	return out;
}

Vector4 vNorm4(Vector4 v) { 
	Vector4 out;
	vNorm4p(&v, &out);
	return out;
}

void vNorm2p(const Vector2* v, Vector2* out) {
	float n;
	n = (v->x * v->x) + (v->y * v->y);
	
	if(n >= 1.0f - FLT_EPSILON && n <= 1.0f + FLT_EPSILON) {
		out->x = v->x;
		out->y = v->y;		
	}
	else if(n == 0.0) {
		out->x = 0.0;
		out->y = 0.0;
	}
	
	n = 1.0f / sqrtf(n);
	out->x = v->x * n;
	out->y = v->y * n;
} 

void  vNorm3p(const Vector3* v, Vector3* out) {
	float n;
	n = (v->x * v->x) + (v->y * v->y) + (v->z * v->z);
	
	if(n >= 1.0f - FLT_EPSILON && n <= 1.0f + FLT_EPSILON) {
		out->x = v->x;
		out->y = v->y;		
		out->z = v->z;		
	}
	else if(n == 0.0) {
		out->x = 0.0;
		out->y = 0.0;
		out->z = 0.0;
	}
	
	n = 1.0f / sqrtf(n);
	out->x = v->x * n;
	out->y = v->y * n;
	out->z = v->z * n;
}

void  vNorm4p(const Vector4* v, Vector4* out) {
	float n;
	n = (v->x * v->x) + (v->y * v->y) + (v->z * v->z) + (v->w * v->w);
	
	if(n >= 1.0f - FLT_EPSILON && n <= 1.0f + FLT_EPSILON) {
		out->x = v->x;
		out->y = v->y;		
		out->z = v->z;		
		out->w = v->w;		
	}
	else if(n == 0.0) {
		out->x = 0.0;
		out->y = 0.0;
		out->z = 0.0;
		out->w = 0.0;
	}
	
	n = 1.0f / sqrtf(n);
	out->x = v->x * n;
	out->y = v->y * n;
	out->z = v->z * n;
	out->w = v->w * n;
}

// alternate name
Vector2 vUnit2(Vector2 v) { return vNorm2(v); }
Vector3 vUnit3(Vector3 v) { return vNorm3(v); }
Vector4 vUnit4(Vector4 v) { return vNorm4(v); }
void  vUnit2p(const Vector2* v, Vector2* out) { return vNorm2p(v, out); } 
void  vUnit3p(const Vector3* v, Vector3* out) { return vNorm3p(v, out); } 
void  vUnit4p(const Vector4* v, Vector4* out) { return vNorm4p(v, out); } 



// Returns the minimum values of each component
Vector2i vMin2i(Vector2i a, Vector2i b) {
	Vector2i out;
	vMin2ip(&a, &b, &out);
	return out;
}

Vector2 vMin2(Vector2 a, Vector2 b) {
	Vector2 out;
	vMin2p(&a, &b, &out);
	return out;
}

Vector3 vMin3(Vector3 a, Vector3 b) {
	Vector3 out;
	vMin3p(&a, &b, &out);
	return out;
}

Vector4 vMin4(Vector4 a, Vector4 b) {
	Vector4 out;
	vMin4p(&a, &b, &out);
	return out;
}

void vMin2ip(Vector2i* a, Vector2i* b, Vector2i* out) {
	out->x = MIN(a->x, b->x);
	out->y = MIN(a->y, b->y);
}
void vMin2p(Vector2* a, Vector2* b, Vector2* out) {
	out->x = fminf(a->x, b->x);
	out->y = fminf(a->y, b->y);
}

void vMin3p(Vector3* a, Vector3* b, Vector3* out) {
	out->x = fminf(a->x, b->x);
	out->y = fminf(a->y, b->y);
	out->z = fminf(a->z, b->z);
}

void vMin4p(Vector4* a, Vector4* b, Vector4* out) {
	out->x = fminf(a->x, b->x);
	out->y = fminf(a->y, b->y);
	out->z = fminf(a->z, b->z);
	out->w = fminf(a->w, b->w);
}



// Returns the maximum values of each component
Vector2i vMax2i(Vector2i a, Vector2i b) {
	Vector2i out;
	vMax2ip(&a, &b, &out);
	return out;
}

Vector2 vMax2(Vector2 a, Vector2 b) {
	Vector2 out;
	vMax2p(&a, &b, &out);
	return out;
}

Vector3 vMax3(Vector3 a, Vector3 b) {
	Vector3 out;
	vMax3p(&a, &b, &out);
	return out;
}

Vector4 vMax4(Vector4 a, Vector4 b) {
	Vector4 out;
	vMax4p(&a, &b, &out);
	return out;
}

void vMax2ip(Vector2i* a, Vector2i* b, Vector2i* out) {
	out->x = MAX(a->x, b->x);
	out->y = MAX(a->y, b->y);
}
void vMax2p(Vector2* a, Vector2* b, Vector2* out) {
	out->x = fmaxf(a->x, b->x);
	out->y = fmaxf(a->y, b->y);
}

void vMax3p(Vector3* a, Vector3* b, Vector3* out) {
	out->x = fmaxf(a->x, b->x);
	out->y = fmaxf(a->y, b->y);
	out->z = fmaxf(a->z, b->z);
}

void vMax4p(Vector4* a, Vector4* b, Vector4* out) {
	out->x = fmaxf(a->x, b->x);
	out->y = fmaxf(a->y, b->y);
	out->z = fmaxf(a->z, b->z);
	out->w = fmaxf(a->w, b->w);
}






// feeding a zero vector into this will cause div/0 and you will deserve it
void  vProject3p(Vector3* what, Vector3* onto, Vector3* out) { // slower; onto may not be normalized
	float wdo = vDot3p(what, onto);
	float odo = vDot3p(onto, onto);
	vScale3p(onto, wdo / odo, out);
}

void  vProjectNorm3p(Vector3* what, Vector3* onto, Vector3* out) { // faster; onto must be normalized
	float wdo = vDot3p(what, onto);
	vScale3p(onto, wdo, out);
}



void vRandom3p(Vector3* end1, Vector3* end2, Vector3* out) {
	out->x = frand(fmin(end1->x, end2->x), fmax(end1->x, end2->x));
	out->y = frand(fmin(end1->y, end2->y), fmax(end1->y, end2->y));
	out->z = frand(fmin(end1->z, end2->z), fmax(end1->z, end2->z));
}

void vRandomNorm3p(Vector3* out) {
	float x = frand(-1, 1);
	float y = frand(-1, 1);
	float z = frand(-1, 1);
	
	float r = sqrt(x*x + y*y + z*z);
	if(r == 0.0) {
		vRandomNorm3p(out); // in the rare case of a zero-length vector, try again
		return;
	}
	
	out->x = x / r;
	out->y = y / r;
	out->z = z / r;
}



// reflects the distance from v to pivot across pivot.
// out, pivot, and v will form a straight line with pivot exactly in the middle.
void vReflectAcross3p(Vector3* v, Vector3* pivot, Vector3* out) {
	Vector3 diff;
	
	vSub3p(pivot, v, &diff);
	vAdd3p(pivot, &diff, out);
}

// calculate a unit vector normal to a triangle's face.
void  vTriFaceNormal3p(Vector3* a, Vector3* b, Vector3* c, Vector3* out) {
	Vector3 b_a, c_a;
	
	vSub3p(b, a, &b_a);
	vSub3p(c, a, &c_a);
	vCross3p(&b_a, &c_a, out);
	vNorm3p(out, out);
}

// calculate a unit vector normal to a triangle's face.
void  vpTriFaceNormal3p(Vector3* tri, Vector3* out) {
	vTriFaceNormal3p(tri+0, tri+1, tri+2, out);
}




void vProjectOntoPlane3p(Vector3* v, Plane* p, Vector3* out) {
	Vector3 v_ortho;
	
	// get the component of v that's perpendicular to the plane
	vProjectNorm3p(v, &p->n, &v_ortho);
	
	// subtract it from v
	vSub3p(v, &v_ortho, out);
}

void vProjectOntoPlaneNormalized3p(Vector3* v, Plane* p, Vector3* out) {
	vProjectOntoPlane3p(v, p, out);
	vNorm3p(out, out);
}


int IntersectPlaneRay3p(Plane* p, Ray3* r, Vector3* out) {
	float d = vDot3p(&p->n, &r->d);
	
	if(fabs(d) < FLT_CMP_EPSILON) return C3DLAS_DISJOINT;
	
	float t = 1.0 - (vDot3(p->n, r->o) + d) / vDot3(p->n, r->d); 
	
	if(t < 0) return C3DLAS_DISJOINT;
	
	*out = vAdd3(r->o, vScale3(r->d, t));
	
	return C3DLAS_INTERSECT;
}


// calculates a plane from a triangle
void planeFromTriangle3p(Vector3* v1, Vector3* v2, Vector3* v3, Plane* out) {
	vTriFaceNormal3p(v1, v2, v3, &out->n);
	out->d = vDot3p(&out->n, v1);
}

// copy a plane
void planeCopy3p(Plane* in, Plane* out) {
	vCopy3p(&in->n, &out->n);
	out->d = in->d;
}

// reverses the plane's direction
void planeInverse3p(Plane* in, Plane* out) {
	vInverse3p(&in->n, &out->n);
	out->d = -in->d;
}

// classifies a point by which side of the plane it's on, default espilon
int planeClassifyPoint3p(Plane* p, Vector3* pt) {
	return planeClassifyPointEps3p(p, pt, FLT_CMP_EPSILON);
}

// classifies a point by which side of the plane it's on, custom espilon
int planeClassifyPointEps3p(Plane* p, Vector3* pt, float epsilon) {
	float dist = vDot3p(&p->n, pt);
	if(fabs(dist - p->d) < epsilon)
		return C3DLAS_COPLANAR;
	else if (dist < p->d)
		return C3DLAS_BACK;
	else
		return C3DLAS_FRONT;
}



// C3DLAS_INTERSECT, _COPLANAR or _DISJOINT
int planeLineFindIntersect3p(Plane* pl, Vector3* la, Vector3* lb, Vector3* out) {
	Vector3 ldir;
	float da, db;
	
	vSub3p(lb, la, &ldir);
	
	da = vDot3p(la, &pl->n) - pl->d;
	
	// bail if the line and plane are parallel
	if(fabs(vDot3p(&pl->n, &ldir)) < FLT_CMP_EPSILON) {
		
		// check coplanarity
		if(fabs(da) < FLT_CMP_EPSILON) {
			return C3DLAS_COPLANAR; // the end is on the plane, so the other is too
		}
		
		return C3DLAS_DISJOINT;
	}
	

	db = vDot3p(lb, &pl->n) - pl->d;
	
	// check if one of the points is on the plane
	if(fabs(da) < FLT_CMP_EPSILON) {
		*out = *la;
		return C3DLAS_INTERSECT;
	}
	if(fabs(db) < FLT_CMP_EPSILON) {
		*out = *lb;
		return C3DLAS_INTERSECT;
	}
	
	Vector3 p0, g, j;
	vScale3p(&pl->n, pl->d, &p0);
	vSub3p(&p0, la, &g);
	float h = vDot3p(&g, &pl->n);
	float i = vDot3p(&ldir, &pl->n);
	float d = i != 0 ? h / i : 0;
	
	// check if the plane intersects outside the two points
	if(d < 0 || d > vDist3p(la, lb)) {
		return C3DLAS_DISJOINT;
	}
	
	vScale3p(&ldir, d, &j);
	vAdd3p(la, &j, out);
	
	return C3DLAS_INTERSECT;
}


// Assumes full proper intersection.
// C3DLAS_INTERSECT
int planeLineFindIntersectFast3p(Plane* pl, Vector3* la, Vector3* lb, Vector3* out) {
	Vector3 ldir, p0, g, j;
	float h, i, d;
	
	vSub3p(lb, la, &ldir);
	
	vScale3p(&pl->n, pl->d, &p0);
	vSub3p(&p0, la, &g);
	h = vDot3p(&g, &pl->n);
	i = vDot3p(&ldir, &pl->n);
	d = i != 0 ? h / i : 0;
	
	vScale3p(&ldir, d, &j);
	vAdd3p(la, &j, out);
	
	return C3DLAS_INTERSECT;
}


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
) {
	Vector3 ldir, c;
	float da, db;
	
	vSub3p(lb, la, &ldir);
	
	da = vDot3p(la, &pl->n) - pl->d;
	
	// bail if the line and plane are parallel
	if(fabs(vDot3p(&pl->n, &ldir)) < FLT_CMP_EPSILON) {
		*aboveCnt = 0;
		*belowCnt = 0;
			
		// check coplanarity
		if(fabs(da) < FLT_CMP_EPSILON) {
			return C3DLAS_COPLANAR; // the end is on the plane, so the other is too
		}
		
		return C3DLAS_PARALLEL;
	}
	

	db = vDot3p(lb, &pl->n) - pl->d;
	
	// check if one of the points is on the plane
	if(fabs(da) < FLT_CMP_EPSILON) {
		if(db > 0) {
			aboveOut[0] = *la; // correct ordering
			aboveOut[1] = *lb;
			*aboveCnt = 1;
			*belowCnt = 0;
		}
		else {
			belowOut[0] = *la; // correct ordering
			belowOut[1] = *lb;
			*aboveCnt = 0;
			*belowCnt = 1;
		}
		
		return C3DLAS_INTERSECT;
	}
	if(fabs(db) < FLT_CMP_EPSILON) {
		if(da > 0) {
			aboveOut[0] = *la; // correct ordering
			aboveOut[1] = *lb;
			*aboveCnt = 1;
			*belowCnt = 0;
		}
		else {
			belowOut[0] = *la; // correct ordering
			belowOut[1] = *lb;
			*aboveCnt = 0;
			*belowCnt = 1;
		}
		
		return C3DLAS_INTERSECT;
	}
	
	// calculate itnersection point, c
	Vector3 p0, g, j;
	vScale3p(&pl->n, pl->d, &p0);
	vSub3p(&p0, la, &g);
	float h = vDot3p(&g, &pl->n);
	float i = vDot3p(&ldir, &pl->n);
	float d = i != 0 ? h / i : 0;
	
	// check if the plane intersects outside the two points
	if(d < 0 || d > vDist3p(la, lb)) {
		if(da > 0) {
			aboveOut[0] = *la; // correct ordering
			aboveOut[1] = *lb;
			*aboveCnt = 1;
			*belowCnt = 0;
		}
		else {
			belowOut[0] = *la; // correct ordering
			belowOut[1] = *lb;
			*aboveCnt = 0;
			*belowCnt = 1;
		}
		
		return C3DLAS_DISJOINT;
	}
	
	vScale3p(&ldir, d, &j);
	vAdd3p(la, &j, &c);
	
	if(da > 0) {
		aboveOut[0] = *la; // correct ordering
		aboveOut[1] = c;
		belowOut[0] = c;
		belowOut[1] = *lb;
	}
	else {
		belowOut[0] = *la; // correct ordering
		belowOut[1] = c;
		aboveOut[0] = c;
		aboveOut[1] = *lb;
	}
	
	*aboveCnt = 1;
	*belowCnt = 1;
	
	return C3DLAS_INTERSECT;
}


// C3DLAS_COPLANAR, _INTERSECT, or _DISJOINT
int triPlaneTestIntersect3p(Vector3* pTri, Plane* pl) {
	Vector3 a, b, c;
	float da, db, dc;
	
	// get distance of each vertex from the plane
	// bail early if any of them are coplanar
	a = pTri[0];
	da = vDot3p(&a, &pl->n) - pl->d;
	if(fabs(da) < FLT_CMP_EPSILON) {
		return C3DLAS_COPLANAR;
	}
	
	b = pTri[1];
	db = vDot3p(&b, &pl->n) - pl->d;
	if(fabs(db) < FLT_CMP_EPSILON) {
		return C3DLAS_COPLANAR;
	}
	
	c = pTri[2];
	dc = vDot3p(&c, &pl->n) - pl->d;
	if(fabs(dc) < FLT_CMP_EPSILON) {
		return C3DLAS_COPLANAR;
	}
	
	// the triangle intersects if the sign of all the distances does not match,
	// ie, on vertex is on the opposite side of the plane from the others
	return (signbit(da) == signbit(db) && signbit(db) == signbit(dc)) ? C3DLAS_DISJOINT : C3DLAS_INTERSECT;
}


// C3DLAS_COPLANAR, _INTERSECT, or _DISJOINT
int triPlaneClip3p(
	Vector3* pTri, 
	Plane* pl, 
	Vector3* aboveOut, 
	Vector3* belowOut, 
	int* aboveCnt,
	int* belowCnt
) {
	
	Vector3 v0, v1, v2;
	float vp_d0, vp_d1, vp_d2;
	
	v0 = pTri[0];
	v1 = pTri[1];
	v2 = pTri[2];

	// get distance of each vertex from the plane
	vp_d0 = vDot3p(&v0, &pl->n) - pl->d;
	vp_d1 = vDot3p(&v1, &pl->n) - pl->d;
	vp_d2 = vDot3p(&v2, &pl->n) - pl->d;
	
	
	// bail early if just one is coplanar
	// split in half with single-edge intersections
	if(fabs(vp_d0) < FLT_CMP_EPSILON) {
		if( // single edge intersection
			signbit(vp_d1) != signbit(vp_d2) && 
			fabs(vp_d1) > FLT_CMP_EPSILON &&
			fabs(vp_d2) > FLT_CMP_EPSILON
		) {
			// get intersection point
			Vector3 c;
			planeLineFindIntersectFast3p(pl, &v1, &v2, &c);
			
			if(vp_d1 > 0) { // v1 is above the plane
				aboveOut[0] = c; // correct winding
				aboveOut[1] = v0;
				aboveOut[2] = v1;
				belowOut[0] = c;
				belowOut[1] = v2;
				belowOut[2] = v0;
			}
			else {
				belowOut[0] = c; // correct winding
				belowOut[1] = v0;
				belowOut[2] = v1;
				aboveOut[0] = c;
				aboveOut[1] = v2;
				aboveOut[2] = v0;
			}
			
			*aboveCnt = 1;
			*belowCnt = 1;
			
			return C3DLAS_INTERSECT;
		}
		
		return C3DLAS_COPLANAR; // one vertex is on the plane, the others all above or below
	}
	
	if(fabs(vp_d1) < FLT_CMP_EPSILON) {
		if( // single edge intersection
			signbit(vp_d0) != signbit(vp_d2) && 
			fabs(vp_d0) > FLT_CMP_EPSILON &&
			fabs(vp_d2) > FLT_CMP_EPSILON
		) {
			// get intersection point
			Vector3 c;
			planeLineFindIntersectFast3p(pl, &v0, &v2, &c);
			
			if(vp_d0 > 0) { // v0 is above the plane
				aboveOut[0] = c; // correct winding
				aboveOut[1] = v0;
				aboveOut[2] = v1;
				belowOut[0] = c;
				belowOut[1] = v1;
				belowOut[2] = v2;
			}
			else {
				belowOut[0] = c; // correct winding
				belowOut[1] = v0;
				belowOut[2] = v1;
				aboveOut[0] = c;
				aboveOut[1] = v1;
				aboveOut[2] = v2;
			}
			
			*aboveCnt = 1;
			*belowCnt = 1;
			
			return C3DLAS_INTERSECT;
		}
		
		return C3DLAS_COPLANAR; // one vertex is on the plane, the others all above or below
	}
	
	if(fabs(vp_d2) < FLT_CMP_EPSILON) {
		if( // single edge intersection
			signbit(vp_d0) != signbit(vp_d1) && 
			fabs(vp_d0) > FLT_CMP_EPSILON &&
			fabs(vp_d1) > FLT_CMP_EPSILON
		) {
			// get intersection point
			Vector3 c;
			planeLineFindIntersectFast3p(pl, &v0, &v1, &c);
			
			if(vp_d0 > 0) { // v0 is above the plane
				aboveOut[0] = c; // correct winding
				aboveOut[1] = v2;
				aboveOut[2] = v0;
				belowOut[0] = c;
				belowOut[1] = v1;
				belowOut[2] = v2;
			}
			else {
				belowOut[0] = c; // correct winding
				belowOut[1] = v2;
				belowOut[2] = v0;
				aboveOut[0] = c;
				aboveOut[1] = v1;
				aboveOut[2] = v2;
			}
			
			*aboveCnt = 1;
			*belowCnt = 1;
			
			return C3DLAS_INTERSECT;
		}
		
		return C3DLAS_COPLANAR; // one vertex is on the plane, the others all above or below
	}
	
	
	// the triangle intersects if the sign of all the distances does not match,
	// ie, on vertex is on the opposite side of the plane from the others
	// bail if disjoint
	if(signbit(vp_d0) == signbit(vp_d1) && signbit(vp_d1) == signbit(vp_d2)) {
		return C3DLAS_DISJOINT;
	}
	
	
	// split on which edges intersect the plane
	if(signbit(vp_d0) == signbit(vp_d1)) {
		// vertex 2 is isolated; edges 0,2 and 1,2 intersect
		
		Vector3 c0, c1;
		planeLineFindIntersectFast3p(pl, &v0, &v2, &c0);
		planeLineFindIntersectFast3p(pl, &v1, &v2, &c1);
		
		if(vp_d2 > 0) { // v2 is above the plane
			aboveOut[0] = v2; // correct winding
			aboveOut[1] = c0;
			aboveOut[2] = c1;
			belowOut[0] = c1;
			belowOut[1] = v0;
			belowOut[2] = v1;
			belowOut[3] = c1;
			belowOut[4] = c0;
			belowOut[5] = v1;
			
			*aboveCnt = 1;
			*belowCnt = 2;
		}
		else {
			belowOut[0] = v2; // correct winding
			belowOut[1] = c0;
			belowOut[2] = c1;
			aboveOut[0] = c1;
			aboveOut[1] = v0;
			aboveOut[2] = v1;
			aboveOut[3] = c1;
			aboveOut[4] = c0;
			aboveOut[5] = v1;
			
			*aboveCnt = 2;
			*belowCnt = 1;
		}
		
	}
	else if(signbit(vp_d1) == signbit(vp_d2)) {
		// vertex 0 is isolated; edges 1,0 and 2,0 intersect
		
		Vector3 c0, c1;
		planeLineFindIntersectFast3p(pl, &v1, &v0, &c0);
		planeLineFindIntersectFast3p(pl, &v2, &v0, &c1);
		
		if(vp_d0 > 0) { // v0 is above the plane
			aboveOut[0] = v0; // correct winding
			aboveOut[1] = c0;
			aboveOut[2] = c1;
			belowOut[0] = c1;
			belowOut[1] = v1;
			belowOut[2] = v2;
			belowOut[3] = c1;
			belowOut[4] = c0;
			belowOut[5] = v1;
			
			*aboveCnt = 1;
			*belowCnt = 2;
		}
		else {
			belowOut[0] = v0; // correct winding
			belowOut[1] = c0;
			belowOut[2] = c1;
			aboveOut[0] = c1;
			aboveOut[1] = v1;
			aboveOut[2] = v2;
			aboveOut[3] = c1;
			aboveOut[4] = c0;
			aboveOut[5] = v1;
			
			*aboveCnt = 2;
			*belowCnt = 1;
		}
		
	}
	else {
		// vertex 1 is isolated; edges 0,1 and 2,1 intersect
		
		Vector3 c0, c1;
		planeLineFindIntersectFast3p(pl, &v0, &v1, &c0);
		planeLineFindIntersectFast3p(pl, &v2, &v1, &c1);
		
		if(vp_d1 > 0) { // v1 is above the plane
			aboveOut[0] = v1; // correct winding
			aboveOut[1] = c1;
			aboveOut[2] = c0;
			belowOut[0] = c1;
			belowOut[1] = v2;
			belowOut[2] = v0;
			belowOut[3] = c0;
			belowOut[4] = c1;
			belowOut[5] = v0;
			
			*aboveCnt = 1;
			*belowCnt = 2;
		}
		else {
			belowOut[0] = v1; // correct winding
			belowOut[1] = c1;
			belowOut[2] = c0;
			aboveOut[0] = c1;
			aboveOut[1] = v2;
			aboveOut[2] = v0;
			aboveOut[3] = c0;
			aboveOut[4] = c1;
			aboveOut[5] = v0;
			
			*aboveCnt = 2;
			*belowCnt = 1;
		}
	}
	
	
	return C3DLAS_INTERSECT;
}

// http://geomalgorithms.com/a07-_distance.html
// _PARALLEL with no output on parallel lines
// _INTERSECT with one point of output on intersection
// _DISJOINT with two outputs otherwise
int shortestLineFromRayToRay3p(Ray3* r1, Ray3* r2, Vector3* pOut) {
	
	Vector3 u, v, w, ps, pt;
	float a, b, c, d, e, s, t;
	
	u = r1->d;
	v = r2->d;
	vSub3p(&r1->o, &r2->o, &w);
	
	a = vDot3p(&u, &u); 
	b = vDot3p(&u, &v); 
	c = vDot3p(&v, &v); 
	d = vDot3p(&u, &w); 
	e = vDot3p(&v, &w); 
	
	float ac_bb = a * c - b * b;
	if(fabs(ac_bb) < FLT_CMP_EPSILON) {
		return C3DLAS_PARALLEL;
	}
	
	s = (b * e - c * d) / ac_bb;
	t = (a * e - b * d) / ac_bb;
	
	vScale3p(&u, s, &ps);
	vScale3p(&v, t, &pt);
	vAdd3p(&r1->o, &ps, &ps);
	vAdd3p(&r2->o, &pt, &pt);
	
	pOut[0] = ps;
	pOut[1] = pt;
	
	if(vDistSq3p(&ps, &pt) < FLT_CMP_EPSILON_SQ) {
		return C3DLAS_INTERSECT;
	}
	
	return C3DLAS_DISJOINT;
}

void frustumFromMatrix(Matrix* m, Frustum* out) {
	
	Matrix inv;
	
	mInverse(m, &inv);
	
	// first the points
	// these MUST be in this order
	// near
	vMatrixMulf3p(-1,-1,-1, &inv, &out->points[0]);
	vMatrixMulf3p(-1, 1,-1, &inv, &out->points[1]);
	vMatrixMulf3p( 1,-1,-1, &inv, &out->points[2]);
	vMatrixMulf3p( 1, 1,-1, &inv, &out->points[3]);
	// far
	vMatrixMulf3p(-1,-1, 1, &inv, &out->points[4]);
	vMatrixMulf3p(-1, 1, 1, &inv, &out->points[5]);
	vMatrixMulf3p( 1,-1, 1, &inv, &out->points[6]);
	vMatrixMulf3p( 1, 1, 1, &inv, &out->points[7]);
	
	// now the planes
	// near and far
	planeFromTriangle3p(&out->points[0], &out->points[1], &out->points[2], &out->planes[0]);
	planeFromTriangle3p(&out->points[4], &out->points[5], &out->points[6], &out->planes[1]);
	// sides
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[1], &out->planes[2]);
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[2], &out->planes[3]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[1], &out->planes[4]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[2], &out->planes[5]);
}


void frustumCenter(Frustum* f, Vector3* out) {
	Vector3 sum = {0.0f,0.0f,0.0f};
	
	for(int i = 0; i < 8; i++) vAdd3p(&f->points[i], &sum, &sum);
	vScale3p(&sum, 1.0f/8.0f, out);
}

// General idea of the algorithm:
// https://lxjk.github.io/2017/04/15/Calculate-Minimal-Bounding-Sphere-of-Frustum.html
// http://archive.is/YACj2
void frustumBoundingSphere(Frustum* f, Sphere* out) {
	Vector3 f0, n0;
	vPointAvg3p(&f->points[0], &f->points[3], &n0);
	vPointAvg3p(&f->points[4], &f->points[7], &f0);
	
	float Dn2 = vDistSq3p(&n0, &f->points[0]); 
	float Df2 = vDistSq3p(&f0, &f->points[4]);
	
	 // check for ortho
	if(Dn2 - Df2 < 0.00001) {
		frustumCenter(f, &out->center);
		out->r = vDist3p(&out->center, &f->points[0]); 
		return;
	}
	
	float Dnf = vDist3p(&f0, &n0);
	float Dnc = (Dn2 - Df2 - Df2) / (2 * Dnf);
	
// 	printf("\n f: %f,%f,%f\n", f->points[4].x,f->points[4].y,f->points[4].z);
// 	printf(" n: %f,%f,%f\n", f->points[0].x,f->points[0].y,f->points[0].z);
// 	printf(" f0: %f,%f,%f\n", f0.x,f0.y,f0.z);
// 	printf(" n0: %f,%f,%f\n", n0.x,n0.y,n0.z);
// 	printf(" dn2, df2, dnf, dnc: %f,%f,%f,%f\n", Dn2, Df2, Dnf, Dnc);
	
	
	if(Dnc > 0 && Dnc < Dnf) {
		vLerp3p(&f0, &n0, Dnc / Dnf, &out->center);
		out->r = sqrt(Dnc * Dnc + Dn2);
	}
	else {
		out->center = f0;
		out->r = sqrt(Df2);
	}
}


void frustumInscribeSphere(Frustum* f, Sphere* out) {
	Vector3 fx, nx;
	vPointAvg3p(&f->points[0], &f->points[3], &nx);
	vPointAvg3p(&f->points[4], &f->points[7], &fx);
	
/*	
	float Dn2 = vDistSq3p(&n0, &f->points[0]); 
	float Df2 = vDistSq3p(&f0, &f->points[4]);
	float Dnf = vDist3p(&f0, n0);
	float Dnc = (Dn2 - Df2 - Df2) / (2 * Dnf);*/

}



void quadCenterp3p(Vector3* a, Vector3* b, Vector3* c, Vector3* d, Vector3* out) {
	Vector3 sum;
	vAdd3p(a, b, &sum);
	vAdd3p(&sum, c, &sum);
	vAdd3p(&sum, d, &sum);
	vScale3p(&sum, 0.25f, out);
}

// closest distance from an arbitrary point to the plane 
float planePointDist3p(Plane* pl, Vector3* p) {
	Vector3 a;
	vScale3p(&pl->n, pl->d, &a);
	return fabs(vDot3p(&a, p));
} 

// signed closest distance from an arbitrary point to the plane 
float planePointDistSigned3p(Plane* pl, Vector3* p) {
	Vector3 a;
	vScale3p(&pl->n, pl->d, &a);
	return vDot3p(&a, p);
} 

void vPointAvg3p(Vector3* a, Vector3* b, Vector3* out) {
	Vector3 sum;
	vAdd3p(a, b, &sum);
	vScale3p(&sum, 0.5f, out);
} 




// reflects the distance from v to pivot across pivot.
// out, pivot, and v will form a straight line with pivot exactly in the middle.
void vReflectAcross2p(Vector2* v, Vector2* pivot, Vector2* out) {
	Vector2 diff;
	
	vSub2p(pivot, v, &diff);
	vAdd2p(pivot, &diff, out);
}


// degenerate cases may not give desired results. GIGO.
void vRoundAway2p(const Vector2* in, const Vector2* center, Vector2i* out) {
	
	if(in->x > center->x) out->x = ceilf(in->x);
	else out->x = floorf(in->x);
	
	if(in->y > center->y) out->y = ceilf(in->y);
	else out->y = floorf(in->y);
}

// degenerate cases may not give desired results. GIGO.
void vRoundToward2p(const Vector2* in, const Vector2* center, Vector2i* out) {
	
	if(in->x > center->x) out->x = floorf(in->x);
	else out->x = ceilf(in->x);
	
	if(in->y > center->y) out->y = floorf(in->y);
	else out->y = ceilf(in->y);
}


// returns the *signed* area of a triangle. useful for determining winding
// positive values mean a clockwise triangle
float triArea2p(Vector2* a, Vector2* b, Vector2* c) {
	return 0.5 * (
		((b->x - a->x) * (b->y + a->y)) +
		((c->x - b->x) * (c->y + b->y)) +
		((a->x - c->x) * (a->y + c->y)));
}


// determines if a point is inside a triangle
int triPointInside2p(Vector2* p, Vector2* a, Vector2* b, Vector2* c) {
	int d = signbit((p->x - b->x) * (a->y - b->y) - (a->x - b->x) * (p->y - b->y));
	int e = signbit((p->x - c->x) * (b->y - c->y) - (b->x - c->x) * (p->y - c->y));
	if(d != e) return 0;
	int f = signbit((p->x - a->x) * (c->y - a->y) - (c->x - a->x) * (p->y - a->y));
	return e == f;
}





// plane-vector operations

// distance from point to plane
float pvDist3p(Plane* p, Vector3* v) {
	return vDot3p(v, &p->n) + p->d;
}



// matrix-vector operations


// multiply a vector by a matrix
void vMatrixMul3p(Vector3* restrict in, Matrix* restrict m, Vector3* restrict out) {
	vMatrixMulf3p(in->x, in->y, in->z, m, out);
}

void vMatrixMulf3p(float x, float y, float z, Matrix* restrict m, Vector3* restrict out) {
	Vector4 v;

	v.x = x * m->m[0+0] + y * m->m[4+0] + z * m->m[8+0] + 1 * m->m[12+0];
	v.y = x * m->m[0+1] + y * m->m[4+1] + z * m->m[8+1] + 1 * m->m[12+1];
	v.z = x * m->m[0+2] + y * m->m[4+2] + z * m->m[8+2] + 1 * m->m[12+2];
	v.w = x * m->m[0+3] + y * m->m[4+3] + z * m->m[8+3] + 1 * m->m[12+3];
	
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
	m.m[15] = 1.0;
	
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
void mPerspSetNF(Matrix* m, float near, float far) {
	m->m[10] = -(far + near) / (far - near);
	m->m[14] = (-2.0 * far * near) / (far - near);
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


void mOrthoFromSphere(Sphere* s, Vector3* eyePos, Matrix* out) {
	Matrix m;
	
	float right = -s->r ;
	float left = s->r;
	float top = s->r;
	float bottom = -s->r;
	float near = -s->r;
	float far = s->r;
	
	// this is the ortho projection matrix
	m = IDENT_MATRIX;
	m.m[0] = 2 / (right - left);
	m.m[5] = 2 / (top - bottom);
	m.m[10] = -2 / (far - near);
	m.m[12] = -(right + left) / (right - left);
	m.m[13] = -(top + bottom) / (top - bottom);
	m.m[14] = -(far + near) / (far - near);
	m.m[15] = 1;
	
	Vector3 d;
	vSub3p(&s->center, eyePos, &d);
	vNorm3p(&d, &d);
	
	Matrix m2;

	
	
	m2 = IDENT_MATRIX;
// 	mRotX(asin(d.z) + (F_PI / 4)  , &m2);
	mRotY(atan2(d.z, d.x) + (F_PI / 2)  , &m2);
	mMul(&m2, &m);

	m2 = IDENT_MATRIX;
	mRotZ(asin(d.y), &m2);
	mMul(&m2, &m);
	

	m2 = IDENT_MATRIX;
	Vector3 ic;
	vScale3p(&s->center, -1, &ic);
	mTransv(&ic, &m2);
	
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


// analgous to gluLookAt
// BUG: very broken apparently
// https://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
void mLookAt(Vector3* eye, Vector3* center, Vector3* up, Matrix* out) {
	
	Vector3 f, upn, s, u, sn;
	Matrix m, m2;
	m2 = IDENT_MATRIX;
	
	vSub3p(center, eye, &f);
	vNorm3p(&f, &f);
	
	vNorm3p(up, &upn);
	
	vCross3p(&f, &upn, &s);
	vNorm3p(&s, &sn);
	
	vCross3p(&sn, &f, &u);
	
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
	mLookAt(eye, center, up, msGetTop(ms));
}



void evalBezier3p(Vector3* e1, Vector3* e2, Vector3* c1, Vector3* c2, float t, Vector3* out) {
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



void evalBezierTangent3p(Vector3* e1, Vector3* e2, Vector3* c1, Vector3* c2, float t, Vector3* out) {
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

void evalBezierNorm3p(Vector3* e1, Vector3* e2, Vector3* c1, Vector3* c2, float t, Vector3* out) {
	out->x = evalBezier1D_ddt(e1->x, e2->x, c1->x, c2->x, t);
	out->y = evalBezier1D_ddt(e1->y, e2->y, c1->y, c2->y, t);
	out->z = evalBezier1D_ddt(e1->z, e2->z, c1->z, c2->z, t);
}


// Quadratic bezier functions
float evalQBezier1D(float e1, float e2, float c1, float t) {
	float mt, mt2;
	mt = 1 - t;
	mt2 = mt * mt;
	
	return (mt2 * e1) + (2 * mt * t * c1) + (t * t * e2);
}

void evalQBezier2D3p(Vector2* e1, Vector2* e2, Vector2* c1, float t, Vector2* out) {
	out->x = evalQBezier1D(e1->x, e2->x, c1->x, t);
	out->y = evalQBezier1D(e1->y, e2->y, c1->y, t);
}

void evalQBezier3p(Vector3* e1, Vector3* e2, Vector3* c1, float t, Vector3* out) {
	out->x = evalQBezier1D(e1->x, e2->x, c1->x, t);
	out->y = evalQBezier1D(e1->y, e2->y, c1->y, t);
	out->z = evalQBezier1D(e1->z, e2->z, c1->z, t);
}




///// bounding box functions


// 3D versions

int boxDisjoint3p(const AABB3* a, const AABB3* b) {
	return a->max.x < b->min.x || b->max.x < a->min.x
		|| a->max.y < b->min.y || b->max.y < a->min.y
		|| a->max.z < b->min.z || b->max.z < a->min.z;
}

int boxOverlaps3p(const AABB3* a, const AABB3* b) {
	return !boxDisjoint3p(a, b);
}



int boxContainsPoint3p(const AABB3* b, const Vector3* p) {
	return b->min.x <= p->x && b->max.x >= p->x
		&& b->min.y <= p->y && b->max.y >= p->y
		&& b->min.z <= p->z && b->max.z >= p->z;
}


void boxCenter3p(const AABB3* b, Vector3* out) {
	out->x = (b->max.x + b->min.x) / 2;
	out->y = (b->max.y + b->min.y) / 2;
	out->z = (b->max.z + b->min.z) / 2;
}

void boxSize3p(const AABB3* b, Vector3* out) {
	out->x = b->max.x - b->min.x;
	out->y = b->max.y - b->min.y;
	out->z = b->max.z - b->min.z;
}



// 2D versions

int boxDisjoint2p(const AABB2* a, const AABB2* b) {

	return a->max.x < b->min.x || b->max.x < a->min.x
		|| a->max.y < b->min.y || b->max.y < a->min.y;
}

int boxOverlaps2p(const AABB2* a, const AABB2* b) {
	return !boxDisjoint2p(a, b);
}



int boxContainsPoint2p(const AABB2* b, const Vector2* p) {
	return b->min.x <= p->x && b->max.x >= p->x
		&& b->min.y <= p->y && b->max.y >= p->y;
}


void boxCenter2p(const AABB2* b, Vector2* out) {
	out->x = (b->max.x + b->min.x) / 2;
	out->y = (b->max.y + b->min.y) / 2;
}

void boxSize2p(const AABB2* b, Vector2* out) {
	out->x = b->max.x - b->min.x;
	out->y = b->max.y - b->min.y;
}

void boxQuadrant2p(const AABB2* in, char ix, char iy, AABB2* out) {
	Vector2 sz, c;
	
	boxCenter2p(in, &c);
	boxSize2p(in, &sz);
	sz.x *= .5;
	sz.y *= .5;
	
	out->min.x = c.x - (ix ? 0.0f : sz.x);
	out->min.y = c.y - (iy ? 0.0f : sz.y);
	out->max.x = c.x + (ix ? sz.x : 0.0f);
	out->max.y = c.y + (iy ? sz.y : 0.0f);
}


// 2D integer versions

int boxDisjoint2ip(const AABB2i* a, const AABB2i* b) {

	return a->max.x < b->min.x || b->max.x < a->min.x
		|| a->max.y < b->min.y || b->max.y < a->min.y;
}

int boxOverlaps2ip(const AABB2i* a, const AABB2i* b) {
	return !boxDisjoint2ip(a, b);
}



int boxContainsPoint2ip(const AABB2i* b, const Vector2i* p) {
	return b->min.x <= p->x && b->max.x >= p->x
		&& b->min.y <= p->y && b->max.y >= p->y;
}


void boxCenter2ip(const AABB2i* b, Vector2* out) {
	out->x = (b->max.x + b->min.x) / 2.0f;
	out->y = (b->max.y + b->min.y) / 2.0f;
}

void boxSize2ip(const AABB2i* b, Vector2* out) {
	out->x = b->max.x - b->min.x;
	out->y = b->max.y - b->min.y;
}

// BUG: needs some fancy math work to keep everything tight. integers don't split nicely
void boxQuadrant2ip(const AABB2i* in, char ix, char iy, AABB2i* out) {
	Vector2 sz, c;
	
	printf("fix me: %s:%d", __FILE__, __LINE__);
	exit(666);
	
	boxCenter2ip(in, &c);
	boxSize2ip(in, &sz);
	sz.x *= .5;
	sz.y *= .5;
	
	out->min.x = c.x - (ix ? 0.0f : sz.x);
	out->min.y = c.y - (iy ? 0.0f : sz.y);
	out->max.x = c.x + (ix ? sz.x : 0.0f);
	out->max.y = c.y + (iy ? sz.y : 0.0f);
}



// find the center of a quad
void quadCenter2p(const Quad2* q, Vector2* out) {
	Vector2 c = {0};
	int i;
	
	for(i = 0; i < 4; i++) {
		c.x += q->v[i].x;
		c.y += q->v[i].y;
	}
	
	out->x = c.x / 4;
	out->y = c.y / 4;
}


void quadRoundOutward2p(const Quad2* in, Quad2i* out) {
	Vector2 c;
	int i;
	
	quadCenter2p(in, &c);
	
	for(i = 0; i < 4; i++)
		vRoundAway2p(&in->v[i], &c, &out->v[i]);
}

void quadRoundInward2p(const Quad2* in, Quad2i* out) {
	Vector2 c;
	int i;
	
	quadCenter2p(in, &c);
	
	for(i = 0; i < 4; i++)
		vRoundToward2p(&in->v[i], &c, &out->v[i]);
}


int quadIsPoint2i(const Quad2i* q) {
	return (
		(q->v[0].x == q->v[1].x) && (q->v[1].x == q->v[2].x) && (q->v[2].x == q->v[3].x) &&
		(q->v[0].y == q->v[1].y) && (q->v[1].y == q->v[2].y) && (q->v[2].y == q->v[3].y)
	);
}

int quadIsAARect2i(const Quad2i* q) {
	return (
		q->v[0].x == q->v[3].x && q->v[1].x == q->v[2].x &&
		q->v[0].y == q->v[1].y && q->v[2].y == q->v[3].y);
}


// ray stuff
void makeRay3p(Vector3* origin, Vector3* direction, Ray3* out) {
	
	out->o.x = origin->x;
	out->o.y = origin->y;
	out->o.z = origin->z;
	
	vNorm3p(direction, &out->d);
}

// ray stuff
void makeRay2(Vector2* origin, Vector2* direction, Ray2* out) {
	
	out->o.x = origin->x;
	out->o.y = origin->y;
	
	vNorm2p(direction, &out->d);
}

// this version has no branching, but only answers yes or no.
// algorithm explanation here. hopefully my extrapolation into 3 dimensions is correct.
// http://tavianator.com/fast-branchless-raybounding-box-intersections/
int boxRayIntersectFast3p(const AABB3* b, const Ray3* r) {
	Vector3 t1, t2;
	float tmin, tmax;
	Vector3 id;
	
	vInverse3p(&r->d, &id);
	
	t1.x = (b->min.x - r->o.x) * id.x;
	t2.x = (b->max.x - r->o.x) * id.x;
	tmin = fmin(t1.x, t2.x);
	tmax = fmax(t1.x, t2.x);
	
	t1.y = (b->min.y - r->o.y) * id.y;
	t2.y = (b->max.y - r->o.y) * id.y;
	tmin = fmax(tmin, fmin(t1.y, t2.y));
	tmax = fmin(tmax, fmax(t1.y, t2.y));
	
	t1.z = (b->min.z - r->o.z) * id.z;
	t2.z = (b->max.z - r->o.z) * id.z;
	tmin = fmax(tmin, fmin(t1.z, t2.z));
	tmax = fmin(tmax, fmax(t1.z, t2.z));

	return tmax >= tmin && tmax > 0.0f;
}

// this version has no branching, but only answers yes or no.
// http://tavianator.com/fast-branchless-raybounding-box-intersections/
int boxRayIntersectFast2(const AABB2* b, const Ray2* r) {
	Vector2 t1, t2;
	float tmin, tmax;
	Vector2 id;
	
	vInverse2p(&r->d, &id);
	
	t1.x = (b->min.x - r->o.x) * id.x;
	t2.x = (b->max.x - r->o.x) * id.x;
	tmin = fmin(t1.x, t2.x);
	tmax = fmax(t1.x, t2.x);
	
	t1.y = (b->min.y - r->o.y) * id.y;
	t2.y = (b->max.y - r->o.y) * id.y;
	tmin = fmax(tmin, fmin(t1.y, t2.y));
	tmax = fmin(tmax, fmax(t1.y, t2.y));
	
	return tmax >= tmin && tmax > 0.0f;
}


// this version gives the point of intersection as well as distance
// algorithm explanation here. hopefully my extrapolation into 3 dimensions is correct.
// http://tavianator.com/fast-branchless-raybounding-box-intersections/
int boxRayIntersect3p(const AABB3* b, const Ray3* r, Vector3* ipoint, float* idist) {
	Vector3 t1, t2, id;
	float tmin, tmax;
	
	vInverse3p(&r->d, &id);
		
	t1.x = (b->min.x - r->o.x) * id.x;
	t2.x = (b->max.x - r->o.x) * id.x;
	tmin = fmin(t1.x, t2.x);
	tmax = fmax(t1.x, t2.x);
	
	t1.y = (b->min.y - r->o.y) * id.y;
	t2.y = (b->max.y - r->o.y) * id.y;
	tmin = fmax(tmin, fmin(t1.y, t2.y));
	tmax = fmin(tmax, fmax(t1.y, t2.y));
	
	t1.z = (b->min.z - r->o.z) * id.z;
	t2.z = (b->max.z - r->o.z) * id.z;
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
	vReflectAcross2p(&n->c, &n->e, &out[2]);
	
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
	/*
	float tIncrement;
	
	
	tIncrement = (float)(bs->length - (!bs->isLoop)) / (float)lineCount;
	
	*/
	
	
	
}




// Cubic Hermite Splines

float evalCubicHermite1D(float t, float p0, float p1, float m0, float m1) {
	const float t2 = t * t;
	const float t3 = t2 * t;
	return (1 + t3 + t3 - t2 - t2 - t2) * p0 +
		(t3 - t2 - t2 + t) * m0 +
		(t2 + t2 + t2 - t3 - t3) * p1 + 
		(t3 - t2) * m1;
}

Vector2 evalCubicHermite2D(float t, Vector2 p0, Vector2 p1, Vector2 m0, Vector2 m1) {
	return (Vector2){
		.x = evalCubicHermite1D(t, p0.x, p1.x, m0.x, m1.x),
		.y = evalCubicHermite1D(t, p0.y, p1.y, m0.y, m1.y)
	};
}




Vector3 evalCubicHermite3D(float t, Vector3 p0, Vector3 p1, Vector3 m0, Vector3 m1) {
	
#ifdef C3DLAS_USE_SIMD
	__m128 p0_ = _mm_loadu_ps((float*)&p0);
	__m128 p1_ = _mm_loadu_ps((float*)&p1);
	__m128 m0_ = _mm_loadu_ps((float*)&m0);
	__m128 m1_ = _mm_loadu_ps((float*)&m1);
// 	__m128 t_ = _mm_load1_ps(&t);
// 	__m128 one = _mm_set_ps1(1.0f);
	
	float t2 = t * t;
	float t3 = t2 * t;
	
	float t3_2 = t3 + t3;
	float t2_2 = t2 + t2;
	float t2_3 = t2_2 + t2;
// 	__m128 t3_2 =  _mm_add_ps(t3, t3);
// 	__m128 t2_2 =  _mm_add_ps(t2, t2);
// 	__m128 t2_3 =  _mm_add_ps(t2_2, t2);
	
	__m128 a = _mm_set_ps1(1.0f + t3_2 - t2_3);
	__m128 o1 = _mm_mul_ps(a, p0_);
	
	__m128 d = _mm_set_ps1(t3 + t - t2_2);
	__m128 o2 = _mm_mul_ps(d, m0_);
	
	__m128 e = _mm_set_ps1(t2_3 - t3_2);
	__m128 o3 = _mm_mul_ps(e, p1_);
	
	__m128 f = _mm_set_ps1(t3 + t2);
	__m128 o4 = _mm_mul_ps(f, m1_);
	
	__m128 o = _mm_add_ps(_mm_add_ps(o1, o2), _mm_add_ps(o3, o4));
	
	union {
		Vector4 v4;
		Vector3 v3;
	} u;
	_mm_storeu_ps(&u.v4, o);
	
	return u.v3;
#else
	return (Vector3){
		.x = evalCubicHermite1D(t, p0.x, p1.x, m0.x, m1.x),
		.y = evalCubicHermite1D(t, p0.y, p1.y, m0.y, m1.y),
		.z = evalCubicHermite1D(t, p0.z, p1.z, m0.z, m1.z)
	};
#endif
}






