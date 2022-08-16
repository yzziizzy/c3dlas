


#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <x86intrin.h>

#include "c3dlas.h"



#ifndef _GNU_SOURCE
static inline sincosf(float x, float* s, float* c) {
	*s = sinf(x);
	*c = cosf(x);
}
#endif


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



float frandPCG(float low, float high, PCG* pcg) {
	return low + ((high - low) * (pcg_f(&pcg->state, pcg->stream) * 0.5 + 0.5));
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


Vector2 vClamp2f(Vector2 in, float min, float max) {
	return (Vector2) {
		.x = fmaxf(min, fminf(in.x, max)), 
		.y = fmaxf(min, fminf(in.y, max)) 
	};
}

Vector2 vClamp2(Vector2 in, Vector2 min, Vector2 max) {
	return (Vector2) {
		.x = fmaxf(min.x, fminf(in.x, max.x)), 
		.y = fmaxf(min.y, fminf(in.y, max.y)) 
	};
}

Vector3 vClamp3f(Vector3 in, float min, float max) {
	return (Vector3) {
		.x = fmaxf(min, fminf(in.x, max)), 
		.y = fmaxf(min, fminf(in.y, max)), 
		.z = fmaxf(min, fminf(in.z, max)) 
	};
}

Vector3 vClamp3(Vector3 in, Vector3 min, Vector3 max) {
	return (Vector3) {
		.x = fmaxf(min.x, fminf(in.x, max.x)), 
		.y = fmaxf(min.y, fminf(in.y, max.y)), 
		.z = fmaxf(min.z, fminf(in.z, max.z)) 
	};
}


// Coordinate system conversions

// Does not check for degenerate vectors
// Cartesian to Spherical
Vector3 vC2S3(Vector3 cart) {
	Vector3 sp;
	sp.rho = vMag3(cart);
	sp.theta = atan2f(cart.x, cart.y);
	sp.phi = acosf(cart.z / sp.rho);
	
	return sp;
}

// Spherical to Cartesian
Vector3 vS2C3(Vector3 s) {
	float st, ct, sp, cp;
	
	// as of July 2022, gcc trunk is smart enough to automatically optimize to sincos, but clang isn't. 
	sincosf(s.phi, &sp, &cp);
	sincosf(s.theta, &st, &ct);

	return (Vector3){
		.x = s.rho * sp * ct,
		.y = s.rho * sp * st,
		.z = s.rho * cp
	};
}



// Muchas gracias, Inigo.  
// https://iquilezles.org/articles/distfunctions2d/
float vDistPointLine2(Vector2 p, Line2 ls) {
	Vector2 pa = vSub2(p, ls.start); // vector from the starting point to p
	Vector2 ba = vSub2(ls.end, ls.start); // vector from the starting point to the ending point
	
	float t = vDot2(pa, ba) / vDot2(ba, ba); // project the pa onto ba, then divide that distance by the length of ba to normalize it
	fclamp(t, 0.0, 1.0); // clamp t to between the endpoints of the line segment
	
	// Consider the starting point to be at the origin, for ease of visualization.
	// ba is the vector from the origin to the endpoint og the line that now passes through the origin.
	// Scaling ba by t gives the intercept point of the line through p that is perpendicular to the test line segment.
	// pa is p if a was the origin. Therefore, pi is the vector from p to the intercept point on the test line segment. 
	Vector2 pi = vSub2(pa, vScale2(ba, t));
	return vMag2(pi); // the answer is the length of pi 
}

float vDistPointLine3(Vector3 p, Line3 ls) {
	Vector3 pa = vSub3(p, ls.start);
	Vector3 ba = vSub3(ls.end, ls.start);
	
	float t = fclamp(vDot3(pa, ba) / vDot3(ba, ba), 0.0, 1.0);
	return vMag3(vSub3(pa, vScale3(ba, t)));
}

// This version also returns the normalized distance along the line of the closest point
float vDistTPointLine2(Vector2 p, Line2 ls, float* T) {
	Vector2 pa = vSub2(p, ls.start);
	Vector2 ba = vSub2(ls.end, ls.start);
	
	float t = fclamp(vDot2(pa, ba) / vDot2(ba, ba), 0.0, 1.0);
	if(T) *T = t;
	return vMag2(vSub2(pa, vScale2(ba, t)));
}

float vDistTPointLine3(Vector3 p, Line3 ls, float* T) {
	Vector3 pa = vSub3(p, ls.start);
	Vector3 ba = vSub3(ls.end, ls.start);
	
	float t = fclamp(vDot3(pa, ba) / vDot3(ba, ba), 0.0, 1.0);
	if(T) *T = t;
	return vMag3(vSub3(pa, vScale3(ba, t)));
}

// ----


int vInsidePolygon(Vector2 p, Polygon* poly) {
	int inside = 0;
	int cnt = poly->pointCount;
	
	if(poly->maxRadiusSq < vDot2(poly->centroid, p)) return 0;
	
	for(int i = 0; i < cnt; i++) {
		Vector2 a = poly->points[i];
		Vector2 b = poly->points[(i + 1) % cnt];
		
		if(a.y == b.y) continue; // horizontal edges are ignored
		
		// we're testing a ray going to the right
		if(a.x < p.x && b.x < p.x) continue; // segment is entirely left of the point
		
		if(a.y >= p.y && b.y >= p.y) continue; // segment entirely above the point
		if(a.y < p.y && b.y < p.y) continue; // segment entirely below the point
		// segment is in the same vertical band as the point
		
		float sx = a.x + (b.x - a.x) * ((p.y - a.y) / (b.y - a.y));
		if(p.x > sx) continue;
		
		inside = !inside;
	}

	return inside;
}


// Muchas gracias, Inigo.  
// https://iquilezles.org/articles/distfunctions2d/
float vDistPolygon(Vector2 p, Polygon* poly) {

	float d = vDot2(vSub2(p, poly->points[0]), vSub2(p, poly->points[0]));
	float s = 1.0;

	for(int i = 0, j = poly->pointCount - 1; i < poly->pointCount; j = i, i++) {
		Vector2 A = poly->points[i];
		Vector2 B = poly->points[j];

		Vector2 e = vSub2(B, A);
		Vector2 w = vSub2(p, A);
		Vector2 b = vSub2(w, vScale2(e, fclamp(vDot2(w, e) / vDot2(e, e), 0.0, 1.0)));

		d = fminf(d, vDot2(b, b));
        
		int c1 = p.y >= A.y;
		int c2 = p.y < B.y;
		int c3 = e.x * w.y > e.y * w.x;
		if((c1 && c2 && c3) || (!c1 && !c2 && !c3)) s *= -1.0;  
	}

	return s * sqrtf(d);
}

// ----


void polyCalcCentroid(Polygon* poly) {
	int cnt = poly->pointCount;
	Vector2 centroid = {0,0};
	
	for(int i = 0; i < cnt; i++) {
		Vector2 a = poly->points[i];
		centroid = vAdd2(centroid, a);
	}
	
	poly->centroid = vScale2(centroid, 1.0 / poly->pointCount);

}


void polyCalcRadiusSq(Polygon* poly) {
	int cnt = poly->pointCount;
	float d = 0;
	
	for(int i = 0; i < cnt; i++) {
		Vector2 a = poly->points[i];
		d = fmaxf(d, vDot2(poly->centroid, a));
	}
	
	poly->maxRadiusSq = d;
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



void vRandomPCG3p(Vector3* end1, Vector3* end2, PCG* pcg, Vector3* out) {
	out->x = frandPCG(fminf(end1->x, end2->x), fmaxf(end1->x, end2->x), pcg);
	out->y = frandPCG(fminf(end1->y, end2->y), fmaxf(end1->y, end2->y), pcg);
	out->z = frandPCG(fminf(end1->z, end2->z), fmaxf(end1->z, end2->z), pcg);
}

Vector3 vRandomPCG3(Vector3 end1, Vector3 end2, PCG* pcg) {
	Vector3 o;
	vRandomPCG3p(&end1, &end2, pcg, &o);
	return o;
}

void vRandom3p(Vector3* end1, Vector3* end2, Vector3* out) {
	out->x = frand(fminf(end1->x, end2->x), fmaxf(end1->x, end2->x));
	out->y = frand(fminf(end1->y, end2->y), fmaxf(end1->y, end2->y));
	out->z = frand(fminf(end1->z, end2->z), fmaxf(end1->z, end2->z));
}

Vector3 vRandom3(Vector3 end1, Vector3 end2) {
	return (Vector3){
		.x = frand(fminf(end1.x, end2.x), fmaxf(end1.x, end2.x)),
		.y = frand(fminf(end1.y, end2.y), fmaxf(end1.y, end2.y)),
		.z = frand(fminf(end1.z, end2.z), fmaxf(end1.z, end2.z))
	};
}

// Uniformly distributed around the unit sphere; ie, no clustering at the poles.
Vector3 vRandomNorm3() {
	Vector3 out;
	vRandomNorm3p(&out);
	return out;
}

void vRandomNorm3p(Vector3* out) {
	float u = frand(-1.0, 1.0);
	float th = frand(0, 2.0 * F_PI);
	float q = sqrtf(1.0 - u * u);
	float sth, cth;
	
	sincosf(th, &sth, &cth);
	out->x = u * cth;
	out->y = u * sth;
	out->z = u;
}



// reflects the distance from v to pivot across pivot.
// out, pivot, and v will all be in the same plane, with pivot half way between v and out
void vReflectAcross3p(Vector3* v, Vector3* pivot, Vector3* out) {
	Vector3 v2 = vScale3(*v, -1);
	float d = vDot3(v2, *pivot) * 2.0;
	*out = vSub3(v2, vScale3(*pivot, d));
}

Vector3 vReflectAcross3(Vector3 v, Vector3 pivot) {
	Vector3 o;
	vReflectAcross3p(&v, &pivot, &o);
	return o;
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

Vector3 vMatrixMulProjectedMagic3(Vector3 in, Matrix* m) {
	Vector4 v;

	v.x = in.x * m->m[0+0] + in.y * m->m[4+0] + in.z * m->m[8+0] + 1.0 * m->m[12+0];
	v.y = in.x * m->m[0+1] + in.y * m->m[4+1] + in.z * m->m[8+1] + 1.0 * m->m[12+1];
	v.z = in.x * m->m[0+2] + in.y * m->m[4+2] + in.z * m->m[8+2] + 1.0 * m->m[12+2];
	v.w = in.x * m->m[0+3] + in.y * m->m[4+3] + in.z * m->m[8+3] + 1.0 * m->m[12+3];
	
	if(v.w == 0) return (Vector3){0,0,0}; 
	if(v.w < 0) v.w = -v.w;
	
	return (Vector3){x: v.x / v.w, y: v.y / v.w, z: v.z / v.w};
}

Vector3 vMatrixMul3(Vector3 in, Matrix* m) {
	Vector4 v;

	v.x = in.x * m->m[0+0] + in.y * m->m[4+0] + in.z * m->m[8+0] + 1.0 * m->m[12+0];
	v.y = in.x * m->m[0+1] + in.y * m->m[4+1] + in.z * m->m[8+1] + 1.0 * m->m[12+1];
	v.z = in.x * m->m[0+2] + in.y * m->m[4+2] + in.z * m->m[8+2] + 1.0 * m->m[12+2];
	v.w = in.x * m->m[0+3] + in.y * m->m[4+3] + in.z * m->m[8+3] + 1.0 * m->m[12+3];
	
	if(v.w == 0) return (Vector3){0,0,0};
	
	return (Vector3){x: v.x / v.w, y: v.y / v.w, z: v.z / v.w};
}

Vector4 vMatrixMul4(Vector4 in, Matrix* m) {
	Vector4 v;

	v.x = in.x * m->m[0+0] + in.y * m->m[4+0] + in.z * m->m[8+0] + in.w * m->m[12+0];
	v.y = in.x * m->m[0+1] + in.y * m->m[4+1] + in.z * m->m[8+1] + in.w * m->m[12+1];
	v.z = in.x * m->m[0+2] + in.y * m->m[4+2] + in.z * m->m[8+2] + in.w * m->m[12+2];
	v.w = in.x * m->m[0+3] + in.y * m->m[4+3] + in.z * m->m[8+3] + in.w * m->m[12+3];
	
	return v;
}

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
	out->x = (b->max.x + b->min.x) / 2.0;
	out->y = (b->max.y + b->min.y) / 2.0;
	out->z = (b->max.z + b->min.z) / 2.0;
}

Vector3 boxCenter3(const AABB3 b) {
	return (Vector3) {
		(b.max.x + b.min.x) / 2.0,
		(b.max.y + b.min.y) / 2.0,
		(b.max.z + b.min.z) / 2.0
	};
}


void boxSize3p(const AABB3* b, Vector3* out) {
	out->x = b->max.x - b->min.x;
	out->y = b->max.y - b->min.y;
	out->z = b->max.z - b->min.z;
}

Vector3 boxSize3(const AABB3 b) {
	return (Vector3){
		b.max.x - b.min.x,
		b.max.y - b.min.y,
		b.max.z - b.min.z
	};
}

void boxExpandTo3p(AABB3* b, Vector3* p) {
	b->min.x = fminf(b->min.x, p->x);
	b->min.y = fminf(b->min.y, p->y);
	b->min.z = fminf(b->min.z, p->z);
	b->max.x = fmaxf(b->max.x, p->x);
	b->max.y = fmaxf(b->max.y, p->y);
	b->max.z = fmaxf(b->max.z, p->z);
}

void boxExpandTo3(AABB3* b, Vector3 p) {
	boxExpandTo3p(b, &p);
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
	out->x = (b->max.x + b->min.x) / 2.0;
	out->y = (b->max.y + b->min.y) / 2.0;
}

Vector2 boxSize2(const AABB2 b) {
	return (Vector2){
		b.max.x - b.min.x,
		b.max.y - b.min.y
	};
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
	tmin = fminf(t1.x, t2.x);
	tmax = fmaxf(t1.x, t2.x);
	
	t1.y = (b->min.y - r->o.y) * id.y;
	t2.y = (b->max.y - r->o.y) * id.y;
	tmin = fmaxf(tmin, fminf(t1.y, t2.y));
	tmax = fminf(tmax, fmaxf(t1.y, t2.y));
	
	t1.z = (b->min.z - r->o.z) * id.z;
	t2.z = (b->max.z - r->o.z) * id.z;
	tmin = fmaxf(tmin, fminf(t1.z, t2.z));
	tmax = fminf(tmax, fmaxf(t1.z, t2.z));

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
	tmin = fminf(t1.x, t2.x);
	tmax = fmaxf(t1.x, t2.x);
	
	t1.y = (b->min.y - r->o.y) * id.y;
	t2.y = (b->max.y - r->o.y) * id.y;
//	tmin = fmaxf(tmin, fminf(t1.y, t2.y));
	tmax = fminf(tmax, fmaxf(t1.y, t2.y));
	
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
	tmin = fminf(t1.x, t2.x);
	tmax = fmaxf(t1.x, t2.x);
	
	t1.y = (b->min.y - r->o.y) * id.y;
	t2.y = (b->max.y - r->o.y) * id.y;
	tmin = fmaxf(tmin, fminf(t1.y, t2.y));
	tmax = fminf(tmax, fmaxf(t1.y, t2.y));
	
	t1.z = (b->min.z - r->o.z) * id.z;
	t2.z = (b->max.z - r->o.z) * id.z;
	tmin = fmaxf(tmin, fminf(t1.z, t2.z));
	tmax = fminf(tmax, fmaxf(t1.z, t2.z));
	
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




#include "matrix4.c"
#include "quaternion.c"



