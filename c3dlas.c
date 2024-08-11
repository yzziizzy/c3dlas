


#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include <x86intrin.h>

#include "c3dlas.h"

#ifdef C3DLAS_USE_BUILTINS

	#define abs_double __builtin_fabs
	#define abs_float  __builtin_fabsf
#else
	#define abs_double fabs
	#define abs_float  fabsf
#endif


#ifdef C3DLAS_NO_TGMATH
	// requires GCC probably
	
		/*
	#define FD_CHOOSE_1(a, b, fn_f, fn_d)\
		__builtin_choose_expr( \
			__builtin_types_compatible_p(__typeof__(a), double), \
				fn_d(a), \
				fn_f(a))
		
	#define FD_CHOOSE_2(a, b, fn_f, fn_d)\
		__builtin_choose_expr( \
			__builtin_types_compatible_p(__typeof__(a), double) || __builtin_types_compatible_p(__typeof__(b), double), \
				fn_d(a, b), \
				fn_f(a, b))
	
	#define fmax(a,b) FD_CHOOSE_2(a, b, fmaxf, fmax)
	#define fmin(a,b) FD_CHOOSE_2(a, b, fminf, fmin)
	#define fabs(a)   FD_CHOOSE_1(a, fabsf, fabs)
	#define sqrt(a)   FD_CHOOSE_1(a, sqrtf, sqrt)
		*/
#else
	#include <tgmath.h>
#endif



#ifndef _GNU_SOURCE
static inline void sincosf(float x, float* s, float* c) {
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
// Thanks to Kaslai (https://github.com/Aslai) for fixing a nasty bug in the previous version
float pcg_f(uint64_t* state, uint64_t stream) {
	union {
		uint32_t fu;
		float ff;
	} u;
	
	
	uint64_t last = *state;
	*state = (last * 6364136223846793005ULL) + (stream | 1);
	uint32_t xs = ((last >> 18) ^ last) >> 27;
	uint32_t rot = last >> 59;
	uint32_t fin = (xs >> rot) | (xs << ((-rot) & 31));
	
	uint32_t exp = (fin & 0x3F800000);
	exp = (0x7F + 33 - __builtin_clzl(exp)) << 23;
	
	u.fu = ((fin) & 0x807fffff) | exp;
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






#define FN(sz, suf, ty, ft, sufft, pref, ...) \
\
int vEqExact##suf(const Vector##suf a, const Vector##suf b) { \
	return vEqExact##suf##p(&a, &b); \
} \
int vEqExact##suf##p(const Vector##suf const * a, const Vector##suf const * b) { \
	int tmp = 0; \
	for(int i = 0; i < sz; i++) \
		tmp |= ((ty*)a)[i] == ((ty*)b)[i]; \
	return tmp; \
} \
\
int vEq##suf(const Vector##suf a, const Vector##suf b) { \
	return vEqEp##suf(a, b, pref##_CMP_EPSILON); \
} \
int vEq##suf##p(const Vector##suf* a, const Vector##suf* b) { \
	return vEqEp##suf(*a, *b, pref##_CMP_EPSILON); \
} \
\
int vEqEp##suf(const Vector##suf a, const Vector##suf b, ft epsilon) { \
	return vEqEp##suf##p(&a, &b, epsilon); \
} \
int vEqEp##suf##p(const Vector##suf* a, const Vector##suf* b, ft epsilon) { \
	return vDistSq##suf(*a, *b) <= epsilon * epsilon; \
} \
\
ft vDistSq##suf(const Vector##suf a, const Vector##suf b) { \
	return vDistSq##suf##p(&a, &b); \
} \
ft vDistSq##suf##p(const Vector##suf* a, const Vector##suf* b) { \
	ft tmp = 0; \
	for(int i = 0; i < sz; i++) { \
		ft q = ((ty*)a)[i] - ((ty*)b)[i]; \
		tmp += q * q; \
	} \
	return tmp;\
} \
ft vDist##suf(const Vector##suf a, const Vector##suf b) { \
	return sqrt(vDistSq##suf##p(&a, &b)); \
} \
ft vDist##suf##p(const Vector##suf* a, const Vector##suf* b) { \
	return sqrt(vDistSq##suf##p(a, b)); \
} \
\
Vector##suf vAdd##suf(const Vector##suf a, const Vector##suf b) { \
	Vector##suf out; \
	vAdd##suf##p(&a, &b, &out); \
	return out; \
} \
void vAdd##suf##p(const Vector##suf* a, const Vector##suf* b, Vector##suf* out) { \
	for(int i = 0; i < sz; i++) { \
		((ty*)out)[i] = ((ty*)a)[i] + ((ty*)b)[i]; \
	} \
} \
\
Vector##suf vSub##suf(const Vector##suf a, const Vector##suf b) { \
	Vector##suf out; \
	vSub##suf##p(&a, &b, &out); \
	return out; \
} \
void vSub##suf##p(const Vector##suf const * a, const Vector##suf const * b, Vector##suf* out) { \
	for(int i = 0; i < sz; i++) { \
		((ty*)out)[i] = ((ty*)a)[i] - ((ty*)b)[i]; \
	} \
} \
\
Vector##suf vMul##suf(const Vector##suf a, const Vector##suf b) { \
	Vector##suf out; \
	vMul##suf##p(&a, &b, &out); \
	return out; \
} \
void vMul##suf##p(const Vector##suf const * a, const Vector##suf const * b, Vector##suf* out) { \
	for(int i = 0; i < sz; i++) { \
		((ty*)out)[i] = ((ty*)a)[i] * ((ty*)b)[i]; \
	} \
} \
\
ft vDot##suf(const Vector##suf a, const Vector##suf b) { \
	return vDot##suf##p(&a, &b); \
} \
ft vDot##suf##p(const Vector##suf* a, const Vector##suf* b) { \
	ft tmp = 0; \
	for(int i = 0; i < sz; i++) { \
		tmp += ((ty*)a)[i] * ((ty*)b)[i]; \
	} \
	return tmp;\
} \
\
Vector##sufft vScale##suf(const Vector##suf v, ft scalar) { \
	Vector##sufft out; \
	vScale##suf##p(&v, scalar, &out); \
	return out; \
} \
void vScale##suf##p(const Vector##suf* v, ft scalar, Vector##sufft* out) { \
	for(int i = 0; i < sz; i++) \
		((ft*)out)[i] = (ft)((ty*)v)[i] * scalar; \
} \
\
\
Vector##sufft vAvg##suf(const Vector##suf a, const Vector##suf b) { \
	Vector##sufft out; \
	vAvg##suf##p(&a, &b, &out); \
	return out; \
} \
void vAvg##suf##p(const Vector##suf* a, const Vector##suf* b, Vector##sufft* out) { \
	for(int i = 0; i < sz; i++) { \
		((ty*)out)[i] = (((ty*)a)[i] + ((ty*)b)[i]) / (ft)2.0; \
	} \
} \
\
Vector##suf vNeg##suf(const Vector##suf v) { \
	Vector##suf out; \
	vNeg##suf##p(&v, &out); \
	return out; \
} \
void vNeg##suf##p(const Vector##suf* v, Vector##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((ty*)out)[i] = -((ty*)v)[i]; \
} \
\
Vector##sufft vLerp##suf(const Vector##suf a, const Vector##suf b, ft t) { \
	Vector##sufft out; \
	vLerp##suf##p(&a, &b, t, &out); \
	return out; \
} \
void vLerp##suf##p(const Vector##suf* a, const Vector##suf* b, ft t, Vector##sufft* out) { \
	for(int i = 0; i < sz; i++) \
		((ft*)out)[i] = (ft)((ty*)a)[i] + ((ft)(((ty*)b)[i] - ((ty*)a)[i]) * t) ; \
} \
\
Vector##sufft vInv##suf(const Vector##suf v) { \
	Vector##sufft out; \
	vInv##suf##p(&v, &out); \
	return out; \
} \
void vInv##suf##p(const Vector##suf* v, Vector##sufft* out) { \
	for(int i = 0; i < sz; i++) \
		((ft*)out)[i] = (((ty*)v)[i] == 0) ? pref##_MAX : ((ft)1.0 /  (ft)((ty*)v)[i]); \
} \
\
ft vLen##suf(const Vector##suf v) { \
	return vLen##suf##p(&v); \
} \
ft vLen##suf##p(const Vector##suf* v) { \
	ft tmp = 0.0; \
	for(int i = 0; i < sz; i++) \
		tmp += (ft)((ty*)v)[i] * (ft)((ty*)v)[i]; \
	return sqrt(tmp); \
} \
\
ft vLenSq##suf(const Vector##suf v) { \
	return vLenSq##suf##p(&v); \
} \
ft vLenSq##suf##p(const Vector##suf* v) { \
	return vDot##suf##p(v, v); \
} \
\
ft vMag##suf(const Vector##suf v) { \
	return vLen##suf##p(&v); \
} \
ft vMag##suf##p(const Vector##suf* v) { \
	return vLen##suf##p(v); \
} \
\
ft vInvLen##suf(const Vector##suf v) { \
	ft tmp = vLen##suf(v); \
	return tmp == 0 ? pref##_MAX : (ft)1.0 / tmp; \
} \
ft vInvLen##suf##p(const Vector##suf* v) { \
	return vInvLen##suf(*v); \
} \
\
Vector##sufft vNorm##suf(const Vector##suf v) { \
	Vector##sufft out; \
	vNorm##suf##p(&v, &out); \
	return out; \
} \
void vNorm##suf##p(const Vector##suf* v, Vector##sufft* out) { \
	ft n = vLenSq##suf(*v); \
	\
	if(n >= (ft)1.0f - pref##_CMP_EPSILON && n <= (ft)1.0 + pref##_CMP_EPSILON) { \
		for(int i = 0; i < sz; i++) \
			((ft*)out)[i] = (ft)((ty*)v)[i]; \
		return;		 \
	} \
	else if(n == 0.0) { \
		for(int i = 0; i < sz; i++) \
			((ft*)out)[i] = 0; \
		return; \
	} \
	 \
	n = (ft)1.0 / sqrt(n); \
	for(int i = 0; i < sz; i++) \
		((ft*)out)[i] = (ft)((ty*)v)[i] * n; \
}  \
\
Vector##sufft vUnit##suf(const Vector##suf v) { \
	return vNorm##suf(v); \
} \
void vUnit##suf##p(const Vector##suf* v, Vector##sufft* out) { \
	return vNorm##suf##p(v, out); \
} \


	C3DLAS_VECTOR_LIST(FN)
#undef FN


 



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





// scalar muliplication



// Dot product (inner product)



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






// Vector Inverse. Returns FLT_MAX on div/0


// Vector magnitude (length)


// Squared distance from one point to another



// Distance from one point to another




// Vector normalize (scale to unit length)




// vMin(a, b)  Returns the minimum values of each component
// vMin(a, b)  Returns the maximum values of each component

#define FN(sz, suf, t, maxval) \
void vMin##sz##suf##p(const Vector##sz##suf* a, const Vector##sz##suf* b, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = fmin(((t*)a)[i], ((t*)b)[i]); \
} \
void vMax##sz##suf##p(const Vector##sz##suf* a, const Vector##sz##suf* b, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = fmax(((t*)a)[i], ((t*)b)[i]); \
} \
Vector##sz##suf vMin##sz##suf(Vector##sz##suf a, Vector##sz##suf b) { \
	Vector##sz##suf out; \
	vMin##sz##suf##p(&a, &b, &out); \
	return out; \
} \
Vector##sz##suf vMax##sz##suf(Vector##sz##suf a, Vector##sz##suf b) { \
	Vector##sz##suf out; \
	vMax##sz##suf##p(&a, &b, &out); \
	return out; \
} \
\
int vMinComp##sz##suf##p(const Vector##sz##suf* a) { \
	t best = ((t*)a)[0]; \
	int best_ind = 0; \
	for(int i = 1; i < sz; i++) { \
		if(((t*)a)[i] < best) { \
			best = ((t*)a)[i]; \
			best_ind = i; \
		} \
	} \
	return best_ind; \
} \
\
int vMaxComp##sz##suf##p(const Vector##sz##suf* a) { \
	t best = ((t*)a)[0]; \
	int best_ind = 0; \
	for(int i = 1; i < sz; i++) { \
		if(((t*)a)[i] > best) { \
			best = ((t*)a)[i]; \
			best_ind = i; \
		} \
	} \
	return best_ind; \
} \
\
int vMinComp##sz##suf(const Vector##sz##suf a) { \
	return vMinComp##sz##suf##p(&a); \
} \
\
int vMaxComp##sz##suf(const Vector##sz##suf a) { \
	return vMaxComp##sz##suf##p(&a); \
} \
\
Vector##sz##suf vClamp##sz##suf(Vector##sz##suf in, Vector##sz##suf min, Vector##sz##suf max) { \
	Vector##sz##suf out; \
	for(int i = 0; i < sz; i++) \
		((t*)&out)[i] = fmax(((t*)&min)[i], fmin(((t*)&in)[i], ((t*)&max)[i])); \
	return out; \
} \
Vector##sz##suf vAbs##sz##suf(const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vAbs##sz##suf##p(&v, &out); \
	return out; \
} \
void vAbs##sz##suf##p(const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = abs_##t( ((t*)v)[i] ); \
} \
Vector##sz##suf vSign##sz##suf(const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vSign##sz##suf##p(&v, &out); \
	return out; \
} \
void vSign##sz##suf##p(const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = copysign((t)1.0, ((t*)v)[i] ); \
} \
Vector##sz##suf vStep##sz##suf(const Vector##sz##suf edge, const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vStep##sz##suf##p(&edge, &v, &out); \
	return out; \
} \
void vStep##sz##suf##p(const Vector##sz##suf* edge, const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = ((t*)v)[i] < ((t*)edge)[i] ? 0.0 : 1.0; \
} \


FN(2,  , float,  FLT_MAX)
FN(3,  , float,  FLT_MAX)
FN(4,  , float,  FLT_MAX)
FN(2, d, double, DBL_MAX)
FN(3, d, double, DBL_MAX)
FN(4, d, double, DBL_MAX)
#undef FN

#define FN(sz, suf, t, maxval) \
void vMin##sz##suf##p(const Vector##sz##suf* a, const Vector##sz##suf* b, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = MIN(((t*)a)[i], ((t*)b)[i]); \
} \
void vMax##sz##suf##p(const Vector##sz##suf* a, const Vector##sz##suf* b, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = MAX(((t*)a)[i], ((t*)b)[i]); \
} \
Vector##sz##suf vMin##sz##suf(Vector##sz##suf a, Vector##sz##suf b) { \
	Vector##sz##suf out; \
	vMin##sz##suf##p(&a, &b, &out); \
	return out; \
} \
Vector##sz##suf vMax##sz##suf(Vector##sz##suf a, Vector##sz##suf b) { \
	Vector##sz##suf out; \
	vMax##sz##suf##p(&a, &b, &out); \
	return out; \
} \
Vector##sz##suf vClamp##sz##suf(Vector##sz##suf in, Vector##sz##suf min, Vector##sz##suf max) { \
	Vector##sz##suf out; \
	for(int i = 0; i < sz; i++) \
		((t*)&out)[i] = MAX(((t*)&min)[i], MIN(((t*)&in)[i], ((t*)&max)[i])); \
	return out; \
} \
Vector##sz##suf vAbs##sz##suf(const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vAbs##sz##suf##p(&v, &out); \
	return out; \
} \
void vAbs##sz##suf##p(const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = labs( ((t*)v)[i] ); \
} \
Vector##sz##suf vSign##sz##suf(const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vSign##sz##suf##p(&v, &out); \
	return out; \
} \
void vSign##sz##suf##p(const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = ((t*)v)[i] < 0 ? -1 : 1; \
} \
Vector##sz##suf vStep##sz##suf(const Vector##sz##suf edge, const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vStep##sz##suf##p(&edge, &v, &out); \
	return out; \
} \
void vStep##sz##suf##p(const Vector##sz##suf* edge, const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = ((t*)v)[i] < ((t*)edge)[i] ? 0.0 : 1.0; \
} \


FN(2, i, int,  DBL_MAX)
FN(3, i, int,  DBL_MAX)
FN(4, i, int,  DBL_MAX)
FN(2, l, long, DBL_MAX)
FN(3, l, long, DBL_MAX)
FN(4, l, long, DBL_MAX)
#undef FN



// Returns an arbitrary unit vector perpendicular to the input
// The input vector does not need to be normalized
void vPerp3p(Vector3* n, Vector3* out) {
	*out = vPerp3(*n);
}

Vector3 vPerp3(Vector3 n) {
	float f, d;
	
	float absx = fabs(n.x);
	float absy = fabs(n.y);
	
	if(absx < absy) {
		if(n.x < n.z) goto MIN_X;
		goto MIN_Z;
	}
	if(absy < fabs(n.z)) goto MIN_Y;
	goto MIN_Z;
	
MIN_X: 
	d = 1.0f / sqrtf(n.z * n.z + n.y * n.y); 
	f = n.z;
	n.z = n.y * d;
	n.y = -f * d;
	n.x = 0;
	return n;
	
MIN_Y: 
	d = 1.0f / sqrtf(n.z * n.z + n.x * n.x); 
	f = n.x;
	n.x = n.z * d;
	n.z = -f * d;
	n.y = 0;
	return n;
	
MIN_Z: 
	d = 1.0f / sqrtf(n.x * n.x + n.y * n.y); 
	f = n.y;
	n.y = n.x * d;
	n.x = -f * d;
	n.z = 0;
	return n;
}


// Returns an arbitrary unit vector perpendicular to the input
// The input vector does not need to be normalized
void vPerp2p(Vector2* n, Vector2* out) {
	*out = vPerp2(*n);
}

Vector2 vPerp2(Vector2 n) {
	return vNorm2((Vector2){.x = -n.y, .y = n.x});
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

Vector3 closestPointToRay3(Vector3 p, Ray3 r) {
	Vector3 po = vSub3(p, r.o); // vector from the starting point to p
	
	float t = vDot3(po, r.d); // project the pa onto the ray direction
	fclamp(t, 0.0, 1.0); // clamp t to between the endpoints of the line segment

	return vSub3(po, vScale3(r.d, t));
}


// completely untested.
// can probably be optimized
// This function is poorly named. It is designed to check if a bounding sphere intersects a cone surrounding a viewing frustum.
int distanceSphereToCone(Vector3 spc, float spr, Vector3 c1, Vector3 c2, float cr1, float cr2) {
	
	Vector3 cnorm = vNorm(vSub(c2, c1)); // normal pointing down the center of the cone
	
	Vector3 sp_c1 = vSub(spc, c1); // vector pointing from c1 to the sphere center
	Vector3 up = vCross3(spc, cnorm); // vector perpendicular to the plane containing the cone's centerline and the sphere center.
	Vector3 perp_toward_sp = vNorm(vCross3(cnorm, up)); // vector perpendicular to the cone's centerline within the plane, towards the sphere
	Vector3 outer_c1 = vAdd(c1, vScale(perp_toward_sp, cr1)); // point in the plane on the outer edge of the cone
	Vector3 outer_c2 = vAdd(c2, vScale(perp_toward_sp, cr2)); // point in the plane on the outer edge of the cone
	Vector3 closest = closestPointToRay3(spc, (Ray3){.o = outer_c1, .d = vNorm(vSub(outer_c2, outer_c1))}); // point on the cone closest to the sphere
	
	// this part is probably wrong
	if(vDot(perp_toward_sp, vSub(spc, closest)) < 0) return 1; // is the sphere center inside the cone?

	return (vDist(closest, spc) - spr) <= 0;
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


float projPointLine2(Vector2 p, Line2 ls) {
	Vector2 pa = vSub2(p, ls.start);
	Vector2 ba = vSub2(ls.end, ls.start);
	
	return vDot2(pa, ba) / vDot2(ba, ba);
}


float distTPointRay3(Vector3 p, Ray3 r, float* T) {
	Vector3 pa = vSub3(p, r.o);
	Vector3 ba = vNeg3(r.d);// vSub3(ls.end, ls.start);
	
	float t = vDot3(pa, ba) / vDot3(ba, ba);
	if(T) *T = t;
	return vLen3(vSub3(pa, vScale3(ba, t)));
}

float dist2TPointRay3(Vector3 p, Ray3 r, float* T) {
	Vector3 pa = vSub3(p, r.o);
	Vector3 ba = vNeg3(r.d);// vSub3(ls.end, ls.start);
	
	float t = vDot3(pa, ba) / vDot3(ba, ba);
	if(T) *T = t;
	return vLenSq3(vSub3(pa, vScale3(ba, t)));
}


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



void vRandomPCG2p(Vector2* end1, Vector2* end2, PCG* pcg, Vector2* out) {
	out->x = frandPCG(fminf(end1->x, end2->x), fmaxf(end1->x, end2->x), pcg);
	out->y = frandPCG(fminf(end1->y, end2->y), fmaxf(end1->y, end2->y), pcg);
}

Vector2 vRandomPCG2(Vector2 end1, Vector2 end2, PCG* pcg) {
	Vector2 o;
	vRandomPCG2p(&end1, &end2, pcg, &o);
	return o;
}

void vRandomNormPCG2p(PCG* pcg, Vector2* out) {
	float th = frandPCG(0, 2.0 * F_PI, pcg);
	float sth, cth;
	
	sincosf(th, &sth, &cth);
	out->x = cth;
	out->y = sth;
}

Vector2 vRandomNormPCG2(PCG* pcg) {
	Vector2 o;
	vRandomNormPCG2p(pcg, &o);
	return o;
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

// This algorithm is uniformly distributed over the surface of a sphere. There is no clustering at the poles.
void vRandomNormPCG3p(PCG* pcg, Vector3* out) {
	float u = frandPCG(-1.f, 1.f, pcg);
	float th = frandPCG(0.f, 2.f * F_PI, pcg);
	float q = sqrtf(1.f - u * u);
	float sth, cth;
	
	sincosf(th, &sth, &cth);
	out->x = u * cth;
	out->y = u * sth;
	out->z = u;
}

Vector3 vRandomNormPCG3(PCG* pcg) {
	Vector3 o;
	vRandomNormPCG3p(pcg, &o);
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


Vector4i vFloor4(const Vector4 v) {
	return (Vector4i){.x = floorf(v.x), .y = floorf(v.y), .z = floorf(v.z), .w = floorf(v.w)};
}

Vector4i vCeil4(const Vector4 v) {
	return (Vector4i){.x = ceilf(v.x), .y = ceilf(v.y), .z = ceilf(v.z), .w = ceilf(v.w)};
}

Vector3i vFloor3(const Vector3 v) {
	return (Vector3i){.x = floorf(v.x), .y = floorf(v.y), .z = floorf(v.z)};
}

Vector3i vCeil3(const Vector3 v) {
	return (Vector3i){.x = ceilf(v.x), .y = ceilf(v.y), .z = ceilf(v.z)};
}

Vector2i vFloor2(const Vector2 v) {
	return (Vector2i){.x = floorf(v.x), .y = floorf(v.y)};
}

Vector2i vCeil2(const Vector2 v) {
	return (Vector2i){.x = ceilf(v.x), .y = ceilf(v.y)};
}



Vector4l vFloor4d(const Vector4d v) {
	return (Vector4l){.x = floor(v.x), .y = floor(v.y), .z = floor(v.z), .w = floor(v.w)};
}

Vector4l vCeil4d(const Vector4d v) {
	return (Vector4l){.x = ceil(v.x), .y = ceil(v.y), .z = ceil(v.z), .w = ceil(v.w)};
}

Vector3l vFloor3d(const Vector3d v) {
	return (Vector3l){.x = floor(v.x), .y = floor(v.y), .z = floor(v.z)};
}

Vector3l vCeil3d(const Vector3d v) {
	return (Vector3l){.x = ceil(v.x), .y = ceil(v.y), .z = ceil(v.z)};
}

Vector2l vFloor2d(const Vector2d v) {
	return (Vector2l){.x = floor(v.x), .y = floor(v.y)};
}

Vector2l vCeil2d(const Vector2d v) {
	return (Vector2l){.x = ceil(v.x), .y = ceil(v.y)};
}



Vector2 vModPositive2(Vector2 v, Vector2 m) {
	return (Vector2){
		.x = fmodf(fmodf(v.x, m.x) + m.x, m.x),
		.y = fmodf(fmodf(v.y, m.y) + m.y, m.y),
	};
}

Vector3 vModPositive3(Vector3 v, Vector3 m) {
	return (Vector3){
		.x = fmodf(fmodf(v.x, m.x) + m.x, m.x),
		.y = fmodf(fmodf(v.y, m.y) + m.y, m.y),
		.z = fmodf(fmodf(v.z, m.z) + m.z, m.z),
	};
}

Vector4 vModPositive4(Vector4 v, Vector4 m) {
	return (Vector4){
		.x = fmodf(fmodf(v.x, m.x) + m.x, m.x),
		.y = fmodf(fmodf(v.y, m.y) + m.y, m.y),
		.z = fmodf(fmodf(v.z, m.z) + m.z, m.z),
		.w = fmodf(fmodf(v.w, m.w) + m.w, m.w),
	};
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
Vector3 vTriFaceNormal3(Vector3 a, Vector3 b, Vector3 c) {
	Vector3 b_a, c_a, out;
	
	b_a = vSub3(b, a);
	c_a = vSub3(c, a);
	return vNorm3(vCross3(b_a, c_a));
}

// calculate a unit vector normal to a triangle's face.
Vector3 vTriFaceNormalArea3(Vector3 a, Vector3 b, Vector3 c, float* area) {
	Vector3 b_a, c_a, out;
	
	b_a = vSub3(b, a);
	c_a = vSub3(c, a);
	
	Vector3 n = vCross3(b_a, c_a);
	
	if(area) *area = vLen(n) * .5f;
	
	return vNorm3(n);
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




void planeFromPointNormal(Vector3* p, Vector3* norm, Plane* out) {
	out->n = *norm;
	out->d = vDot3p(norm, p);
}

// calculates a plane from a triangle
void planeFromTriangle3p(Vector3* v1, Vector3* v2, Vector3* v3, Plane* out) {
	vTriFaceNormal3p(v1, v2, v3, &out->n);
	out->d = vDot3p(&out->n, v1);
}

// copy a plane
void planeCopy3p(Plane* in, Plane* out) {
	out->n = in->n;
	out->d = in->d;
}

// reverses the plane's direction
void planeInverse3p(Plane* in, Plane* out) {
	vInv3p(&in->n, &out->n);
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


// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
// returns _INTERSECT or _DISJOINT
int rayTriangleIntersect(
	Vector3* a, Vector3* b, Vector3* c, // triangle
	Vector3* ray_origin, Vector3* ray_dir, // ray
	float* u, float* v, float* t // barycentric out coords, t of intersection point along ray 
) {

	Vector3 ab = vSub3(*b, *a);
	Vector3 ac = vSub3(*c, *a);
	
	Vector3 n = vCross3(ab, ac);
	
	float det = -vDot3(*ray_dir, n);
	if(fabsf(det) <= FLT_CMP_EPSILON) {
		return C3DLAS_DISJOINT; // the ray is parallel to the triangle
	}
	
	float idet = 1.0f / det;
	
	Vector3 ao = vSub3(*ray_origin, *a);
	Vector3 dao = vCross3(ao, *ray_dir);
	
	*u = vDot3(ac, dao) * idet;
	if(*u < 0.f) return C3DLAS_DISJOINT; // barycentric coord is outside the triangle
	
	*v = -vDot3(ab, dao) * idet;
	if(*v < 0.f || *u + *v > 1.f) return C3DLAS_DISJOINT; // barycentric coord is outside the triangle
	
	
	*t = vDot3(ao, n) * idet;
//	if(*t < 0.0f) return C3DLAS_DISJOINT; // the ray intersects the triangle behind the origin
	
	return C3DLAS_INTERSECT;
}




// returns _INTERSECT or _DISJOINT
Vector3 triangleClosestPoint_Reference(
	Vector3* a, Vector3* b, Vector3* c, // triangle
	Vector3* p, // test point
	float* out_u, float* out_v // barycentric out coords of closest point 
) {

	Vector3 ab = vSub3(*b, *a);
	Vector3 ac = vSub3(*c, *a);
	
	Vector3 n = vCross3(ab, ac);
	Vector3 ray_dir = vNeg3(vNorm(n));
	
	float idet = 1.0f / -vDot3(ray_dir, n);
	
	Vector3 ao = vSub3(*p, *a);
	Vector3 dao = vCross3(ao, ray_dir);
	
//	printf("idet = %f, n = %f,%f,%f\n", idet, n.x, n.y, n.z);
	
	float u = vDot3(ac, dao) * idet;
	float v = -vDot3(ab, dao) * idet;
//	printf("u,v = %f, %f\n", u, v);
	if(u >= 0 && v >= 0.f && u + v <= 1.f) {
		float nt = vDot3(ao, n);
		Vector3 planep = vAdd3(*p, vScale3(vNeg3(n), nt));
		return planep; // the ray intersects the triangle
	}
	
	float t_ab, t_bc, t_ca;
	
	// collect all the possible locations
	float dist[6];
	
	dist[0] = vDistTPointLine3(*p, (Line3){*a, *b}, &t_ab);
	dist[1] = vDistTPointLine3(*p, (Line3){*b, *c}, &t_bc);
	dist[2] = vDistTPointLine3(*p, (Line3){*c, *a}, &t_ca);

	dist[3] = vDist(*a, *p);
	dist[4] = vDist(*b, *p);
	dist[5] = vDist(*c, *p);

	// find the smallest distance
	float min = dist[0];
	int mini = 0;
	
	for(int i = 1; i < 6; i++) {
		if(dist[i] < min) {
			min = dist[i];
			mini = i;
		}
	}
	
	switch(mini) {
		case 0: return vLerp(*a, *b, t_ab);
		case 1: return vLerp(*b, *c, t_bc);
		case 2: return vLerp(*c, *a, t_ca);
		case 3: return *a;
		case 4: return *b;
		case 5: return *c;
	}
	
	return (Vector3){0,0,0}; // HACK just return something
}


/*
// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
// returns _INTERSECT or _DISJOINT
vec3 triangleClosestPoint(
	Vector3* a, Vector3* b, Vector3* c, // triangle
	Vector3* p, // test point
	float* out_u, float* out_v // barycentric out coords of closest point 
) {

	Vector3 ab = vSub3(*b, *a); 
	Vector3 ac = vSub3(*c, *a);
	
	Vector3 n = vCross3(ab, ac); // triangle plane normal
	

	
	Vector3 ap = vSub3(*p, *a);
	Vector3 dap = vCross3(ap, vNeg(n)); // p projected onto the triangle's plane, relative to a
	
//	p = w*a + u*b + v*c;
	
	float u = vDot3(ac, dap); // inversely proportional to distance from _b_, aka "beta"
	// u < 0 means outside the triangle past the a/c edge
	// u > 1 means outside the triangle past b
	// u == 0 means somewhere on the a/c edge

//	if(*u < 0.f) return C3DLAS_DISJOINT; // barycentric coord is outside the triangle
	
	float v = -vDot3(ab, dap); // inversely proportional to distance from _c_, aka "gamma"
	// v < 0 means outside the triangle past the a/b edge
	// v > 1 means outside the triangle past c
	// v == 0 means somewhere on the a/b edge

//	if(*u < 0.f || *u + *v < 1.0f) return C3DLAS_DISJOINT; // barycentric coord is outside the triangle
	
	float w = 1.0f - u - v; // inversely proportional to distance from _a_, aka "alpha"
	// w < 0 means outside the triangle past the b/c edge
	// w > 1 means outside the triangle past a
	// w == 0 means somewhere on the b/c edge
	
	
	if(u > 0 && v > 0 && w > 0) // point inside triangle
		float t = vDot3(ap, n);
		vec3 closest = vAdd3(*p, vScale3(vNeg(n), t)); 
		return closest;
	}
	
	float new_u = 0, new_v = 0, new_w = 0;
	
	
	if(w < 0) {
		float t = fclamp(0f, 1f, vDot3() / vDot3());
		new_v = 1.f - t;
		new_w = t;
	}
	else if(v < 0) {
		float t = fclamp(0f, 1f, vDot3() / vDot3());
		new_u = t;
		new_w = 1.f - t;
	}
	else if(u < 0) {
		float t = fclamp(0f, 1f, vDot3() / vDot3());
		new_u = 1.f - t;
		new_v = t;
	}
	
	
//	if(*t < 0.0f) return C3DLAS_DISJOINT; // the ray intersects the triangle behind the origin
	
	return C3DLAS_INTERSECT;
}


*/

float distLineLine3(Line3* a, Line3* b) {
	
	Vector3 ea = vSub3(a->end, a->start);
	Vector3 eb = vSub3(b->end, b->start);
	Vector3 q = vSub(b->start, a->start);
	
	
	float vaa = vLenSq3(ea);
	float vbb = vLenSq3(eb);
	float vba = vDot3(ea, eb);
	float vba_a = vDot3(q, ea);
	
	float den = vba * vba - vbb * vaa;
	
	float ta, tb;
	if(fabs(den) < 1e-6) {
		ta = 0;
		tb = -vba_a / vba; // vba can never be zero here
	}
	else {
		float vba_b = vDot3(q, eb);
		
		ta = (vba_b * vba - vbb * vba_a) / den;
		tb = (-vba_a * vba + vaa * vba_b) / den;
	}
	
	ta = fclamp(0, 1, ta);
	tb = fclamp(0, 1, tb);
	
	Vector3 pa = vAdd(a->start, vScale(ea, ta));
	Vector3 pb = vAdd(b->start, vScale(eb, tb));

	return vDist(pa, pb);
}


/*
float distLineLine3Slow(Line3* a, Line3* b) {
	
	Vector3 ea = vSub3(a->end, a->start);
	Vector3 eb = vSub3(b->end, b->start);
	
	Vector3 n = vCross3(ea, eb);
	float nsq = vLenSq3(n);
	
	// TODO: handle parallel lines
	
	vec3 b1ma1 = vSub(b->start, a->start);
	
	float ta = vDot3(vCross3(eb, n), b1ma1) / nsq;
	float tb = vDot3(vCross3(ea, n), b1ma1) / nsq;
	
	ta = fclamp(0, 1, ta);
	tb = fclamp(0, 1, tb);
	
	vec3 pa = vAdd(a->start, vScale(ea, ta));
	vec3 pb = vAdd(b->start, vScale(eb, tb));
	
	return vDist3(pa, pb);
}
*/


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

//
Line3 shortestLineFromLineToLine(Line3* a, Line3* b) {
	
	Vector3 ea = vSub3(a->end, a->start);
	Vector3 eb = vSub3(b->end, b->start);
	Vector3 q = vSub(b->start, a->start);
	
	float vaa = vLenSq3(ea);
	float vbb = vLenSq3(eb);
	
	float vba = vDot3(ea, eb);
	float vba_a = vDot3(q, ea);
	
	float den = vba * vba - vbb * vaa;
	
	float ta, tb;
	if(fabs(den) < 1e-6) {
		ta = 0;
		tb = -vba_a / vba; // vba can never be zero here
	}
	else {
		float vba_b = vDot3(q, eb);
		
		ta = (vba_b * vba - vbb * vba_a) / den;
		tb = (-vba_a * vba + vaa * vba_b) / den;
	}
	
	ta = fclamp(ta, 0, 1);
	tb = fclamp(tb, 0, 1);
	
	Vector3 pa = vAdd(a->start, vScale(ea, ta));
	Vector3 pb = vAdd(b->start, vScale(eb, tb));

	return (Line3){pa, pb};
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

void frustumFromMatrixVK(Matrix* m, Frustum* out) {
	
	Matrix inv;
	
	mInverse(m, &inv);
	
	// first the points
	// these MUST be in this order
	// near
	vMatrixMulf3p(-1,-1, 0, &inv, &out->points[0]);
	vMatrixMulf3p(-1, 1, 0, &inv, &out->points[1]);
	vMatrixMulf3p( 1,-1, 0, &inv, &out->points[2]);
	vMatrixMulf3p( 1, 1, 0, &inv, &out->points[3]);
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

void frustumFromMatrixVK_ZUP(Matrix* m, Frustum* out) {
	
	Matrix inv;
	
	mInverse(m, &inv);
	
	*out = (Frustum){0};
	// first the points
	// these MUST be in this order
	// near
	vMatrixMulf3p(-1,-1, 0, &inv, &out->points[0]); // BUG this order is likely wrong for the planes but results in a sane wireframe.
	vMatrixMulf3p( 1,-1, 0, &inv, &out->points[1]);
	vMatrixMulf3p(-1, 1, 0, &inv, &out->points[2]);
	vMatrixMulf3p( 1, 1, 0, &inv, &out->points[3]);
	// far
	vMatrixMulf3p(-1,-1, 1, &inv, &out->points[4]);
	vMatrixMulf3p(1, -1, 1, &inv, &out->points[5]);
	vMatrixMulf3p( -1,1, 1, &inv, &out->points[6]);
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
	
	return (Vector3){.x = v.x / v.w, .y = v.y / v.w, .z = v.z / v.w};
}

Vector3 vMatrixMul3(Vector3 in, Matrix* m) {
	Vector4 v;

	v.x = in.x * m->m[0+0] + in.y * m->m[4+0] + in.z * m->m[8+0] + 1.0 * m->m[12+0];
	v.y = in.x * m->m[0+1] + in.y * m->m[4+1] + in.z * m->m[8+1] + 1.0 * m->m[12+1];
	v.z = in.x * m->m[0+2] + in.y * m->m[4+2] + in.z * m->m[8+2] + 1.0 * m->m[12+2];
	v.w = in.x * m->m[0+3] + in.y * m->m[4+3] + in.z * m->m[8+3] + 1.0 * m->m[12+3];
	
	if(v.w == 0) return (Vector3){0,0,0};
	
	return (Vector3){.x = v.x / v.w, .y = v.y / v.w, .z = v.z / v.w};
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



// Catmull-Rom Spline

float evalCatmullRom1D(float t, float a, float b, float c, float d) {
	
	float t2 = t * t;
	float t3 = t2 * t;

	return 0.5f * (   	
		(2.f * b) + 
		(-a + c) * t +
		(2.f * a - 5.f * b + 4.f * c - d) * t2 +
		(-a + 3.f * b - 3.f * c + d) * t3
	);

}

Vector2 evalCatmullRom2D(float t, Vector2 a, Vector2 b, Vector2 c, Vector2 d) {
	
	float t2 = t * t;
	float t3 = t2 * t;
	
	float q0 = -t3 + 2.f * t2 - t;
	float q1 = 3.f * t3 - 5.f * t2 + 2.f;
	float q2 = -3.f * t3 + 4.f * t2 + t;
	float q3 = t3 - t2;
	
	return (Vector2){
		.x = 0.5f * (a.x * q0 + b.x * q1 + c.x * q2 + d.x * q3),
		.y = 0.5f * (a.y * q0 + b.y * q1 + c.y * q2 + d.y * q3)
	};
}

Vector3 evalCatmullRom3D(float t, Vector3 a, Vector3 b, Vector3 c, Vector3 d) {
	
	float t2 = t * t;
	float t3 = t2 * t;
	
	float q0 = -t3 + 2.f * t2 - t;
	float q1 = 3.f * t3 - 5.f * t2 + 2.f;
	float q2 = -3.f * t3 + 4.f * t2 + t;
	float q3 = t3 - t2;
	
	return (Vector3){
		.x = 0.5f * (a.x * q0 + b.x * q1 + c.x * q2 + d.x * q3),
		.y = 0.5f * (a.y * q0 + b.y * q1 + c.y * q2 + d.y * q3),
		.z = 0.5f * (a.z * q0 + b.z * q1 + c.z * q2 + d.z * q3)
	};
}


float evalCatmullRom1D_dt(float t, float a, float b, float c, float d) {
	
	float t2 = t * t;
	
	float q0 = -3.f * t2 + 4.f * t - 1.f;
	float q1 = 9.f * t2 - 10.f * t;
	float q2 = -9.f * t2 + 8.f * t + 1.f;
	float q3 = 3.f * t2 - 2.f * t;
	
	return 0.5f * (a * q0 + b * q1 + c * q2 + d * q3);
}

Vector2 evalCatmullRom2D_dt(float t, Vector2 a, Vector2 b, Vector2 c, Vector2 d) {
	
	float t2 = t * t;
	
	float q0 = -3.f * t2 + 4.f * t - 1.f;
	float q1 = 9.f * t2 - 10.f * t;
	float q2 = -9.f * t2 + 8.f * t + 1.f;
	float q3 = 3.f * t2 - 2.f * t;
	
	return (Vector2){
		.x = 0.5f * (a.x * q0 + b.x * q1 + c.x * q2 + d.x * q3),
		.y = 0.5f * (a.y * q0 + b.y * q1 + c.y * q2 + d.y * q3)
	};
}

Vector3 evalCatmullRom3D_dt(float t, Vector3 a, Vector3 b, Vector3 c, Vector3 d) {
	
	float t2 = t * t;
	
	float q0 = -3.f * t2 + 4.f * t - 1.f;
	float q1 = 9.f * t2 - 10.f * t;
	float q2 = -9.f * t2 + 8.f * t + 1.f;
	float q3 = 3.f * t2 - 2.f * t;
	
	return (Vector3){
		.x = 0.5f * (a.x * q0 + b.x * q1 + c.x * q2 + d.x * q3),
		.y = 0.5f * (a.y * q0 + b.y * q1 + c.y * q2 + d.y * q3),
		.z = 0.5f * (a.z * q0 + b.z * q1 + c.z * q2 + d.z * q3)
	};
}


float evalCatmullRom1D_both(float t, float a, float b, float c, float d, float* dt) {
	
	float t2 = t * t;
	float t3 = t2 * t;
	
	float q0 = -t3 + 2.f * t2 - t;
	float q1 = 3.f * t3 - 5.f * t2 + 2.f;
	float q2 = -3.f * t3 + 4.f * t2 + t;
	float q3 = t3 - t2;
	
	float dq0 = -3.f * t2 + 4.f * t - 1.f;
	float dq1 = 9.f * t2 - 10.f * t;
	float dq2 = -9.f * t2 + 8.f * t + 1.f;
	float dq3 = 3.f * t2 - 2.f * t;
	
	*dt = 0.5f * (a * dq0 + b * dq1 + c * dq2 + d * dq3);
	
	return 0.5f * (a * q0 + b * q1 + c * q2 + d * q3);
}


Vector2 evalCatmullRom2D_both(float t, Vector2 a, Vector2 b, Vector2 c, Vector2 d, Vector2* dt) {
	
	float t2 = t * t;
	float t3 = t2 * t;
	
	float q0 = -t3 + 2.f * t2 - t;
	float q1 = 3.f * t3 - 5.f * t2 + 2.f;
	float q2 = -3.f * t3 + 4.f * t2 + t;
	float q3 = t3 - t2;
	
	float dq0 = -3.f * t2 + 4.f * t - 1.f;
	float dq1 = 9.f * t2 - 10.f * t;
	float dq2 = -9.f * t2 + 8.f * t + 1.f;
	float dq3 = 3.f * t2 - 2.f * t;
	
	*dt = (Vector2){
		.x = 0.5f * (a.x * dq0 + b.x * dq1 + c.x * dq2 + d.x * dq3),
		.y = 0.5f * (a.y * dq0 + b.y * dq1 + c.y * dq2 + d.y * dq3)
	};
	
	return (Vector2){
		.x = 0.5f * (a.x * q0 + b.x * q1 + c.x * q2 + d.x * q3),
		.y = 0.5f * (a.y * q0 + b.y * q1 + c.y * q2 + d.y * q3)
	};
}

Vector3 evalCatmullRom3D_both(float t, Vector3 a, Vector3 b, Vector3 c, Vector3 d, Vector3* dt) {
	
	float t2 = t * t;
	float t3 = t2 * t;
	
	float q0 = -t3 + 2.f * t2 - t;
	float q1 = 3.f * t3 - 5.f * t2 + 2.f;
	float q2 = -3.f * t3 + 4.f * t2 + t;
	float q3 = t3 - t2;
	
	float dq0 = -3.f * t2 + 4.f * t - 1.f;
	float dq1 = 9.f * t2 - 10.f * t;
	float dq2 = -9.f * t2 + 8.f * t + 1.f;
	float dq3 = 3.f * t2 - 2.f * t;
	
	*dt = (Vector3){
		.x = 0.5f * (a.x * dq0 + b.x * dq1 + c.x * dq2 + d.x * dq3),
		.y = 0.5f * (a.y * dq0 + b.y * dq1 + c.y * dq2 + d.y * dq3),
		.z = 0.5f * (a.z * dq0 + b.z * dq1 + c.z * dq2 + d.z * dq3)
	};
	
	return (Vector3){
		.x = 0.5f * (a.x * q0 + b.x * q1 + c.x * q2 + d.x * q3),
		.y = 0.5f * (a.y * q0 + b.y * q1 + c.y * q2 + d.y * q3),
		.z = 0.5f * (a.z * q0 + b.z * q1 + c.z * q2 + d.z * q3)
	};
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




#include "matrix3.c"
#include "matrix4.c"
#include "quaternion.c"
#include "intersect/plane.c"
#include "intersect/box.c"



