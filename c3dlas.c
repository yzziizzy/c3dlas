


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




char const * const C3DLAS_retval__names[] = {
#define X(name, ...) [C3DLAS_##name] = "C3DLAS_" #name,
	C3DLAS_RETURN_VALUE_LIST(X)
#undef X
	0
};





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

// prepare the pcg for use
void pcg_init(PCG* pcg, uint64_t seed) {
	pcg->stream = seed;
	pcg->state = (seed * 6364136223846793005ULL) + (seed | 1);
}


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


// returns a random number in [0, UINT32_MAX] inclusive
uint32_t pcg_u32(uint64_t* state, uint64_t stream) {
	
	uint64_t last = *state;
	*state = (last * 6364136223846793005ULL) + (stream | 1);
	uint32_t xs = ((last >> 18) ^ last) >> 27;
	uint32_t rot = last >> 59;
	uint32_t fin = (xs >> rot) | (xs << ((-rot) & 31));
	
	return fin;
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


uint32_t urandPCG(uint32_t low, uint32_t high, PCG* pcg) {
	return pcg_u32(&pcg->state, pcg->stream) % (high - low) + low;
}

int32_t irandPCG(int32_t low, int32_t high, PCG* pcg) {
	return pcg_u32(&pcg->state, pcg->stream) % (high - low) + low;
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
	out->d = -vDot3p(p, norm);
	
	// negation is suspicious but makes these planes work with intersectPlaneRay3p
}

// calculates a plane from a triangle
void planeFromTriangle3p(Vector3* v1, Vector3* v2, Vector3* v3, Plane* out) {
	vTriFaceNormal3p(v1, v2, v3, &out->n);
	out->d = -vDot3p(&out->n, v1);
	
	// negation is suspicious but makes these planes work with intersectPlaneRay3p
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







#include "vector_fns.c"
#include "matrix3.c"
#include "matrix4.c"
#include "frustum.c"
#include "quaternion.c"
#include "quad.c"
#include "line.c"
#include "polygon.c"
#include "splines.c"
#include "intersect/circle.c"
#include "intersect/plane.c"
#include "intersect/box.c"
#include "intersect/line.c"
#include "intersect/triangle.c"



