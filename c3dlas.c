

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
	uint64_t last = *state;
	*state = (last * 6364136223846793005ULL) + (stream | 1);
	uint32_t xs = ((last >> 18) ^ last) >> 27;
	uint32_t rot = last >> 59;
	uint32_t f = (((xs >> rot) | (xs << ((-rot) & 31))) & 0x807fffff) | 0x3f000000;
	return *(float*)&f;
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





// vector operations

int vEq(Vector* a, Vector* b) {
	return vEqEp(a, b, FLT_CMP_EPSILON);
}

int vEqEp(Vector* a, Vector* b, float epsilon) {
	float x, y, z, n;
	
	x = a->x - b->x;
	y = a->y - b->y;
	z = a->z - b->z;
	
	n = fabs(x * x + y * y + z * z);
	
	return n <= epsilon * epsilon;
}


void vCopy(const Vector* src, Vector* dst) {
	*dst = *src;
}

void vSwap(Vector* a, Vector* b) { // swap two vectors
	Vector t;
	t = *b;
	*b = *a;
	*a = t;
}

void vAdd(Vector* a, Vector* b, Vector* out) {
// #ifdef C3DLAS_USE_SIMD
// 	__m128 a_ = _mm_loadu_ps((float*)a);
// 	__m128 b_ = _mm_loadu_ps((float*)b);
// 	b_ = _mm_add_ps(a_, b_);
	// requires avx
// 	_mm_maskstore_ps((float*)out, _mm_set_epi32(0, -1, -1, -1), b_);
// #else
	out->x = a->x + b->x;
	out->y = a->y + b->y;
	out->z = a->z + b->z;
// #endif
}

void vSub(Vector* from, Vector* what, Vector* diff) { // diff = from - what
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

void vScale(Vector* v, float scalar, Vector* out) {
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

void  vLerp(Vector* a, Vector* b, float t, Vector* out) { // Linear interpolation between two vectors
#ifdef C3DLAS_USE_SIMD
	__m128 a_ = _mm_loadu_ps((float*)&a);
	__m128 b_ = _mm_loadu_ps((float*)&b);
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

void vInverse(Vector* v, Vector* out) {
	// yeah yeah yeah shut up. in games, pesky details are just annoying. this function does what you mean rather than sucking Gauss and Euclid.
	out->x = v->x == 0.0f ? FLT_MAX : 1.0f / v->x;
	out->y = v->y == 0.0f ? FLT_MAX : 1.0f / v->y;
	out->z = v->z == 0.0f ? FLT_MAX : 1.0f / v->z;
}

float vMag(Vector* v) {
#ifdef C3DLAS_USE_SIMD
	// needs sse4.1
	__m128 a = _mm_loadu_ps((float*)&v);
	__m128 b = _mm_dp_ps(a, a, (_1110b << 4) & _0001b); // BUG: check mask
	__m128 c = _mm_sqrt_ss(b);
	return c[0];
#else
	return sqrt((float)((v->x * v->x) + (v->y * v->y) + (v->z * v->z)));
#endif
}

float vDot(Vector* a, Vector* b) {
#ifdef broken_C3DLAS_USE_SIMD
	// needs sse4.1
	// BUG does not work
	__m128 a_ = _mm_loadu_ps((float*)&a);
	__m128 b_ = _mm_loadu_ps((float*)&b);
	__m128 c = _mm_dp_ps(a_, b_, (_1110b << 4) & _0001b); // BUG: check mask
	
	return c[0];
#else
	return (float)((a->x * b->x) + (a->y * b->y) + (a->z * b->z));
#endif
}

// distance from one point to another
float vDist(Vector* from, Vector* to) {
// #ifdef C3DLAS_USE_SIMD
// 	// needs sse4.1
// 	__m128 f = _mm_loadu_ps((float*)&from);
// 	__m128 t = _mm_loadu_ps((float*)&to);
// 	__m128 a = _mm_sub_ps(f, t);
// 	
// 	float b = a[1] + a[2] + a[3];
// 	return sqrt(b);
// // 	__m128 c = _mm_sqrt_ss(b);
// // 	return c[0];
// #else
	float dx = (from->x - to->x);
	float dy = (from->y - to->y);
	float dz = (from->z - to->z);
	
	return sqrt(dx*dx + dy*dy + dz*dz);
// #endif
}


// squared distance from one point to another
float vDistSq(Vector* from, Vector* to) {
// #ifdef C3DLAS_USE_SIMD
// 	// needs sse4.1
// 	__m128 f = _mm_loadu_ps((float*)&from);
// 	__m128 t = _mm_loadu_ps((float*)&to);
// 	__m128 a = _mm_sub_ps(f, t);
// 	__m128 b = _mm_dp_ps(a, a, (_1110b << 4) & _0001b); // BUG: check mask
// 	return b[0];
// #else
	float dx = (from->x - to->x);
	float dy = (from->y - to->y);
	float dz = (from->z - to->z);
	
	return dx*dx + dy*dy + dz*dz;
// #endif
}


void vNorm(Vector* v, Vector* out) {
	vUnit(v, out);
}

void vUnit(Vector* v, Vector* out) {
	float n;
	n = (v->x * v->x) + (v->y * v->y) + (v->z * v->z);
	
	if(n >= 1.0f - FLT_EPSILON && n <= 1.0f + FLT_EPSILON) return; // very exact here
	
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
#ifdef C3DLAS_USE_SIMD
	__m128 a_ = _mm_loadu_ps((float*)&a);
	__m128 b_ = _mm_loadu_ps((float*)&b);
	__m128 c = _mm_min_ps(a_, b_);
	_mm_maskstore_ps((float*)out, _mm_set_epi32(0, -1, -1, -1), c);
#else
	out->x = fmin(a->x, b->x);
	out->y = fmin(a->y, b->y);
	out->z = fmin(a->z, b->z);
#endif
}

// returns the maximum values of each component
void  vMax(Vector* a, Vector* b, Vector* out) {
#ifdef C3DLAS_USE_SIMD
	__m128 a_ = _mm_loadu_ps((float*)&a);
	__m128 b_ = _mm_loadu_ps((float*)&b);
	__m128 c = _mm_max_ps(a_, b_);
	
	_mm_maskstore_ps((float*)out, _mm_set_epi32(0, -1, -1, -1), c);
#else
	out->x = fmax(a->x, b->x);
	out->y = fmax(a->y, b->y);
	out->z = fmax(a->z, b->z);
#endif
}

void inline vSet(float x, float y, float z, Vector* out) {
	out->x = x;
	out->y = y;
	out->z = z;
}

void vRandom(Vector* end1, Vector* end2, Vector* out) {
	out->x = frand(fmin(end1->x, end2->x), fmax(end1->x, end2->x));
	out->y = frand(fmin(end1->y, end2->y), fmax(end1->y, end2->y));
	out->z = frand(fmin(end1->z, end2->z), fmax(end1->z, end2->z));
}

void vRandomNorm(Vector* out) {
	float x = frand(-1, 1);
	float y = frand(-1, 1);
	float z = frand(-1, 1);
	
	float r = sqrt(x*x + y*y + z*z);
	if(r == 0.0) {
		vRandomNorm(out); // in the rare case of a zero-length vector, try again
		return;
	}
	
	out->x = x / r;
	out->y = y / r;
	out->z = z / r;
}


void  vLerp4(Vector4* a, Vector4* b, float t, Vector4* out) { // Linear interpolation between two vectors
#ifdef C3DLAS_USE_SIMD
	__m128 a_ = _mm_loadu_ps((float*)&a);
	__m128 b_ = _mm_loadu_ps((float*)&b);
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




void vProjectOntoPlane(Vector* v, Plane* p, Vector* out) {
	Vector v_ortho;
	
	// get the component of v that's perpendicular to the plane
	vProjectNorm(v, &p->n, &v_ortho);
	
	// subtract it from v
	vSub(v, &v_ortho, out);
}

void vProjectOntoPlaneNormalized(Vector* v, Plane* p, Vector* out) {
	vProjectOntoPlane(v, p, out);
	vNorm(out, out);
}



// calculates a plane from a triangle
void planeFromTriangle(Vector* v1, Vector* v2, Vector* v3, Plane* out) {
	vTriFaceNormal(v1, v2, v3, &out->n);
	out->d = vDot(&out->n, v1);
}

// copy a plane
void planeCopy(Plane* in, Plane* out) {
	vCopy(&in->n, &out->n);
	out->d = in->d;
}

// reverses the plane's direction
void planeInverse(Plane* in, Plane* out) {
	vInverse(&in->n, &out->n);
	out->d = -in->d;
}

// classifies a point by which side of the plane it's on, default espilon
int planeClassifyPoint(Plane* p, Vector* pt) {
	return planeClassifyPointEps(p, pt, FLT_CMP_EPSILON);
}

// classifies a point by which side of the plane it's on, custom espilon
int planeClassifyPointEps(Plane* p, Vector* pt, float epsilon) {
	float dist = vDot(&p->n, pt);
	if(fabs(dist - p->d) < epsilon)
		return C3DLAS_COPLANAR;
	else if (dist < p->d)
		return C3DLAS_BACK;
	else
		return C3DLAS_FRONT;
}



// C3DLAS_INTERSECT, _COPLANAR or _DISJOINT
int planeLineFindIntersect(Plane* pl, Vector* la, Vector* lb, Vector* out) {
	Vector ldir;
	float da, db;
	
	vSub(lb, la, &ldir);
	
	da = vDot(&la, &pl->n) - pl->d;
	
	// bail if the line and plane are parallel
	if(fabs(vDot(&pl->n, &ldir)) < FLT_CMP_EPSILON) {
		
		// check coplanarity
		if(fabs(da) < FLT_CMP_EPSILON) {
			return C3DLAS_COPLANAR; // the end is on the plane, so the other is too
		}
		
		return C3DLAS_DISJOINT;
	}
	

	db = vDot(&lb, &pl->n) - pl->d;
	
	// check if one of the points if on the plane
	if(fabs(da) < FLT_CMP_EPSILON) {
		*out = *la;
		return C3DLAS_INTERSECT; // the end is on the plane, so the other is too
	}
	if(fabs(db) < FLT_CMP_EPSILON) {
		*out = *lb;
		return C3DLAS_INTERSECT;
	}
	
	Vector p0, g, j;
	vScale(&pl->n, pl->d, &p0);
	vSub(&p0, la, &g);
	float h = vDot(&g, &pl->n);
	float i = vDot(&ldir, &pl->n);
	float d = i != 0 ? h / i : 0;
	
	// check if the plane intersects outside the two points
	if(d < 0 || d > vDist(la, lb)) {
		return C3DLAS_DISJOINT;
	}
	
	vScale(&ldir, d, &j);
	vAdd(la, &j, out);
	
	return C3DLAS_INTERSECT;
}



// C3DLAS_COPLANAR, _INTERSECT, or _DISJOINT
int triPlaneTestIntersect(Vector* pTri, Plane* pl) {
	Vector a, b, c;
	float da, db, dc;
	
	// get distance of each vertex from the plane
	// bail early if any of them are coplanar
	a = pTri[0];
	da = vDot(&a, &pl->n) - pl->d;
	if(fabs(da) < FLT_CMP_EPSILON) {
		return C3DLAS_COPLANAR;
	}
	
	b = pTri[1];
	db = vDot(&b, &pl->n) - pl->d;
	if(fabs(db) < FLT_CMP_EPSILON) {
		return C3DLAS_COPLANAR;
	}
	
	c = pTri[2];
	dc = vDot(&c, &pl->n) - pl->d;
	if(fabs(dc) < FLT_CMP_EPSILON) {
		return C3DLAS_COPLANAR;
	}
	
	// the triangle intersects if the sign of all the distances does not match,
	// ie, on vertex is on the opposite side of the plane from the others
	return (signbit(da) == signbit(db) && signbit(db) == signbit(dc)) ? C3DLAS_DISJOINT : C3DLAS_INTERSECT;
}


// C3DLAS_COPLANAR, _INTERSECT, or _DISJOINT
int triPlaneFindIntersect(Vector* pTri, Plane* pl, Vector* pOut, unsigned char* lineMask) {
	Vector a, b, c;
	float da, db, dc;
	
	// get distance of each vertex from the plane
	// bail early if any of them are coplanar
	a = pTri[0];
	da = vDot(&a, &pl->n) - pl->d;
	if(fabs(da) < FLT_CMP_EPSILON) {
		return C3DLAS_COPLANAR;
	}
	
	b = pTri[1];
	db = vDot(&b, &pl->n) - pl->d;
	if(fabs(db) < FLT_CMP_EPSILON) {
		return C3DLAS_COPLANAR;
	}
	
	c = pTri[2];
	dc = vDot(&c, &pl->n) - pl->d;
	if(fabs(dc) < FLT_CMP_EPSILON) {
		return C3DLAS_COPLANAR;
	}
	
	// the triangle intersects if the sign of all the distances does not match,
	// ie, on vertex is on the opposite side of the plane from the others
	// bail if disjoint
	if(signbit(da) == signbit(db) && signbit(db) == signbit(dc)) { 
		return C3DLAS_DISJOINT;
	}
	
	
	
	
	
	return C3DLAS_INTERSECT;
}


void frustumFromMatrix(Matrix* m, Frustum* out) {
	
	Matrix inv;
	
	mInverse(m, &inv);
	
	// first the points
	// these MUST be in this order
	// near
	vMatrixMulf(-1,-1,-1, &inv, &out->points[0]);
	vMatrixMulf(-1, 1,-1, &inv, &out->points[1]);
	vMatrixMulf( 1,-1,-1, &inv, &out->points[2]);
	vMatrixMulf( 1, 1,-1, &inv, &out->points[3]);
	// far
	vMatrixMulf(-1,-1, 1, &inv, &out->points[4]);
	vMatrixMulf(-1, 1, 1, &inv, &out->points[5]);
	vMatrixMulf( 1,-1, 1, &inv, &out->points[6]);
	vMatrixMulf( 1, 1, 1, &inv, &out->points[7]);
	
	// now the planes
	// near and far
	planeFromTriangle(&out->points[0], &out->points[1], &out->points[2], &out->planes[0]);
	planeFromTriangle(&out->points[4], &out->points[5], &out->points[6], &out->planes[1]);
	// sides
	planeFromTriangle(&out->points[0], &out->points[4], &out->points[1], &out->planes[2]);
	planeFromTriangle(&out->points[0], &out->points[4], &out->points[2], &out->planes[3]);
	planeFromTriangle(&out->points[3], &out->points[7], &out->points[1], &out->planes[4]);
	planeFromTriangle(&out->points[3], &out->points[7], &out->points[2], &out->planes[5]);
}


void frustumCenter(Frustum* f, Vector* out) {
	Vector sum = {0.0f,0.0f,0.0f};
	
	for(int i = 0; i < 8; i++) vAdd(&f->points[i], &sum, &sum);
	vScale(&sum, 1.0f/8.0f, out);
}

// General idea of the algorithm:
// https://lxjk.github.io/2017/04/15/Calculate-Minimal-Bounding-Sphere-of-Frustum.html
// http://archive.is/YACj2
void frustumBoundingSphere(Frustum* f, Sphere* out) {
	Vector f0, n0;
	vPointAvg(&f->points[0], &f->points[3], &n0);
	vPointAvg(&f->points[4], &f->points[7], &f0);
	
	float Dn2 = vDistSq(&n0, &f->points[0]); 
	float Df2 = vDistSq(&f0, &f->points[4]);
	
	 // check for ortho
	if(Dn2 - Df2 < 0.00001) {
		frustumCenter(f, &out->center);
		out->r = vDist(&out->center, &f->points[0]); 
		return;
	}
	
	float Dnf = vDist(&f0, &n0);
	float Dnc = (Dn2 - Df2 - Df2) / (2 * Dnf);
	
// 	printf("\n f: %f,%f,%f\n", f->points[4].x,f->points[4].y,f->points[4].z);
// 	printf(" n: %f,%f,%f\n", f->points[0].x,f->points[0].y,f->points[0].z);
// 	printf(" f0: %f,%f,%f\n", f0.x,f0.y,f0.z);
// 	printf(" n0: %f,%f,%f\n", n0.x,n0.y,n0.z);
// 	printf(" dn2, df2, dnf, dnc: %f,%f,%f,%f\n", Dn2, Df2, Dnf, Dnc);
	
	
	if(Dnc > 0 && Dnc < Dnf) {
		vLerp(&f0, &n0, Dnc / Dnf, &out->center);
		out->r = sqrt(Dnc * Dnc + Dn2);
	}
	else {
		out->center = f0;
		out->r = sqrt(Df2);
	}
}


void frustumInscribeSphere(Frustum* f, Sphere* out) {
	Vector fx, nx;
	vPointAvg(&f->points[0], &f->points[3], &nx);
	vPointAvg(&f->points[4], &f->points[7], &fx);
	
/*	
	float Dn2 = vDistSq(&n0, &f->points[0]); 
	float Df2 = vDistSq(&f0, &f->points[4]);
	float Dnf = vDist(&f0, n0);
	float Dnc = (Dn2 - Df2 - Df2) / (2 * Dnf);*/

}



void quadCenterp(Vector* a, Vector* b, Vector* c, Vector* d, Vector* out) {
	Vector sum;
	vAdd(a, b, &sum);
	vAdd(&sum, c, &sum);
	vAdd(&sum, d, &sum);
	vScale(&sum, 0.25f, out);
}

// closest distance from an arbitrary point to the plane 
float planePointDist(Plane* pl, Vector* p) {
	Vector a;
	vScale(&pl->n, pl->d, &a);
	return fabs(vDot(&a, p));
} 

// signed closest distance from an arbitrary point to the plane 
float planePointDistSigned(Plane* pl, Vector* p) {
	Vector a;
	vScale(&pl->n, pl->d, &a);
	return vDot(&a, p);
} 

void vPointAvg(Vector* a, Vector* b, Vector* out) {
	Vector sum;
	vAdd(a, b, &sum);
	vScale(&sum, 0.5f, out);
} 

// 2d vector stuff

int vEq2(Vector2* a, Vector2* b) {
	return vEqEp2(a, b, FLT_CMP_EPSILON);
}

int vEqEp2(Vector2* a, Vector2* b, float epsilon) {
	float x, y, n;
	
	x = a->x - b->x;
	y = a->y - b->y;
	
	n = fabs(x * x + y * y);
	
	return n <= epsilon * epsilon;
}

void vCopy2(const Vector2* src, Vector2* dst) {
	dst->x = src->x;
	dst->y = src->y;
}

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

float vDist2(Vector2* a, Vector2* b) {
	float x = a->x - b->x;
	float y = a->y - b->y;
	return sqrt(x * x + y * y);
}

void  vLerp2(Vector2* a, Vector2* b, float t, Vector2* out) {
	out->x = a->x + ((b->x - a->x) * t);
	out->y = a->y + ((b->y - a->y) * t);
}

void vInverse2(Vector2* v, Vector2* out) {
	// see vInverse for snark
	out->x = v->x == 0.0f ? FLT_MAX : 1.0f / v->x;
	out->y = v->y == 0.0f ? FLT_MAX : 1.0f / v->y;
}

float vMag2(Vector2* v) {
	return sqrt((float)((v->x * v->x) + (v->y * v->y)));
}

float vDot2(Vector2* a, Vector2* b) {
	return (float)((a->x * b->x) + (a->y * b->y));
}

void vNorm2(Vector2* v, Vector2* out) {
	vUnit2(v, out);
}

void vUnit2(Vector2* v, Vector2* out) {
	float n;
	n = (v->x * v->x) + (v->y * v->y);
	
	if(n >= 1.0f - FLT_EPSILON && n <= 1.0f + FLT_EPSILON) return; // very exact here
	
	n = 1.0f / sqrtf(n);
	out->x = v->x * n;
	out->y = v->y * n;
}

// returns the minimum values of each component
void  vMin2(Vector2* a, Vector2* b, Vector2* out) {
	out->x = fmin(a->x, b->x);
	out->y = fmin(a->y, b->y);
}

// returns the maximum values of each component
void  vMax2(Vector2* a, Vector2* b, Vector2* out) {
	out->x = fmax(a->x, b->x);
	out->y = fmax(a->y, b->y);
}

void inline vSet2(float x, float y, Vector2* out) {
	out->x = x;
	out->y = y;
}


// reflects the distance from v to pivot across pivot.
// out, pivot, and v will form a straight line with pivot exactly in the middle.
void vReflectAcross2(Vector2* v, Vector2* pivot, Vector2* out) {
	Vector2 diff;
	
	vSub2(pivot, v, &diff);
	vAdd2(pivot, &diff, out);
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


// returns the *signed* area of a triangle. useful for determining winding
// positive values mean a clockwise triangle
float triArea2(Vector2* a, Vector2* b, Vector2* c) {
	return 0.5 * (
		((b->x - a->x) * (b->y + a->y)) +
		((c->x - b->x) * (c->y + b->y)) +
		((a->x - c->x) * (a->y + c->y)));
}


// determines if a point is inside a triangle
int triPointInside2(Vector2* p, Vector2* a, Vector2* b, Vector2* c) {
	int d = signbit((p->x - b->x) * (a->y - b->y) - (a->x - b->x) * (p->y - b->y));
	int e = signbit((p->x - c->x) * (b->y - c->y) - (b->x - c->x) * (p->y - c->y));
	if(d != e) return 0;
	int f = signbit((p->x - a->x) * (c->y - a->y) - (c->x - a->x) * (p->y - a->y));
	return e == f;
}




// 2d integer vector stuff

int vEq2i(Vector2i* a, Vector2i* b) {
	return a->x == b->x && a->y == b->y;
}

void vCopy2i(const Vector2i* src, Vector2i* dst) {
	dst->x = src->x;
	dst->y = src->y;
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

void vAdd2i(Vector2i* a, Vector2i* b, Vector2i* out) {
	out->x = a->x + b->x;
	out->y = a->y + b->y;
}

void vSub2i(Vector2i* from, Vector2i* what, Vector2i* diff) { // diff = from - what
	diff->x = from->x - what->x;
	diff->y = from->y - what->y;
}

void vScale2i(Vector2i* v, int scalar, Vector2i* out) {
	out->x = v->x * scalar;
	out->y = v->y * scalar;
}

int vDot2i(Vector2i* a, Vector2i* b) {
	return ((a->x * b->x) + (a->y * b->y));
}

// returns the minimum values of each component
void  vMin2i(Vector2i* a, Vector2i* b, Vector2i* out) {
	out->x = MIN(a->x, b->x);
	out->y = MIN(a->y, b->y);
}

// returns the maximum values of each component
void  vMax2i(Vector2i* a, Vector2i* b, Vector2i* out) {
	out->x = MAX(a->x, b->x);
	out->y = MAX(a->y, b->y);
}

void inline vSet2i(int x, int y, Vector2i* out) {
	out->x = x;
	out->y = y;
}

// returns the absolute distance between two vectors
float vDist2i(Vector2i* a, Vector2i* b) {
	float x = (float)a->x - (float)b->x;
	float y = (float)a->y - (float)b->y;
	return sqrt(x * x + y * y);
}





// plane-vector operations

// distance from point to plane
float pvDist(Plane* p, Vector* v) {
	return vDot(v, &p->n) + p->d;
}



// matrix-vector operations


// multiply a vector by a matrix
void vMatrixMul(Vector* restrict in, Matrix* restrict m, Vector* restrict out) {
	vMatrixMulf(in->x, in->y, in->z, m, out);
}

void vMatrixMulf(float x, float y, float z, Matrix* restrict m, Vector* restrict out) {
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
	int i;
	
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


void mOrthoFromSphere(Sphere* s, Vector* eyePos, Matrix* out) {
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
	
	Vector d;
	vSub(&s->center, eyePos, &d);
	vNorm(&d, &d);
	
	Matrix m2;

	
	
	m2 = IDENT_MATRIX;
// 	mRotX(asin(d.z) + (F_PI / 4)  , &m2);
	mRotY(atan2(d.z, d.x) + (F_PI / 2)  , &m2);
	mMul(&m2, &m);

	m2 = IDENT_MATRIX;
	mRotZ(asin(d.y), &m2);
	mMul(&m2, &m);
	

	m2 = IDENT_MATRIX;
	Vector ic;
	vScale(&s->center, -1, &ic);
	mTransv(&ic, &m2);
	
	mFastMul(&m2, &m, out);
}

// analgous to gluLookAt
// BUG: very broken apparently
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


// Quadratic bezier functions
float evalQBezier1D(float e1, float e2, float c1, float t) {
	float mt, mt2;
	mt = 1 - t;
	mt2 = mt * mt;
	
	return (mt2 * e1) + (2 * mt * t * c1) + (t * t * e2);
}

void evalQBezier2D(Vector2* e1, Vector2* e2, Vector2* c1, float t, Vector2* out) {
	out->x = evalQBezier1D(e1->x, e2->x, c1->x, t);
	out->y = evalQBezier1D(e1->y, e2->y, c1->y, t);
}

void evalQBezier(Vector* e1, Vector* e2, Vector* c1, float t, Vector* out) {
	out->x = evalQBezier1D(e1->x, e2->x, c1->x, t);
	out->y = evalQBezier1D(e1->y, e2->y, c1->y, t);
	out->z = evalQBezier1D(e1->z, e2->z, c1->z, t);
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

void boxSize2i(const AABB2i* b, Vector2* out) {
	out->x = b->max.x - b->min.x;
	out->y = b->max.y - b->min.y;
}

// BUG: needs some fancy math work to keep everything tight. integers don't split nicely
void boxQuadrant2i(const AABB2i* in, char ix, char iy, AABB2i* out) {
	Vector2 sz, c;
	
	printf("fix me: %s:%d", __FILE__, __LINE__);
	exit(666);
	
	boxCenter2i(in, &c);
	boxSize2i(in, &sz);
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
}

// ray stuff
void makeRay2(Vector2* origin, Vector2* direction, Ray2* out) {
	
	out->o.x = origin->x;
	out->o.y = origin->y;
	
	vNorm2(direction, &out->d);
}

// this version has no branching, but only answers yes or no.
// algorithm explanation here. hopefully my extrapolation into 3 dimensions is correct.
// http://tavianator.com/fast-branchless-raybounding-box-intersections/
int boxRayIntersectFast(const AABB* b, const Ray* r) {
	Vector t1, t2;
	float tmin, tmax;
	Vector id;
	
	vInverse(&r->d, &id);
	
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
	
	vInverse(&r->d, &id);
	
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
int boxRayIntersect(const AABB* b, const Ray* r, Vector* ipoint, float* idist) {
	Vector t1, t2, id;
	float tmin, tmax;
	
	vInverse(&r->d, &id);
		
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
	vReflectAcross2(&n->c, &n->e, &out[2]);
	
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




Vector evalCubicHermite3D(float t, Vector p0, Vector p1, Vector m0, Vector m1) {
	
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
		Vector v3;
	} u;
	_mm_storeu_ps(&u.v4, o);
	
	return u.v3;
#else
	return (Vector){
		.x = evalCubicHermite1D(t, p0.x, p1.x, m0.x, m1.x),
		.y = evalCubicHermite1D(t, p0.y, p1.y, m0.y, m1.y),
		.z = evalCubicHermite1D(t, p0.z, p1.z, m0.z, m1.z)
	};
#endif
}






