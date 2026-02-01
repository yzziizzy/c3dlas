#ifndef __c3dlas__simd_h__
#define __c3dlas__simd_h__


#ifndef C3DLAS_USE_AVX  


typedef struct {
	union { float xf[4]; };
	union { float yf[4]; };
} mm4_vec2;

typedef struct {
	union { float xf[8]; };
	union { float yf[8]; };
} mm8_vec2;


typedef struct {
	union { mm8_vec2 start, a; };
	union { mm8_vec2 end, b; };
} mm8_Line2;



#else 

#include <immintrin.h>






////////////////////////
//
//  Missing Functions
//
////////////////////////


static inline __m128 mm_not_ps(__m128 a) {
	return _mm_xor_ps(a, _mm_cmp_ps(a, a, _CMP_EQ_OQ));
}
static inline __m128d mm_not_pd(__m128d a) {
	return _mm_xor_pd(a, _mm_cmp_pd(a, a, _CMP_EQ_OQ));
}

static inline __m256 mm256_not_ps(__m256 a) {
	return _mm256_andnot_ps(a, _mm256_cmp_ps(a, a, _CMP_EQ_OQ));
}
static inline __m256d mm256_not_pd(__m256d a) {
	return _mm256_xor_pd(a, _mm256_cmp_pd(a, a, _CMP_EQ_OQ));
}


static inline __m128 mm_abs_ps(__m128 a) {
    static const __m128 nz = {-0.f, -0.f, -0.f, -0.f};
    return _mm_andnot_ps(nz, a);
}

static inline __m128d mm_abs_pd(__m128d a) {
    static const __m128d nz = {-0., -0.};
    return _mm_andnot_pd(nz, a);
}


static inline __m256 mm256_abs_ps(__m256 a) {
    static const __m256 nz = {-0.f, -0.f, -0.f, -0.f, -0.f, -0.f, -0.f, -0.f};
    return _mm256_andnot_ps(nz, a);
}

static inline __m256d mm256_abs_pd(__m256d a) {
    static const __m256d nz = {-0., -0., -0., -0.};
    return _mm256_andnot_pd(nz, a);
}


static inline __m128 mm_clamp_ps(__m128 a, __m128 vmin, __m128 vmax) {
    return _mm_max_ps(vmin, _mm_min_ps(a, vmax));
}

static inline __m128d mm_clamp_pd(__m128d a, __m128d vmin, __m128d vmax) {
    return _mm_max_pd(vmin, _mm_min_pd(a, vmax));
}

static inline __m256 mm256_clamp_ps(__m256 a, __m256 vmin, __m256 vmax) {
    return _mm256_max_ps(vmin, _mm256_min_ps(a, vmax));
}

static inline __m256d mm256_clamp_pd(__m256d a, __m256d vmin, __m256d vmax) {
    return _mm256_max_pd(vmin, _mm256_min_pd(a, vmax));
}



static inline __m128 mm_clamp_scalar_ps(__m128 a, float min, float max) {
	__m128 vmin = _mm_set1_ps(min);
	__m128 vmax = _mm_set1_ps(max);
    return _mm_max_ps(vmin, _mm_min_ps(a, vmax));
}

static inline __m128d mm_clamp_scalar_pd(__m128d a, double min, double max) {
	__m128d vmin = _mm_set1_pd(min);
	__m128d vmax = _mm_set1_pd(max);
    return _mm_max_pd(vmin, _mm_min_pd(a, vmax));
}

static inline __m256 mm256_clamp_scalar_ps(__m256 a, float min, float max) {
	__m256 vmin = _mm256_set1_ps(min);
	__m256 vmax = _mm256_set1_ps(max);
    return _mm256_max_ps(vmin, _mm256_min_ps(a, vmax));
}

static inline __m256d mm256_clamp_scalar_pd(__m256d a, double min, double max) {
	__m256d vmin = _mm256_set1_pd(min);
	__m256d vmax = _mm256_set1_pd(max);
    return _mm256_max_pd(vmin, _mm256_min_pd(a, vmax));
}



////////////////////////
//
//  2D Vectors
//
////////////////////////



typedef struct {
	union { __m128 x; float xf[4]; };
	union { __m128 y; float yf[4]; };
} mm4_vec2;

typedef struct {
	union { __m256 x; float xf[8]; };
	union { __m256 y; float yf[8]; };
} mm8_vec2;


typedef struct {
	union { mm8_vec2 start, a; };
	union { mm8_vec2 end, b; };
} mm8_Line2;



// __m128


// __m256



static inline mm8_vec2 mm8_load_arrays_v2(float x[8], float y[8]) {
	return (mm8_vec2) {
		.x = _mm256_loadu_ps(x),
		.y = _mm256_loadu_ps(y)
	};
}

static inline mm8_vec2 mm8_set1_v2(Vector2 v) {
    return (mm8_vec2) {
		.x = _mm256_set1_ps(v.x),
		.y = _mm256_set1_ps(v.y)
	};
}


static inline mm8_vec2 mm8_sub_v2(mm8_vec2 a, mm8_vec2 b) {
	return (mm8_vec2) {
		.x = _mm256_sub_ps(a.x, b.x),
		.y = _mm256_sub_ps(a.y, b.y)
	};
}

static inline mm8_vec2 mm8_add_v2(mm8_vec2 a, mm8_vec2 b) {
	return (mm8_vec2) {
		.x = _mm256_add_ps(a.x, b.x),
		.y = _mm256_add_ps(a.y, b.y)
	};
}

static inline mm8_vec2 mm8_mul_v2(mm8_vec2 a, mm8_vec2 b) {
	return (mm8_vec2) {
		.x = _mm256_mul_ps(a.x, b.x),
		.y = _mm256_mul_ps(a.y, b.y)
	};
}

static inline mm8_vec2 mm8_div_v2(mm8_vec2 a, mm8_vec2 b) {
	return (mm8_vec2) {
		.x = _mm256_mul_ps(a.x, b.x),
		.y = _mm256_mul_ps(a.y, b.y)
	};
}

static inline mm8_vec2 mm8_scale_v2(mm8_vec2 v, __m256 s) {
	 return (mm8_vec2) {
		.x = _mm256_mul_ps(v.x, s),
		.y = _mm256_mul_ps(v.y, s)
	};
}



static inline __m256 mm8_dot_v2(mm8_vec2 a, mm8_vec2 b) {
	return _mm256_add_ps(
		_mm256_mul_ps(a.x, b.x), 
		_mm256_mul_ps(a.y, b.y)
	); 
}

static inline __m256 mm8_len_v2(mm8_vec2 v) {
	return _mm256_sqrt_ps(mm8_dot_v2(v, v));
}

static inline mm8_vec2 mm8_norm_v2(mm8_vec2 v) {
    __m256 len = mm8_len_v2(v);
    return (mm8_vec2) {
		.x = _mm256_div_ps(v.x, len),
	    .y = _mm256_div_ps(v.y, len)
	};
}

static inline __m256 mm8_dist_v2(mm8_vec2 a, mm8_vec2 b) {
	return mm8_len_v2(mm8_sub_v2(a, b));
}

static inline mm8_vec2 mm256_normDir_v2(mm8_vec2 to, mm8_vec2 from) {
	return mm8_norm_v2(mm8_sub_v2(to, from));
}





#endif // C3DLAS_USE_AVX




#endif // __c3dlas__simd_h__
