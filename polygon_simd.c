




#ifdef C3DLAS_USE_AVX2

// like the scalar one, but logical-OR's all the points' results
int mm8_polyCheckMaybeWithinRadius_group(Polygon* poly, float px[8], float py[8], float radius) {
	mm8_vec2 p_1 = mm8_load_arrays_v2(px, py);
	
	mm8_vec2 centroid = mm8_set1_v2(poly->centroid);
	
	__m256 dist_1 = mm8_dist_v2(p_1, centroid);
	
	__m256 rad = _mm256_set1_ps(sqrtf(poly->maxRadiusSq) + radius);
	
	__m256 res = _mm256_cmp_ps(dist_1, rad, _CMP_GT_OQ);
	
	unsigned int r = *(unsigned int*)&res;
	for(int i = 1; i < 8; i++) {
		r &= ((unsigned int*)&res)[i];
	}
	
	return !!r;
}


// like the scalar one, but logical-OR's all the points' results
int mm16_polyCheckMaybeWithinRadius_group(Polygon* poly, float px[16], float py[16], float radius) {
	mm8_vec2 p_1 = mm8_load_arrays_v2(px, py);
	mm8_vec2 p_2 = mm8_load_arrays_v2(px + 8, py + 8);
	
	mm8_vec2 centroid = mm8_set1_v2(poly->centroid);
	
	__m256 dist_1 = mm8_dist_v2(p_1, centroid);
	__m256 dist_2 = mm8_dist_v2(p_2, centroid);
	
	__m256 rad = _mm256_set1_ps(sqrtf(poly->maxRadiusSq) + radius);
	
	__m256 res_1 = _mm256_cmp_ps(dist_1, rad, _CMP_GT_OQ);
	__m256 res_2 = _mm256_cmp_ps(dist_2, rad, _CMP_GT_OQ);
	
	__m256 res = _mm256_and_ps(res_1, res_2);
	
	
	unsigned int r = *(unsigned int*)&res;
	for(int i = 1; i < 8; i++) {
		r &= ((unsigned int*)&res)[i];
	}
	
	return !!r;
}


// like the scalar one, but logical-OR's all the points' results
int mm32_polyCheckMaybeWithinRadius_group(Polygon* poly, float px[32], float py[32], float radius) {
	mm8_vec2 p_1 = mm8_load_arrays_v2(px, py);
	mm8_vec2 p_2 = mm8_load_arrays_v2(px + 8, py + 8);
	mm8_vec2 p_3 = mm8_load_arrays_v2(px + 16, py + 16);
	mm8_vec2 p_4 = mm8_load_arrays_v2(px + 24, py + 24);
	
	mm8_vec2 centroid = mm8_set1_v2(poly->centroid);
	
	__m256 dist_1 = mm8_dist_v2(p_1, centroid);
	__m256 dist_2 = mm8_dist_v2(p_2, centroid);
	__m256 dist_3 = mm8_dist_v2(p_3, centroid);
	__m256 dist_4 = mm8_dist_v2(p_4, centroid);
	
	__m256 rad = _mm256_set1_ps(sqrtf(poly->maxRadiusSq) + radius);
	
	__m256 res_1 = _mm256_cmp_ps(dist_1, rad, _CMP_GT_OQ);
	__m256 res_2 = _mm256_cmp_ps(dist_2, rad, _CMP_GT_OQ);
	__m256 res_3 = _mm256_cmp_ps(dist_3, rad, _CMP_GT_OQ);
	__m256 res_4 = _mm256_cmp_ps(dist_4, rad, _CMP_GT_OQ);
	
	__m256 res = _mm256_and_ps(_mm256_and_ps(res_1, res_2), _mm256_and_ps(res_1, res_2));
	
	
	unsigned int r = *(unsigned int*)&res;
	for(int i = 1; i < 8; i++) {
		r &= ((unsigned int*)&res)[i];
	}
	
	return !!r;
}
#else //#ifdef C3DLAS_USE_AVX2

int mm8_polyCheckMaybeWithinRadius_group(Polygon* poly, float px[8], float py[8], float radius) {
	int r = 0;
	for(int i = 0; i < 8; i++) {
		r &= polyCheckMaybeWithinRadius(poly, (Vector2){px[i], py[i]}, radius);
	}
	return r;
}
int mm16_polyCheckMaybeWithinRadius_group(Polygon* poly, float px[16], float py[16], float radius) {
	int r = 0;
	for(int i = 0; i < 16; i++) {
		r &= polyCheckMaybeWithinRadius(poly, (Vector2){px[i], py[i]}, radius);
	}
	return r;
}
int mm32_polyCheckMaybeWithinRadius_group(Polygon* poly, float px[32], float py[32], float radius) {
	int r = 0;
	for(int i = 0; i < 32; i++) {
		r &= polyCheckMaybeWithinRadius(poly, (Vector2){px[i], py[i]}, radius);
	}
	return r;
}

#endif //#else #ifndef C3DLAS_USE_AVX2








#ifdef C3DLAS_USE_AVX2

//                                                  (percent slower)
// mapgen.c:670 avx32: 57.774661 ms, avx16: 63.600085 ms (10.1%), avx8: 76.575981 ms (32.5%), scalar: 552.618354 ms (856.5%)

void mm8_polyDistToPoint(Polygon* poly, f32 px[8], f32 py[8], f32 out[8]) {

	mm8_vec2 p = mm8_load_arrays_v2(px, py);
	mm8_vec2 pp0 = mm8_set1_v2(poly->points[0]);
	
	__m256 d = mm8_dot_v2(mm8_sub_v2(p, pp0), mm8_sub_v2(p, pp0));
	__m256 s = _mm256_set1_ps(1.0);
	__m256 one = _mm256_set1_ps(1.0);
	__m256 zero = _mm256_set1_ps(0.0);
	__m256 negOne = _mm256_set1_ps(-1.0);
	__m256 negZero = _mm256_set1_ps(-.0f);

	for(long i = 0, j = poly->pointCount - 1; i < poly->pointCount; j = i, i++) {
		
		mm8_vec2 A = mm8_set1_v2(poly->points[i]);
		mm8_vec2 B = mm8_set1_v2(poly->points[j]);

		mm8_vec2 e = mm8_sub_v2(B, A);
		mm8_vec2 w = mm8_sub_v2(p, A);
		mm8_vec2 b = mm8_sub_v2(w, mm8_scale_v2(e, mm256_clamp_scalar_ps(mm8_dot_v2(w, e) / mm8_dot_v2(e, e), 0.0, 1.0)));

		d = _mm256_min_ps(d, mm8_dot_v2(b, b));

		__m256 c1 = _mm256_cmp_ps(p.y, A.y, _CMP_GE_OQ);
		__m256 c2 = _mm256_cmp_ps(p.y, B.y, _CMP_LT_OQ);
		__m256 c3 = _mm256_cmp_ps(_mm256_mul_ps(e.x, w.y), _mm256_mul_ps(e.y, w.x), _CMP_GT_OQ);
		
		__m256 r1 = _mm256_and_ps(c1, _mm256_and_ps(c2, c3));
		__m256 r2 = _mm256_andnot_ps(c1, _mm256_andnot_ps(c2, mm256_not_ps(c3)));
		__m256 r3 = _mm256_or_ps(r1, r2);		
		
		s = _mm256_mul_ps(s, _mm256_blendv_ps(one, negOne, r3));
	}

	_mm256_storeu_ps(
		out,
		_mm256_blendv_ps(
			_mm256_mul_ps(s, _mm256_sqrt_ps(d)),
			zero,
			_mm256_cmp_ps(d, negZero, _CMP_EQ_OQ)
		)
	);
}


void mm16_polyDistToPoint(Polygon* poly, f32 px[16], f32 py[16], f32 out[16]) {

	mm8_vec2 p_1 = mm8_load_arrays_v2(px, py);
	mm8_vec2 p_2 = mm8_load_arrays_v2(px + 8, py + 8);
	mm8_vec2 pp0 = mm8_set1_v2(poly->points[0]);
	
	__m256 d_1 = mm8_dot_v2(mm8_sub_v2(p_1, pp0), mm8_sub_v2(p_1, pp0));
	__m256 d_2 = mm8_dot_v2(mm8_sub_v2(p_2, pp0), mm8_sub_v2(p_2, pp0));
	__m256 s_1 = _mm256_set1_ps(1.0);
	__m256 s_2 = _mm256_set1_ps(1.0);
	__m256 one = _mm256_set1_ps(1.0);
	__m256 zero = _mm256_set1_ps(0.0);
	__m256 negOne = _mm256_set1_ps(-1.0);
	__m256 negZero = _mm256_set1_ps(-.0f);

	for(long i = 0, j = poly->pointCount - 1; i < poly->pointCount; j = i, i++) {
		
		mm8_vec2 A = mm8_set1_v2(poly->points[i]);
		mm8_vec2 B = mm8_set1_v2(poly->points[j]);

		mm8_vec2 e = mm8_sub_v2(B, A);
		__m256 dote = mm8_dot_v2(e, e);
		mm8_vec2 w_1 = mm8_sub_v2(p_1, A);
		mm8_vec2 w_2 = mm8_sub_v2(p_2, A);
		
		mm8_vec2 b_1 = mm8_sub_v2(w_1, mm8_scale_v2(e, mm256_clamp_scalar_ps(mm8_dot_v2(w_1, e) / dote, 0.0, 1.0)));
		mm8_vec2 b_2 = mm8_sub_v2(w_2, mm8_scale_v2(e, mm256_clamp_scalar_ps(mm8_dot_v2(w_2, e) / dote, 0.0, 1.0)));

		d_1 = _mm256_min_ps(d_1, mm8_dot_v2(b_1, b_1));
		d_2 = _mm256_min_ps(d_2, mm8_dot_v2(b_2, b_2));

		__m256 c1_1 = _mm256_cmp_ps(p_1.y, A.y, _CMP_GE_OQ);
		__m256 c1_2 = _mm256_cmp_ps(p_2.y, A.y, _CMP_GE_OQ);
		__m256 c2_1 = _mm256_cmp_ps(p_1.y, B.y, _CMP_LT_OQ);
		__m256 c2_2 = _mm256_cmp_ps(p_2.y, B.y, _CMP_LT_OQ);
		__m256 c3_1 = _mm256_cmp_ps(_mm256_mul_ps(e.x, w_1.y), _mm256_mul_ps(e.y, w_1.x), _CMP_GT_OQ);
		__m256 c3_2 = _mm256_cmp_ps(_mm256_mul_ps(e.x, w_2.y), _mm256_mul_ps(e.y, w_2.x), _CMP_GT_OQ);
		
		__m256 r1_1 = _mm256_and_ps(c1_1, _mm256_and_ps(c2_1, c3_1));
		__m256 r1_2 = _mm256_and_ps(c1_2, _mm256_and_ps(c2_2, c3_2));
		__m256 r2_1 = _mm256_andnot_ps(c1_1, _mm256_andnot_ps(c2_1, mm256_not_ps(c3_1)));
		__m256 r2_2 = _mm256_andnot_ps(c1_2, _mm256_andnot_ps(c2_2, mm256_not_ps(c3_2)));
		__m256 r3_1 = _mm256_or_ps(r1_1, r2_1);		
		__m256 r3_2 = _mm256_or_ps(r1_2, r2_2);		
		
		
		s_1 = _mm256_mul_ps(s_1, _mm256_blendv_ps(one, negOne, r3_1));
		s_2 = _mm256_mul_ps(s_2, _mm256_blendv_ps(one, negOne, r3_2));
	}

	_mm256_storeu_ps(
		out,
		_mm256_blendv_ps(
			_mm256_mul_ps(s_1, _mm256_sqrt_ps(d_1)),
			zero,
			_mm256_cmp_ps(d_1, negZero, _CMP_EQ_OQ)
		)
	);
	_mm256_storeu_ps(
		out + 8,
		_mm256_blendv_ps(
			_mm256_mul_ps(s_2, _mm256_sqrt_ps(d_2)),
			zero,
			_mm256_cmp_ps(d_2, negZero, _CMP_EQ_OQ)
		)
	);
}

void mm16_2_polyDistToPoint(Polygon* poly, f32 px[16], f32 py[16], f32 out[16]) {

	mm8_vec2 p_1 = mm8_load_arrays_v2(px, py);
	mm8_vec2 p_2 = mm8_load_arrays_v2(px + 8, py + 8);
	mm8_vec2 pp0 = mm8_set1_v2(poly->points[0]);
	
	__m256 d_1 = mm8_dot_v2(mm8_sub_v2(p_1, pp0), mm8_sub_v2(p_1, pp0));
	__m256 d_2 = mm8_dot_v2(mm8_sub_v2(p_2, pp0), mm8_sub_v2(p_2, pp0));
	__m256 s_1 = _mm256_set1_ps(1.0);
	__m256 s_2 = _mm256_set1_ps(1.0);
	__m256 one = _mm256_set1_ps(1.0);
	__m256 zero = _mm256_set1_ps(0.0);
	__m256 negOne = _mm256_set1_ps(-1.0);
	__m256 negZero = _mm256_set1_ps(-.0f);

	for(long i = 0, j = poly->pointCount - 1; i < poly->pointCount; j = i, i++) {
		
		mm8_vec2 A = mm8_set1_v2(poly->points[i]);
		mm8_vec2 B = mm8_set1_v2(poly->points[j]);

		mm8_vec2 e = mm8_sub_v2(B, A);
		__m256 dote = mm8_dot_v2(e, e);
		
		
		//------
		
		mm8_vec2 w_1 = mm8_sub_v2(p_1, A);		
		mm8_vec2 b_1 = mm8_sub_v2(w_1, mm8_scale_v2(e, mm256_clamp_scalar_ps(mm8_dot_v2(w_1, e) / dote, 0.0, 1.0)));
		d_1 = _mm256_min_ps(d_1, mm8_dot_v2(b_1, b_1));

		__m256 c1_1 = _mm256_cmp_ps(p_1.y, A.y, _CMP_GE_OQ);
		__m256 c2_1 = _mm256_cmp_ps(p_1.y, B.y, _CMP_LT_OQ);
		__m256 c3_1 = _mm256_cmp_ps(_mm256_mul_ps(e.x, w_1.y), _mm256_mul_ps(e.y, w_1.x), _CMP_GT_OQ);

		__m256 r1_1 = _mm256_and_ps(c1_1, _mm256_and_ps(c2_1, c3_1));
		__m256 r2_1 = _mm256_andnot_ps(c1_1, _mm256_andnot_ps(c2_1, mm256_not_ps(c3_1)));
		__m256 r3_1 = _mm256_or_ps(r1_1, r2_1);		

		s_1 = _mm256_mul_ps(s_1, _mm256_blendv_ps(one, negOne, r3_1));
	
		//------

		mm8_vec2 w_2 = mm8_sub_v2(p_2, A);
		mm8_vec2 b_2 = mm8_sub_v2(w_2, mm8_scale_v2(e, mm256_clamp_scalar_ps(mm8_dot_v2(w_2, e) / dote, 0.0, 1.0)));
		d_2 = _mm256_min_ps(d_2, mm8_dot_v2(b_2, b_2));

		__m256 c1_2 = _mm256_cmp_ps(p_2.y, A.y, _CMP_GE_OQ);
		__m256 c2_2 = _mm256_cmp_ps(p_2.y, B.y, _CMP_LT_OQ);
		__m256 c3_2 = _mm256_cmp_ps(_mm256_mul_ps(e.x, w_2.y), _mm256_mul_ps(e.y, w_2.x), _CMP_GT_OQ);
		
		__m256 r1_2 = _mm256_and_ps(c1_2, _mm256_and_ps(c2_2, c3_2));
		__m256 r2_2 = _mm256_andnot_ps(c1_2, _mm256_andnot_ps(c2_2, mm256_not_ps(c3_2)));
		__m256 r3_2 = _mm256_or_ps(r1_2, r2_2);		
		
		s_2 = _mm256_mul_ps(s_2, _mm256_blendv_ps(one, negOne, r3_2));
	}

	_mm256_storeu_ps(
		out,
		_mm256_blendv_ps(
			_mm256_mul_ps(s_1, _mm256_sqrt_ps(d_1)),
			zero,
			_mm256_cmp_ps(d_1, negZero, _CMP_EQ_OQ)
		)
	);
	_mm256_storeu_ps(
		out + 8,
		_mm256_blendv_ps(
			_mm256_mul_ps(s_2, _mm256_sqrt_ps(d_2)),
			zero,
			_mm256_cmp_ps(d_2, negZero, _CMP_EQ_OQ)
		)
	);
}

void mm32_polyDistToPoint(Polygon* poly, f32 px[32], f32 py[32], f32 out[32]) {

	mm8_vec2 p_1 = mm8_load_arrays_v2(px, py);
	mm8_vec2 p_2 = mm8_load_arrays_v2(px + 8, py + 8);
	mm8_vec2 p_3 = mm8_load_arrays_v2(px + 16, py + 16);
	mm8_vec2 p_4 = mm8_load_arrays_v2(px + 24, py + 24);
	mm8_vec2 pp0 = mm8_set1_v2(poly->points[0]);
	
	__m256 d_1 = mm8_dot_v2(mm8_sub_v2(p_1, pp0), mm8_sub_v2(p_1, pp0));
	__m256 d_2 = mm8_dot_v2(mm8_sub_v2(p_2, pp0), mm8_sub_v2(p_2, pp0));
	__m256 d_3 = mm8_dot_v2(mm8_sub_v2(p_3, pp0), mm8_sub_v2(p_3, pp0));
	__m256 d_4 = mm8_dot_v2(mm8_sub_v2(p_4, pp0), mm8_sub_v2(p_4, pp0));
	__m256 s_1 = _mm256_set1_ps(1.0);
	__m256 s_2 = _mm256_set1_ps(1.0);
	__m256 s_3 = _mm256_set1_ps(1.0);
	__m256 s_4 = _mm256_set1_ps(1.0);
	__m256 one = _mm256_set1_ps(1.0);
	__m256 zero = _mm256_set1_ps(0.0);
	__m256 negOne = _mm256_set1_ps(-1.0);
	__m256 negZero = _mm256_set1_ps(-.0f);

	for(long i = 0, j = poly->pointCount - 1; i < poly->pointCount; j = i, i++) {
		
		mm8_vec2 A = mm8_set1_v2(poly->points[i]);
		mm8_vec2 B = mm8_set1_v2(poly->points[j]);

		mm8_vec2 e = mm8_sub_v2(B, A);
		__m256 dote = mm8_dot_v2(e, e);
		mm8_vec2 w_1 = mm8_sub_v2(p_1, A);
		mm8_vec2 w_2 = mm8_sub_v2(p_2, A);
		mm8_vec2 w_3 = mm8_sub_v2(p_3, A);
		mm8_vec2 w_4 = mm8_sub_v2(p_4, A);
		
		mm8_vec2 b_1 = mm8_sub_v2(w_1, mm8_scale_v2(e, mm256_clamp_scalar_ps(mm8_dot_v2(w_1, e) / dote, 0.0, 1.0)));
		mm8_vec2 b_2 = mm8_sub_v2(w_2, mm8_scale_v2(e, mm256_clamp_scalar_ps(mm8_dot_v2(w_2, e) / dote, 0.0, 1.0)));
		mm8_vec2 b_3 = mm8_sub_v2(w_3, mm8_scale_v2(e, mm256_clamp_scalar_ps(mm8_dot_v2(w_3, e) / dote, 0.0, 1.0)));
		mm8_vec2 b_4 = mm8_sub_v2(w_4, mm8_scale_v2(e, mm256_clamp_scalar_ps(mm8_dot_v2(w_4, e) / dote, 0.0, 1.0)));

		d_1 = _mm256_min_ps(d_1, mm8_dot_v2(b_1, b_1));
		d_2 = _mm256_min_ps(d_2, mm8_dot_v2(b_2, b_2));
		d_3 = _mm256_min_ps(d_3, mm8_dot_v2(b_3, b_3));
		d_4 = _mm256_min_ps(d_4, mm8_dot_v2(b_4, b_4));

		__m256 c1_1 = _mm256_cmp_ps(p_1.y, A.y, _CMP_GE_OQ);
		__m256 c1_2 = _mm256_cmp_ps(p_2.y, A.y, _CMP_GE_OQ);
		__m256 c1_3 = _mm256_cmp_ps(p_3.y, A.y, _CMP_GE_OQ);
		__m256 c1_4 = _mm256_cmp_ps(p_4.y, A.y, _CMP_GE_OQ);
		__m256 c2_1 = _mm256_cmp_ps(p_1.y, B.y, _CMP_LT_OQ);
		__m256 c2_2 = _mm256_cmp_ps(p_2.y, B.y, _CMP_LT_OQ);
		__m256 c2_3 = _mm256_cmp_ps(p_3.y, B.y, _CMP_LT_OQ);
		__m256 c2_4 = _mm256_cmp_ps(p_4.y, B.y, _CMP_LT_OQ);
		__m256 c3_1 = _mm256_cmp_ps(_mm256_mul_ps(e.x, w_1.y), _mm256_mul_ps(e.y, w_1.x), _CMP_GT_OQ);
		__m256 c3_2 = _mm256_cmp_ps(_mm256_mul_ps(e.x, w_2.y), _mm256_mul_ps(e.y, w_2.x), _CMP_GT_OQ);
		__m256 c3_3 = _mm256_cmp_ps(_mm256_mul_ps(e.x, w_3.y), _mm256_mul_ps(e.y, w_3.x), _CMP_GT_OQ);
		__m256 c3_4 = _mm256_cmp_ps(_mm256_mul_ps(e.x, w_4.y), _mm256_mul_ps(e.y, w_4.x), _CMP_GT_OQ);
		
		__m256 r1_1 = _mm256_and_ps(c1_1, _mm256_and_ps(c2_1, c3_1));
		__m256 r1_2 = _mm256_and_ps(c1_2, _mm256_and_ps(c2_2, c3_2));
		__m256 r1_3 = _mm256_and_ps(c1_3, _mm256_and_ps(c2_3, c3_3));
		__m256 r1_4 = _mm256_and_ps(c1_4, _mm256_and_ps(c2_4, c3_4));
		__m256 r2_1 = _mm256_andnot_ps(c1_1, _mm256_andnot_ps(c2_1, mm256_not_ps(c3_1)));
		__m256 r2_2 = _mm256_andnot_ps(c1_2, _mm256_andnot_ps(c2_2, mm256_not_ps(c3_2)));
		__m256 r2_3 = _mm256_andnot_ps(c1_3, _mm256_andnot_ps(c2_3, mm256_not_ps(c3_3)));
		__m256 r2_4 = _mm256_andnot_ps(c1_4, _mm256_andnot_ps(c2_4, mm256_not_ps(c3_4)));
		__m256 r3_1 = _mm256_or_ps(r1_1, r2_1);
		__m256 r3_2 = _mm256_or_ps(r1_2, r2_2);
		__m256 r3_3 = _mm256_or_ps(r1_3, r2_3);
		__m256 r3_4 = _mm256_or_ps(r1_4, r2_4);
		
		
		s_1 = _mm256_mul_ps(s_1, _mm256_blendv_ps(one, negOne, r3_1));
		s_2 = _mm256_mul_ps(s_2, _mm256_blendv_ps(one, negOne, r3_2));
		s_3 = _mm256_mul_ps(s_3, _mm256_blendv_ps(one, negOne, r3_3));
		s_4 = _mm256_mul_ps(s_4, _mm256_blendv_ps(one, negOne, r3_4));
	}

	_mm256_storeu_ps(
		out,
		_mm256_blendv_ps(
			_mm256_mul_ps(s_1, _mm256_sqrt_ps(d_1)),
			zero,
			_mm256_cmp_ps(d_1, negZero, _CMP_EQ_OQ)
		)
	);
	_mm256_storeu_ps(
		out + 8,
		_mm256_blendv_ps(
			_mm256_mul_ps(s_2, _mm256_sqrt_ps(d_2)),
			zero,
			_mm256_cmp_ps(d_2, negZero, _CMP_EQ_OQ)
		)
	);
	_mm256_storeu_ps(
		out + 16,
		_mm256_blendv_ps(
			_mm256_mul_ps(s_3, _mm256_sqrt_ps(d_3)),
			zero,
			_mm256_cmp_ps(d_3, negZero, _CMP_EQ_OQ)
		)
	);
	_mm256_storeu_ps(
		out + 24,
		_mm256_blendv_ps(
			_mm256_mul_ps(s_4, _mm256_sqrt_ps(d_4)),
			zero,
			_mm256_cmp_ps(d_4, negZero, _CMP_EQ_OQ)
		)
	);
}


#else // #ifdef C3DLAS_USE_AVX2


void mm8_polyDistToPoint(Polygon* poly, f32 px[8], f32 py[8], f32 out[8]) {
	for(int i = 0; i < 8; i++) {
		out[i] = polyDistToPoint(poly, (Vector2){px[i], py[i]});
	}
}
void mm16_polyDistToPoint(Polygon* poly, f32 px[16], f32 py[16], f32 out[16]) {
	for(int i = 0; i < 15; i++) {
		out[i] = polyDistToPoint(poly, (Vector2){px[i], py[i]});
	}
}
void mm32_polyDistToPoint(Polygon* poly, f32 px[32], f32 py[32], f32 out[32]) {
	for(int i = 0; i < 32; i++) {
		out[i] = polyDistToPoint(poly, (Vector2){px[i], py[i]});
	}
}


#endif //#else #ifdef C3DLAS_USE_AVX2


