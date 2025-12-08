


void evalBezier3p(Vector3* e1, Vector3* e2, Vector3* c1, Vector3* c2, float t, Vector3* out) {
	out->x = evalBezier1D(e1->x, e2->x, c1->x, c2->x, t);
	out->y = evalBezier1D(e1->y, e2->y, c1->y, c2->y, t);
	out->z = evalBezier1D(e1->z, e2->z, c1->z, c2->z, t);
}

Vector2 evalBezier2D(Vector2 e1, Vector2 e2, Vector2 c1, Vector2 c2, float t) {
	return (Vector2){
		.x = evalBezier1D(e1.x, e2.x, c1.x, c2.x, t),
		.y = evalBezier1D(e1.y, e2.y, c1.y, c2.y, t)
	};
}
Vector2 evalBezier2D_dt(Vector2 e1, Vector2 e2, Vector2 c1, Vector2 c2, float t) {
	return (Vector2){
		.x = evalBezier1D_dt(e1.x, e2.x, c1.x, c2.x, t),
		.y = evalBezier1D_dt(e1.y, e2.y, c1.y, c2.y, t)
	};
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

void evalQBezier2Dp(Vector2* e1, Vector2* e2, Vector2* c1, float t, Vector2* out) {
	out->x = evalQBezier1D(e1->x, e2->x, c1->x, t);
	out->y = evalQBezier1D(e1->y, e2->y, c1->y, t);
}

Vector2 evalQBezier2D(Vector2 e1, Vector2 e2, Vector2 c1, float t) {
	return (Vector2){
		.x = evalQBezier1D(e1.x, e2.x, c1.x, t),
		.y = evalQBezier1D(e1.y, e2.y, c1.y, t)
	};
}

void evalQBezier3Dp(Vector3* e1, Vector3* e2, Vector3* c1, float t, Vector3* out) {
	out->x = evalQBezier1D(e1->x, e2->x, c1->x, t);
	out->y = evalQBezier1D(e1->y, e2->y, c1->y, t);
	out->z = evalQBezier1D(e1->z, e2->z, c1->z, t);
}

Vector3 evalQBezier3D(Vector3 e1, Vector3 e2, Vector3 c1, float t) {
	return (Vector3){
		.x = evalQBezier1D(e1.x, e2.x, c1.x, t),
		.y = evalQBezier1D(e1.y, e2.y, c1.y, t),
		.z = evalQBezier1D(e1.z, e2.z, c1.z, t)
	};
}


float evalQBezier1D_dt(float e1, float e2, float c1, float t) {
	return 2.f * (-e1 + t*e1 + 2.f*t*c1 - 3*t*t*c1 + t*e2);	
}

Vector2 evalQBezier2D_dt(Vector2 e1, Vector2 e2, Vector2 c1, float t) {
	return (Vector2){
		.x = evalQBezier1D_dt(e1.x, e2.x, c1.x, t),
		.y = evalQBezier1D_dt(e1.y, e2.y, c1.y, t)
	};
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
	
	*out = evalBezier2D(cp[0], cp[3], cp[1], cp[2], localT);
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




