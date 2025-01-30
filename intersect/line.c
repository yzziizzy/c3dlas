


int intersectLine2Line2(Line2 a, Line2 b) {
	
	vec2 s1 = vSub(a.end, a.start);
	vec2 s2 = vSub(b.end, b.start);
	
	float det = vCross2(s1, s2);
	if(fabsf(det) < 1e-5f) {
		return C3DLAS_PARALLEL;
	}
	
	float invdet = 1.f / det;
    
	float x02 = a.start.x - b.start.x;
	float y02 = a.start.y - b.start.y;
	
	float s = (s1.x * y02 - s1.y * x02) * invdet;
	if(s >= 0.f && s <= 1.f) {
		float t = (s2.x * y02 - s2.y * x02) * invdet;
		if(t >= 0.f && t <= 1.f) {
			return C3DLAS_INTERSECT;
		}
	}
	
	return C3DLAS_DISJOINT;
}


float distLine2Line2(Line2 a, Line2 b) {
	
	vec2 s1 = vSub(a.end, a.start);
	vec2 s2 = vSub(b.end, b.start);
	
	float det = vCross2(s1, s2);
	if(fabsf(det) < 1e-5f) { // parallel lines
		
		if(s1.x == 0) return fabs(a.start.y - b.start.y); // vertical lines
		
		float m = s1.y / s1.x;
		float c1 = a.start.y - (m * a.start.x);
		float c2 = b.start.y - (m * b.start.x);
		
		return fabs((c1 - c2) / sqrt(1 + m * m));
	}
	
	float invdet = 1.f / det;
    
	float x02 = a.start.x - b.start.x;
	float y02 = a.start.y - b.start.y;
	
	float s = (s1.x * y02 - s1.y * x02) * invdet;
	float t = (s2.x * y02 - s2.y * x02) * invdet;
	
	if(s >= 0.f && s <= 1.f && t >= 0.f && t <= 1.f) {
		return 0; // intersecting lines
	}
	
	// one or both of the ends is the closest point on its line 
	s = fclamp(s, 0, 1);
	t = fclamp(t, 0, 1);
	
	vec2 pa = vLerp(a.start, a.end, s); 
	vec2 pb = vLerp(b.start, b.end, t); 
	
	return vDist(pa, pb);
}


// Quad *must* be a rectangle
float distLine2Rect2(Line2 a, Quad2 q) {
	
	// check the lines around the outside
	float d[4];
	
	for(int i = 0; i < 4; i++) {
		d[i] = distLine2Line2(a, (Line2){q.v[i], q.v[(i + 1) % 4]});
		if(d[i] == 0) return 0;
	}
	
	// check if the line is entirely inside the rectangle,
	//   if one point (and therefore both points) is inside the quad 
	vec2 ab = vSub(q.v[1], q.v[0]);
	vec2 am = vSub(a.start, q.v[0]);
	
	vec2 bc = vSub(q.v[2], q.v[1]);
	vec2 bm = vSub(a.start, q.v[1]);
	float dabam = vDot(ab, am);
	if(dabam <= 0) {
		float dabab = vDot(ab, ab);
		if(dabam <= dabab) {
			float dbcbm = vDot(bc, bm);
			if(0 <= dbcbm) {
				float dbcbc = vDot(bc, bc);
				if(dbcbm <= dbcbc) return 0;
			}
		}
	}
	
	// the line is outside; one of the corners is closest, find which one
	return MIN(MIN(d[0], d[1]), MIN(d[2], d[3]));
}


// Quad *must* be a rectangle
int intersectLine2Rect2(Line2 a, Quad2 q) {
	
	// check if it intersects the boundaries
	if(intersectLine2Line2(a, (Line2){q.v[0], q.v[1]}) == C3DLAS_INTERSECT) return C3DLAS_INTERSECT;
	if(intersectLine2Line2(a, (Line2){q.v[1], q.v[2]}) == C3DLAS_INTERSECT) return C3DLAS_INTERSECT;
	if(intersectLine2Line2(a, (Line2){q.v[2], q.v[3]}) == C3DLAS_INTERSECT) return C3DLAS_INTERSECT;
	if(intersectLine2Line2(a, (Line2){q.v[3], q.v[0]}) == C3DLAS_INTERSECT) return C3DLAS_INTERSECT;
	
	// check if one point (and therefore both points) is inside the quad 
	vec2 ab = vSub(q.v[1], q.v[0]);
	vec2 bc = vSub(q.v[2], q.v[1]);
	vec2 am = vSub(a.start, q.v[0]);
	vec2 bm = vSub(a.start, q.v[1]);
	float dabam = vDot(ab, am);
	float dbcbm = vDot(bc, bm);
	float dabab = vDot(ab, ab);
	float dbcbc = vDot(bc, bc);
	
	if(0 <= dabam && dabam <= dabab && 0 <= dbcbm && dbcbm <= dbcbc) return C3DLAS_INTERSECT;
	
	return C3DLAS_DISJOINT; 
}



