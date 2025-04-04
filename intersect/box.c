
AABB3 boxUnion(AABB3 a, AABB3 b) {
	return (AABB3) {
		{fmin(a.min.x, b.min.x), fmin(a.min.y, b.min.y), fmin(a.min.z, b.min.z)},
		{fmax(a.max.x, b.max.x), fmax(a.max.y, b.max.y), fmax(a.max.z, b.max.z)}
	};
}

// this version has no branching, but only answers yes or no.
// algorithm explanation here. hopefully my extrapolation into 3 dimensions is correct.
// http://tavianator.com/fast-branchless-raybounding-box-intersections/
int boxRayIntersectFast3p(const AABB3* b, const Ray3* r) {
	Vector3 t1, t2;
	float tmin, tmax;
	Vector3 id;
	
	vInv3p(&r->d, &id);
	
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
	
	vInv2p(&r->d, &id);
	
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
	
	vInv3p(&r->d, &id);
		
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
	
	if(tmax < tmin) return C3DLAS_DISJOINT;
	
	if(idist) *idist = tmin;
	
	if(ipoint) {
		ipoint->x = r->o.x + (r->d.x * tmin);
		ipoint->y = r->o.y + (r->d.y * tmin);
		ipoint->z = r->o.z + (r->d.z * tmin);
	}
	
	return C3DLAS_INTERSECT;
}



// this version gives the point of intersection as well as distance
// algorithm explanation here. hopefully my extrapolation into 3 dimensions is correct.
// http://tavianator.com/fast-branchless-raybounding-box-intersections/
int intersectBoxLine3p(const AABB3* b, const Line3* l, Vector3* ipoint, float* idist) {
	Vector3 t1, t2, id;
	float tmin, tmax;
	
	Vector3 o = l->start;
	Vector3 d = vNorm3(vSub3(l->end, l->start));  
	id = vInv3(d);
		
	t1.x = (b->min.x - o.x) * id.x;
	t2.x = (b->max.x - o.x) * id.x;
	tmin = fminf(t1.x, t2.x);
	tmax = fmaxf(t1.x, t2.x);
	
	t1.y = (b->min.y - o.y) * id.y;
	t2.y = (b->max.y - o.y) * id.y;
	tmin = fmaxf(tmin, fminf(t1.y, t2.y));
	tmax = fminf(tmax, fmaxf(t1.y, t2.y));
	
	t1.z = (b->min.z - o.z) * id.z;
	t2.z = (b->max.z - o.z) * id.z;
	tmin = fmaxf(tmin, fminf(t1.z, t2.z));
	tmax = fminf(tmax, fmaxf(t1.z, t2.z));
	
	if(tmax < tmin) return C3DLAS_DISJOINT;
	
	if(idist) *idist = tmin;
	
	if(ipoint) {
		ipoint->x = o.x + (d.x * tmin);
		ipoint->y = o.y + (d.y * tmin);
		ipoint->z = o.z + (d.z * tmin);
	}
	
	return C3DLAS_INTERSECT;
}

int intersectBoxLine3(AABB3 b, Line3 l, Vector3* ipoint, float* idist) {
	return intersectBoxLine3p(&b, &l, ipoint, idist);
}


bool boxContainsPoint3(AABB3 b, Vector3 p) {
	return b.min.x <= p.x && b.max.x >= p.x
		&& b.min.y <= p.y && b.max.y >= p.y
		&& b.min.z <= p.z && b.max.z >= p.z;
}


int boxClipRay(AABB3* b, Ray3 r, Line3* out) {
	
	vec3 norms[] = {
		{-1,0,0},
		{0,-1,0},
		{0,0,-1},
		{1,0,0},
		{0,1,0},
		{0,0,1},
	};
	
	f32 ds[] = {
		b->max.x,
		b->max.y,
		b->max.z,
		-b->min.x,
		-b->min.y,
		-b->min.z,
	}; 
	
	f32 mint = FLT_MAX;
	f32 maxt = 0;
	
	for(int i = 0; i < 6; i++) {
	
		float d = vDot3(norms[i], r.d);
	
		if(fabs(d) < FLT_CMP_EPSILON) continue; // parallel
		
		float t = (vDot3(norms[i], r.o) + ds[i]) / -d; 
		
		vec3 pt = vAdd3(r.o, vScale3(r.d, t));
		switch(i) { // clamp to the other planes
			case 0:
			case 3:
				if(pt.y > b->max.y || pt.y < b->min.y) continue;
				if(pt.z > b->max.z || pt.z < b->min.z) continue;
				break;
				
			case 1:
			case 4:
				if(pt.x > b->max.x || pt.x < b->min.x) continue;
				if(pt.z > b->max.z || pt.z < b->min.z) continue;
				break;
				
			case 2:
			case 5:
				if(pt.x > b->max.x || pt.x < b->min.x) continue;
				if(pt.y > b->max.y || pt.y < b->min.y) continue;
				break;
		}
		
		mint = fminf(mint, t);
		maxt = fmaxf(maxt, t);
	}
	
	mint = fmaxf(0, mint); // clamp to ray origin
	
	
	if(maxt != 0) {
		out->a = vAdd3(r.o, vScale3(r.d, mint));
		out->b = vAdd3(r.o, vScale3(r.d, maxt));
	
		return C3DLAS_INTERSECT;
	}
	else {
		return C3DLAS_DISJOINT;
	}
}


float distPoint2Rect2(Vector2 a, Quad2 q) {
	
	// check the lines around the outside
	float d[4];
	
	for(int i = 0; i < 4; i++) {
		d[i] = vDistPointLine2(a, (Line2){q.v[i], q.v[(i + 1) % 4]});
		if(d[i] == 0) return 0;
	}
	
	// check if the line is entirely inside the rectangle,
	//   if one point (and therefore both points) is inside the quad 
	vec2 ab = vSub(q.v[1], q.v[0]);
	vec2 am = vSub(a, q.v[0]);
	
	vec2 bc = vSub(q.v[2], q.v[1]);
	vec2 bm = vSub(a, q.v[1]);
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


int intersectPoint2Rect2(vec2 a, Quad2 q) {
	
	vec2 ab = vSub(q.v[1], q.v[0]);
	vec2 am = vSub(a, q.v[0]);
	
	vec2 bc = vSub(q.v[2], q.v[1]);
	vec2 bm = vSub(a, q.v[1]);
	float dabam = vDot(ab, am);
	if(dabam <= 0) {
		float dabab = vDot(ab, ab);
		if(dabam <= dabab) {
			float dbcbm = vDot(bc, bm);
			if(0 <= dbcbm) {
				float dbcbc = vDot(bc, bc);
				if(dbcbm <= dbcbc) return C3DLAS_INTERSECT;
			}
		}
	}

	return C3DLAS_DISJOINT;
}


int intersectRect2Rect2(Quad2 a, Quad2 b) {
	
	// early rejection test based on distance between centers
	vec2 ca = {0, 0};
	vec2 cb = {0, 0};
	
	for(int i = 0; i < 4; i++) {
		ca = vAdd(ca, a.v[i]);
		cb = vAdd(cb, b.v[i]);
	}
	ca = vScale(ca, .25f);
	cb = vScale(cb, .25f);
	
	float radSq = vDistSq(ca, a.v[0]) + vDistSq(cb, b.v[0]);
	if(radSq > vDistSq(ca, cb)) return C3DLAS_DISJOINT;
	
	
	// Check if any of the points are inside the other rectangle
	// This handles situations where one rect is fully contained in the other
	for(int i = 0; i < 4; i++) {
		if(C3DLAS_INTERSECT == intersectPoint2Rect2(a.v[i], b)) return C3DLAS_INTERSECT;
		if(C3DLAS_INTERSECT == intersectPoint2Rect2(b.v[i], a)) return C3DLAS_INTERSECT;
	}
	
	
	// check all the lines against all the other lines
	for(int i = 0; i < 4; i++)
	for(int j = 0; j < 4; j++) {
		if(intersectLine2Line2((Line2){a.v[i], a.v[(i + 1) % 4]}, (Line2){b.v[j], b.v[(j + 1) % 4]}) == C3DLAS_INTERSECT) return C3DLAS_INTERSECT;
	}
	
	return C3DLAS_DISJOINT;
}



