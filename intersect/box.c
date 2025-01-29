
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
	
	// intercept distances
	f32 tx_1 = (b->min.x - r.o.x) / r.d.x;
	f32 ty_1 = (b->min.y - r.o.y) / r.d.y;
	f32 tz_1 = (b->min.z - r.o.z) / r.d.z;
	f32 tx_2 = (b->max.x - r.o.x) / r.d.x;
	f32 ty_2 = (b->max.y - r.o.y) / r.d.y;
	f32 tz_2 = (b->max.z - r.o.z) / r.d.z;
	
	f32 tx_min = MIN(tx_1, tx_2); 
	f32 ty_min = MIN(ty_1, ty_2); 
	f32 tz_min = MIN(tz_1, tz_2); 
	f32 tx_max = MIN(tx_1, tx_2); 
	f32 ty_max = MIN(ty_1, ty_2); 
	f32 tz_max = MIN(tz_1, tz_2); 
	
	// ALL HORRIBLY WRONG
	
	f32 min_t = -FLT_MAX;
	
	// check if the ray intersects on various axes
	if(tx_min < ty_max && tx_min < tz_max) {
		min_t = tx_min;
	}
	
	if(ty_min < tx_max && ty_min < tz_max) {
//		if(
	}
//		&& (tz_min < tx_max && tz_min < ty_max)
//	) {
//		return C3DLAS_DISJOINT;
//	}
	
//	f32 entry = MAX(tx_min, ty_min, tz_min);
//	f32 exit = MIN(tx_max, ty_max, tz_max);
	
//	if(out) {
//		out->start = vAdd(r.o, vScale(r.d, entry));
//		out->end = vAdd(r.o, vScale(r.d, exit));
//	}
	
	return C3DLAS_INTERSECT;
}



