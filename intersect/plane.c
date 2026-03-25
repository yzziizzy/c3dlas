


// closest distance from an arbitrary point to the plane 
float planePointDist3p(Plane* pl, Vector3* p) {
	Vector3 a;
	vScale3p(&pl->n, pl->d, &a);
	return fabs(vDot3(vSub3(a, *p), pl->n));
} 


// signed closest distance from an arbitrary point to the plane 
float planePointDistSigned3p(Plane* pl, Vector3* p) {
	Vector3 a;
	vScale3p(&pl->n, pl->d, &a);
	return vDot3(vSub3(a, *p), pl->n);
}


// C3DLAS_INTERSECT, _COPLANAR or _DISJOINT
int findIntersectLinePlane(Line3 l, Plane* pl, Vector3* out) {
//int planeLineFindIntersect3p(Plane* pl, Vector3* la, Vector3* lb, Vector3* out) {
	Vector3 ldir;
	float da, db;
	
	ldir = vSub3(l.b, l.a);
	
	da = vDot3(l.a, pl->n) - pl->d;
	
	// bail if the line and plane are parallel
	if(fabs(vDot3(pl->n, ldir)) < FLT_CMP_EPSILON) {
		
		// check coplanarity
		if(fabs(da) < FLT_CMP_EPSILON) {
			return C3DLAS_COPLANAR; // the end is on the plane, so the other is too
		}
		
		return C3DLAS_DISJOINT;
	}
	

	db = vDot3(l.b, pl->n) - pl->d;
	
	// check if one of the points is on the plane
	if(fabs(da) < FLT_CMP_EPSILON) {
		*out = l.a;
		return C3DLAS_INTERSECT;
	}
	if(fabs(db) < FLT_CMP_EPSILON) {
		*out = l.b;
		return C3DLAS_INTERSECT;
	}
	
	float h = vDot3(vSub3(vScale3(pl->n, pl->d), l.a), pl->n);
	float i = vDot3(ldir, pl->n);
	float d = i != 0 ? h / i : 0;
	
	// check if the plane intersects outside the two points
	if(d < 0 || d > vDist3(l.a, l.b)) {
		return C3DLAS_DISJOINT;
	}
	
	*out = vAdd3(l.a, vScale3(ldir, d));
	
	return C3DLAS_INTERSECT;
}




// Assumes full proper intersection.
// C3DLAS_INTERSECT
int findIntersectFastLinePlane(Line3 l, Plane* pl, Vector3* out) {
//int planeLineFindIntersectFast3p(Plane* pl, Vector3* la, Vector3* lb, Vector3* out) {
	Vector3 ldir;
	float h, i, d;
	
	ldir = vSub3(l.b, l.a);
	
	h = vDot3(vSub3(vScale3(pl->n, pl->d), l.a), pl->n);
	i = vDot3(ldir, pl->n);
	d = i != 0 ? h / i : 0;
	
	*out = vAdd3(l.a, vScale3(ldir, d));
	
	return C3DLAS_INTERSECT;
}



// C3DLAS_INTERSECT, _PARALLEL or _DISJOINT
// negative values of idist are "behind" ray->o
int findIntersectTRayPlane(Ray3* ray, Plane* pl, Vector3* ipoint, float* T) {
	float d = vDot3p(&pl->n, &ray->d);
	
	if(fabs(d) < FLT_CMP_EPSILON) return C3DLAS_PARALLEL; // TODO: check for coplanarity?
	
	// negation seems suspicious, may be causing or relying on normal/distance towards origin instead of away from origin
	float t = -(vDot3(pl->n, ray->o) + pl->d) / d;
	
	*ipoint = vAdd3(ray->o, vScale3(ray->d, t));
	*T = t;
	
	return t >= 0 ? C3DLAS_INTERSECT : C3DLAS_DISJOINT;
}










