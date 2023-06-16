


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
int planeLineFindIntersect3p(Plane* pl, Vector3* la, Vector3* lb, Vector3* out) {
	Vector3 ldir;
	float da, db;
	
	vSub3p(lb, la, &ldir);
	
	da = vDot3p(la, &pl->n) - pl->d;
	
	// bail if the line and plane are parallel
	if(fabs(vDot3p(&pl->n, &ldir)) < FLT_CMP_EPSILON) {
		
		// check coplanarity
		if(fabs(da) < FLT_CMP_EPSILON) {
			return C3DLAS_COPLANAR; // the end is on the plane, so the other is too
		}
		
		return C3DLAS_DISJOINT;
	}
	

	db = vDot3p(lb, &pl->n) - pl->d;
	
	// check if one of the points is on the plane
	if(fabs(da) < FLT_CMP_EPSILON) {
		*out = *la;
		return C3DLAS_INTERSECT;
	}
	if(fabs(db) < FLT_CMP_EPSILON) {
		*out = *lb;
		return C3DLAS_INTERSECT;
	}
	
	Vector3 p0, g, j;
	vScale3p(&pl->n, pl->d, &p0);
	vSub3p(&p0, la, &g);
	float h = vDot3p(&g, &pl->n);
	float i = vDot3p(&ldir, &pl->n);
	float d = i != 0 ? h / i : 0;
	
	// check if the plane intersects outside the two points
	if(d < 0 || d > vDist3p(la, lb)) {
		return C3DLAS_DISJOINT;
	}
	
	vScale3p(&ldir, d, &j);
	vAdd3p(la, &j, out);
	
	return C3DLAS_INTERSECT;
}




// Assumes full proper intersection.
// C3DLAS_INTERSECT
int planeLineFindIntersectFast3p(Plane* pl, Vector3* la, Vector3* lb, Vector3* out) {
	Vector3 ldir, p0, g, j;
	float h, i, d;
	
	vSub3p(lb, la, &ldir);
	
	vScale3p(&pl->n, pl->d, &p0);
	vSub3p(&p0, la, &g);
	h = vDot3p(&g, &pl->n);
	i = vDot3p(&ldir, &pl->n);
	d = i != 0 ? h / i : 0;
	
	vScale3p(&ldir, d, &j);
	vAdd3p(la, &j, out);
	
	return C3DLAS_INTERSECT;
}




