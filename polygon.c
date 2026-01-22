


int polyContainsPoint(Polygon* poly, Vector2 p) {
	int inside = 0;
	int cnt = poly->pointCount;
	
	if(isfinite(poly->maxRadiusSq) && poly->maxRadiusSq < vDot2(poly->centroid, p)) return 0;
	
	for(int i = 0; i < cnt; i++) {
		Vector2 a = poly->points[i];
		Vector2 b = poly->points[(i + 1) % cnt];
		
		if(a.y == b.y) continue; // horizontal edges are ignored
		
		// we're testing a ray going to the right
		if(a.x < p.x && b.x < p.x) continue; // segment is entirely left of the point
		
		if(a.y >= p.y && b.y >= p.y) continue; // segment entirely above the point
		if(a.y < p.y && b.y < p.y) continue; // segment entirely below the point
		// segment is in the same vertical band as the point
		
		float sx = a.x + (b.x - a.x) * ((p.y - a.y) / (b.y - a.y));
		if(p.x > sx) continue;
		
		inside = !inside;
	}

	return inside;
}


// Muchas gracias, Inigo.  
// https://iquilezles.org/articles/distfunctions2d/
float polyDistToPoint(Polygon* poly, Vector2 p) {

	float d = vDot2(vSub2(p, poly->points[0]), vSub2(p, poly->points[0]));
	float s = 1.0;

	for(int i = 0, j = poly->pointCount - 1; i < poly->pointCount; j = i, i++) {
		Vector2 A = poly->points[i];
		Vector2 B = poly->points[j];

		Vector2 e = vSub2(B, A);
		Vector2 w = vSub2(p, A);
		Vector2 b = vSub2(w, vScale2(e, fclamp(vDot2(w, e) / vDot2(e, e), 0.0, 1.0)));

		d = fminf(d, vDot2(b, b));
        
		int c1 = p.y >= A.y;
		int c2 = p.y < B.y;
		int c3 = e.x * w.y > e.y * w.x;
		if((c1 && c2 && c3) || (!c1 && !c2 && !c3)) s *= -1.0;  
	}

	return s * sqrtf(d);
}

// ----


// adds it to the end of the list
// stats need to be recalculated afterward
void polyPushPoint(Polygon* poly, Vector2 p) {
	if(poly->pointCount >= poly->pointAlloc) {
		if(poly->pointAlloc == 0) {
			poly->pointAlloc = C3DLAS_INITIAL_POLY_POINT_ALLOC;
			poly->points = C3DLAS_malloc(sizeof(*poly->points) * poly->pointAlloc);
		}
		else {
			poly->pointAlloc *= 2;
			poly->points = C3DLAS_realloc(poly->points, sizeof(*poly->points) * poly->pointAlloc);
		}
	}
	
	poly->points[poly->pointCount++] = p;
	poly->centroid = (Vector2){NAN, NAN};
	poly->maxRadiusSq = NAN;
}

// the point is inserted the specified index, pushing all further ones forward
// index can be negative, but should be in the range [-pointCount+1, pointCount] inclusive
// stats need to be recalculated afterward
void polyInsertPoint(Polygon* poly, int index, Vector2 p) {
	if(poly->pointCount >= poly->pointAlloc) {
		if(poly->pointAlloc == 0) {
			poly->pointAlloc = C3DLAS_INITIAL_POLY_POINT_ALLOC;
			poly->points = C3DLAS_malloc(sizeof(*poly->points) * poly->pointAlloc);
		}
		else {
			poly->pointAlloc *= 2;
			poly->points = C3DLAS_realloc(poly->points, sizeof(*poly->points) * poly->pointAlloc);
		}
	}
	
	// make it positive
	index = (index + poly->pointCount) % poly->pointCount;
	
	if(index < poly->pointCount) {
		// TODO: verify this
		memmove(poly->points + index + 1, poly->points + index, sizeof(poly->points) * poly->pointCount - index);
	}
	
	poly->points[index] = p;
	poly->pointCount++;
	poly->centroid = (Vector2){NAN, NAN};
	poly->maxRadiusSq = NAN;
}



void polyExtrude(Polygon* poly, f32 dist, Polygon* out) {
	
	long chlen = poly->pointCount;
	long chlenm1 = chlen - 1;
	
	if(out->pointAlloc != poly->pointAlloc) {
		out->points = realloc(out->points, sizeof(out->points) * poly->pointAlloc);
		out->pointAlloc = poly->pointAlloc;
	}
	out->pointCount = poly->pointCount;

	Vector2 centroid = {0,0};
	for(long i = 0; i < poly->pointCount; i++) {
//	VEC_EACH(&z->convexHull, i, pi) {
//		vec3 p3 = VEC_item(&z->points, pi);
		vec2 p2 = poly->points[i];//V2(p3);
		
		vec2 pa =  poly->points[(i + 1) % chlen];
		vec2 pb =  poly->points[(i + chlenm1) % chlen];
		
		vec2 A = vNorm2(vSub2(p2, pa));
		vec2 B = vNorm2(vSub2(p2, pb));
		
		
		f32 sinth = vCross2(A, B);
		f32 q = dist / sinth;
		
		vec2 expt2 = vAdd2(p2, vAdd2(vScale2(A, q), vScale2(B, q)) );
		
		out->points[i] = expt2;
		centroid = vAdd2(centroid, expt2);
	}
	
	// update stats
	out->centroid = vScale2(centroid, 1.0f / out->pointCount); // should be the same(?), but i don't care to prove that mathematically atm
	
	float d = 0;
	for(int i = 0; i < out->pointCount; i++) {
		d = fmaxf(d, vDot2(out->centroid, out->points[i]));
	}
	
	out->maxRadiusSq = d;
}



// breaks each segment into 'degree' number of segments
void polySubdivide(Polygon* poly, int degree) {
	if(degree <= 1) return;
	
	long newCount = poly->pointCount * (degree + 1);
	if(newCount > poly->pointAlloc) {
		poly->pointAlloc = newCount;
		poly->points = realloc(poly->points, sizeof(poly->points) * poly->pointAlloc);
	}
	
	vec2 prev = poly->points[0];
	
	long n = newCount - 1;
	for(long i = poly->pointCount - 1; i >= 0; i--) {
		
		
		for(int d = 0; d < degree; d++) {
			float t = (float)d / (float)degree;
			poly->points[n--] = vLerp2(prev, poly->points[i], t);
		}
		
		poly->points[n--] = poly->points[i];
		prev = poly->points[i];
	}

	poly->pointCount = newCount;
}



// if out_fn returns nonzero, polyIntersect returns that value immediately, otherwise zero
// seg_a and seg_b are the index of the point starting the segment which intersects
int polyIntersect(Polygon* a, Polygon* b, int (*out_fn)(Vector2 p, long seg_a, long seg_b, void* data), void* out_data) {
	int ret;
	
	for(long ai = 0; ai < a->pointCount; ai++) {
		Line2 aline = {a->points[ai], a->points[(ai + 1) % a->pointCount]};
		
		for(long bi = 0; bi < b->pointCount; bi++) {
			
			Vector2 p;
			if(C3DLAS_INTERSECT == findIntersectLine2Line2(aline, (Line2){b->points[bi], b->points[(bi + 1) % b->pointCount]}, &p)) {
				if(ret = out_fn(p, ai, bi, out_data)) {
					return ret;
				}
			}
		}
	}
	
	return 0;
}


// eliminates holes between the two
// both must be CCW sorted
// returns C3DLAS_INTERSECT if the union was successful and C3DLAS_DISJOINT if they don't intersect
// if it returns disjoint then out will be a copy of either a or b
int polyExteriorUnion(Polygon* a, Polygon* b, Polygon* out) {
	int ret = C3DLAS_DISJOINT;
	
	
	// find an extreme point, guaranteed to be on the union. it could be on either polygon
	int is_a = 1;
	long start_i;
	vec2 start = {FLT_MAX, FLT_MAX};
	for(long ai = 0; ai < a->pointCount; ai++) {
		if(a->points[ai].x <= start.x) {
			if(a->points[ai].y < start.y) {
				start = a->points[ai];
				start_i = ai;
			}
		}
	}
	for(long bi = 0; bi < b->pointCount; bi++) {
		if(b->points[bi].x < start.x) {
			if(b->points[bi].y < start.y) {
				start = b->points[bi];
				start_i = bi;
				is_a = 0;
			}
		}
	}
	
	// the algorithm flips back and forth walking on a and b
	// c is the first polygon being walked, with the extreme point, and d is the other one initially being searched for intersections
	Polygon* c = is_a ? a : b;
	Polygon* d = is_a ? b : a;
	long end_ci = is_a ? (start_i - 1 + a->pointCount) % a->pointCount : (start_i - 1 + b->pointCount) % b->pointCount;
	
	
	for(long ci = start_i; ci != end_ci; ci = (ci + 1) % c->pointCount) {
//		long ci = (nci + start_i) % c->pointCount;
		Line2 cline = {c->points[ci], c->points[(ci + 1) % c->pointCount]};
		
		// find the first intersection on line c
		float best_tc = FLT_MAX;
		long best_di = -1;
		Vector2 best_p;
		
		for(long di = 0; di <= d->pointCount; di++) {
			Line2 dline = {d->points[di], d->points[(di + 1) % d->pointCount]};
			
			float tc, td;
			Vector2 p;
			if(C3DLAS_INTERSECT == findIntersectLine2Line2T(cline, dline, &p, &tc, &td)) {
				
				if(tc < best_tc) {
					best_tc = tc;
					best_di = di;
					best_p = p;
					ret = C3DLAS_INTERSECT;
				}
			}
		}

		polyPushPoint(out, c->points[ci]);
		
		if(best_tc == FLT_MAX) { // no intersection, keep walking c
			continue;
		}
		
		// push the new point best_p then start walking on d until we get back to c		
		polyPushPoint(out, best_p);
		
		long end_di = best_di;
		for(long di = (best_di + 1) % d->pointCount; di != end_di; di = (di + 1) % d->pointCount) {
			Line2 dline = {d->points[di], d->points[(di + 1) % d->pointCount]};
			
			// find the first intersection on line d
			float best_td = FLT_MAX;
			long best_cci = -1;
			Vector2 best_p;
			
			for(long cci = 0; cci <= c->pointCount; cci++) {
				Line2 ccline = {c->points[cci % c->pointCount], c->points[(cci + 1) % c->pointCount]};

				float td, tcc;
				Vector2 p;
				if(C3DLAS_INTERSECT == findIntersectLine2Line2T(dline, ccline, &p, &td, &tcc)) {
					if(td < best_td) {
						best_td = td;
						best_cci = cci % c->pointCount;
						best_p = p;
					}
				}
			}
			
			polyPushPoint(out, d->points[di]);

			if(best_td == FLT_MAX) { // no intersection, keep walking c
				continue;
			}
			
			polyPushPoint(out, best_p);
			
			// go back to polygon c
			ci = (best_cci + 1) % c->pointCount;
			break;
		}
		
	}
	
	return ret;
}


void polyCalcCentroid(Polygon* poly) {
	int cnt = poly->pointCount;
	Vector2 centroid = {0,0};
	
	for(int i = 0; i < cnt; i++) {
		Vector2 a = poly->points[i];
		centroid = vAdd2(centroid, a);
	}
	
	poly->centroid = vScale2(centroid, 1.0 / poly->pointCount);

}


void polyCalcRadiusSq(Polygon* poly) {
	int cnt = poly->pointCount;
	float d = 0;
	
	for(int i = 0; i < cnt; i++) {
		Vector2 a = poly->points[i];
		d = fmaxf(d, vDot2(poly->centroid, a));
	}
	
	poly->maxRadiusSq = d;
}


// centroid and radius squared
void polyCalcStats(Polygon* poly) {
	polyCalcCentroid(poly);
	polyCalcRadiusSq(poly); // depends on the centroid
}



static int poly_point_sort_ccw_fn(Vector2* a, Vector2* b, Vector2* center) {
	Vector2 da = vSub2(*a, *center);
	Vector2 db = vSub2(*b, *center);
	
	if(da.x >= 0 && db.x < 0) return 1;
	if(da.x < 0 && db.x >= 0) return -1;
	if(da.x == 0 && db.x == 0) return 0;
	
	float ret = vCross2(da, db);
	return ret == 0 ? 0 : (ret > 0 ? -1 : 1);
}

void polySortCCW(Polygon* poly) {
	if(!isfinite(poly->centroid.x)) {
		polyCalcCentroid(poly);
	}
	
	C3DLAS_sort_r(poly->points, poly->pointCount, sizeof(*poly->points), (void*)poly_point_sort_ccw_fn, &poly->centroid);
}




// a full deep copy
void polyCopy(Polygon* dst, Polygon* src) {
	if(dst->pointAlloc != src->pointAlloc) {
		dst->points = realloc(dst->points, sizeof(dst->points) * src->pointAlloc);
		dst->pointAlloc = src->pointAlloc;
	}
	
	if(dst->pointCount) {
		memcpy(dst->points, src->points, sizeof(dst->points) * src->pointCount);
	}
	
	dst->pointCount = src->pointCount;
	dst->centroid = src->centroid;
	dst->maxRadiusSq = src->maxRadiusSq;
}



void polyFreeInternals(Polygon* poly) {
	if(poly->points) C3DLAS_free(poly->points);
	poly->points = NULL;
	poly->pointCount = 0;
	poly->pointAlloc = 0;
}




