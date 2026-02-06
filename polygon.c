

#include "internal_array.h"


int polyContainsPoint(Polygon* poly, Vector2 p) {
	int inside = 0;
	int cnt = poly->pointCount;
	
	if(isfinite(poly->maxRadiusSq) && poly->maxRadiusSq < vDistSq(poly->centroid, p)) return 0;
	
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


// returns 1 if the polygon needs a more precise check, and 0 if the polygon is definitely farther away than the radius
int polyCheckMaybeWithinRadius(Polygon* poly, vec2 p, float radius) {
	return vDist(p, poly->centroid) > sqrtf(poly->maxRadiusSq) + radius;
}





// Muchas gracias, Inigo.  
// https://iquilezles.org/articles/distfunctions2d/
float polyDistToPoint(Polygon* poly, Vector2 p) {

	float d = vDot2(vSub2(p, poly->points[0]), vSub2(p, poly->points[0]));
	float s = 1.0;

	for(long i = 0, j = poly->pointCount - 1; i < poly->pointCount; j = i, i++) {
		Vector2 A = poly->points[i];
		Vector2 B = poly->points[j];

		Vector2 e = vSub2(B, A);
		Vector2 w = vSub2(p, A);
		Vector2 b = vSub2(w, vScale2(e, fclamp(vDot2(w, e) / vDot2(e, e), 0.0, 1.0)));

		d = fminf(d, vDot2(b, b));
//		if(fabsf(d) == 0) return 0;         

		int c1 = p.y >= A.y;
		int c2 = p.y < B.y;
		int c3 = e.x * w.y > e.y * w.x;

		if((c1 && c2 && c3) || (!c1 && !c2 && !c3)) s = -s;
	}

	return fabsf(d) == 0 ? 0 : s * sqrtf(d);
}



///
// there are three scenarios:
//	intersect; return 0
//	one fully inside the other; negative distance being closest to edge
//	both fully disjoint; positive dist being the closest approach
float polyMinDistToPoly(Polygon* a, Polygon* b, int* a_is_inside_b) {
	
	// supposedly this O(n^2) algorithm is the best there is
	float d = polyDistToPoint(b, a->points[0]);
	float d2 = polyDistToPoint(a, b->points[0]);
	if(d == 0 || d2 == 0) return 0; // intersecting 
	
	bool was_swapped = 0;
	if(d2 < 0) { // swap
		Polygon* tmp = a;
		a = b;
		b = tmp;
		d = d2;
		was_swapped = 1;
	} 	
	
	if(d > 0) { // the polygons are either intersecting or fully disjoint; return 0 or a positive distance
		
		for(long i = 1; i < a->pointCount; i++) {
			float da = polyDistToPoint(b, a->points[i]);
			if(da <= 0) return 0; // intersecting 
			if(da < d) {
				d = da;
			}
		}
		for(long i = 0; i < b->pointCount; i++) {
			float db = polyDistToPoint(a, b->points[i]);
			if(db <= 0) return 0; // intersecting 
			if(db < d) {
				d = db;
			}
		}
		
	}
	else { // the polygons are either intersecting or nested; return 0 or a negative distance
		
		for(long i = 1; i < a->pointCount; i++) {
			float da = polyDistToPoint(b, a->points[i]);
			if(da >= 0) {
				return 0; // intersecting 
			}
			
			if(da > d) {
				d = da;
			}
		}
		for(long i = 0; i < b->pointCount; i++) {
			float db = -polyDistToPoint(a, b->points[i]); // negative because the perspective is reversed; containing points are always outside
			if(db >= 0) {
				return 0; // intersecting 
			}
			if(db > d) {
				d = db;
			}
		}	
	}

	if(d < 0 && a_is_inside_b) *a_is_inside_b = !was_swapped; 
	
	return d;
}

// ----

static inline void poly_checkAlloc(Polygon* poly) {
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
}

// should be called before p is added to poly->points
static inline void poly_checkMinExtrema(Polygon* poly, vec2 p) {
	if(poly->pointCount == 0 || poly->points[poly->minExtrema].x < p.x) {
		poly->minExtrema = poly->pointCount;
	}
	else if(poly->points[poly->minExtrema].x == p.x && poly->points[poly->minExtrema].y < p.y) {
		poly->minExtrema = poly->pointCount;
	}
}

// adds the point to centroid as an accumulator
// poly->maxRadiusSq is left to the caller
static void polyPushPoint_Internal(Polygon* poly, Vector2 p) {
	poly_checkAlloc(poly);
	poly_checkMinExtrema(poly, p);
	
	poly->points[poly->pointCount++] = p;
	
	poly->centroid = vAdd2(poly->centroid, p);
}



// adds it to the end of the list
// stats need to be recalculated afterward
void polyPushPoint(Polygon* poly, Vector2 p) {
	poly_checkAlloc(poly);
	poly_checkMinExtrema(poly, p);
	
	poly->points[poly->pointCount++] = p;
	poly->centroid = (Vector2){NAN, NAN};
	poly->maxRadiusSq = NAN;
}

// the point is inserted the specified index, pushing all further ones forward
// index can be negative, but should be in the range [-pointCount+1, pointCount] inclusive
// stats need to be recalculated afterward
void polyInsertPoint(Polygon* poly, int index, Vector2 p) {
	poly_checkAlloc(poly);
	
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
	
	if(poly->pointCount == 1 || poly->points[poly->minExtrema].x < p.x) {
		poly->minExtrema = index;
	}
	else if(poly->points[poly->minExtrema].x == p.x && poly->points[poly->minExtrema].y < p.y) {
		poly->minExtrema = index;
	}
}




// only works on convex polygons
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
		
		if(i == 0 || expt2.x < out->points[out->minExtrema].x) {
			out->minExtrema = i;
		}
		else if(out->points[out->minExtrema].x == expt2.x && expt2.y < out->points[out->minExtrema].y) {
			out->minExtrema = i;
		}
	}
	
	// update stats
	out->centroid = vScale2(centroid, 1.0f / out->pointCount); // should be the same(?), but i don't care to prove that mathematically atm
	
	float d = 0;
	for(int i = 0; i < out->pointCount; i++) {
		d = fmaxf(d, vDot2(out->centroid, out->points[i]));
	}
	
	out->maxRadiusSq = d;
}


// requires CCW wound polys
void polyExtrude2(Polygon* poly, f32 dist, Polygon* out) {
	
	long plen = poly->pointCount;
	
	if(out->pointAlloc != poly->pointAlloc) {
		out->points = realloc(out->points, sizeof(out->points) * poly->pointAlloc);
		out->pointAlloc = poly->pointAlloc;
	}

	out->centroid = (Vector2){0,0};
	
	for(long i = 0; i < poly->pointCount; i++) {
		
		Vector2 pa =  poly->points[i];
		Vector2 pb =  poly->points[(i + 1) % plen];
		
		Vector2 A = vNorm2(vSub2(pb, pa));
		Vector2 perp = vScale(vPerpCW2(A), dist);
		
		// TODO: calculate exterior cusps and add extra points
		
		polyPushPoint_Internal(out, vAdd(pa, perp));
		polyPushPoint_Internal(out, vAdd(pb, perp));
	}
	
	// update stats
	out->centroid = vScale2(out->centroid, 1.0f / out->pointCount); // should be the same(?), but i don't care to prove that mathematically atm
	
	float d = 0;
	for(int i = 0; i < out->pointCount; i++) {
		d = fmaxf(d, vDot2(out->centroid, out->points[i]));
	}
	
	out->maxRadiusSq = d;
}

// technically yes, but actually no. It works ok-sh for a certain class of self-intersecting polygons oftencreated by polyExtrude2.
void polyUnSelfIntersect(Polygon* poly) {
	
	
	long totalCopied = 0;
	long plen = poly->pointCount;
	
	long i = poly->minExtrema;
	long writeHead = 1;//(i + 1) % plen; // the first point is always copied in, as the minExtrema can't change
	
	Vector2* tmp = C3DLAS_malloc(sizeof(*tmp) * poly->pointAlloc);
	tmp[0] = poly->points[i];
	
	float minTa = 0; // when continuing on a line; should be >, not >=
	
//	long newMinExtrema = i;
//	Vector2 newMinExtrema_p = poly->points[i];
	Vector2 centroid = {};
	
	
	bool first = 1;
	do {

		Line2 l1 = {.a = poly->points[i], .b = poly->points[(i + 1) % plen]};
		
		
		long besti;
		float bestTa = FLT_MAX;
		float bestTb;
		Vector2 bestp = {1000, 1000};
		
		
		for(long jj = 0, j = (i + 1) % plen; jj < plen; jj++, j = (j + 1) % plen) { // BUG: loops too far
		
			Line2 l2 = {.a = poly->points[j], .b = poly->points[(j + 1) % plen]};
			
			f32 Ta, Tb;
			Vector2 p;
			
			if(C3DLAS_INTERSECT == findIntersectLine2Line2T(l1, l2, &p, &Ta, &Tb)) {
				if(Ta > bestTa) continue; 
				if(Ta <= minTa) continue; // also should catch the endpoints intersecting
				if(i == j) continue;
				
				bestTa = Ta;
				bestTb = Tb;
				bestp = p;
				besti = j;
			}
		}
		
		if(bestTa == FLT_MAX) { // no intersection
			
			totalCopied++;
			
			i = (i + 1) % plen;
						
			tmp[writeHead] = poly->points[i];
			centroid = vAdd2(centroid, poly->points[i]);

			writeHead = (writeHead + 1) % plen;
			

			if(!first && i == poly->minExtrema) goto DONE;
			if(i != poly->minExtrema) first = 0;
			
			minTa = 0;
			
			continue;
		}
		
		// self-intersect.
		
		// push the intersected point
		tmp[writeHead] = bestp;
		writeHead = (writeHead + 1) % plen;
		
		// stats
		centroid = vAdd2(centroid, bestp);
		totalCopied++;
		
		// move on to the intersecting line, starting at the bestTb
		minTa = bestTb;
		i = besti;
		
		if(!first && i == poly->minExtrema) goto DONE;
		if(i != poly->minExtrema) first = 0;
		
	} while(1);
	

	
DONE:
		
	poly->minExtrema = 0; // it was were we started
	memcpy(poly->points, tmp, totalCopied * sizeof(*poly->points));
	
	C3DLAS_free(tmp);
	
	poly->pointCount = totalCopied;
	
	// update stats
	poly->centroid = vScale2(centroid, 1.0f / totalCopied); // should be the same(?), but i don't care to prove that mathematically atm
	
	float d = 0;
	for(long i = 0; i < poly->pointCount; i++) {
		d = fmaxf(d, vDot2(poly->centroid, poly->points[i]));
	}
	
	poly->maxRadiusSq = d;
}



void polyCombineClosePoints(Polygon* poly, float dist) {
	
	long plen = poly->pointCount;
	float distSq = dist * dist;
	
	long writeHead = 1;
	
	long a = 0, b = 1;
	while(b < poly->pointCount) {
		float d2 = vDistSq(poly->points[a], poly->points[b]);
		if(d2 > distSq) {
			poly->points[writeHead] = poly->points[b];
			writeHead++;
			a++; 
			b++;
			continue;
		}
		
		poly->points[writeHead] = vAvg(poly->points[a], poly->points[b]);
		b++;
	}
	
	poly->pointCount = writeHead;
	
	polyCalcStats(poly);
	polyCalcMinExtrema(poly);
}




// breaks each segment into 'degree' number of segments
void polySubdivide(Polygon* poly, int degree) {
	if(degree <= 1) return;
	
	long newCount = poly->pointCount * (degree + 1);
	if(newCount > poly->pointAlloc) {
		poly->pointAlloc = newCount;
		poly->points = realloc(poly->points, sizeof(poly->points) * poly->pointAlloc);
	}
	
	Vector2 prev = poly->points[0];
	
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

	if(
		isfinite(a->maxRadiusSq) && isfinite(b->maxRadiusSq) && 
		vDistSq(a->centroid, b->centroid) > sqrtf(a->maxRadiusSq) + sqrtf(a->maxRadiusSq)
	) return C3DLAS_DISJOINT;
	
	
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


// returns C3DLAS_INTERSECT or C3DLAS_DISJOINT
int polyCheckIntersect(Polygon* a, Polygon* b) {
	int ret;
	
	if(
		isfinite(a->maxRadiusSq) && isfinite(b->maxRadiusSq) && 
		vDist(a->centroid, b->centroid) > sqrtf(a->maxRadiusSq) + sqrtf(a->maxRadiusSq)
	) return C3DLAS_DISJOINT;
	
	
	for(long ai = 0; ai < a->pointCount; ai++) {
		Line2 aline = {a->points[ai], a->points[(ai + 1) % a->pointCount]};
		
		for(long bi = 0; bi < b->pointCount; bi++) {
			
			Vector2 p;
			if(C3DLAS_INTERSECT == intersectLine2Line2(aline, (Line2){b->points[bi], b->points[(bi + 1) % b->pointCount]})) {
				return C3DLAS_INTERSECT;
			}
		}
	}
	
	return C3DLAS_DISJOINT;
}



// atrocious algorithm, but i don't have time to waste on a fancier one
// returns 0 or 1
int polyIsSelfIntersecting(Polygon* poly) {

	for(long i = 0; i < poly->pointCount; i++) {
		for(long j = 0; j < poly->pointCount; j++) {
			Line2 line1 = {poly->points[i % poly->pointCount], poly->points[(i + 1) % poly->pointCount]};
			Line2 line2 = {poly->points[j % poly->pointCount], poly->points[(j + 1) % poly->pointCount]};
			
			bool ba = vEqExact(line1.b, line2.a);
			bool ab = vEqExact(line1.a, line2.b);
			if(ab && ba) {
				return 1;
			}
			
			if(C3DLAS_INTERSECT == intersectLine2Line2(line1, line2)) {
				bool aa = vEqExact(line1.a, line2.a);
				bool bb = vEqExact(line1.b, line2.b);
				
				if(!aa && !ba && !ab && !bb) {
//					printf("%f,%f -> %f,%f | %f,%f -> %f,%f\n", line1.a.x,line1.a.y, line1.b.x,line1.b.y, line2.a.x,line2.a.y, line2.b.x,line2.b.y);
					return 1;
				}
			}
		}
	}

	return 0;
}


// unsigned area
float polyCalcArea(Polygon* poly) {
	return fabs(polyCalcSignedArea(poly));
}


void polyCalcMinExtrema(Polygon* poly) {
	poly->minExtrema = 0;
	for(long i = 0; i < poly->pointCount; i++) {
		if(poly->points[poly->minExtrema].x < poly->points[i].x) {
			poly->minExtrema = poly->pointCount;
		}
		else if(poly->points[poly->minExtrema].x == poly->points[i].x && poly->points[poly->minExtrema].y < poly->points[i].y) {
			poly->minExtrema = poly->pointCount;
		}
	}
}


AABB2 polyCalcBounds(Polygon* poly) {
	AABB2 b = {{FLT_MAX, FLT_MAX}, {-FLT_MAX, -FLT_MAX}};
	
	for(long i = 0; i < poly->pointCount; i++) {
		Vector2 p = poly->points[i];
	
		b.min.x = fminf(b.min.x, p.x);
		b.min.y = fminf(b.min.y, p.y);
		b.max.x = fmaxf(b.max.x, p.x);
		b.max.y = fmaxf(b.max.y, p.y);
	}
	
	return b;
}


float polyCalcSignedArea(Polygon* poly) {
	f32 area = 0;
	
	for(long i = 0; i < poly->pointCount; i++) {
		Vector2 p0 = poly->points[i];
		Vector2 p1 = poly->points[(i + 1) % poly->pointCount];
		
		area += p0.x * p1.y;
		area -= p0.y * p1.x;
	}
	
	return 0.5f * area;
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
		d = fmaxf(d, vDistSq(poly->centroid, a));
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

// very niche use; doesn't do what you think it does
void polySortCCW(Polygon* poly) {
	if(!isfinite(poly->centroid.x)) {
		polyCalcCentroid(poly);
	}
	
	C3DLAS_sort_r(poly->points, poly->pointCount, sizeof(*poly->points), (void*)poly_point_sort_ccw_fn, &poly->centroid);
}


// returns C3DLAS_CCW or C3DLAS_CW;
int polyGetWinding(Polygon* poly) {
	Vector2 p = poly->points[poly->minExtrema];
	Vector2 a = vSub2(poly->points[(poly->minExtrema + 1) % poly->pointCount], p);
	Vector2 b = vSub2(poly->points[(poly->minExtrema - 1 + poly->pointCount) % poly->pointCount], p);
	
	return vCross2(a, b) > 0 ? C3DLAS_CCW : C3DLAS_CW;
}

void polyEnsureCCW(Polygon* poly) {
	if(polyGetWinding(poly) == C3DLAS_CW) polyReverse(poly);
}

void polyEnsureCW(Polygon* poly) {
	if(polyGetWinding(poly) == C3DLAS_CCW) polyReverse(poly);
}

void polyReverse(Polygon* poly) {
	size_t i, j;
	Vector2 tmp;
	
	for(i = 0, j = poly->pointCount - 1; i < j; i++, j--) {
		tmp = poly->points[i];
		poly->points[i] = poly->points[j];
		poly->points[j] = tmp;
	}
}




// a full deep copy
void polyCopy(Polygon* dst, Polygon* src) {
	if(dst->pointAlloc != src->pointAlloc) {
		dst->points = realloc(dst->points, sizeof(dst->points) * src->pointAlloc);
		dst->pointAlloc = src->pointAlloc;
	}
	
	if(src->pointCount) {
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





// Polygon Union

struct poly_CrossingInfo {
	long i;
	long crossing_count;
};

struct poly_CrossingPoint {
	long ci, di; // every crossing has a unique index set; a straight line can't cross twice
	long pi; // index of the crossing point
	float tc, td; // t values for the c and d line respectively
	Vector2 p;
};

static int poly_point_sort_c_fn(struct poly_CrossingPoint** a, struct poly_CrossingPoint** b, void* meh) {
	if((**a).tc < (**b).tc) return -1;
	if((**a).tc > (**b).tc) return 1;
	return 0;
}
static int poly_point_sort_d_fn(struct poly_CrossingPoint** a, struct poly_CrossingPoint** b, void* meh) {
	if((**a).td < (**b).td) return -1;
	if((**a).td > (**b).td) return 1;
	return 0;
}


// eliminates holes between the two
// both must be CCW sorted
// returns C3DLAS_INTERSECT if the union was successful and C3DLAS_DISJOINT if they don't intersect
// if it returns disjoint then out will be a copy of either a or b
int polyExteriorUnion(Polygon* a, Polygon* b, Polygon* out) {
	
	if(b->pointCount == 0 || a->pointCount == 0) {
		return C3DLAS_DISJOINT;
	}
	
	long failsafe_limit;
	
	// figure out all the segments that have crossings, and how many they have.
	long* crossinglist_a = C3DLAS_calloc(1, sizeof(*crossinglist_a) * a->pointCount); // all from c
	long* crossinglist_b = C3DLAS_calloc(1, sizeof(*crossinglist_b) * b->pointCount); // all from d
	long total_crossings = 0;
	long num_multiples = 0;
	
	for(long ai = 0; ai < a->pointCount; ai++) {
		Line2 aline = {a->points[ai], a->points[(ai + 1) % a->pointCount]};
		
		for(long bi = 0; bi < b->pointCount; bi++) {
			Line2 bline = {b->points[bi], b->points[(bi + 1) % b->pointCount]};
			
			if(C3DLAS_INTERSECT == intersectLine2Line2(aline, bline)) {
				crossinglist_a[ai]++;
				crossinglist_b[bi]++;
				
				if(crossinglist_b[bi] == 2 || crossinglist_a[ai] == 2) {
					num_multiples++;
				} 
				
				total_crossings++;
			}
		}
	}
	
	failsafe_limit = total_crossings + 5 + a->pointCount + b->pointCount;
	
	// some trivial cases
	if(total_crossings == 0) {
		// check if one is inside the other
		
		C3DLAS_free(crossinglist_a);
		C3DLAS_free(crossinglist_b);
		
		if(polyContainsPoint(a, b->points[0])) { // b is inside a, return a
			polyCopy(out, a);
			return C3DLAS_INTERSECT;
		}
		
		if(polyContainsPoint(b, a->points[0])) { // a is inside b, return b
			polyCopy(out, b);
			return C3DLAS_INTERSECT;
		}
		
		return C3DLAS_DISJOINT;
	}
	
	
	// TODO: combine this part into the complicated one's setup?
	// find an extreme point, guaranteed to be on the union. it could be on either polygon
	int is_a = 1;
	long start_i;
	vec2 start = {FLT_MAX, FLT_MAX};
	for(long ai = 0; ai < a->pointCount; ai++) {
		if(a->points[ai].x < start.x) {
			start = a->points[ai];
			start_i = ai;
		}
		else if(a->points[ai].x == start.x && a->points[ai].y <= start.y) {
			start = a->points[ai];
			start_i = ai;
		}
	}
	for(long bi = 0; bi < b->pointCount; bi++) {
		if(b->points[bi].x < start.x) {
			start = b->points[bi];
			start_i = bi;
			is_a = 0;
		}
		else if(b->points[bi].x == start.x && b->points[bi].y <= start.y) {
			start = b->points[bi];
			start_i = bi;
			is_a = 0;
		}
	}
	
	
	
	
	// the algorithm flips back and forth walking on a and b
	// c is the first polygon being walked, with the extreme point, and d is the other one initially being searched for intersections
	Polygon* c = is_a ? a : b;
	Polygon* d = is_a ? b : a;
	long* crossinglist_c = is_a ? crossinglist_a : crossinglist_b;
	long* crossinglist_d = is_a ? crossinglist_b : crossinglist_a;
	long end_ci = is_a ? (start_i - 1 + a->pointCount) % a->pointCount : (start_i - 1 + b->pointCount) % b->pointCount;
	
	
	
	if(num_multiples == 0) { // somewhat simple algorithm
		
		// walk from the start point, adding the obvious segment if there's only one but going into
		//  a more complicated algorithm if there's crossings
		
		long ci = start_i;
		do {
			
			if(out->pointCount > failsafe_limit) {
//				printf("failsafe simpler \n");
//				exit(1);
				return C3DLAS_INTERSECT;
			}
			
//			if(debug1++ >= 12) return C3DLAS_INTERSECT;
			
			long nc = crossinglist_c[ci];
			if(nc == 0) { // easy case: one edge, push it and carry on
				polyPushPoint(out, c->points[ci % c->pointCount]);
				ci = (ci + 1) % c->pointCount;
				continue;
			}
		
			
			Line2 cline = {c->points[ci], c->points[(ci + 1) % c->pointCount]};
			
			long best_di;
			for(long di = 0; di < d->pointCount; di++) { // maybe could be cached
				Line2 dline = {d->points[di % d->pointCount], d->points[(di + 1) % d->pointCount]};
			
				Vector2 p;
				if(C3DLAS_INTERSECT == findIntersectLine2Line2(cline, dline, &p)) {
					polyPushPoint(out, c->points[ci % c->pointCount]);
					polyPushPoint(out, p);
					best_di = di;
					break;
				}
			}
			
			// the entire d-search
			long end_di = best_di;
			for(long di = (best_di + 1) % d->pointCount; di != end_di; di = (di + 1) % d->pointCount) {
		
		
//				if(debug2++ >= 12) return C3DLAS_INTERSECT;
		
				long nc = crossinglist_d[di];
				if(nc == 0) { // easy case: one edge, push it and carry on
					polyPushPoint(out, d->points[di % d->pointCount]);
					continue;
				}
			
				Line2 dline = {d->points[di], d->points[(di + 1) % d->pointCount]};
				
				long cci = 0;
				for(; cci <= c->pointCount; cci++) { // maybe could be cached
					Line2 ccline = {c->points[cci % c->pointCount], c->points[(cci + 1) % c->pointCount]};
				
					Vector2 p;
					if(C3DLAS_INTERSECT == findIntersectLine2Line2(ccline, dline, &p)) {
						polyPushPoint(out, d->points[di % d->pointCount]);
						polyPushPoint(out, p);
						break;
					}
				}
				
				// back to the c-search
				ci = (cci + 1) % c->pointCount;
				
				break;
			}
				
		} while(ci != start_i);
		
		
		// returns below
	}
	else { // there are multiple crossings of the same edge
		
		// complicated fucked up graph-based algorithm
		
		ARR(struct poly_CrossingInfo) cr_infos_c = {};
		ARR(struct poly_CrossingInfo) cr_infos_d = {};
		ARR(struct poly_CrossingPoint) cr_points = {};
		
		long max_crossing_count = 0; 
		
		// make a list of all the crossing points and a list of metadata fro both c and d
		for(long ci = 0; ci < c->pointCount; ci++) {
			if(0 == crossinglist_c[ci]) continue;
			
			
			struct poly_CrossingInfo* cri_c = NULL;
			
			Line2 cline = {c->points[ci], c->points[(ci + 1) % c->pointCount]};
			
			for(long di = 0; di < d->pointCount; di++) {
				if(0 == crossinglist_d[di]) continue;
				
				Line2 dline = {d->points[di], d->points[(di + 1) % d->pointCount]};
			
				float tc, td;
				Vector2 p;
				if(C3DLAS_INTERSECT == findIntersectLine2Line2T(cline, dline, &p, &tc, &td)) {
					
					struct poly_CrossingPoint* crp;
					ARR_inc(&cr_points, crp);
					crp->ci = ci;
					crp->di = di;
					crp->pi = cr_points.count - 1;
					crp->tc = tc;
					crp->td = td;
					crp->p = p;
					
					// Just checking for uniqueness; linear search should be fine considering how few crossings there should be in most
					//   use cases. If you put ridiculously bad polygons in you get ridiculously bad performance out.
					struct poly_CrossingInfo* cri_d = NULL;
					ARR_EACH(&cr_infos_d, i) {
						if(cr_infos_d.data[i].i != di) continue;
						cri_d = &cr_infos_d.data[i];
						break;
					}
					if(!cri_d) {
						ARR_inc(&cr_infos_d, cri_d);
						cri_d->i = di;
						cri_d->crossing_count = 0;
					}
					
					if(!cri_c) {
						ARR_EACH(&cr_infos_c, i) {
							if(cr_infos_c.data[i].i != ci) continue;
							cri_c = &cr_infos_c.data[i];
							break;
						}
						if(!cri_c) {
							ARR_inc(&cr_infos_c, cri_c);
							cri_c->i = ci;
							cri_c->crossing_count = 0;
						}	
					}
					
					cri_c->crossing_count++;
					cri_d->crossing_count++;
					
					max_crossing_count = MAX(max_crossing_count, MAX(cri_c->crossing_count, cri_d->crossing_count));
				}
			}
		}
		
		struct poly_CrossingPoint** point_list = C3DLAS_malloc(sizeof(*point_list) * max_crossing_count);
		
		struct Edge {
			long a, b; // [c->points | d->points | cr_points]
		};
		
		ARR(struct Edge) edges = {}; 
		
		#define add_edge(x, y) \
			ARR_push(&edges, (struct Edge){.a = x, .b = y})
			
		
		// add all the edges
		
		long coff = c->pointCount;
		long cdoff = c->pointCount + d->pointCount;
		ARR_EACH(&cr_infos_c, i) {
			struct poly_CrossingInfo* cri = &cr_infos_c.data[i];
			long pt_cnt = 0;
			
			// collect c crossings
			ARR_EACH(&cr_points, j) {
				if(cri->i != cr_points.data[j].ci) continue;
				point_list[pt_cnt++] = &cr_points.data[j];
			}
			
			//sort by tc
			C3DLAS_sort_r(point_list, pt_cnt, sizeof(*point_list), (void*)poly_point_sort_c_fn, NULL);
			
			// add edges
			add_edge(cri->i, point_list[0]->pi + cdoff);
			for(int j = 1; j < pt_cnt; j++) {
				add_edge(point_list[j - 1]->pi + cdoff, point_list[j]->pi + cdoff);
			}
			add_edge(point_list[pt_cnt - 1]->pi + cdoff, (cri->i + 1) % c->pointCount);
		}
		
		ARR_EACH(&cr_infos_d, i) {
			struct poly_CrossingInfo* cri = &cr_infos_d.data[i];
			long pt_cnt = 0;
			
			// collect d crossings
			ARR_EACH(&cr_points, j) {
				if(cri->i != cr_points.data[j].di) continue;
				point_list[pt_cnt++] = &cr_points.data[j];
			}
			
			//sort by td
			C3DLAS_sort_r(point_list, pt_cnt, sizeof(*point_list), (void*)poly_point_sort_d_fn, NULL);
			
			// add edges
			add_edge(cri->i + coff, point_list[0]->pi + cdoff);
			for(int j = 1; j < pt_cnt; j++) {
				add_edge(point_list[j - 1]->pi + cdoff, point_list[j]->pi + cdoff);
			}
			add_edge(point_list[pt_cnt - 1]->pi + cdoff, ((cri->i + 1) % d->pointCount) + coff);
		}
		
	
		// now for the actual edge-following algorithm
		
		long prevprev = (start_i - 1 + c->pointCount) % c->pointCount;
		long prev = start_i;
		do {
			
			if(out->pointCount > failsafe_limit) {
//				printf("\n\n\nvec2 p_points[] = {\n");
//				FOR(k, a->pointCount) {
//					printf("\t{%f, %f},\n", a->points[k].x, a->points[k].y);
//				}
//				printf("};\n\n\n");
//				
//				printf("\n\n\nvec2 p2_points[] = {\n");
//				FOR(k, b->pointCount) {
//					printf("\t{%f, %f},\n", b->points[k].x, b->points[k].y);
//				}
//				printf("};\n\n\n");
//				*(int*)0 = 0;
				return C3DLAS_INTERSECT;
			}
			
			// TODO: some sort of absolute limit, like croff+vec_len(cr_points)+1
			
			if(prev < coff) {
				polyPushPoint(out, c->points[prev]);

				if(crossinglist_c[prev] == 0) {
					prevprev = prev;
					prev = (prev + 1) % c->pointCount;
				}
				else {
					
					// search the edges
					float best_angle = FLT_MAX;
					struct Edge* best_edge = NULL;
					
					vec2 p3 = (prevprev < coff) ? c->points[prevprev] : ((prevprev < cdoff) ? d->points[prevprev - coff] : VEC_item(&cr_points, prevprev - cdoff).p);
					vec2 p1 = c->points[prev];
					vec2 best_p = {};
					
					ARR_EACH(&edges, j) {
						struct Edge* e = &edges.data[j];
						if(e->a != prev) continue;
						
						// should only be two outgoing edges; one on a c line and one on a d line
						// p1 is the pivot // positive = CCW, negative = CW
						
						vec2 p2 = (e->b < coff) ? c->points[e->b] : ((e->b < cdoff) ? d->points[e->b - coff] : VEC_item(&cr_points, e->b - cdoff).p);
						float angle = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x); 
						
						if(angle < best_angle) {
							best_angle = angle;
							best_edge = e;
							best_p = p2;
						}
					}
					
					prevprev = prev;
					prev = best_edge->b;
				}
			
			}
			else if(prev < cdoff) {
				int di = prev - coff;
				polyPushPoint(out, d->points[prev - coff]);

				if(crossinglist_d[di] == 0) {
					prevprev = prev;
					prev = coff + ((di + 1) % d->pointCount);
				}
				else {
					// search the edges
					
					float best_angle = FLT_MAX;
					struct Edge* best_edge = NULL;
					
					vec2 p3 = (prevprev < coff) ? c->points[prevprev] : ((prevprev < cdoff) ? d->points[prevprev - coff] : VEC_item(&cr_points, prevprev - cdoff).p);
					vec2 p1 = d->points[prev - coff];
					vec2 best_p = {};
					
					ARR_EACH(&edges, j) {
						struct Edge* e = &edges.data[j];
						if(e->a != prev) continue;
						
						// should only be two outgoing edges; one on a c line and one on a d line
						// p1 is the pivot // positive = CCW, negative = CW
						
						vec2 p2 = (e->b < coff) ? c->points[e->b] : ((e->b < cdoff) ? d->points[e->b - coff] : VEC_item(&cr_points, e->b - cdoff).p);
						float angle = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x); 
						
						if(angle < best_angle) {
							best_angle = angle;
							best_edge = e;
							best_p = p2;
						}
					}
					
					prevprev = prev;
					prev = best_edge->b;
				}
			
			}
			else { // crossing point

				// search the edges				
				
				float best_angle = FLT_MAX;
				struct Edge* best_edge = NULL;
				
				vec2 p3;
				if(prevprev < coff) p3 = c->points[prevprev];
				else if(prevprev < cdoff) p3 = d->points[prevprev - coff];
				else p3 = VEC_item(&cr_points, prevprev - cdoff).p;
				
				
				vec2 p1 = VEC_item(&cr_points, prev - cdoff).p;
				vec2 best_p = {};

				polyPushPoint(out, p1);

				ARR_EACH(&edges, j) {
					struct Edge* e = &edges.data[j];
					if(e->a != prev) continue;
					
					// should only be two outgoing edges; one on a c line and one on a d line
					// p1 is the pivot // positive = CCW, negative = CW
					
					vec2 p2 = (e->b < coff) ? c->points[e->b] : ((e->b < cdoff) ? d->points[e->b - coff] : VEC_item(&cr_points, e->b - cdoff).p);
					float angle = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x); 
					
					if(angle < best_angle) {
						best_angle = angle;
						best_edge = e;
						best_p = p2;
					}
				}
				
				prevprev = prev;
				prev = best_edge->b;
			}
			
		} while(prev != start_i);
		
		
		ARR_free(&cr_infos_c);
		ARR_free(&cr_infos_d);
		ARR_free(&cr_points);
		ARR_free(&edges);
		C3DLAS_free(point_list);
	}
	
	
	C3DLAS_free(crossinglist_a);
	C3DLAS_free(crossinglist_b);
	
	return C3DLAS_INTERSECT;

}

