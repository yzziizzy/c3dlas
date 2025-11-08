
// Muchas gracias, Inigo.  
// https://iquilezles.org/articles/distfunctions2d/
float vDistPointLine2(Vector2 p, Line2 ls) {
	Vector2 pa = vSub2(p, ls.start); // vector from the starting point to p
	Vector2 ba = vSub2(ls.end, ls.start); // vector from the starting point to the ending point
	
	float t = vDot2(pa, ba) / vDot2(ba, ba); // project the pa onto ba, then divide that distance by the length of ba to normalize it
	t = fclamp(t, 0.0, 1.0); // clamp t to between the endpoints of the line segment
	
	// Consider the starting point to be at the origin, for ease of visualization.
	// ba is the vector from the origin to the endpoint of the line that now passes through the origin.
	// Scaling ba by t gives the intercept point of the line through p that is perpendicular to the test line segment.
	// pa is p if a was the origin. Therefore, pi is the vector from p to the intercept point on the test line segment. 
	Vector2 pi = vSub2(pa, vScale2(ba, t));
	return vLen2(pi); // the answer is the length of pi 
}

float vDistPointLine3(Vector3 p, Line3 ls) {
	Vector3 pa = vSub3(p, ls.start);
	Vector3 ba = vSub3(ls.end, ls.start);
	
	float t = fclamp(vDot3(pa, ba) / vDot3(ba, ba), 0.0, 1.0);
	return vMag3(vSub3(pa, vScale3(ba, t)));
}

// This version also returns the normalized distance along the line of the closest point
float vDistTPointLine2(Vector2 p, Line2 ls, float* T) {
	Vector2 pa = vSub2(p, ls.start);
	Vector2 ba = vSub2(ls.end, ls.start);
	
	float t = fclamp(vDot2(pa, ba) / vDot2(ba, ba), 0.0, 1.0);
	if(T) *T = t;
	return vMag2(vSub2(pa, vScale2(ba, t)));
}

float vDistTPointLine3(Vector3 p, Line3 ls, float* T) {
	Vector3 pa = vSub3(p, ls.start);
	Vector3 ba = vSub3(ls.end, ls.start);
	
	float t = fclamp(vDot3(pa, ba) / vDot3(ba, ba), 0.0, 1.0);
	if(T) *T = t;
	return vMag3(vSub3(pa, vScale3(ba, t)));
}

// ----




float projPointLine2(Vector2 p, Line2 ls) {
	Vector2 pa = vSub2(p, ls.start);
	Vector2 ba = vSub2(ls.end, ls.start);
	
	return vDot2(pa, ba) / vDot2(ba, ba);
}





float distLineLine3(Line3* a, Line3* b) {
	
	Vector3 ea = vSub3(a->end, a->start);
	Vector3 eb = vSub3(b->end, b->start);
	Vector3 q = vSub(b->start, a->start);
	
	
	float vaa = vLenSq3(ea);
	float vbb = vLenSq3(eb);
	float vba = vDot3(ea, eb);
	float vba_a = vDot3(q, ea);
	
	float den = vba * vba - vbb * vaa;
	
	float ta, tb;
	if(fabs(den) < 1e-6) {
		ta = 0;
		tb = -vba_a / vba; // vba can never be zero here
	}
	else {
		float vba_b = vDot3(q, eb);
		
		ta = (vba_b * vba - vbb * vba_a) / den;
		tb = (-vba_a * vba + vaa * vba_b) / den;
	}
	
	ta = fclamp(0, 1, ta);
	tb = fclamp(0, 1, tb);
	
	Vector3 pa = vAdd(a->start, vScale(ea, ta));
	Vector3 pb = vAdd(b->start, vScale(eb, tb));

	return vDist(pa, pb);
}


/*
float distLineLine3Slow(Line3* a, Line3* b) {
	
	Vector3 ea = vSub3(a->end, a->start);
	Vector3 eb = vSub3(b->end, b->start);
	
	Vector3 n = vCross3(ea, eb);
	float nsq = vLenSq3(n);
	
	// TODO: handle parallel lines
	
	vec3 b1ma1 = vSub(b->start, a->start);
	
	float ta = vDot3(vCross3(eb, n), b1ma1) / nsq;
	float tb = vDot3(vCross3(ea, n), b1ma1) / nsq;
	
	ta = fclamp(0, 1, ta);
	tb = fclamp(0, 1, tb);
	
	vec3 pa = vAdd(a->start, vScale(ea, ta));
	vec3 pb = vAdd(b->start, vScale(eb, tb));
	
	return vDist3(pa, pb);
}
*/



//
Line2 shortestLineFromLineToLine2(Line2* a, Line2* b) {
	
	float s, t;
	vec2 da = vSub(a->end, a->start);
	vec2 db = vSub(b->end, b->start);
	// BUG: if a or b has zero length, it returns nans
	
	// first check if they intersect
	float det = vCross2(da, db);
	if(fabsf(det) > 1e-5f) {
	
		float invdet = 1.f / det;
	    
		float x02 = a->start.x - b->start.x;
		float y02 = a->start.y - b->start.y;
		
		s = (da.x * y02 - da.y * x02) * invdet;
		if(s >= 0.f && s <= 1.f) {
			t = (db.x * y02 - db.y * x02) * invdet;
			if(t >= 0.f && t <= 1.f) {
				vec2 pa = vAdd(a->start, vScale(da, t));
				
				return (Line2){pa, pa};
			}
		}
	}	
	
	
	// otherwise get the closest point
	vec2 ab = vSub(b->start, a->start);
	
	float d_abda = vDot(ab, da);
	float d_abdb = vDot(ab, db);
	float d_dadb = vDot(da, db);
	float l2_da = vLenSq(da);
	float l2_db = vLenSq(db);
	float invl2_da = 1.f / l2_da;
	
	float den = l2_da * l2_db - d_dadb * d_dadb;
	
	if(den < 1e-6f * l2_db * l2_da) {
		s = fclamp(d_abda * invl2_da, 0, 1);
		t = 0;
	}
	else {
		float invden = 1.f / den;
		s = fclamp((d_abda * l2_db - d_abdb * d_dadb) * invden, 0, 1);
		t = fclamp((d_abda * d_dadb - d_abdb * l2_da) * invden, 0, 1);
	}
	
	return (Line2){
		vAdd(a->start, vScale(da, fclamp((t * d_dadb + d_abda) * invl2_da, 0, 1))),
		vAdd(b->start, vScale(db, fclamp((s * d_dadb - d_abdb) / l2_db, 0, 1)))
	};
}



////
// BUG, maybe: The 2d port of this gave wrong results 
Line3 shortestLineFromLineToLine(Line3* a, Line3* b) {
	
	Vector3 ea = vSub3(a->end, a->start);
	Vector3 eb = vSub3(b->end, b->start);
	Vector3 q = vSub(b->start, a->start);
	
	float vaa = vLenSq3(ea);
	float vbb = vLenSq3(eb);
	
	float vba = vDot3(ea, eb);
	float vba_a = vDot3(q, ea);
	
	float den = vba * vba - vbb * vaa;
	
	float ta, tb;
	if(fabs(den) < 1e-6) {
		ta = 0;
		tb = -vba_a / vba; // vba can never be zero here
	}
	else {
		float vba_b = vDot3(q, eb);
		
		ta = (vba_b * vba - vbb * vba_a) / den;
		tb = (-vba_a * vba + vaa * vba_b) / den;
	}
	
	ta = fclamp(ta, 0, 1);
	tb = fclamp(tb, 0, 1);
	
	Vector3 pa = vAdd(a->start, vScale(ea, ta));
	Vector3 pb = vAdd(b->start, vScale(eb, tb));

	return (Line3){pa, pb};
}




// C3DLAS_COPLANAR, _PARALLEL, _INTERSECT, or _DISJOINT
// aboveCnt and belowCnt are always set.
int linePlaneClip3p(
	Vector3* la, 
	Vector3* lb, 
	Plane* pl, 
	Vector3* aboveOut, 
	Vector3* belowOut,
	int* aboveCnt,
	int* belowCnt
) {
	Vector3 ldir, c;
	float da, db;
	
	vSub3p(lb, la, &ldir);
	
	da = vDot3p(la, &pl->n) - pl->d;
	
	// bail if the line and plane are parallel
	if(fabs(vDot3p(&pl->n, &ldir)) < FLT_CMP_EPSILON) {
		*aboveCnt = 0;
		*belowCnt = 0;
			
		// check coplanarity
		if(fabs(da) < FLT_CMP_EPSILON) {
			return C3DLAS_COPLANAR; // the end is on the plane, so the other is too
		}
		
		return C3DLAS_PARALLEL;
	}
	

	db = vDot3p(lb, &pl->n) - pl->d;
	
	// check if one of the points is on the plane
	if(fabs(da) < FLT_CMP_EPSILON) {
		if(db > 0) {
			aboveOut[0] = *la; // correct ordering
			aboveOut[1] = *lb;
			*aboveCnt = 1;
			*belowCnt = 0;
		}
		else {
			belowOut[0] = *la; // correct ordering
			belowOut[1] = *lb;
			*aboveCnt = 0;
			*belowCnt = 1;
		}
		
		return C3DLAS_INTERSECT;
	}
	if(fabs(db) < FLT_CMP_EPSILON) {
		if(da > 0) {
			aboveOut[0] = *la; // correct ordering
			aboveOut[1] = *lb;
			*aboveCnt = 1;
			*belowCnt = 0;
		}
		else {
			belowOut[0] = *la; // correct ordering
			belowOut[1] = *lb;
			*aboveCnt = 0;
			*belowCnt = 1;
		}
		
		return C3DLAS_INTERSECT;
	}
	
	// calculate itnersection point, c
	Vector3 p0, g, j;
	vScale3p(&pl->n, pl->d, &p0);
	vSub3p(&p0, la, &g);
	float h = vDot3p(&g, &pl->n);
	float i = vDot3p(&ldir, &pl->n);
	float d = i != 0 ? h / i : 0;
	
	// check if the plane intersects outside the two points
	if(d < 0 || d > vDist3p(la, lb)) {
		if(da > 0) {
			aboveOut[0] = *la; // correct ordering
			aboveOut[1] = *lb;
			*aboveCnt = 1;
			*belowCnt = 0;
		}
		else {
			belowOut[0] = *la; // correct ordering
			belowOut[1] = *lb;
			*aboveCnt = 0;
			*belowCnt = 1;
		}
		
		return C3DLAS_DISJOINT;
	}
	
	vScale3p(&ldir, d, &j);
	vAdd3p(la, &j, &c);
	
	if(da > 0) {
		aboveOut[0] = *la; // correct ordering
		aboveOut[1] = c;
		belowOut[0] = c;
		belowOut[1] = *lb;
	}
	else {
		belowOut[0] = *la; // correct ordering
		belowOut[1] = c;
		aboveOut[0] = c;
		aboveOut[1] = *lb;
	}
	
	*aboveCnt = 1;
	*belowCnt = 1;
	
	return C3DLAS_INTERSECT;
}


