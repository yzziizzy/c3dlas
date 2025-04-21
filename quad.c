

#include <assert.h>






int intersectQuadPoint2(vec2 p, vec2 q[4]) {
	Vector2 outp;
	Vector2 outd;
	
	// two barycentric coordinate calculations
	vec2 a = q[1];
	vec2 b = q[0];
	vec2 c = q[2];
	vec2 d = q[3];
	
	vec2 ab = vSub(b, a);
	vec2 ac = vSub(c, a);
	vec2 ap = vSub(p, a);
	
	float invden_a = 1.f / (ab.x * ac.y - ac.x * ab.y);
	outp.x = (ap.x * ac.y - ac.x * ap.y) * invden_a;
	outp.y = (ab.x * ap.y - ap.x * ab.y) * invden_a;
	
	// see if p's coords are inside the triangle and less than d's
	if(outp.x < 0 || outp.y < 0) return C3DLAS_DISJOINT;
	
	vec2 db = vSub(b, d);
	vec2 dc = vSub(c, d);
	vec2 dp = vSub(p, d);

	float invden_d = 1.f / (db.x * dc.y - dc.x * db.y);
	outd.x = (dp.x * dc.y - dc.x * dp.y) * invden_d;
	outd.y = (db.x * dp.y - dp.x * db.y) * invden_d;
	
	// see if p's coords are inside the triangle and less than d's
	if(outd.x < 0 || outd.y < 0) return C3DLAS_DISJOINT;
	
	return C3DLAS_INTERSECT;
}




void intersectQuadPoint2_precalc(vec2 q[4], struct intersectQuadPoint2_precalc* pc) {
	// two barycentric coordinate calculations
	vec2 a = q[1];
	vec2 b = q[0];
	vec2 c = q[2];
	vec2 d = q[3];
	pc->a = a;
	pc->d = d;
	
	pc->ab = vSub(b, a);
	pc->ac = vSub(c, a);
	pc->db = vSub(b, d);
	pc->dc = vSub(c, d);
	pc->invden_a = 1.f / (pc->ab.x * pc->ac.y - pc->ac.x * pc->ab.y);
	pc->invden_d = 1.f / (pc->db.x * pc->dc.y - pc->dc.x * pc->db.y);
}


int intersectQuadPoint2_withprecalc(vec2 p, struct intersectQuadPoint2_precalc* pc) {
	
	vec2 ap = vSub(p, pc->a);	
	if((ap.x * pc->ac.y - pc->ac.x * ap.y) * pc->invden_a < 0) return C3DLAS_DISJOINT;
	if((pc->ab.x * ap.y - ap.x * pc->ab.y) * pc->invden_a < 0) return C3DLAS_DISJOINT;
	
	vec2 dp = vSub(p, pc->d);
	if((dp.x * pc->dc.y - pc->dc.x * dp.y) * pc->invden_d < 0) return C3DLAS_DISJOINT;
	if((pc->db.x * dp.y - dp.x * pc->db.y) * pc->invden_d < 0) return C3DLAS_DISJOINT;
	
	return C3DLAS_INTERSECT;
}









