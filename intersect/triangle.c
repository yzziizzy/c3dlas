

Vector3 baryCoords2(Vector2 p, Vector2 a, Vector2 b, Vector2 c) {

	vec2 ab = vSub(b, a);
	vec2 ac = vSub(c, a);
	vec2 ap = vSub(p, a);
	float invden = 1.f / (ab.x * ac.y - ac.x * ab.y);
	
	Vector3 out;
	
	out.y = (ap.x * ac.y - ac.x * ap.y) * invden;
	out.z = (ab.x * ap.y - ap.x * ab.y) * invden;
	out.x = 1.0f - out.y - out.z;
	
	return out;
}


// returns the *signed* area of a triangle. useful for determining winding
// positive values mean a clockwise triangle
float triArea2p(Vector2* a, Vector2* b, Vector2* c) {
	return 0.5 * (
		((b->x - a->x) * (b->y + a->y)) +
		((c->x - b->x) * (c->y + b->y)) +
		((a->x - c->x) * (a->y + c->y)));
}


// determines if a point is inside a triangle
int triPointInside2p(Vector2* p, Vector2* a, Vector2* b, Vector2* c) {
	int d = signbit((p->x - b->x) * (a->y - b->y) - (a->x - b->x) * (p->y - b->y));
	int e = signbit((p->x - c->x) * (b->y - c->y) - (b->x - c->x) * (p->y - c->y));
	if(d != e) return 0;
	int f = signbit((p->x - a->x) * (c->y - a->y) - (c->x - a->x) * (p->y - a->y));
	return e == f;
}




float distLine2Triangle2(Line2 a, Vector2 tri[3]) {
	
	// check the lines around the outside
	float d[3];
	
	for(int i = 0; i < 3; i++) {
		d[i] = distLine2Line2(a, (Line2){tri[i], tri[(i + 1) % 3]});
		if(d[i] == 0) return 0;
	}
	
	// check if the line is entirely inside the triangle,
	//   if one point (and therefore both points) is inside 
	
	if(triPointInside2p(&a.start, tri, tri + 1, tri + 2)) {
		return 0;
	}
	
	// the line is outside; one of the corners is closest, find which one
	return MIN(MIN(d[0], d[1]), d[2]);
}


float distPoint2Triangle2(Vector2 a, Vector2 tri[3]) {
	
	// check if the point is inside the triangle,
	if(triPointInside2p(&a, tri, tri + 1, tri + 2)) {
		return 0;
	}
	
	// check the lines around the outside
	float d[3];
	
	for(int i = 0; i < 3; i++) {
		d[i] = vDistPointLine2(a, (Line2){tri[i], tri[(i + 1) % 3]});
		if(d[i] == 0) return 0;
	}
	
	// the point is outside; one of the corners is closest, find which one
	return MIN(MIN(d[0], d[1]), d[2]);
}


