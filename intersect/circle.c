


// out vectors are normilzed directions to the tangent point from the center of the circle
// returns number of answers present
int findCircleTangents(Vector2 c, float r, Vector2 p, Vector2 out[2]) {
	Vector2 q = vSub(p, c);
	
	float r2 = r * r;
	float d2 = vDot(q, q) - r2;
	if(d2 < 0) return 0;  // p is inside the circle
	
	float d = sqrtf(d2); 
	
	out[0] = V2(
		(q.x * r + q.y * d) / (r2 + d2),
		(q.y * r - q.x * d) / (r2 + d2)
	);
	
	out[1] = V2(
		(q.x * r + q.y * -d) / (r2 + d2),
		(q.y * r - q.x * -d) / (r2 + d2)
	);
	
	return 2;
}


// returns the radius
float circleFromPoints(vec2 a, vec2 b, vec2 c, vec2* center) {

	float ycma = c.y - a.y;
	float xamb = a.x - b.x;
	float ybma = b.y - a.y;
	float xamc = a.x - c.x;

	float yamb = a.y - b.y;
	float yamc = a.y - c.y;
	float xcma = c.x - a.x;
	float xbma = b.x - a.x;

	float qdenom = ycma*xamb - ybma*xamc;
	if(fabsf(qdenom) < FLT_EPSILON) {
		// colinear points
		goto COLINEAR;
	}
	
	float rdenom = xcma*yamb - xbma*yamc;
	if(fabsf(rdenom) < FLT_EPSILON) {
		// colinear points
		goto COLINEAR;
	}
	
	float ax2 = a.x * a.x;
	float ay2 = a.y * a.y;
	float bx2 = b.x * b.x;
	float by2 = b.y * b.y;
	float cx2 = c.x * c.x;
	float cy2 = c.y * c.y;
	float x2amc = ax2 - cx2;
	float y2amc = ay2 - cy2;
	float x2bma = bx2 - ax2;
	float y2bma = by2 - ay2;
	
	float cy = (x2amc*xamb + y2amc*xamb + x2bma*xamc + y2bma*xamc) / (-2.f * qdenom);
	float cx = (x2amc*yamb + y2amc*yamb + x2bma*yamc + y2bma*yamc) / (-2.f * rdenom);
	float q = -ax2 - ay2 + 2.f*cx*a.x + 2.f*cy*a.y;
	
	center->x = cx;
	center->y = cy;
	return sqrtf(cx*cx + cy*cy - q);

COLINEAR:
	float xbmc = b.x - c.x;
	float ybmc = b.y - c.y;
	
	float dab2 = xamb*xamb + yamb*yamb; 
	float dbc2 = xbmc*xbmc + ybmc*ybmc;
	float dac2 = xamc*xamc + yamc*yamc;
	
	if(dab2 > dbc2) {
		if(dac2 > dab2) { // dac
			*center = vAvg(a, c);
			return sqrtf(dac2) * .5f;
		}
		else { // dab
			*center = vAvg(a, b);
			return sqrtf(dab2) * .5f;
		}
	}
	else {
		if(dac2 > dbc2) {// dac
			*center = vAvg(a, c);
			return sqrtf(dac2) * .5f;
		}
		else { // dbc
			*center = vAvg(b, c);
			return sqrtf(dbc2) * .5f;
		}
	}
}



// returns C3DLAS_INTERSECT or C3DLAS_DISJOINT
int pointInsideCircleFromPoints(vec2 p1, vec2 p2, vec2 p3, vec2 test) {
	
	float a = p1.x - test.x;
	float b = p1.y - test.y;
	float c = a * a + b * b;
	float d = p2.x - test.x;
	float e = p2.y - test.y;
	float f = d * d + e * e;
	float g = p3.x - test.x;
	float h = p3.y - test.y;
	float i = g * g + h * h;

	float det = a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;
	
	return det < 0 ? C3DLAS_INTERSECT : C3DLAS_DISJOINT;
}
	
	
	

