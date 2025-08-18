


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




