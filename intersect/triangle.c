



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

