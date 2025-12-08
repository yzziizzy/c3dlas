
void frustumFromMatrix(Matrix* m, Frustum* out) {
	
	Matrix inv;
	
	mInverse(m, &inv);
	
	// first the points
	// these MUST be in this order
	// near
	vMatrixMulf3p(-1,-1,-1, &inv, &out->points[0]);
	vMatrixMulf3p(-1, 1,-1, &inv, &out->points[1]);
	vMatrixMulf3p( 1,-1,-1, &inv, &out->points[2]);
	vMatrixMulf3p( 1, 1,-1, &inv, &out->points[3]);
	// far
	vMatrixMulf3p(-1,-1, 1, &inv, &out->points[4]);
	vMatrixMulf3p(-1, 1, 1, &inv, &out->points[5]);
	vMatrixMulf3p( 1,-1, 1, &inv, &out->points[6]);
	vMatrixMulf3p( 1, 1, 1, &inv, &out->points[7]);
	
	// now the planes
	// near and far
	planeFromTriangle3p(&out->points[0], &out->points[1], &out->points[2], &out->planes[0]);
	planeFromTriangle3p(&out->points[4], &out->points[5], &out->points[6], &out->planes[1]);
	// sides
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[1], &out->planes[2]);
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[2], &out->planes[3]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[1], &out->planes[4]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[2], &out->planes[5]);
}

void frustumFromMatrixVK(Matrix* m, Frustum* out) {
	
	Matrix inv;
	
	mInverse(m, &inv);
	
	// first the points
	// these MUST be in this order
	// near
	vMatrixMulf3p(-1,-1, 0, &inv, &out->points[0]);
	vMatrixMulf3p(-1, 1, 0, &inv, &out->points[1]);
	vMatrixMulf3p( 1,-1, 0, &inv, &out->points[2]);
	vMatrixMulf3p( 1, 1, 0, &inv, &out->points[3]);
	// far
	vMatrixMulf3p(-1,-1, 1, &inv, &out->points[4]);
	vMatrixMulf3p(-1, 1, 1, &inv, &out->points[5]);
	vMatrixMulf3p( 1,-1, 1, &inv, &out->points[6]);
	vMatrixMulf3p( 1, 1, 1, &inv, &out->points[7]);
	
	// now the planes
	// near and far
	planeFromTriangle3p(&out->points[0], &out->points[1], &out->points[2], &out->planes[0]);
	planeFromTriangle3p(&out->points[4], &out->points[5], &out->points[6], &out->planes[1]);
	// sides
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[1], &out->planes[2]);
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[2], &out->planes[3]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[1], &out->planes[4]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[2], &out->planes[5]);
}

void frustumFromMatrixVK_ZUP(Matrix* m, Frustum* out) {
	
	Matrix inv;
	
	mInverse(m, &inv);
	
	*out = (Frustum){0};
	// first the points
	// these MUST be in this order
	// near
	vMatrixMulf3p(-1,-1, 0, &inv, &out->points[0]); // BUG this order is likely wrong for the planes but results in a sane wireframe.
	vMatrixMulf3p( 1,-1, 0, &inv, &out->points[1]);
	vMatrixMulf3p(-1, 1, 0, &inv, &out->points[2]);
	vMatrixMulf3p( 1, 1, 0, &inv, &out->points[3]);
	// far
	vMatrixMulf3p(-1,-1, 1, &inv, &out->points[4]);
	vMatrixMulf3p(1, -1, 1, &inv, &out->points[5]);
	vMatrixMulf3p( -1,1, 1, &inv, &out->points[6]);
	vMatrixMulf3p( 1, 1, 1, &inv, &out->points[7]);
	
	
	// now the planes
	// near and far
	planeFromTriangle3p(&out->points[0], &out->points[1], &out->points[2], &out->planes[0]);
	planeFromTriangle3p(&out->points[4], &out->points[5], &out->points[6], &out->planes[1]);
	// sides
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[1], &out->planes[2]);
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[2], &out->planes[3]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[1], &out->planes[4]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[2], &out->planes[5]);
}



void frustumFromMatrixVK_RDepth(Matrix* m, Frustum* out) {
	
	Matrix inv;
	
	mInverse(m, &inv);
	
	// first the points
	// these MUST be in this order
	// near
	vMatrixMulf3p(-1,-1, 1, &inv, &out->points[0]);
	vMatrixMulf3p(-1, 1, 1, &inv, &out->points[1]);
	vMatrixMulf3p( 1,-1, 1, &inv, &out->points[2]);
	vMatrixMulf3p( 1, 1, 1, &inv, &out->points[3]);
	// far
	vMatrixMulf3p(-1,-1, 0, &inv, &out->points[4]);
	vMatrixMulf3p(-1, 1, 0, &inv, &out->points[5]);
	vMatrixMulf3p( 1,-1, 0, &inv, &out->points[6]);
	vMatrixMulf3p( 1, 1, 0, &inv, &out->points[7]);
	
	// now the planes
	// near and far
	planeFromTriangle3p(&out->points[0], &out->points[1], &out->points[2], &out->planes[0]);
	planeFromTriangle3p(&out->points[4], &out->points[5], &out->points[6], &out->planes[1]);
	// sides
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[1], &out->planes[2]);
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[2], &out->planes[3]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[1], &out->planes[4]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[2], &out->planes[5]);
}

void frustumFromMatrixVK_ZUP_RDepth(Matrix* m, Frustum* out) {
	
	Matrix inv;
	
	mInverse(m, &inv);
	
	*out = (Frustum){0};
	// first the points
	// these MUST be in this order
	// near
	vMatrixMulf3p(-1,-1, 1, &inv, &out->points[0]); // BUG this order is likely wrong for the planes but results in a sane wireframe.
	vMatrixMulf3p( 1,-1, 1, &inv, &out->points[1]);
	vMatrixMulf3p(-1, 1, 1, &inv, &out->points[2]);
	vMatrixMulf3p( 1, 1, 1, &inv, &out->points[3]);
	// far
	vMatrixMulf3p(-1,-1, 0, &inv, &out->points[4]);
	vMatrixMulf3p(1, -1, 0, &inv, &out->points[5]);
	vMatrixMulf3p( -1,1, 0, &inv, &out->points[6]);
	vMatrixMulf3p( 1, 1, 0, &inv, &out->points[7]);
	
	
	// now the planes
	// near and far
	planeFromTriangle3p(&out->points[0], &out->points[1], &out->points[2], &out->planes[0]);
	planeFromTriangle3p(&out->points[4], &out->points[5], &out->points[6], &out->planes[1]);
	// sides
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[1], &out->planes[2]);
	planeFromTriangle3p(&out->points[0], &out->points[4], &out->points[2], &out->planes[3]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[1], &out->planes[4]);
	planeFromTriangle3p(&out->points[3], &out->points[7], &out->points[2], &out->planes[5]);
}



void frustumCenter(Frustum* f, Vector3* out) {
	Vector3 sum = {0.0f,0.0f,0.0f};
	
	for(int i = 0; i < 8; i++) vAdd3p(&f->points[i], &sum, &sum);
	vScale3p(&sum, 1.0f/8.0f, out);
}

// General idea of the algorithm:
// https://lxjk.github.io/2017/04/15/Calculate-Minimal-Bounding-Sphere-of-Frustum.html
// http://archive.is/YACj2
void frustumBoundingSphere(Frustum* f, Sphere* out) {
	Vector3 f0, n0;
	vPointAvg3p(&f->points[0], &f->points[3], &n0);
	vPointAvg3p(&f->points[4], &f->points[7], &f0);
	
	float Dn2 = vDistSq3p(&n0, &f->points[0]); 
	float Df2 = vDistSq3p(&f0, &f->points[4]);
	
	 // check for ortho
	if(Dn2 - Df2 < 0.00001) {
		frustumCenter(f, &out->center);
		out->r = vDist3p(&out->center, &f->points[0]); 
		return;
	}
	
	float Dnf = vDist3p(&f0, &n0);
	float Dnc = (Dn2 - Df2 - Df2) / (2 * Dnf);
	
// 	printf("\n f: %f,%f,%f\n", f->points[4].x,f->points[4].y,f->points[4].z);
// 	printf(" n: %f,%f,%f\n", f->points[0].x,f->points[0].y,f->points[0].z);
// 	printf(" f0: %f,%f,%f\n", f0.x,f0.y,f0.z);
// 	printf(" n0: %f,%f,%f\n", n0.x,n0.y,n0.z);
// 	printf(" dn2, df2, dnf, dnc: %f,%f,%f,%f\n", Dn2, Df2, Dnf, Dnc);
	
	
	if(Dnc > 0 && Dnc < Dnf) {
		vLerp3p(&f0, &n0, Dnc / Dnf, &out->center);
		out->r = sqrt(Dnc * Dnc + Dn2);
	}
	else {
		out->center = f0;
		out->r = sqrt(Df2);
	}
}


void frustumInscribeSphere(Frustum* f, Sphere* out) {
	Vector3 fx, nx;
	vPointAvg3p(&f->points[0], &f->points[3], &nx);
	vPointAvg3p(&f->points[4], &f->points[7], &fx);
	
/*	
	float Dn2 = vDistSq3p(&n0, &f->points[0]); 
	float Df2 = vDistSq3p(&f0, &f->points[4]);
	float Dnf = vDist3p(&f0, n0);
	float Dnc = (Dn2 - Df2 - Df2) / (2 * Dnf);*/

}


