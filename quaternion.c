



Quaternion qAdd(Quaternion l, Quaternion r) {
	return vAdd4(l, r);
}

Quaternion qSub(Quaternion l, Quaternion r) {
	return vSub4(l, r);
}

Quaternion qScale(Quaternion q, float s) {
	return (Quaternion){.i = q.i*s, .j = q.j*s, .k = q.k*s, .real = q.real*s }; 
}

Quaternion qMul(Quaternion l, Quaternion r) {
	return (Quaternion){
		.i    = l.real*r.i    + l.i*r.real + l.j*r.k    - l.k*r.j,
		.j    = l.real*r.j    - l.i*r.k    + l.j*r.real + l.k*r.i,
		.k    = l.real*r.k    + l.i*r.j    - l.j*r.i    + l.k*r.real,
		.real = l.real*r.real - l.i*r.i    - l.j*r.j    - l.k*r.k
	};
}

Quaternion qDiv(Quaternion n, Quaternion d) {
	float m = vDot4(d, d);

	return (Quaternion){
		.i = (d.real*n.i - d.i*n.real - d.k*n.k    + d.k*n.j   ) / m,
		.j = (d.real*n.j + d.i*n.k    - d.k*n.real - d.k*n.i   ) / m,
		.k = (d.real*n.k - d.i*n.j    - d.k*n.i    - d.k*n.real) / m,
		.real = (n.real*d.real + n.i*d.i + n.j*d.j + n.k*d.k) / m
	};
}


// rotate one quaternion by another
// aka "conjugation"
Quaternion qRot(Quaternion a, Quaternion r) {
	return qMul(qMul(r, a), qConj(r));
}

Vector3 qRot3(Vector3 a, Quaternion r) {
	Vector4 a4 = {a.x, a.y, a.z, 0.0};
	a4 = qMul(qMul(r, a4), qConj(r));
	return (Vector3){a4.x, a4.y, a4.z};
}

/*
vec3 quaternion_apply_to_vec3(const Quaternion* q, vec3 v) {
	r64 ix = q->w * v.x + q->y * v.z - q->z * v.y;
	r64 iy = q->w * v.y + q->z * v.x - q->x * v.z;
	r64 iz = q->w * v.z + q->x * v.y - q->y * v.x;
	r64 iw = - q->x * v.x - q->y * v.y - q->z * v.z;

	return (vec3) {
		(ix * q->w) + (iw * -q->x) + (iy * -q->z) - (iz * -q->y),
		(iy * q->w) + (iw * -q->y) + (iz * -q->x) - (ix * -q->z),
		(iz * q->w) + (iw * -q->z) + (ix * -q->y) - (iy * -q->x)
	};
}
*/


Quaternion qConjugation(Quaternion a, Quaternion r) {
	return qMul(qMul(r, a), qInv(r));
}

// create a quaternion representing a rotation of theta around vector r
Quaternion qFromRTheta(Vector3 r, float theta) {
	float t_2 = theta / 2;
	float st_2, ct_2;
	sincosf(t_2, &st_2, &ct_2);
	
	return (Quaternion){
		.real = ct_2,
		.i = st_2 * r.x,
		.j = st_2 * r.y,
		.k = st_2 * r.z
	};
}


// conjugate
Quaternion qConj(Quaternion q) {
	return (Quaternion){
		.i = -q.i,
		.j = -q.j,
		.k = -q.k,
		.real = q.real
	};
}

// inverse
Quaternion qInv(Quaternion q) {
	float m = vDot4(q, q);
	return (Quaternion){
		.i = -q.i / m,
		.j = -q.j / m,
		.k = -q.k / m,
		.real = q.real / m
	};
}

Quaternion qUnit(Quaternion q) {
	return vNorm4(q);
}

// aliases of the same thing for conceptual convenience
float qMod(Quaternion q) {
	return vMag4(q);
}

float qMag(Quaternion q) {
	return vMag4(q);
}

float qLen(Quaternion q) {
	return vMag4(q);
}

Quaternion qNorm(Quaternion q) {
	return vNorm4(q);
}

Quaternion qSlerp(Quaternion a, Quaternion b, float t) {
	float d = vDot4(a, b);
	if(d >= 1.0) return a;
	
	float th = acos(d);
	float inv_s_th = 1.0 / sin(th);
	
	float ca = sin((1.0 - t) * th) * inv_s_th;
	float cb = sin(t * th) * inv_s_th;
	
	return (Quaternion){
		.x = a.x * ca + b.x * cb,
		.y = a.y * ca + b.y * cb,
		.z = a.z * ca + b.z * cb,
		.w = a.w * ca + b.w * cb,
	};
}



// Applies the full conjugate multiplication qvq*
void qNonUnitToMatrix(Quaternion q, Matrix* out) {
	float s2 = q.real * q.real;
	float x2 = q.i * q.i;
	float y2 = q.j * q.j;
	float z2 = q.k * q.k;
	
	float sx = 2.0 * q.real * q.i;
	float sy = 2.0 * q.real * q.j;
	float sz = 2.0 * q.real * q.k;
	float xy = 2.0 * q.i * q.j;
	float xz = 2.0 * q.i * q.k;
	float yz = 2.0 * q.j * q.k;

	out->m[0] = s2 + x2 - y2 - z2;
	out->m[1] = xy + sz;
	out->m[2] = xz - sy;
	out->m[3] = 0;
	
	out->m[4] = xy - sz;
	out->m[5] = s2 - x2 + y2 - z2;
	out->m[6] = yz + sx;
	out->m[7] = 0;
	
	out->m[8] = xz + sy;
	out->m[9] = yz - sx;
	out->m[10] = s2 - x2 - y2 + z2;
	out->m[11] = 0;
	
	out->m[12] = 0;
	out->m[13] = 0;
	out->m[14] = 0;
	out->m[15] = 1;
}


void qUnitToMatrix(Quaternion q, Matrix* out) {
	float x2 = 2.0 * q.i * q.i;
	float y2 = 2.0 * q.j * q.j;
	float z2 = 2.0 * q.k * q.k;
	
	float sx = 2.0 * q.real * q.i;
	float sy = 2.0 * q.real * q.j;
	float sz = 2.0 * q.real * q.k;
	float xy = 2.0 * q.i * q.j;
	float xz = 2.0 * q.i * q.k;
	float yz = 2.0 * q.j * q.k;
	
	out->m[0] = 1.0 - y2 - z2;
	out->m[1] = xy + sz;
	out->m[2] = xz - sy;
	out->m[3] = 0;
	
	out->m[4] = xy - sz;
	out->m[5] = 1.0 - x2 - z2;
	out->m[6] = yz + sx;
	out->m[7] = 0;
	
	out->m[8] = xz + sy;
	out->m[9] = yz - sx;
	out->m[10] = 1.0 - x2 - y2;
	out->m[11] = 0;
	
	out->m[12] = 0;
	out->m[13] = 0;
	out->m[14] = 0;
	out->m[15] = 1;
}


void qUnitToMatrix3(Quaternion q, Matrix3* out) {
	float x2 = 2.0 * q.i * q.i;
	float y2 = 2.0 * q.j * q.j;
	float z2 = 2.0 * q.k * q.k;
	
	float sx = 2.0 * q.real * q.i;
	float sy = 2.0 * q.real * q.j;
	float sz = 2.0 * q.real * q.k;
	float xy = 2.0 * q.i * q.j;
	float xz = 2.0 * q.i * q.k;
	float yz = 2.0 * q.j * q.k;
	
	out->m[0] = 1.0 - y2 - z2;
	out->m[1] = xy + sz;
	out->m[2] = xz - sy;
	
	out->m[3] = xy - sz;
	out->m[4] = 1.0 - x2 - z2;
	out->m[5] = yz + sx;
	
	out->m[6] = xz + sy;
	out->m[7] = yz - sx;
	out->m[8] = 1.0 - x2 - y2;
}


