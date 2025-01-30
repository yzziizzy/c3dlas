



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
		.i    = r.real*l.i    + r.i*l.real + l.j*r.k    - l.k*r.j,
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
Quaternion qRot(Quaternion l, Quaternion r) {
	return qMul(qConj(l), qMul(r, l));
}

Vector3 qRot3(Vector3 a, Quaternion r) {
	Vector4 a4 = {a.x, a.y, a.z, 0.0};
	a4 = qMul(qMul(r, a4), qConj(r));
	return (Vector3){a4.x, a4.y, a4.z};
}

Vector2 qRot2(Vector2 a, Quaternion r) {
	Vector4 a4 = {a.x, a.y, 0.f, 0.f};
	a4 = qMul(qMul(r, a4), qConj(r));
	return (Vector2){a4.x, a4.y};
}


/*
Quaternion qConjugation(Quaternion a, Quaternion r) {
	return qMul(qMul(r, a), qInv(r));
}
*/

// create a quaternion representing a rotation of theta around vector r
Quaternion qFromRTheta(Vector3 r, float theta) {
	
	float tm2p = fmod(fmod(theta, F_2PI) + F_2PI, F_2PI);
	if(tm2p < 0.0001 || tm2p > 6.2831) {
		return (Quaternion){
			.real = 1,
			.i = 0,
			.j = 0,
			.k = 0
		};
	}
	
	float t_2 = tm2p / 2.0f;
	float st_2, ct_2;
	sincosf(t_2, &st_2, &ct_2);
	
	return (Quaternion){
		.real = ct_2,
		.i = st_2 * r.x,
		.j = st_2 * r.y,
		.k = st_2 * r.z
	};
}

// transform a quaternion back into axis-angle representation
// the axis may not be numerically stable if the rotation represented by the quaternion is very small
//    this is an inherent property of quaterions, not a deficiency in the function below
void qToRTheta(Quaternion q, Vector3* r, float* theta) {
//	float l2 = q.x * q.x + q.y * q.y + q.z * q.z;
	
	float l = sqrtf(1.f - (q.w * q.w));

	*theta = 2.f * acos(q.w); // original algorithm from the internet
	
	if(fabs(l) > 0.0001) {
		
//		*theta = -2.f * acos(q.w); // algorithm that works; maybe left-handed vs right-handed?
		r->x = q.x / l;
		r->y = q.y / l;
		r->z = q.z / l;
	}
	else {
//		*theta = 0;
		r->x = q.x;
		r->y = q.y;
		r->z = q.z;	
	}
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
	
	if (d < 0.0f) { // reverse direction for the shorter way
		d = -d;
		b.x = -b.x;
		b.y = -b.y;
		b.z = -b.z;
		b.w = -b.w;
	}

	float th = acosf(d);
	float inv_s_th = 1.0 / sinf(th);
	
	float ca = sin((1.0 - t) * th) * inv_s_th;
	float cb = sin(t * th) * inv_s_th;
	
	return (Quaternion){
		.x = a.x * ca + b.x * cb,
		.y = a.y * ca + b.y * cb,
		.z = a.z * ca + b.z * cb,
		.w = a.w * ca + b.w * cb,
	};
}


Quaternion qNlerp(Quaternion a, Quaternion b, float t) {
	return qNorm(vLerp4(a, b, t));
}



float qAngleBetween(Quaternion a, Quaternion b) {
	float g = vDot4(a, b);
	return acos(2.f * g * g - 1.f);
}


// WARNING: these vectors cannot represent a reflection about any axis.
// det(M4(bx,by,bz)) cannot equal -1, even if the vectors are all orthogonal.
// Very few places on the internet talk about this, and never with sufficient
//   emphasis, espectially not SO. The matrix formed by the vectors must be 
//   "special orthogonal", not simply orthogonal to each other.
//
// If your quaternions are coming out invalid, especially after some cross 
//   products or 90-degree rotations through swapping, try inverting one axis.
//
Quaternion qFromBasis(Vector3 bx, Vector3 by, Vector3 bz) {

	Quaternion q;
	
	float trace = bx.x + by.y + bz.z;
	
	if(trace > 0) { 
		float s = sqrtf(trace + 1.f) * 2.f;
		float invs = 1.f / s;
		q.w = 0.25f * s;
		q.x = (by.z - bz.y) * invs;
		q.y = (bz.x - bx.z) * invs; 
		q.z = (bx.y - by.x) * invs; 
	} else if((bx.x > by.y) && (bx.x > bz.z)) { 
		float s = sqrtf(1.0f + bx.x - by.y - bz.z) * 2.f;
		float invs = 1.f / s;
		q.w = (by.z - bz.y) * invs;
		q.x = 0.25f * s;
		q.y = (by.x + bx.y) * invs; 
		q.z = (bz.x + bx.z) * invs; 
	} else if(by.y > bz.z) { 
		float s = sqrtf(1.0f + by.y - bx.x - bz.z) * 2.f;
		float invs = 1.f / s;
		q.w = (bz.x - bx.z) * invs;
		q.x = (by.x + bx.y) * invs; 
		q.y = 0.25f * s;
		q.z = (bz.y + by.z) * invs; 
	} else { 
		float s = sqrtf(1.0f + bz.z - bx.x - by.y) * 2.f;
		float invs = 1.f / s;
		q.w = (bx.y - by.x) * invs;
		q.x = (bz.x + bx.z) * invs;
		q.y = (bz.y + by.z) * invs;
		q.z = 0.25f * s;
	}
	
	return q;
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
	
	float x2 = q.i + q.i;
	float y2 = q.j + q.j;
	float z2 = q.k + q.k;
	
	float xx = q.i * x2;
	float yy = q.j * y2;
	float zz = q.k * z2;
	
	float sx = q.real * x2;
	float sy = q.real * y2;
	float sz = q.real * z2;
	
	float yx = q.j * x2;
	float zx = q.k * x2;
	float zy = q.k * y2;
	
	
	out->m[0] = 1.0f - yy - zz;
	out->m[1] = yx + sz;
	out->m[2] = zx - sy;
	out->m[3] = 0;
	
	out->m[4] = yx - sz;
	out->m[5] = 1.0f - xx - zz;
	out->m[6] = zy + sx;
	out->m[7] = 0;
	
	out->m[8] = zx + sy;
	out->m[9] = zy - sx;
	out->m[10] = 1.0f - xx - yy;
	out->m[11] = 0;
	
	out->m[12] = 0;
	out->m[13] = 0;
	out->m[14] = 0;
	out->m[15] = 1.0f;
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


