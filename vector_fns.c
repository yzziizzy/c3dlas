






// Coordinate system conversions

// Does not check for degenerate vectors
// Cartesian to Spherical
Vector3 vC2S3(Vector3 cart) {
	Vector3 sp;
	sp.rho = vMag3(cart);
	sp.theta = atan2f(cart.x, cart.y);
	sp.phi = acosf(cart.z / sp.rho);
	
	return sp;
}

// Spherical to Cartesian
Vector3 vS2C3(Vector3 s) {
	float st, ct, sp, cp;
	
	// as of July 2022, gcc trunk is smart enough to automatically optimize to sincos, but clang isn't. 
	sincosf(s.phi, &sp, &cp);
	sincosf(s.theta, &st, &ct);

	return (Vector3){
		.x = s.rho * sp * ct,
		.y = s.rho * sp * st,
		.z = s.rho * cp
	};
}



// Some functions only exist in certain dimensions and base types:



Vector4i vFloor4(const Vector4 v) {
	return (Vector4i){.x = floorf(v.x), .y = floorf(v.y), .z = floorf(v.z), .w = floorf(v.w)};
}

Vector4i vCeil4(const Vector4 v) {
	return (Vector4i){.x = ceilf(v.x), .y = ceilf(v.y), .z = ceilf(v.z), .w = ceilf(v.w)};
}

Vector3i vFloor3(const Vector3 v) {
	return (Vector3i){.x = floorf(v.x), .y = floorf(v.y), .z = floorf(v.z)};
}

Vector3i vCeil3(const Vector3 v) {
	return (Vector3i){.x = ceilf(v.x), .y = ceilf(v.y), .z = ceilf(v.z)};
}

Vector2i vFloor2(const Vector2 v) {
	return (Vector2i){.x = floorf(v.x), .y = floorf(v.y)};
}

Vector2i vCeil2(const Vector2 v) {
	return (Vector2i){.x = ceilf(v.x), .y = ceilf(v.y)};
}



Vector4l vFloor4d(const Vector4d v) {
	return (Vector4l){.x = floor(v.x), .y = floor(v.y), .z = floor(v.z), .w = floor(v.w)};
}

Vector4l vCeil4d(const Vector4d v) {
	return (Vector4l){.x = ceil(v.x), .y = ceil(v.y), .z = ceil(v.z), .w = ceil(v.w)};
}

Vector3l vFloor3d(const Vector3d v) {
	return (Vector3l){.x = floor(v.x), .y = floor(v.y), .z = floor(v.z)};
}

Vector3l vCeil3d(const Vector3d v) {
	return (Vector3l){.x = ceil(v.x), .y = ceil(v.y), .z = ceil(v.z)};
}

Vector2l vFloor2d(const Vector2d v) {
	return (Vector2l){.x = floor(v.x), .y = floor(v.y)};
}

Vector2l vCeil2d(const Vector2d v) {
	return (Vector2l){.x = ceil(v.x), .y = ceil(v.y)};
}




// Cross product: out = a x b
// Cross products *proper* only exist in 3 and 7 dimensions
Vector3 vCross3(Vector3 a, Vector3 b) {
	return (Vector3) {
		.x = (a.y * b.z) - (a.z * b.y),
		.y = (a.z * b.x) - (a.x * b.z),
		.z = (a.x * b.y) - (a.y * b.x)
	};
}

void vCross3p(Vector3* a, Vector3* b, Vector3* out) {
	out->x = (a->y * b->z) - (a->z * b->y);
	out->y = (a->z * b->x) - (a->x * b->z);
	out->z = (a->x * b->y) - (a->y * b->x);
}

Vector3d vCross3d(Vector3d a, Vector3d b) {
	return (Vector3d) {
		.x = (a.y * b.z) - (a.z * b.y),
		.y = (a.z * b.x) - (a.x * b.z),
		.z = (a.x * b.y) - (a.y * b.x)
	};
}

void vCross3dp(Vector3d* a, Vector3d* b, Vector3d* out) {
	out->x = (a->y * b->z) - (a->z * b->y);
	out->y = (a->z * b->x) - (a->x * b->z);
	out->z = (a->x * b->y) - (a->y * b->x);
}


// ... however, if you apply it to two dimensions it yields
//  the sine of the angle between the two vectors, and the sign
//  of the value determines which side of a that b is on. If
//  the value is zero, the vectors are parallel.

float vCross2(Vector2 a, Vector2 b) {
	return (a.x * b.y) - (a.y * b.x);
}

float vCross2p(Vector2* a, Vector2* b) {
	return (a->x * b->y) - (a->y * b->x);
}

double vCross2d(Vector2d a, Vector2d b) {
	return (a.x * b.y) - (a.y * b.x);
}

double vCross2dp(Vector2d* a, Vector2d* b) {
	return (a->x * b->y) - (a->y * b->x);
}





// Scalar triple product: a . (b x c)
float vScalarTriple3(Vector3 a, Vector3 b, Vector3 c) {
	return (float)((a.x * b.y * c.z) - (a.x * b.z * c.y) - (a.y * b.x * c.z)
				 + (a.z * b.x * c.y) + (a.y * b.z * c.x) - (a.z * b.y * c.x));
}

float vScalarTriple3p(Vector3* a, Vector3* b, Vector3* c) {
	return (float)((a->x * b->y * c->z) - (a->x * b->z * c->y) - (a->y * b->x * c->z)
				 + (a->z * b->x * c->y) + (a->y * b->z * c->x) - (a->z * b->y * c->x));
}

double vScalarTriple3d(Vector3d a, Vector3d b, Vector3d c) {
	return (double)((a.x * b.y * c.z) - (a.x * b.z * c.y) - (a.y * b.x * c.z)
				 + (a.z * b.x * c.y) + (a.y * b.z * c.x) - (a.z * b.y * c.x));
}

double vScalarTriple3dp(Vector3d* a, Vector3d* b, Vector3d* c) {
	return (double)((a->x * b->y * c->z) - (a->x * b->z * c->y) - (a->y * b->x * c->z)
				 + (a->z * b->x * c->y) + (a->y * b->z * c->x) - (a->z * b->y * c->x));
}





// Returns an arbitrary unit vector perpendicular to the input
// The input vector does not need to be normalized
void vPerp3p(Vector3* n, Vector3* out) {
	*out = vPerp3(*n);
}

Vector3 vPerp3(Vector3 n) {
	float f, d;
	
	float absx = fabs(n.x);
	float absy = fabs(n.y);
	
	if(absx < absy) {
		if(n.x < n.z) goto MIN_X;
		goto MIN_Z;
	}
	if(absy < fabs(n.z)) goto MIN_Y;
	goto MIN_Z;
	
MIN_X: 
	d = 1.0f / sqrtf(n.z * n.z + n.y * n.y); 
	f = n.z;
	n.z = n.y * d;
	n.y = -f * d;
	n.x = 0;
	return n;
	
MIN_Y: 
	d = 1.0f / sqrtf(n.z * n.z + n.x * n.x); 
	f = n.x;
	n.x = n.z * d;
	n.z = -f * d;
	n.y = 0;
	return n;
	
MIN_Z: 
	d = 1.0f / sqrtf(n.x * n.x + n.y * n.y); 
	f = n.y;
	n.y = n.x * d;
	n.x = -f * d;
	n.z = 0;
	return n;
}


// Returns an arbitrary unit vector perpendicular to the input
// The input vector does not need to be normalized
void vPerp2p(Vector2* n, Vector2* out) {
	*out = vPerp2(*n);
}

Vector2 vPerp2(Vector2 n) {
	return vNorm2((Vector2){.x = -n.y, .y = n.x});
}

Vector2 vPerpCW2(Vector2 n) {
	return vNorm2((Vector2){.x = n.y, .y = -n.x});
}
Vector2 vPerpCCW2(Vector2 n) {
	return vNorm2((Vector2){.x = -n.y, .y = n.x});
}





// a and b must be normalized
float vAngleBetween2(Vector2 a, Vector2 b) {
	return atan2f(vCross2(a, b), vDot2(a, b));
}

// a and b must be normalized
float vAngleBetween3(Vector3 a, Vector3 b) {
	float theta = acosf(vDot3(a, b));
	Vector3 cross = vCross3(a, b);
	
	Quaternion q = qFromRTheta(cross, theta);
	
	Vector3 check = qRot3(a, q);
	
	return vEq3(check, b) ? theta : -theta;
}





Vector2 vModPositive2(Vector2 v, Vector2 m) {
	return (Vector2){
		.x = fmodf(fmodf(v.x, m.x) + m.x, m.x),
		.y = fmodf(fmodf(v.y, m.y) + m.y, m.y),
	};
}

Vector3 vModPositive3(Vector3 v, Vector3 m) {
	return (Vector3){
		.x = fmodf(fmodf(v.x, m.x) + m.x, m.x),
		.y = fmodf(fmodf(v.y, m.y) + m.y, m.y),
		.z = fmodf(fmodf(v.z, m.z) + m.z, m.z),
	};
}

Vector4 vModPositive4(Vector4 v, Vector4 m) {
	return (Vector4){
		.x = fmodf(fmodf(v.x, m.x) + m.x, m.x),
		.y = fmodf(fmodf(v.y, m.y) + m.y, m.y),
		.z = fmodf(fmodf(v.z, m.z) + m.z, m.z),
		.w = fmodf(fmodf(v.w, m.w) + m.w, m.w),
	};
}



// All the basics, implemented as macros out of temporary convenience


#define FN(sz, suf, ty, ft, sufft, pref, ...) \
\
int vEqExact##suf(const Vector##suf a, const Vector##suf b) { \
	return vEqExact##suf##p(&a, &b); \
} \
int vEqExact##suf##p(const Vector##suf const * a, const Vector##suf const * b) { \
	int tmp = 0; \
	for(int i = 0; i < sz; i++) \
		tmp += ((ty*)a)[i] == ((ty*)b)[i]; \
	return tmp == sz; \
} \
\
int vEq##suf(const Vector##suf a, const Vector##suf b) { \
	return vEqEp##suf(a, b, pref##_CMP_EPSILON); \
} \
int vEq##suf##p(const Vector##suf* a, const Vector##suf* b) { \
	return vEqEp##suf(*a, *b, pref##_CMP_EPSILON); \
} \
\
int vEqEp##suf(const Vector##suf a, const Vector##suf b, ft epsilon) { \
	return vEqEp##suf##p(&a, &b, epsilon); \
} \
int vEqEp##suf##p(const Vector##suf* a, const Vector##suf* b, ft epsilon) { \
	return vDistSq##suf(*a, *b) <= epsilon * epsilon; \
} \
\
ft vDistSq##suf(const Vector##suf a, const Vector##suf b) { \
	return vDistSq##suf##p(&a, &b); \
} \
ft vDistSq##suf##p(const Vector##suf* a, const Vector##suf* b) { \
	ft tmp = 0; \
	for(int i = 0; i < sz; i++) { \
		ft q = ((ty*)a)[i] - ((ty*)b)[i]; \
		tmp += q * q; \
	} \
	return tmp;\
} \
ft vDist##suf(const Vector##suf a, const Vector##suf b) { \
	return sqrt(vDistSq##suf##p(&a, &b)); \
} \
ft vDist##suf##p(const Vector##suf* a, const Vector##suf* b) { \
	return sqrt(vDistSq##suf##p(a, b)); \
} \
\
Vector##suf vAdd##suf(const Vector##suf a, const Vector##suf b) { \
	Vector##suf out; \
	vAdd##suf##p(&a, &b, &out); \
	return out; \
} \
void vAdd##suf##p(const Vector##suf* a, const Vector##suf* b, Vector##suf* out) { \
	for(int i = 0; i < sz; i++) { \
		((ty*)out)[i] = ((ty*)a)[i] + ((ty*)b)[i]; \
	} \
} \
\
Vector##suf vSub##suf(const Vector##suf a, const Vector##suf b) { \
	Vector##suf out; \
	vSub##suf##p(&a, &b, &out); \
	return out; \
} \
void vSub##suf##p(const Vector##suf const * a, const Vector##suf const * b, Vector##suf* out) { \
	for(int i = 0; i < sz; i++) { \
		((ty*)out)[i] = ((ty*)a)[i] - ((ty*)b)[i]; \
	} \
} \
\
Vector##suf vMul##suf(const Vector##suf a, const Vector##suf b) { \
	Vector##suf out; \
	vMul##suf##p(&a, &b, &out); \
	return out; \
} \
void vMul##suf##p(const Vector##suf const * a, const Vector##suf const * b, Vector##suf* out) { \
	for(int i = 0; i < sz; i++) { \
		((ty*)out)[i] = ((ty*)a)[i] * ((ty*)b)[i]; \
	} \
} \
\
Vector##suf vDiv##suf(const Vector##suf top, const Vector##suf bot) { \
	Vector##suf out; \
	vDiv##suf##p(&top, &bot, &out); \
	return out; \
} \
void vDiv##suf##p(const Vector##suf const * top, const Vector##suf const * bot, Vector##suf* out) { \
	for(int i = 0; i < sz; i++) { \
		((ty*)out)[i] = ((ty*)top)[i] / ((ty*)bot)[i]; \
	} \
} \
\
ft vDot##suf(const Vector##suf a, const Vector##suf b) { \
	ft tmp = 0; \
	for(int i = 0; i < sz; i++) { \
		tmp += a.f[i] * b.f[i]; \
	} \
	return tmp;\
} \
ft vDot##suf##p(const Vector##suf* a, const Vector##suf* b) { \
	ft tmp = 0; \
	for(int i = 0; i < sz; i++) { \
		tmp += a->f[i] * b->f[i]; \
	} \
	return tmp;\
} \
\
Vector##sufft vScale##suf(const Vector##suf v, ft scalar) { \
	Vector##sufft out; \
	vScale##suf##p(&v, scalar, &out); \
	return out; \
} \
void vScale##suf##p(const Vector##suf* v, ft scalar, Vector##sufft* out) { \
	for(int i = 0; i < sz; i++) \
		((ft*)out)[i] = (ft)((ty*)v)[i] * scalar; \
} \
\
\
Vector##sufft vFMA##suf(const Vector##suf a, const Vector##suf m, ft scalar) { \
	Vector##sufft out; \
	vFMA##suf##p(&a, &m, scalar, &out); \
	return out; \
} \
void vFMA##suf##p(const Vector##suf* a, const Vector##suf* m, ft scalar, Vector##sufft* out) { \
	for(int i = 0; i < sz; i++) \
		((ft*)out)[i] = (ft)((ty*)a)[i] + (ft)((ty*)m)[i] * scalar; \
} \
\
\
Vector##sufft vAvg##suf(const Vector##suf a, const Vector##suf b) { \
	Vector##sufft out; \
	vAvg##suf##p(&a, &b, &out); \
	return out; \
} \
void vAvg##suf##p(const Vector##suf* a, const Vector##suf* b, Vector##sufft* out) { \
	for(int i = 0; i < sz; i++) { \
		((ty*)out)[i] = (((ty*)a)[i] + ((ty*)b)[i]) / (ft)2.0; \
	} \
} \
\
Vector##suf vNeg##suf(const Vector##suf v) { \
	Vector##suf out; \
	vNeg##suf##p(&v, &out); \
	return out; \
} \
void vNeg##suf##p(const Vector##suf* v, Vector##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((ty*)out)[i] = -((ty*)v)[i]; \
} \
\
Vector##sufft vLerp##suf(const Vector##suf a, const Vector##suf b, ft t) { \
	Vector##sufft out; \
	vLerp##suf##p(&a, &b, t, &out); \
	return out; \
} \
void vLerp##suf##p(const Vector##suf* a, const Vector##suf* b, ft t, Vector##sufft* out) { \
	for(int i = 0; i < sz; i++) \
		((ft*)out)[i] = (ft)((ty*)a)[i] + ((ft)(((ty*)b)[i] - ((ty*)a)[i]) * t) ; \
} \
\
Vector##sufft vInv##suf(const Vector##suf v) { \
	Vector##sufft out; \
	vInv##suf##p(&v, &out); \
	return out; \
} \
void vInv##suf##p(const Vector##suf* v, Vector##sufft* out) { \
	for(int i = 0; i < sz; i++) \
		((ft*)out)[i] = (((ty*)v)[i] == 0) ? pref##_MAX : ((ft)1.0 /  (ft)((ty*)v)[i]); \
} \
\
ft vLen##suf(const Vector##suf v) { \
	return vLen##suf##p(&v); \
} \
ft vLen##suf##p(const Vector##suf* v) { \
	ft tmp = 0.0; \
	for(int i = 0; i < sz; i++) \
		tmp += (ft)((ty*)v)[i] * (ft)((ty*)v)[i]; \
	return sqrt(tmp); \
} \
\
ft vLenSq##suf(const Vector##suf v) { \
	return vLenSq##suf##p(&v); \
} \
ft vLenSq##suf##p(const Vector##suf* v) { \
	return vDot##suf##p(v, v); \
} \
\
ft vMag##suf(const Vector##suf v) { \
	return vLen##suf##p(&v); \
} \
ft vMag##suf##p(const Vector##suf* v) { \
	return vLen##suf##p(v); \
} \
\
ft vInvLen##suf(const Vector##suf v) { \
	ft tmp = vLen##suf(v); \
	return tmp == 0 ? pref##_MAX : (ft)1.0 / tmp; \
} \
ft vInvLen##suf##p(const Vector##suf* v) { \
	return vInvLen##suf(*v); \
} \
\
Vector##sufft vNorm##suf(const Vector##suf v) { \
	Vector##sufft out; \
	vNorm##suf##p(&v, &out); \
	return out; \
} \
void vNorm##suf##p(const Vector##suf* v, Vector##sufft* out) { \
	ft n = vLenSq##suf(*v); \
	\
	if(n >= (ft)1.0f - pref##_CMP_EPSILON && n <= (ft)1.0 + pref##_CMP_EPSILON) { \
		for(int i = 0; i < sz; i++) \
			((ft*)out)[i] = (ft)((ty*)v)[i]; \
		return;		 \
	} \
	else if(n == 0.0) { \
		for(int i = 0; i < sz; i++) \
			((ft*)out)[i] = 0; \
		return; \
	} \
	 \
	n = (ft)1.0 / sqrt(n); \
	for(int i = 0; i < sz; i++) \
		((ft*)out)[i] = (ft)((ty*)v)[i] * n; \
}  \
\
Vector##sufft vUnit##suf(const Vector##suf v) { \
	return vNorm##suf(v); \
} \
void vUnit##suf##p(const Vector##suf* v, Vector##sufft* out) { \
	return vNorm##suf##p(v, out); \
} \


	C3DLAS_VECTOR_LIST(FN)
#undef FN


 


#define FN(sz, suf, t, maxval) \
void vMin##sz##suf##p(const Vector##sz##suf* a, const Vector##sz##suf* b, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = fmin(((t*)a)[i], ((t*)b)[i]); \
} \
void vMax##sz##suf##p(const Vector##sz##suf* a, const Vector##sz##suf* b, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = fmax(((t*)a)[i], ((t*)b)[i]); \
} \
Vector##sz##suf vMin##sz##suf(Vector##sz##suf a, Vector##sz##suf b) { \
	Vector##sz##suf out; \
	vMin##sz##suf##p(&a, &b, &out); \
	return out; \
} \
Vector##sz##suf vMax##sz##suf(Vector##sz##suf a, Vector##sz##suf b) { \
	Vector##sz##suf out; \
	vMax##sz##suf##p(&a, &b, &out); \
	return out; \
} \
\
int vMinComp##sz##suf##p(const Vector##sz##suf* a) { \
	t best = ((t*)a)[0]; \
	int best_ind = 0; \
	for(int i = 1; i < sz; i++) { \
		if(((t*)a)[i] < best) { \
			best = ((t*)a)[i]; \
			best_ind = i; \
		} \
	} \
	return best_ind; \
} \
\
int vMaxComp##sz##suf##p(const Vector##sz##suf* a) { \
	t best = ((t*)a)[0]; \
	int best_ind = 0; \
	for(int i = 1; i < sz; i++) { \
		if(((t*)a)[i] > best) { \
			best = ((t*)a)[i]; \
			best_ind = i; \
		} \
	} \
	return best_ind; \
} \
\
int vMinComp##sz##suf(const Vector##sz##suf a) { \
	return vMinComp##sz##suf##p(&a); \
} \
\
int vMaxComp##sz##suf(const Vector##sz##suf a) { \
	return vMaxComp##sz##suf##p(&a); \
} \
\
Vector##sz##suf vClamp##sz##suf(Vector##sz##suf in, Vector##sz##suf min, Vector##sz##suf max) { \
	Vector##sz##suf out; \
	for(int i = 0; i < sz; i++) \
		((t*)&out)[i] = fmax(((t*)&min)[i], fmin(((t*)&in)[i], ((t*)&max)[i])); \
	return out; \
} \
Vector##sz##suf vAbs##sz##suf(const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vAbs##sz##suf##p(&v, &out); \
	return out; \
} \
void vAbs##sz##suf##p(const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = abs_##t( ((t*)v)[i] ); \
} \
Vector##sz##suf vSign##sz##suf(const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vSign##sz##suf##p(&v, &out); \
	return out; \
} \
void vSign##sz##suf##p(const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = copysign((t)1.0, ((t*)v)[i] ); \
} \
Vector##sz##suf vStep##sz##suf(const Vector##sz##suf edge, const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vStep##sz##suf##p(&edge, &v, &out); \
	return out; \
} \
void vStep##sz##suf##p(const Vector##sz##suf* edge, const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = ((t*)v)[i] < ((t*)edge)[i] ? 0.0 : 1.0; \
} \


FN(2,  , float,  FLT_MAX)
FN(3,  , float,  FLT_MAX)
FN(4,  , float,  FLT_MAX)
FN(2, d, double, DBL_MAX)
FN(3, d, double, DBL_MAX)
FN(4, d, double, DBL_MAX)
#undef FN

#define FN(sz, suf, t, maxval) \
void vMin##sz##suf##p(const Vector##sz##suf* a, const Vector##sz##suf* b, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = MIN(((t*)a)[i], ((t*)b)[i]); \
} \
void vMax##sz##suf##p(const Vector##sz##suf* a, const Vector##sz##suf* b, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = MAX(((t*)a)[i], ((t*)b)[i]); \
} \
Vector##sz##suf vMin##sz##suf(Vector##sz##suf a, Vector##sz##suf b) { \
	Vector##sz##suf out; \
	vMin##sz##suf##p(&a, &b, &out); \
	return out; \
} \
Vector##sz##suf vMax##sz##suf(Vector##sz##suf a, Vector##sz##suf b) { \
	Vector##sz##suf out; \
	vMax##sz##suf##p(&a, &b, &out); \
	return out; \
} \
Vector##sz##suf vClamp##sz##suf(Vector##sz##suf in, Vector##sz##suf min, Vector##sz##suf max) { \
	Vector##sz##suf out; \
	for(int i = 0; i < sz; i++) \
		((t*)&out)[i] = MAX(((t*)&min)[i], MIN(((t*)&in)[i], ((t*)&max)[i])); \
	return out; \
} \
Vector##sz##suf vAbs##sz##suf(const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vAbs##sz##suf##p(&v, &out); \
	return out; \
} \
void vAbs##sz##suf##p(const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = labs( ((t*)v)[i] ); \
} \
Vector##sz##suf vSign##sz##suf(const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vSign##sz##suf##p(&v, &out); \
	return out; \
} \
void vSign##sz##suf##p(const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = ((t*)v)[i] < 0 ? -1 : 1; \
} \
Vector##sz##suf vStep##sz##suf(const Vector##sz##suf edge, const Vector##sz##suf v) { \
	Vector##sz##suf out; \
	vStep##sz##suf##p(&edge, &v, &out); \
	return out; \
} \
void vStep##sz##suf##p(const Vector##sz##suf* edge, const Vector##sz##suf* v, Vector##sz##suf* out) { \
	for(int i = 0; i < sz; i++) \
		((t*)out)[i] = ((t*)v)[i] < ((t*)edge)[i] ? 0.0 : 1.0; \
} \


FN(2, i, int,  DBL_MAX)
FN(3, i, int,  DBL_MAX)
FN(4, i, int,  DBL_MAX)
FN(2, l, long, DBL_MAX)
FN(3, l, long, DBL_MAX)
FN(4, l, long, DBL_MAX)
#undef FN




