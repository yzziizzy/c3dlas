
#include "inertia.h"




void inertiaSphere(Matrix3* mat, float mass, float radius) {
	float s = (2.0 / 5.0) * mass * radius * radius;
	*mat = (Matrix3){0};
	mat->m[0] = s;
	mat->m[4] = s;
	mat->m[8] = s;
}


// axis aligned box
void inertiaBox(Matrix3* mat, float mass, vec3 dims) {
	float otm = mass / 12.f;
	float x2 = dims.x * dims.x;
	float y2 = dims.y * dims.y;
	float z2 = dims.z * dims.z;
	
	*mat = (Matrix3){0};
	mat->m[0] = otm * (y2 + z2);
	mat->m[4] = otm * (x2 + z2);
	mat->m[8] = otm * (y2 + x2);
}



