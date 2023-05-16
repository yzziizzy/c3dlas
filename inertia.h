#ifndef __c3dlas__inertia_h__
#define __c3dlas__inertia_h__

#include "c3dlas.h"


void inertiaSphere(Matrix3* mat, float mass, float radius);

// axis aligned box
void inertiaBox(Matrix3* mat, float mass, vec3 dims);


#endif // __c3dlas__inertia_h__
