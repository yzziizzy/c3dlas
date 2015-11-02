
#ifndef __c3dlas_h__
#define __c3dlas_h__


#define F_PI ((float)3.1415926535897932384626433832795028841971693993751)
#define D_PI ((double)3.1415926535897932384626433832795028841971693993751)
#define F_2PI ((float)6.2831853071795864769252867665590057683943387987502)
#define D_2PI ((double)6.2831853071795864769252867665590057683943387987502)
#define F_1_PI ((float)0.3183098861837906715377675267450287240689192914809)
#define D_1_PI ((double)0.3183098861837906715377675267450287240689192914809)
#define F_PI_2 ((float)1.5707963267948966192313216916397514420985846996875)
#define D_PI_2 ((double)1.5707963267948966192313216916397514420985846996875)
#define F_3PI_2 ((float)4.7123889803846898576939650749192543262957540990626)
#define D_3PI_2 ((double)4.7123889803846898576939650749192543262957540990626)

#define RAD2DEG (57.29577951308232087679815481410517033240547246656432154916024386)
#define DEG2RAD (0.0174532925199432957692369076848861271344287188854172545609719144)



typedef struct {
	float x,y;
} Vector2;

typedef struct {
	float x,y,z;
} Vector;

typedef struct {
	float x,y,z,w;
} Vector4;


typedef struct {
	Vector o; // origin
	Vector d; // normalized direction
	Vector id; // inverse normalized direction (handy enough to keep around)
} Ray;


typedef struct {
	Vector n; // normal
	float d; // distance along normal to the origin
} Plane;

typedef struct { // does not have to be coplanar
	Vector v[4];
} Quad;

typedef struct {
	float m[16];
} Matrix;


typedef struct MatrixStack {
	short size;
	short top;
	Matrix* stack;
} MatrixStack;


// axis-aligned bounding box
typedef struct AABB {
	Vector min;
	Vector max;
} AABB;

typedef struct AABB2 {
	Vector min;
	Vector max;
} AABB2;


extern const Matrix IDENT_MATRIX;





void  vAdd(Vector* a, Vector* b, Vector* out); // add two vectors
void  vSub(Vector* from, Vector* what, Vector* diff); // diff = from - what
void  vScale(Vector* v, float scalar, Vector* out); // scalar muliplication 
void  vInverse(Vector* v, Vector* out); // inverse
float vMag(Vector* v); // return the magnitude
float vDot(Vector* a, Vector* b); // dot product
void  vNorm(Vector* v, Vector* out); // normalize the vector
void  vUnit(Vector* v, Vector* out); // normalise the vector, alternate name
void  vCross(Vector* a, Vector* b, Vector* out); // cross product: out = a x b
float vScalarTriple(Vector* a, Vector* b, Vector* c); // scalar triple product: a . (b x c)
void  vProject(Vector* what, Vector* onto, Vector* out); // slower; onto may not be normalized
void  vProjectNorm(Vector* what, Vector* onto, Vector* out); // faster; onto must be normalized

float pvDist(Plane* p, Vector* v);

void vMatrixMul(Vector* in, Matrix* m, Vector* out); // multiply a vector by a matrix
void vMatrixMulf(float x, float y, float z, Matrix* m, Vector* out); // multiply a vector by a matrix


void mIdent(Matrix* m); // set m to the identity matrix
void mCopy(Matrix* in, Matrix* out); 
void mFastMul(Matrix* a, Matrix* b, Matrix* out); // a and b cannot also be out. mostly internal use.
void mMul(Matrix* a, Matrix* out); // makes a copy of out before multiplying over it
void mTransv(Vector* v, Matrix* out); // translation
void mTrans3f(float x, float y, float z, Matrix* out); // translation
void mScalev(Vector* v, Matrix* out);
void mScale3f(float x, float y, float z, Matrix* out);
void mRotv(Vector* v, float theta, Matrix* out); // rotate about a vector
void mRot3f(float x, float y, float z, float theta, Matrix* out); // rotate about a vector
void mRotX(float theta, Matrix* out); // 
void mRotY(float theta, Matrix* out); // rotate about axes
void mRotZ(float theta, Matrix* out); //
void mTranspose(Matrix* in, Matrix* out);
void mTransposeFast(Matrix* in, Matrix* out); // in cannot be out
float mDeterminate(Matrix* m);
int mInverse(Matrix* in, Matrix* out); // returns 0 on success, 1 if there is no inverse; out remains unchanged


// analogous to glFrustum
// no div/0 checking here for right == left etc. just don't be an idiot.
void mFrustum(float left, float right, float top, float bottom, float near, float far, Matrix* out);

// analogous to gluPerspective
// same div/0 warnings apply. if you get an FP exception you deserve it. 
// use a double for fov; the precision matters often.
// https://www.opengl.org/archives/resources/faq/technical/transformations.htm
// https://www.opengl.org/sdk/docs/man2/xhtml/gluPerspective.xml
void mPerspective(double fov, float aspect, float near, float far, Matrix* out);

// orthographic projection. use this for a "2D" look.
// same div/0 warnings. 
void mOrtho(float left, float right, float top, float bottom, float near, float far, Matrix* out);

// analgous to gluLookAt
// https://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
void mLookAt(Vector* eye, Vector* center, Vector* up, Matrix* out);

void mPrint(Matrix* m, FILE* f);



// matrix stack functions

// make sure you allocate enough. when it's out, it's out. no surprise mallocs later on. (yet)
void msAlloc(int size, MatrixStack* ms);
void msFree(MatrixStack* ms);

int msPush(MatrixStack* ms); 
void msPop(MatrixStack* ms);
Matrix* msGetTop(MatrixStack* ms);

void msPrintAll(MatrixStack* ms, FILE* f);

// these are all wrappers around the functions listed above
void msIdent(MatrixStack* ms); // set to the identity matrix
void msCopy(Matrix* in, MatrixStack* ms); 
void msMul(Matrix* a, MatrixStack* ms); // makes a copy of out before multiplying over it
void msTransv(Vector* v, MatrixStack* ms); // translation
void msTrans3f(float x, float y, float z, MatrixStack* ms); // translation
void msScalev(Vector* v, MatrixStack* ms);
void msScale3f(float x, float y, float z, MatrixStack* ms);
void msRotv(Vector* v, float theta, MatrixStack* ms); // rotate about a vector
void msRot3f(float x, float y, float z, float theta, MatrixStack* ms); // rotate about a vector
void msFrustum(float left, float right, float top, float bottom, float near, float far, MatrixStack* ms);
void msPerspective(double fov, float aspect, float near, float far, MatrixStack* ms);
void msOrtho(float left, float right, float top, float bottom, float near, float far, MatrixStack* ms);
void msLookAt(Vector* eye, Vector* center, Vector* up, MatrixStack* ms);


void evalBezier(Vector* e1, Vector* e2, Vector* c1, Vector* c2, float t, Vector* out);
void evalBezierTangent(Vector* e1, Vector* e2, Vector* c1, Vector* c2, float t, Vector* out); // tangent vector; not normalized 
void evalBezierNorm(Vector* e1, Vector* e2, Vector* c1, Vector* c2, float t, Vector* out); // normal vector; not normalized
float evalBezier1D(float e1, float e2, float c1, float c2, float t);
float evalBezier1D_dt(float e1, float e2, float c1, float c2, float t); // first derivative with respect to t
float evalBezier1D_ddt(float e1, float e2, float c1, float c2, float t); // second derivative with respect to t

///// bounding box functions

// 3D versions
int boxDisjoint(const AABB* a, const AABB* b);
int boxOverlaps(const AABB* a, const AABB* b);
int boxContainsPoint(const AABB* b, const Vector* p);

void boxCenter(const AABB* b, Vector* out); // calculates the center of the box
void boxSize(const AABB* b, Vector* out); // calculates the size of the box


void makeRay(Vector* origin, Vector* direction, Ray* out);
int boxRayIntersectFast(const AABB* b, const Ray* r);
int boxRayIntersect(const AABB* b, const Ray* r, Vector* ipoint, float* idist);


// 2D versions
int boxDisjoint2(const AABB2* a, const AABB2* b);
int boxOverlaps2(const AABB2* a, const AABB2* b);
int boxContainsPoint2(const AABB2* b, const Vector2* p);



void boxCenter2(const AABB2* b, Vector2* out); // calcuates the center of the box
void boxSize2(const AABB2* b, Vector2* out); // calculates the size of the box
void boxQuadrant2(const AABB2* in, char ix, char iy, AABB2* out);



#endif // __c3dlas_h__

