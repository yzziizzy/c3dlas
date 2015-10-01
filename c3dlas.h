

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
	float x,y,z;
} Vector;

typedef struct {
	union {
		struct { float x,y,z,w; };
		float v[4];
	};
} Vector4;




typedef struct {
	float m[16];
} Matrix;



extern const Matrix IDENT_MATRIX;





void  vAdd(Vector* a, Vector* b, Vector* out); // add two vectors
void  vSub(Vector* from, Vector* what, Vector* diff); // diff = from - what
void  vScale(Vector* v, float scalar, Vector* out); // scalar muliplication 
float vMag(Vector* v); // return the magnitude
float vDot(Vector* a, Vector* b); // dot product
void  vNorm(Vector* v, Vector* out); // normalize the vector
void  vUnit(Vector* v, Vector* out); // normalise the vector, alternate name
void  vCross(Vector* a, Vector* b, Vector* out); // cross product: out = a x b
float vScalarTriple(Vector* a, Vector* b, Vector* c); // scalar triple product: a . (b x c)

void vMatrixMul(Vector* in, Matrix* m, Vector* out); // multiply a vector by a matrix
void vMatrixMulf(float x, float y, float z, Matrix* m, Vector* out); // multiply a vector by a matrix


void mIdent(Matrix* m); // set m to the identity matrix
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






