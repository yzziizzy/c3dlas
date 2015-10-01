



typedef struct {
	float x,y,z;
} Vector;




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






