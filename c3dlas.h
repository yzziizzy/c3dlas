



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




void mIdent(Matrix* m);
void mFastMul(Matrix* a, Matrix* b, Matrix* out); // a and b cannot also be out. mostly internal use.
void mMul(Matrix* a, Matrix* out); // makes a copy of out before multiplying over it
void mTransv(Vector* v, Matrix* out);
void mTrans3f(float x, float y, float z, Matrix* out);
void mScalev(Vector* v, Matrix* out);
void mScale3f(float x, float y, float z, Matrix* out);
void mRotv(Vector* v, float theta, Matrix* out);
void mRot3f(float x, float y, float z, float theta, Matrix* out);
void mRotX(float theta, Matrix* out);
void mRotY(float theta, Matrix* out);
void mRotZ(float theta, Matrix* out);









