






const Matrix3 IDENT_MATRIX3 = { { 1, 0, 0,
                                 0, 1, 0,
                                 0, 0, 1,
                                  } };







void mIdent3(Matrix3* m) {
	*m = IDENT_MATRIX3;
}

// out cannot overlap with a or b
// with restrict and -O2, this vectorizes nicely.
void mFastMul3(Matrix3* restrict a, Matrix3* restrict b, Matrix3* restrict out) {
	int r, c;
	
	for(r = 0; r < 3; r++) {
		for(c = 0; c < 3; c++) {
			out->m[c + r * 3] =
				(a->m[r * 3 + 0] * b->m[c + 0]) +
				(a->m[r * 3 + 1] * b->m[c + 3]) +
				(a->m[r * 3 + 2] * b->m[c + 6]);
		}
	}
}


void mMul3(Matrix3* a, Matrix3* out) {
	Matrix3 b = *out;
	mFastMul3(a, &b, out);
}


Vector3 vMatrix3Mul(Vector3 v, Matrix3* restrict m) {
	return (Vector3){
		.x = v.x * m->m[0+0] + v.y * m->m[3+0] + v.z * m->m[6+0],
		.y = v.x * m->m[0+1] + v.y * m->m[3+1] + v.z * m->m[6+1],
		.z = v.x * m->m[0+2] + v.y * m->m[3+2] + v.z * m->m[6+2]
	};
}


float mDeterminate3(Matrix3* m) {
// computes the inverse of a matrix m
	return m->m[0+0] * (m->m[3+1] * m->m[6+2] - m->m[6+1] * m->m[3+2]) -
	       m->m[0+1] * (m->m[3+0] * m->m[6+2] - m->m[3+2] * m->m[6+0]) +
	       m->m[0+2] * (m->m[3+0] * m->m[6+1] - m->m[3+1] * m->m[6+0]);
}


void mInverse3(Matrix3* m, Matrix3* out) {
	double id = 1.0 / mDeterminate3(m);

	out->m[0+0] = (m->m[3+1] * m->m[6+2] - m->m[6+1] * m->m[3+2]) * id;
	out->m[0+1] = (m->m[0+2] * m->m[6+1] - m->m[0+1] * m->m[6+2]) * id;
	out->m[0+2] = (m->m[0+1] * m->m[3+2] - m->m[0+2] * m->m[3+1]) * id;
	out->m[3+0] = (m->m[3+2] * m->m[6+0] - m->m[3+0] * m->m[6+2]) * id;
	out->m[3+1] = (m->m[0+0] * m->m[6+2] - m->m[0+2] * m->m[6+0]) * id;
	out->m[3+2] = (m->m[3+0] * m->m[0+2] - m->m[0+0] * m->m[3+2]) * id;
	out->m[6+0] = (m->m[3+0] * m->m[6+1] - m->m[6+0] * m->m[3+1]) * id;
	out->m[6+1] = (m->m[6+0] * m->m[0+1] - m->m[0+0] * m->m[6+1]) * id;
	out->m[6+2] = (m->m[0+0] * m->m[3+1] - m->m[3+0] * m->m[0+1]) * id;
}



void mScalarMul3(Matrix3* m, float scalar, Matrix3* out) {
	for(int i = 0; i < 9; i++) out->m[i] = m->m[i] * scalar;
}





