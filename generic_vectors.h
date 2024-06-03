#ifdef __c3dlas__generic_vectors_h__
#define __c3dlas__generic_vectors_h__




#define vAdd(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vAdd) default: ((void)0))(a, __VA_ARGS__)
#define vSub(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vSub) default: ((void)0))(a, __VA_ARGS__)
#define vMul(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMul) default: ((void)0))(a, __VA_ARGS__)
#define vScale(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vScale) default: ((void)0))(a, __VA_ARGS__)
#define vDot(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vDot) default: ((void)0))(a, __VA_ARGS__)
#define vCross(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vCross) default: ((void)0))(a, __VA_ARGS__)
#define vScalarTriple(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vScalarTriple) default: ((void)0))(a, __VA_ARGS__)
#define vLerp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vLerp) default: ((void)0))(a, __VA_ARGS__)
#define vInverse(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vInverse) default: ((void)0))(a, __VA_ARGS__)
#define vNorm(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vNorm) default: ((void)0))(a, __VA_ARGS__)










#endif // __c3dlas__generic_vectors_h__
