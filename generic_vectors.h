#ifndef __c3dlas__generic_vectors_h__
#define __c3dlas__generic_vectors_h__




#define vAdd(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vAdd) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vSub(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vSub) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMul(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMul) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vDiv(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vDiv) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vScale(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vScale) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vFMA(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vFMA) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vLerp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vLerp) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vRecip(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vRecip) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vInv(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vInv) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMag(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMag) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vInvLen(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vInvLen) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vLen(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vLen) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vLenSq(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vLenSq) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vDot(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vDot) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vDist(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vDist) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vDistSq(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vDistSq) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vNorm(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vNorm) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vNeg(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vNeg) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vSign(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vSign) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vStep(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vStep) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMin(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMin) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMax(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMax) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMinComp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMinComp) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vMaxComp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vMaxComp) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vAbs(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vAbs) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vAvg(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vAvg) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vClamp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vClamp) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vEq(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vEq) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vEqEp(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vEqEp) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)
#define vEqExact(a, ...) _Generic((a), C3DLAS_VECTOR_LIST(C3DLAS_GEN_HELPER, vEqExact) default: ((void)0))(a __VA_OPT__(,) __VA_ARGS__)

#define vCross(a, ...) _Generic((a), \
	Vector2:  vCross2, \
	Vector2*: vCross2p, \
	Vector2d: vCross2d, \
	Vector2d*: vCross2dp, \
	Vector3:  vCross3, \
	Vector3*: vCross3p, \
	Vector3d: vCross3d, \
	Vector3d*: vCross3dp, \
	default: ((void)0)) \
(a __VA_OPT__(,) __VA_ARGS__) \



#define vPerp(a, ...) _Generic((a), \
	Vector2:  vPerp2, \
	Vector3:  vPerp3, \
	default: ((void)0)) \
(a __VA_OPT__(,) __VA_ARGS__) \







#endif // __c3dlas__generic_vectors_h__
