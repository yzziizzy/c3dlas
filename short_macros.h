#ifndef __c3dlas__short_macros_h__
#define __c3dlas__short_macros_h__


#ifndef PP_NARG
	#define PP_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37, _38, _39, _40, _41, _42, _43, _44, _45, _46, _47, _48, _49, _50, _51, _52, _53, _54, _55, _56, _57, _58, _59, _60, _61, _62, _63, N, ...) N
	#define PP_RSEQ_N() 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
	#define PP_NARG_(...) PP_ARG_N(__VA_ARGS__)
	#define PP_NARG(...)  PP_NARG_(__VA_ARGS__, PP_RSEQ_N())
#endif

typedef uint8_t bool;

typedef int8_t s8;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef float f32;
typedef double f64;

typedef Vector4 quat;
typedef Matrix Matrix4;
typedef Matrix4 mat4;
typedef Matrix3 mat3;

typedef Vector2 vec2;
typedef Vector3 vec3;
typedef Vector4 vec4;
typedef Vector2d vec2d;
typedef Vector3d vec3d;
typedef Vector4d vec4d;
typedef Vector2i vec2i;
typedef Vector3i vec3i;
typedef Vector4i vec4i;
typedef Vector2l vec2l;
typedef Vector3l vec3l;
typedef Vector4l vec4l;


// __attribute__((always_inline))  forces gcc to inline these trivial functions even during -O0
#define _GETTER_FN(type1, type2, contentsx, contentsy, contentsz, contentsw) \
	static inline type2 __attribute__((always_inline)) _##type1##_getx(type1 a, type2 c) { return contentsx; } \
	static inline type2 __attribute__((always_inline)) _##type1##_gety(type1 a, type2 c) { return contentsy; } \
	static inline type2 __attribute__((always_inline)) _##type1##_getz(type1 a, type2 c) { return contentsz; } \
	static inline type2 __attribute__((always_inline)) _##type1##_getw(type1 a, type2 c) { return contentsw; } \

#define _VECTOR_GETTERS(X, ...) \
	X(__VA_ARGS__ __VA_OPT__(,) int,    int,    a,   a,   a,   a  ) \
	X(__VA_ARGS__ __VA_OPT__(,) float,  float,  a,   a,   a,   a  ) \
	X(__VA_ARGS__ __VA_OPT__(,) long,   long,   a,   a,   a,   a  ) \
	X(__VA_ARGS__ __VA_OPT__(,) double, double, a,   a,   a,   a  ) \
	X(__VA_ARGS__ __VA_OPT__(,) vec2,   float,  a.x, a.y, c,   c  ) \
	X(__VA_ARGS__ __VA_OPT__(,) vec3,   float,  a.x, a.y, a.z, c  ) \
	X(__VA_ARGS__ __VA_OPT__(,) vec4,   float,  a.x, a.y, a.z, a.w) \
	X(__VA_ARGS__ __VA_OPT__(,) vec2d,  double, a.x, a.y, c,   c  ) \
	X(__VA_ARGS__ __VA_OPT__(,) vec3d,  double, a.x, a.y, a.z, c  ) \
	X(__VA_ARGS__ __VA_OPT__(,) vec4d,  double, a.x, a.y, a.z, a.w) \
	X(__VA_ARGS__ __VA_OPT__(,) vec2l,  s64,    a.x, a.y, c,   c  ) \
	X(__VA_ARGS__ __VA_OPT__(,) vec3l,  s64,    a.x, a.y, a.z, c  ) \
	X(__VA_ARGS__ __VA_OPT__(,) vec4l,  s64,    a.x, a.y, a.z, a.w) \
	X(__VA_ARGS__ __VA_OPT__(,) vec2i,  s32,    a.x, a.y, c,   c  ) \
	X(__VA_ARGS__ __VA_OPT__(,) vec3i,  s32,    a.x, a.y, a.z, c  ) \
	X(__VA_ARGS__ __VA_OPT__(,) vec4i,  s32,    a.x, a.y, a.z, a.w) 


_VECTOR_GETTERS(_GETTER_FN)
#ifndef CAT
	#define CAT(a,b) a##b
#endif

#define CAT4(a,b,c,d) a##b##c##d

#define V__N(a, narg, ...) CAT(a, narg)(__VA_ARGS__)

#define _V_HELPER(a, b, ...) b: CAT4(_, b, _get, a), 

#define _V_GET(a, b, c) _Generic(a, \
	_VECTOR_GETTERS(_V_HELPER, b) \
	default: _float_getx \
)(a, c)


#define V2(...) V__N(V2_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V2_1(a)  ((Vector2){.x=_V_GET(a, x, 0),.y=_V_GET(a, y, 0)})
#define V2_2(_x,_y) ((Vector2){.x=(_x),.y=(_y)})

#define V3(...) V__N(V3_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V3_1(a)  ((Vector3){.x=_V_GET(a, x, 0),.y=_V_GET(a, y, 0),.z=_V_GET(a, z, 0)})
#define V3_2(a, _z)  ((Vector3){.x=_V_GET(a, x, 0),.y=_V_GET(a, y, 0),.z=_z})
#define V3_3(_x,_y,_z) ((Vector3){.x=(_x),.y=(_y),.z=(_z)})

#define V4(...) V__N(V4_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V4_1(a) ((Vector4){.x=_V_GET(a, x, 0),.y=_V_GET(a, y, 0),.z=_V_GET(a, z, 0),.w=_V_GET(a, w, 1)})
#define V4_2(a, b) ((Vector4){.x=_V_GET(a, x, 0),.y=_V_GET(a, y, 0),.z=_V_GET(a, z, 0),b})
#define V4_3(_x,_y,_z) ((Vector4){.x=(_x),.y=(_y),.z=(_z),.w=1})
#define V4_4(_x,_y,_z,_w) ((Vector4){.x=(_x),.y=(_y),.z=(_z),.w=(_w)})



#define _Vi_GET(a, b, c) _Generic(a, \
	_VECTOR_GETTERS(_V_HELPER, b) \
	default: _int_getx \
)(a, c)

#define V2i(...) V__N(V2i_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V2i_1(a)  ((Vector2i){.x=_Vi_GET(a, x, 0),.y=_Vi_GET(a, y, 0)})
#define V2i_2(_x,_y) ((Vector2i){.x=(_x),.y=(_y)})

#define V3i(...) V__N(V3i_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V3i_1(a)  ((Vector3i){.x=_Vi_GET(a, x, 0),.y=_Vi_GET(a, y, 0),.z=_Vi_GET(a, z, 0)})
#define V3i_2(a, _z)  ((Vector3i){.x=_Vi_GET(a, x, 0),.y=_Vi_GET(a, y, 0),.z=_z})
#define V3i_3(_x,_y,_z) ((Vector3i){.x=(_x),.y=(_y),.z=(_z)})

#define V4i(...) V__N(V4i_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V4i_1(a) ((Vector4i){.x=_Vi_GET(a, x, 0),.y=_Vi_GET(a, y, 0),.z=_Vi_GET(a, z, 0),.w=_Vi_GET(a, w, 1)})
#define V4i_2(a, b) ((Vector4i){.x=_Vi_GET(a, x, 0),.y=_Vi_GET(a, y, 0),.z=_Vi_GET(a, z, 0),b})
#define V4i_3(_x,_y,_z) ((Vector4i){.x=(_x),.y=(_y),.z=(_z),.w=1})
#define V4i_4(_x,_y,_z,_w) ((Vector4i){.x=(_x),.y=(_y),.z=(_z),.w=(_w)})



#define _Vl_GET(a, b, c) _Generic(a, \
	_VECTOR_GETTERS(_V_HELPER, b) \
	default: _long_getx \
)(a, c)

#define V2l(...) V__N(V2l_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V2l_1(a)  ((Vector2l){.x=_Vl_GET(a, x, 0),.y=_Vl_GET(a, y, 0)})
#define V2l_2(_x,_y) ((Vector2l){.x=(_x),.y=(_y)})

#define V3l(...) V__N(V3l_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V3l_1(a)  ((Vector3l){.x=_Vl_GET(a, x, 0),.y=_Vl_GET(a, y, 0),.z=_Vl_GET(a, z, 0)})
#define V3l_2(a, _z)  ((Vector3l){.x=_Vl_GET(a, x, 0),.y=_Vl_GET(a, y, 0),.z=_z})c,
#define V3l_3(_x,_y,_z) ((Vector3l){.x=(_x),.y=(_y),.z=(_z)})

#define V4l(...) V__N(V4l_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V4l_1(a) ((Vector4l){.x=_Vl_GET(a, x, 0),.y=_Vl_GET(a, y, 0),.z=_Vl_GET(a, z, 0),.w=_Vl_GET(a, w, 1)})
#define V4l_2(a, b) ((Vector4l){.x=_Vl_GET(a, x, 0),.y=_Vl_GET(a, y, 0),.z=_Vl_GET(a, z, 0),b})
#define V4l_3(_x,_y,_z) ((Vector4l){.x=(_x),.y=(_y),.z=(_z),.w=1})
#define V4l_4(_x,_y,_z,_w) ((Vector4l){.x=(_x),.y=(_y),.z=(_z),.w=(_w)})



#define _Vd_GET(a, b, c) _Generic(a, \
	_VECTOR_GETTERS(_V_HELPER, b) \
	default: _double_getx \
)(a, c)

#define V2d(...) V__N(V2d_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V2d_1(a)  ((Vector2d){.x=_Vd_GET(a, x, 0),.y=_Vd_GET(a, y, 0)})
#define V2d_2(_x,_y) ((Vector2d){.x=(_x),.y=(_y)})

#define V3d(...) V__N(V3d_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V3d_1(a)  ((Vector3d){.x=_Vd_GET(a, x, 0),.y=_Vd_GET(a, y, 0),.z=_Vd_GET(a, z, 0)})
#define V3d_2(a, _z)  ((Vector3d){.x=_Vd_GET(a, x, 0),.y=_Vd_GET(a, y, 0),.z=_z})c,
#define V3d_3(_x,_y,_z) ((Vector3d){.x=(_x),.y=(_y),.z=(_z)})

#define V4d(...) V__N(V4d_, PP_NARG(__VA_ARGS__), __VA_ARGS__)
#define V4d_1(a) ((Vector4d){.x=_Vd_GET(a, x, 0),.y=_Vd_GET(a, y, 0),.z=_Vd_GET(a, z, 0),.w=_Vd_GET(a, w, 1)})
#define V4d_2(a, b) ((Vector4d){.x=_Vd_GET(a, x, 0),.y=_Vd_GET(a, y, 0),.z=_Vd_GET(a, z, 0),b})
#define V4d_3(_x,_y,_z) ((Vector4d){.x=(_x),.y=(_y),.z=(_z),.w=1})
#define V4d_4(_x,_y,_z,_w) ((Vector4d){.x=(_x),.y=(_y),.z=(_z),.w=(_w)})



#endif // __c3dlas__short_macros_h__

