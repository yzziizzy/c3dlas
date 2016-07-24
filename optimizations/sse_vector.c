


#include <stdio.h>

#include "../c3dlas.h"



#ifndef __SSE__
	#define c3dlas_sse3_vAdd c3dlas_pure_vAdd
#else
void c3dlas_sse3_vAdd(Vector* a, Vector* b, Vector* out) {
	/*
	out->x = a->x + b->x;
	out->y = a->y + b->y;
	out->z = a->z + b->z;
	*/
	asm (
		// pointer dereference
		"mov (%0), %%rax;"
		"mov (%1), %%rbx;"
		"mov (%2), %%rcx;"
		
		// move 3 floats into the xmm regs
		"movq (%%rax), %%xmm1;"
		"movlhps %%xmm1,%%xmm1;"
		"movd 0x8(%%rax),%%xmm1;"
		
		"movq (%%rbx), %%xmm2;"
		"movlhps %%xmm2,%%xmm2;"
		"movd 0x8(%%rbx),%%xmm2;"
		
		// THE MATH
		"addps %%xmm2, %%xmm1;"
		
		// move 3 floats to output
		"movd %%xmm1, 0x8(%%rcx);"
		"movhlps %%xmm1,%%xmm1;"
		"movq %%xmm1, (%%rax);"
		: // no outputs
		: "r" (a), "r" (b), "r" (out)
		: "rax", "rbx", "rcx", "memory"
		);
}
#endif


#ifndef __SSE__
	#define c3dlas_sse3_vAdd4 c3dlas_pure_vAdd
#else
void c3dlas_sse3_vAdd4(Vector* a, Vector* b, Vector* out) {
	/*
	out->x = a->x + b->x;
	out->y = a->y + b->y;
	out->z = a->z + b->z;
	out->w = a->w + b->w;
	*/
	asm (
		"mov (%0), %%rax;"
		"mov (%1), %%rbx;"
		"mov (%2), %%rcx;"
		"movups (%%rax), %%xmm1;"
		"movups (%%rbx), %%xmm2;"
		"addps %%xmm2, %%xmm1;"
		"movups %%xmm1, (%%rcx);"
		: // no outputs
		: "r" (a), "r" (b), "r" (out)
		: "rax", "rbx", "rcx", "memory"
		);
}




