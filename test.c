#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "c3dlas.h"









void checklinelineline() {

	
	Ray3 r1 = {
		.o = { 0, 0, 3},
		.d = { 0, 2, 0},
	};
	
	Ray3 r2 = {
		.o = { 0, 0, 1},
		.d = { 0, 2, 0},
	};

	
	Vector3 o1[] = {
		{99,99,99},
		{99,99,99},
	};
	
	int ret;
	
	ret = shortestLineFromRayToRay3p(&r1, &r2, o1);
	
	printf("%s\n", c3dlas_EnumString(ret));
	printf("  out: % .1f,% .1f,% .1f - % .1f,% .1f,% .1f\n", o1[0].x, o1[0].y, o1[0].z, o1[1].x, o1[1].y, o1[1].z);
	
}


void checklinePlaneClip() {

	Plane pl = {
		{0,0,1}, 1,
	};
	
	Vector3 line[] = {
		{0,0, -3},
		{0,0, -4},
	};
	
	Vector3 o1[] = {
		{99,99,99},
		{99,99,99},
	};
	Vector3 o2[] = {
		{99,99,99},
		{99,99,99},
	};
	
	int ret, i1, i2;
	
	ret = linePlaneClip3p(line, line+1, &pl, o1, o2, &i1, &i2);
	
	printf("%s\n", c3dlas_EnumString(ret));
	printf("  above1 %d: % .1f,% .1f,% .1f - % .1f,% .1f,% .1f\n", i1, o1[0].x, o1[0].y, o1[0].z, o1[1].x, o1[1].y, o1[1].z);
	printf("  below1 %d: % .1f,% .1f,% .1f - % .1f,% .1f,% .1f\n", i2, o2[0].x, o2[0].y, o2[0].z, o2[1].x, o2[1].y, o2[1].z);
	
}



void checkTriPlaneClip() {
	Plane pl = {
		{0,0,1}, 1,
	};
	
	Vector3 tri1[] = { // intersecting
		{ 0,0, 5},
		{-3,0,-3}, // l,h
		{ 3,0,-3},
	};
	
	Vector3 tri11[] = { // intersecting
		{-3,0, 3}, // l,h
		{ 3,0, 3},
		{ 0,0,-5},
	};
	
	Vector3 tri2[] = { // entirely above
		{3,4,2},
		{3,5,3},
		{3,6,7},
	};
	Vector3 tri3[] = {
		{-2,0, 1},
		{ 3,0, 2},
		{ 3,0,-2},
	};
	
	Vector3 l[] = {
		{3, 3, -4},
		{4, 4, 5},
		{4, 4, -9},
	};
	
	Vector3 o1[] = {
		{999,999,999},{999,999,999},{999,999,999},
		{999,999,999},{999,999,999},{999,999,999},
		{888,888,888},{888,888,888},{888,888,888}
	};
	Vector3 o2[] = {
		{999,999,999},{999,999,999},{999,999,999},
		{999,999,999},{999,999,999},{999,999,999},
		{888,888,888},{888,888,888},{888,888,888}
	};
	
	int ret, i1, i2;
	
	ret = triPlaneClip3p(tri3, &pl, o1, o2, &i1, &i2);
	
	printf("%s\n", c3dlas_EnumString(ret));
	printf("  above1 %d: % .1f,% .1f,% .1f - % .1f,% .1f,% .1f - % .1f,% .1f,% .1f \n", i1, o1[0].x, o1[0].y, o1[0].z, o1[1].x, o1[1].y, o1[1].z, o1[2].x, o1[2].y, o1[2].z);
	printf("  above2 %d: % .1f,% .1f,% .1f - % .1f,% .1f,% .1f - % .1f,% .1f,% .1f - %f \n", i1, o1[3].x, o1[3].y, o1[3].z, o1[4].x, o1[4].y, o1[4].z, o1[5].x, o1[5].y, o1[5].z, o1[6].x );
	printf("  below1 %d: % .1f,% .1f,% .1f - % .1f,% .1f,% .1f - % .1f,% .1f,% .1f \n", i2, o2[0].x, o2[0].y, o2[0].z, o2[1].x, o2[1].y, o2[1].z, o2[2].x, o2[2].y, o2[2].z);
	printf("  below2 %d: % .1f,% .1f,% .1f - % .1f,% .1f,% .1f - % .1f,% .1f,% .1f - %f \n", i2, o2[3].x, o2[3].y, o2[3].z, o2[4].x, o2[4].y, o2[4].z, o2[5].x, o2[5].y, o2[5].z, o2[6].x );
	
	
	Vector3 n[] = {{999,999,999},{999,999,999},{999,999,999},{999,999,999},{999,999,999},};
	
	vpTriFaceNormal3p(o1, &n[0]);
	vpTriFaceNormal3p(o1+3, &n[1]);
	vpTriFaceNormal3p(o2, &n[2]);
	vpTriFaceNormal3p(o2+3, &n[3]);
	vpTriFaceNormal3p(tri3, &n[4]);
	
	printf("  nin: % .2f,% .2f,% .2f\n", n[4].x, n[4].y, n[4].z);
	printf("  na1: % .2f,% .2f,% .2f\n", n[0].x, n[0].y, n[0].z);
	printf("  na2: % .2f,% .2f,% .2f\n", n[1].x, n[1].y, n[1].z);
	printf("  nb1: % .2f,% .2f,% .2f\n", n[2].x, n[2].y, n[2].z);
	printf("  nb2: % .2f,% .2f,% .2f\n", n[3].x, n[3].y, n[3].z);
	printf("  %f, %f, %f\n", vDot3p(n, n+2), vDot3p(n, n+1), vDot3p(n+2, n+3));
	
}





int main(int argc, char* argv) {
	
	Matrix m = IDENT_MATRIX;
	
	mOrtho(-2.2, 3.3, 4.4, -5.5, 0.66, 7.7, &m);
	
	float a,b,c,d,e,f;
	
	mOrthoExtractPlanes(&m, &a, &b, &c, &d, &e, &f);
	
	printf("%f,%f,%f,%f,%f,%f\n", a, b,c,d,e,f);
	
	
	return 0;
}
