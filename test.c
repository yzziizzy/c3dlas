#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "c3dlas.h"




int main(int argc, char* argv) {
	
	Plane pl = {
		{0,0,1}, 1,
	};
	
	Vector tri[] = {
		{3,4,2},
		{3,5,3},
		{3,6,-5},
	};
	
	Vector tri2[] = {
		{3,4,2},
		{3,5,3},
		{3,6,7},
	};
	Vector tri3[] = {
		{3,4,2},
		{3,5,3},
		{3,6,1},
	};
	
	Vector l[] = {
		{3, 3, -4},
		{4, 4, 5},
		{4, 4, 9},
	};
	
	Vector o = {999,999, 999};
	
	printf("%s, %f,%f,%f\n", c3dlas_EnumString(planeLineFindIntersect(&pl, &l[0], &l[1], &o)), o.x, o.y, o.z);
	
	o = (Vector){999,999, 999};

	printf("%s, %f,%f,%f\n", c3dlas_EnumString(planeLineFindIntersect(&pl, &l[1], &l[2], &o)), o.x, o.y, o.z);

	return 0;
}

