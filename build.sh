#!/bin/bash

	
# -ffast-math but without reciprocal approximations
MATH_CFLAGS="\
	-fno-math-errno \
	-fexcess-precision=fast \
	-fno-signed-zeros -fno-trapping-math -fassociative-math \
	-ffinite-math-only -fno-rounding-math \
	-fno-signaling-nans
"

# flags for first-gen Ryzen
# change if your cpu is different
CPU_CFLAGS="-march=native -mtune=native"

ERR_CFLAGS="\
	-Wall \
	-Wextra \
	-Wno-unused-result \
	-Wno-unused-function \
	-Werror-implicit-function-declaration \
	-Wno-discarded-qualifiers
"

gcc -o test_c3dlas test.c c3dlas.c -O0 -ggdb -lm  $ERR_CFLAGS $CPU_CFLAGS $MATH_CFLAGS

