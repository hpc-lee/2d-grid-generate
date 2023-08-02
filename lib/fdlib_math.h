#ifndef FDLIB_MATH_H
#define FDLIB_MATH_H

void
fdlib_math_invert2x2(float matrix[2][2]);

void
fdlib_math_invert3x3(float m[][3]);

void
fdlib_math_matmul2x2(float A[][2], float B[][2], float C[][2]);

void
fdlib_math_matmul3x3(float A[][3], float B[][3], float C[][3]);

void
fdlib_math_cross_product(float *A, float *B, float *C);

float
fdlib_math_dot_product(float *A, float *B);

#endif
