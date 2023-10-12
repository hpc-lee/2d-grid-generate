#ifndef FDLIB_MATH_H
#define FDLIB_MATH_H

int
mat_invert2x2(double matrix[][2]);

int
mat_mul2x2(double A[][2], double B[][2], double C[][2]);

int
mat_mul2x1(double A[][2], double *B, double *C);

int
mat_add2x2(double A[][2], double B[][2], double C[][2]);

int
vec_add2x1(double *A, double *B, double *C);

int
vec_sub2x1(double *A, double *B, double *C);

int
mat_sub2x2(double A[][2], double B[][2], double C[][2]);

int
mat_copy2x2(double A[][2], double B[][2]);

int
mat_iden2x2(double A[][2]);

#endif
