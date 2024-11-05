#include <stdio.h>
#include <math.h>

#include "lib_math.h"

int
mat_invert3x3(float m[][3])
{
  float inv[3][3];
  double det;

  inv[0][0] = m[1][1]*m[2][2] - m[2][1]*m[1][2];
  inv[0][1] = m[2][1]*m[0][2] - m[0][1]*m[2][2];
  inv[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
  inv[1][0] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
  inv[1][1] = m[0][0]*m[2][2] - m[2][0]*m[0][2];
  inv[1][2] = m[1][0]*m[0][2] - m[0][0]*m[1][2];
  inv[2][0] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
  inv[2][1] = m[2][0]*m[0][1] - m[0][0]*m[2][1];
  inv[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];

  det = inv[0][0] * m[0][0] 
      + inv[0][1] * m[1][0] 
      + inv[0][2] * m[2][0];

  det = 1.0f / det;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) 
    {
      m[i][j] = inv[i][j] * det;
    }
  }

  return 0;
}

int
mat_mul3x3(float A[][3], float B[][3], float C[][3])
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C[i][j] = 0.0;
      for (int k = 0; k < 3; k++)
      {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  return 0;
}

int
cross_product(float *A, float *B, float *C)
{
  C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0];

  return 0;
}

float
dot_product(float *A, float *B)
{
  float result = 0.0;
  for(int i=0; i<3; i++)
  {
    result += A[i]*B[i];
  }

  return result;
}

int 
mat_mul3x1(float A[][3], float *B, float *C)
{
  for(int i=0; i<3; i++)
  {
    C[i] = A[i][0]*B[0] + A[i][1]*B[1] + A[i][2]*B[2];
  }

  return 0;
}

int
mat_add3x3(float A[][3], float B[][3], float C[][3])
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
        C[i][j] = A[i][j] + B[i][j];
    }
  }

  return 0;
}

int
vec_add3x1(float *A, float *B, float *C)
{
  for (int i=0; i<3; i++){
    C[i] = A[i] + B[i];
  }

  return 0;
}

int
vec_sub3x1(float *A, float *B, float *C)
{
  for (int i=0; i<3; i++){
    C[i] = A[i] - B[i];
  }

  return 0;
}

int
mat_sub3x3(float A[][3], float B[][3], float C[][3])
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
        C[i][j] = A[i][j] - B[i][j];
    }
  }

  return 0;
}

int
mat_copy3x3(float A[][3], float B[][3])
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
        B[i][j] = A[i][j];
    }
  }

  return 0;
}

int
mat_iden3x3(float A[][3])
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      A[i][j] = 0.0;
      if(i==j) {
        A[i][j] = 1;
      }
    }
  }

  return 0;
}

float
dist_point2line(float p0[2], float p1[2], float p2[2])
{
  float L;
  float A,B,C;

  A=p2[1]-p1[1];
  B=p1[0]-p2[0];
  C=-p1[1]*B-p1[0]*A;

  L=fabs((A*p0[0]+B*p0[1]+C)/sqrt(A*A+B*B));

  return L;
}
