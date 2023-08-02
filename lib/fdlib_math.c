#include <stdio.h>
#include <math.h>

#include "fdlib_math.h"

void
fdlib_math_invert2x2(float matrix[2][2])
{
  for (int k=0; k<2; k++)
  {
     float con = matrix[k][k];
     matrix[k][k] = 1.0;

     for (int i=0; i<2; i++) {
       matrix[k][i] = matrix[k][i]/con;
     }

     for (int i=0; i<2; i++)
     {
        if (i!=k) {
           con = matrix[i][k];
           matrix[i][k] = 0.0;
           for (int j=0; j<2; j++) {
            matrix[i][j] = matrix[i][j] - matrix[k][j] * con;
           }
        }
     }
  }

  return ;
}

void
fdlib_math_invert3x3(float m[][3])
{
  float inv[3][3];
  float det;
  int i, j;

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

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      m[i][j] = inv[i][j] * det;

  return ;
}

void
fdlib_math_matmul2x2(float A[][2], float B[][2], float C[][2])
{
  int i, j, k;
  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++){
      C[i][j] = 0.0;
      for (k = 0; k < 2; k++)
        C[i][j] += A[i][k] * B[k][j];
    }

  return ;
}

void
fdlib_math_matmul3x3(float A[][3], float B[][3], float C[][3])
{
  int i, j, k;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
      C[i][j] = 0.0;
      for (k = 0; k < 3; k++)
        C[i][j] += A[i][k] * B[k][j];
    }

  return ;
}

void
fdlib_math_cross_product(float *A, float *B, float *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1];
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
  return;
}

float
fdlib_math_dot_product(float *A, float *B)
{
  int i;
  float result = 0.0;
  for (i = 0; i < 3; i++)
    result += A[i] * B[i];

  return result;
}
