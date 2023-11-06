#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "solver.h"
#include "lib_math.h"
#include "lib_mem.h"
#include "constants.h"

// solve tridiagonal linear system equaltion
// using thomas method
// Due to the sparse coefficient matrix of the tridiagonal equation,
// the computational complexity is proportional to n,
// rather than the n^3 of Gaussian elimination
//
// [b1 c1        ]
// |a2 b2 c2     |
// |   a3 b3 c3  |
// [             ]
//
// a = element of lower diagonal
// b = element of main diagonal
// c = element of up diagonal
// d = right hand item
// u is unknow vector
int
thomas(int n, float *a, float *b, float *c, float *d_x, 
       float *d_z, float *u_x, float *u_z)
{
  float factor;
  for(int i=1; i<n; i++)
  {
    factor = a[i]/b[i-1];
    b[i] = b[i] - factor*c[i-1];
    d_x[i] = d_x[i] - factor*d_x[i-1];
    d_z[i] = d_z[i] - factor*d_z[i-1];
  }

  u_x[n-1] = d_x[n-1]/b[n-1];
  u_z[n-1] = d_z[n-1]/b[n-1];
  for(int i=n-2; i>=0; i--)
  {
    u_x[i] = (d_x[i]-c[i]*u_x[i+1])/b[i];
    u_z[i] = (d_z[i]-c[i]*u_z[i+1])/b[i];
  }

  return 0;
}

/*
  solve block tridiagonal linear system equaltion
  using thomas method
  Due to the sparse coefficient matrix of the tridiagonal equation,
  the computational complexity is proportional to n,
  rather than the n^3 of Gaussian elimination 

  [b1 c1        ]
  |a2 b2 c2     |
  |   a3 b3 c3  |
  [             ]

  a = square matrix element of lower diagonal
  b = square matrix element of main diagonal
  c = square matrix element of up diagonal
  d = right hand item, each element size is n*1

    [G1 0  0       ]
    |a2 G2 0       |
    |   a3 G3      |
  L=|              |
    |              |
    |              |
    [         an Gn]
                    
    [I1 D1 0            ]
    |0  I2 D2           |
    |   0  I3           |
  U=|                   |
    |                   |
    |          In-1 Dn-1|
    [          0    In  ]

  LU*x = d
  Ly=d
  Ux=y
*/
int
thomas_block(int n, double *a, double *b, double *c, double *d,
             double *xz, double *D, double *y)
{
  double mat_a[2][2], mat_b[2][2], mat_c[2][2], vec_d[2];
  double mat_G[2][2], mat_D[2][2], vec_y[2], vec_xz[2];
  double det, mat1[2][2], vec1[2], vec2[2];
  size_t iptr1,iptr2,iptr3,iptr4;

  // i=0
  for(int ii=0; ii<2; ii++) {
    for(int jj=0; jj<2; jj++) {
      mat_G[ii][jj] = b[ii*2+jj];
      mat_c[ii][jj] = c[ii*2+jj];
    }
    vec_d[ii] = d[ii];
  }

  det = mat_G[0][0]*mat_G[1][1] - mat_G[0][1]*mat_G[1][0];
  if(fabs(det)<0.000001)
  {
    fprintf(stderr,"b{1} is singular, stop calculate");
    exit(1);
  } 
  mat_invert2x2(mat_G);
  mat_mul2x2(mat_G,mat_c,mat_D);
  mat_mul2x1(mat_G,vec_d,vec_y);
  for(int ii=0; ii<2; ii++) { 
    for(int jj=0; jj<2; jj++) {
      D[ii*2+jj] = mat_D[ii][jj];
    }
    y[ii] = vec_y[ii];
  }

  // i=1~n-1
  for(int i=1; i<n; i++)
  {
    iptr1 = i*CONST_NDIM*CONST_NDIM;
    iptr2 = i*CONST_NDIM;
    iptr3 = (i-1)*CONST_NDIM*CONST_NDIM;
    iptr4 = (i-1)*CONST_NDIM;
    for(int ii=0; ii<2; ii++) { 
      for(int jj=0; jj<2; jj++) {
        mat_a[ii][jj] = a[iptr1+ii*2+jj];
        mat_b[ii][jj] = b[iptr1+ii*2+jj];
        mat_c[ii][jj] = c[iptr1+ii*2+jj];
        mat_D[ii][jj] = D[iptr3+ii*2+jj];
      }
      vec_d[ii] = d[iptr2+ii];
      vec_y[ii] = y[iptr4+ii];
    }
    mat_mul2x2(mat_a,mat_D,mat1);  // a(i)*D(i-1)
    mat_sub2x2(mat_b,mat1,mat_G);  // G(i) = b(i)-a(i)*D(i-1)
    mat_invert2x2(mat_G);          // inv(G(i))
    mat_mul2x2(mat_G,mat_c,mat_D); // D(i) = inv(G(i))*c(i)
    mat_mul2x1(mat_a,vec_y,vec1);  // a(i) * y(i-1)
    vec_sub2x1(vec_d,vec1,vec2);   // d(i) - a(i) * y(i-1)
    mat_mul2x1(mat_G,vec2,vec_y);  // y(i) = inv(G(i))*(d(i) - a(i) * y(i-1))

    for(int ii=0; ii<2; ii++) { 
      for(int jj=0; jj<2; jj++) {
        D[iptr1+ii*2+jj] = mat_D[ii][jj]; 
      }
      y[iptr2+ii] = vec_y[ii];
    }
  }
  
  // i=n-1
  iptr2 = (n-1)*CONST_NDIM;
  for(int ii=0; ii<2; ii++) { 
    xz[iptr2+ii] = y[iptr2+ii];
  }

  for(int i=n-2; i>=0; i--)
  {
    iptr1 = i*CONST_NDIM*CONST_NDIM;
    iptr2 = i*CONST_NDIM;
    iptr4 = (i+1)*CONST_NDIM;
    
    for(int ii=0; ii<2; ii++) { 
      for(int jj=0; jj<2; jj++) {
        mat_D[ii][jj] = D[iptr1+ii*2+jj]; 
      }
      vec_xz[ii] = xz[iptr4+ii];
    }
    mat_mul2x1(mat_D,vec_xz,vec1);  // D(i)*xz(i+1)
    for(int ii=0; ii<2; ii++) 
    { 
      xz[iptr2+ii] = y[iptr2+ii] - vec1[ii]; // xz(i) = y(i) - D(i)*xz(i+1) 
    }
  }
  
  return 0;
}

