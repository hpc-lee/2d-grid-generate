#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "solver.h"

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

int
update_SOR(float *x2d, float *z2d, float *x2d_tmp, float *z2d_tmp,
           int nx, int nz, float *P, float *Q, float w)
{
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float x_xi,z_xi,x_zt,z_zt,x_xizt,z_xizt;
  float g11,g22,g12,coef;
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 =  k*nx + (i+1); // (i+1,k)
      iptr2 =  k*nx + (i-1); // (i-1,k)
      iptr3 = (k+1)*nx + i;  // (i,k+1)
      iptr4 = (k-1)*nx + i;  // (i,k-1)

      x_xi = 0.5*(x2d[iptr1] - x2d_tmp[iptr2]);
      z_xi = 0.5*(z2d[iptr1] - z2d_tmp[iptr2]);

      x_zt = 0.5*(x2d[iptr3] - x2d_tmp[iptr4]);
      z_zt = 0.5*(z2d[iptr3] - z2d_tmp[iptr4]);

      iptr1 = (k-1)*nx + (i-1);   // (i-1,k-1)
      iptr2 = (k-1)*nx + (i+1);   // (i+1,k-1)
      iptr3 = (k+1)*nx + (i+1);   // (i+1,k+1)
      iptr4 = (k+1)*nx + (i-1);   // (i-1,k+1)

      x_xizt = 0.25*(x2d[iptr3] + x2d_tmp[iptr1] - x2d_tmp[iptr2] - x2d[iptr4]);
      z_xizt = 0.25*(z2d[iptr3] + z2d_tmp[iptr1] - z2d_tmp[iptr2] - z2d[iptr4]);

      g11 = x_xi*x_xi + z_xi*z_xi;
      g22 = x_zt*x_zt + z_zt*z_zt;
      g12 = x_xi*x_zt + z_xi*z_zt;

      coef = 0.5/(g22+g11);
      
      iptr  =  k*nx + i;     // (i,k)
      iptr1 =  k*nx + (i+1); // (i+1,k)
      iptr2 =  k*nx + (i-1); // (i-1,k)
      iptr3 = (k+1)*nx + i;  // (i,k+1)
      iptr4 = (k-1)*nx + i;  // (i,k-1)

      x2d_tmp[iptr] = coef*(g22*(x2d[iptr1]+x2d_tmp[iptr2]) + g11*(x2d[iptr3]
                      + x2d_tmp[iptr4]) - 2*g12*x_xizt + g22*P[iptr]*x_xi
                      + g11*Q[iptr]*x_zt);

      x2d_tmp[iptr] = w*x2d_tmp[iptr] + (1-w)*x2d[iptr];

      z2d_tmp[iptr] = coef*(g22*(z2d[iptr1]+z2d_tmp[iptr2]) + g11*(z2d[iptr3]
                      + z2d_tmp[iptr4]) - 2*g12*z_xizt + g22*P[iptr]*z_xi
                      + g11*Q[iptr]*z_zt);

      z2d_tmp[iptr] = w*z2d_tmp[iptr] + (1-w)*z2d[iptr];

    }
  }

  return 0;
}
