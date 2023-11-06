#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "solver.h"
#include "lib_math.h"
#include "lib_mem.h"
#include "constants.h"

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
