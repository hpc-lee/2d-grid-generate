#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "algebra.h"
#include "lib_mem.h"

// linear tfi interpolation
// U means x-direction
// W means z-direction
int linear_tfi(gd_t *gdcurv)
{
  int nx = gdcurv->nx; 
  int nz = gdcurv->nz; 
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float xi,zeta;
  float a0,a1,c0,c1;
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  size_t iptr5,iptr6,iptr7,iptr8;
  float U_x,W_x,UW_x;
  float U_z,W_z,UW_z;
  
  for (int k=1; k<nz-1; k++) {
    for (int i=1; i<nx-1; i++)
    {
      xi   = (1.0*i)/(nx-1); // 1.0* equal int to float
      zeta = (1.0*k)/(nz-1);
      a0 = 1-xi;
      a1 = xi;
      c0 = 1-zeta;
      c1 = zeta;

      iptr = k*nx + i;  //(i,k) 
      // 4 boudary points
      iptr1 = k*nx;     //(0,k)
      iptr2 = k*nx + nx-1; //(nx-1,k)
      iptr3 = i;         //(i,0)
      iptr4 = (nz-1)*nx + i; //(i,nz-1)
      // 4 corner points
      iptr5 = 0;   //(0,0)
      iptr6 = (nz-1)*nx;   //(0,nz-1)
      iptr7 = nx-1;   //(nx-1,0)
      iptr8 = (nz-1)*nx + nx-1;   //(nx-1,nz-1)

      U_x = a0*x2d[iptr1] + a1*x2d[iptr2];
      W_x = c0*x2d[iptr3] + c1*x2d[iptr4];
      UW_x = a0*c0*x2d[iptr5] + a0*c1*x2d[iptr6] \
           + a1*c0*x2d[iptr7] + a1*c1*x2d[iptr8];

      U_z = a0*z2d[iptr1] + a1*z2d[iptr2];
      W_z = c0*z2d[iptr3] + c1*z2d[iptr4];
      UW_z = a0*c0*z2d[iptr5] + a0*c1*z2d[iptr6] \
           + a1*c0*z2d[iptr7] + a1*c1*z2d[iptr8];

      x2d[iptr] = U_x + W_x - UW_x;
      z2d[iptr] = U_z + W_z - UW_z;
    }
  }
  
  return 0;
}

