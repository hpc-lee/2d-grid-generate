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

// hermite interpolation
// four boundary x1(left), x2(right), z1(bottom) ,z2(top)
int
one_hermite(gd_t *gdcurv, float coef)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float zeta,c0,c1,c10,c11;
  float tan_z1_x,tan_z1_z,tan_z2_x,tan_z2_z;
  float x_len,z_len;
  float abs_z1, abs_z2;
  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float *dh_len;
  float *nor_z1_x;
  float *nor_z1_z;
  float *nor_z2_x;
  float *nor_z2_z;

  dh_len = (float *)mem_calloc_1d_float(nx, 0.0, "init");
  nor_z1_x = (float *)mem_calloc_1d_float(nx, 0.0, "init");
  nor_z1_z = (float *)mem_calloc_1d_float(nx, 0.0, "init");
  nor_z2_x = (float *)mem_calloc_1d_float(nx, 0.0, "init");
  nor_z2_z = (float *)mem_calloc_1d_float(nx, 0.0, "init");

  for (int i=0; i<nx; i++)
  {
    iptr1 = (nz-1)*nx + i; //(i,nz-1)
    iptr2 = i;             //(i,1)
    x_len = x2d[iptr1] - x2d[iptr2];
    z_len = z2d[iptr1] - z2d[iptr2];
    dh_len[i] = sqrt(pow(x_len,2) + pow(z_len,2));
  }

  for (int i=1; i<nx-1; i++) 
  {
    iptr1 = i+1;   //(i+1,1)
    iptr2 = i-1;   //(i-1,1)
    iptr3 = (nz-1)*nx + (i+1);   //(i+1,nz-1)
    iptr4 = (nz-1)*nx + (i-1);   //(i-1,nz-1)
    tan_z1_x = (x2d[iptr1] - x2d[iptr2])/2;
    tan_z1_z = (z2d[iptr1] - z2d[iptr2])/2;
    tan_z2_x = (x2d[iptr3] - x2d[iptr4])/2;
    tan_z2_z = (z2d[iptr3] - z2d[iptr4])/2;
    nor_z1_x[i] = -tan_z1_z;
    nor_z1_z[i] = tan_z1_x;
    nor_z2_x[i] = -tan_z2_z;
    nor_z2_z[i] = tan_z2_x;
  }
  
  // i=0 
  iptr1 = 1;   //(1,1)
  iptr2 = 0;   //(0,1)
  iptr3 = (nz-1)*nx + 1;   //(1,nz-1)
  iptr4 = (nz-1)*nx + 0;   //(0,nz-1)
  nor_z1_x[0] = -(z2d[iptr1]-z2d[iptr2]);
  nor_z1_z[0] = x2d[iptr1]-x2d[iptr2];
  nor_z2_x[0] = -(z2d[iptr3]-z2d[iptr4]);
  nor_z2_z[0] = x2d[iptr3]-x2d[iptr4];
  // i=nx-1
  iptr1 = nx-1;   //(nx-1,1)
  iptr2 = nx-2;   //(nx-2,1)
  iptr3 = (nz-1)*nx + nx-1;   //(nx-1,nz-1)
  iptr4 = (nz-1)*nx + nx-2;   //(nx-2,nz-1)
  nor_z1_x[nx-1] = -(z2d[iptr1]-z2d[iptr2]);
  nor_z1_z[nx-1] = x2d[iptr1]-x2d[iptr2];
  nor_z2_x[nx-1] = -(z2d[iptr3]-z2d[iptr4]);
  nor_z2_z[nx-1] = x2d[iptr3]-x2d[iptr4];

  // cal 1st order derivative coef, this effect 
  // boundary orth length
  for (int i=0; i<nx; i++) 
  {
    abs_z1 = sqrt(pow(nor_z1_x[i],2) + pow(nor_z1_z[i],2));
    abs_z2 = sqrt(pow(nor_z2_x[i],2) + pow(nor_z2_z[i],2));
    nor_z1_x[i] = coef*dh_len[i]*(nor_z1_x[i]/abs_z1);
    nor_z1_z[i] = coef*dh_len[i]*(nor_z1_z[i]/abs_z1);
    nor_z2_x[i] = coef*dh_len[i]*(nor_z2_x[i]/abs_z2);
    nor_z2_z[i] = coef*dh_len[i]*(nor_z2_z[i]/abs_z2);
  }

  for (int i=0; i<nx; i++) {
    for (int k=1; k<nz-1; k++)
    {
      zeta = (1.0*k)/(nz-1);
      c0 =  2*pow(zeta,3) - 3*pow(zeta,2) + 1;
      c1 = -2*pow(zeta,3) + 3*pow(zeta,2);
      c10 = pow(zeta,3)-2*pow(zeta,2)+zeta;
      c11 = pow(zeta,3)-pow(zeta,2);
      
      iptr = k*nx + i;
      iptr1 = i;
      iptr2 = (nz-1)*nx+i;
      x2d[iptr] = c0*x2d[iptr1] + c1*x2d[iptr2] + c10*nor_z1_x[i] + c11*nor_z2_x[i];
      z2d[iptr] = c0*z2d[iptr1] + c1*z2d[iptr2] + c10*nor_z1_z[i] + c11*nor_z2_z[i];
    }
  }

  free(dh_len);
  free(nor_z1_x); 
  free(nor_z1_z);
  free(nor_z2_x);
  free(nor_z2_z);

  return 0;
}

