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
int linear_tfi(gd_t *gdcurv, bdry_t *bdry, mympi_t *mympi)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int gni1 = gdcurv->gni1;
  int gnk1 = gdcurv->gnk1;
  int gni2 = gdcurv->gni2;
  int gnk2 = gdcurv->gnk2;
  size_t siz_iz = gdcurv->siz_iz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;

  int nx_all = bdry->number_of_grid_points_x;
  int nz_all = bdry->number_of_grid_points_z;
  float *x1 = bdry->x1;
  float *x2 = bdry->x2;
  float *z1 = bdry->z1;
  float *z2 = bdry->z2;


  int gni, gnk;
  float xi,zt;
  float a0,a1,c0,c1;
  size_t iptr;
  float U_x,W_x,UW_x;
  float U_z,W_z,UW_z;
  
  for (int k=nk1; k<=nk2; k++) {
    for (int i=ni1; i<=ni2; i++)
    {
      gni = gni1 + i; 
      gnk = gnk1 + k; 
      xi = (1.0*gni)/(nx_all-1); // 1.0* equal int to float
      zt = (1.0*gnk)/(nz_all-1);
      a0 = 1-xi;
      a1 = xi;
      c0 = 1-zt;
      c1 = zt;

      iptr = k*siz_iz + i;  //(i,k) 
      // 4 boudary points
      // 4 corner points

      U_x = a0*x1[gnk] + a1*x2[gnk];
      W_x = c0*z1[gni] + c1*z2[gni];
      UW_x = a0*c0*x1[0] + a0*c1*x1[nz_all-1] 
           + a1*c0*x2[0] + a1*c1*x2[nz_all-1];

      U_z = a0*x1[gnk+nz_all] + a1*x2[gnk+nz_all];
      W_z = c0*z1[gni+nx_all] + c1*z2[gni+nx_all];
      UW_z = a0*c0*x1[0+nz_all] + a0*c1*x1[nz_all-1+nz_all] 
           + a1*c0*x2[0+nz_all] + a1*c1*x2[nz_all-1+nz_all];

      x2d[iptr] = U_x + W_x - UW_x;
      z2d[iptr] = U_z + W_z - UW_z;
    }
  }

  // assgn four bdry point to ghost
  // bdry x1
  if(mympi->neighid[0] == MPI_PROC_NULL)
  {
    for (int k=0; k<nz; k++)
    {
      iptr = k*siz_iz + 0;  //(0,k) 
      gnk = gnk1 + k; 

      x2d[iptr] = x1[gnk];
      z2d[iptr] = x1[gnk+nz_all];
    }
  }
  // bdry x2
  if(mympi->neighid[1] == MPI_PROC_NULL)
  {
    for (int k=0; k<nz; k++)
    {
      iptr = k*siz_iz + nx-1;  //(nx-1,k) 
      gnk = gnk1 + k; 

      x2d[iptr] = x2[gnk];
      z2d[iptr] = x2[gnk+nz_all];
    }
  }
  // bdry z1
  if(mympi->neighid[2] == MPI_PROC_NULL)
  {
    for (int i=0; i<nx; i++)
    {
      iptr = 0*siz_iz + i;  //(i,0) 
      gni = gni1 + i; 

      x2d[iptr] = z1[gni];
      z2d[iptr] = z1[gni+nx_all];
    }
  }
  // bdry z2
  if(mympi->neighid[3] == MPI_PROC_NULL)
  {
    for (int i=0; i<nx; i++)
    {
      iptr = (nz-1)*siz_iz + i;  //(i,nz-1) 
      gni = gni1 + i; 

      x2d[iptr] = z2[gni];
      z2d[iptr] = z2[gni+nx_all];
    }
  }
  
  return 0;
}

