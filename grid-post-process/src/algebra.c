#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "algebra.h"
#include "lib_mem.h"

// strech grid base on arc length 
// use exponential function
int 
zt_arc_stretch(gd_t *gdcurv, float *arc_len)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t siz_icmp = gdcurv->siz_icmp;
  size_t iptr,iptr1,iptr2;
  float x_len,z_len,dh_len;
  float r, ratio, zeta;
  int n;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *x2d_temp; 
  float *z2d_temp;
  float *s;
  float *u;

  x2d_temp = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");
  z2d_temp = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");
  s = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");
  u = (float *)mem_calloc_1d_float(
              nz, 0.0, "init");

 
  // line by line. i=0 -> i=nx-1 
  for(int i=0; i<nx; i++)
  {
    // copy old coords to temp space
    for(int k=0; k<nz; k++)
    {
      iptr1 = k*nx + i;     //(i,k)
      x2d_temp[k] = x2d[iptr1];
      z2d_temp[k] = z2d[iptr1];
    }
    // cal arc length
    for(int k=1; k<nz; k++)
    {
      x_len = x2d_temp[k] - x2d_temp[k-1];
      z_len = z2d_temp[k] - z2d_temp[k-1];
      dh_len = sqrt(pow(x_len,2) + pow(z_len,2));
      s[k] = s[k-1] + dh_len;
    }
    // arc length normalized
    for(int k=0; k<nz; k++)
    {
      u[k] = s[k]/s[nz-1];
    }
    for(int k=1; k<nz-1; k++)
    {
      r = arc_len[k];
      for(int m=0; m<nz-1; m++)
      {
        if(r>=u[m] && r<u[m+1]) {
          n=m; 
          break;
        }
      }

      // linear interp
      iptr = k*nx + i;
      x_len = x2d_temp[n+1] - x2d_temp[n];
      z_len = z2d_temp[n+1] - z2d_temp[n];
      ratio = (r - u[n])/(u[n+1]-u[n]);
      x2d[iptr] = x2d_temp[n] + x_len*ratio;
      z2d[iptr] = z2d_temp[n] + z_len*ratio;
    }
  }

  free(x2d_temp);
  free(z2d_temp);
  free(s);
  free(u);

  return 0;
}

// strech grid base on arc length 
// use exponential function
int 
xi_arc_stretch(gd_t *gdcurv, float *arc_len)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t siz_icmp = gdcurv->siz_icmp;
  size_t iptr,iptr1;
  float x_len,z_len,dh_len;
  float r, ratio, xi;
  int n;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *x2d_temp; 
  float *z2d_temp;
  float *s;
  float *u;

  x2d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  z2d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  s = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  u = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
 
  // line by line. k=0 -> k=nz-1 
  for(int k=0; k<nz; k++)
  {
    // copy old coords to temp space
    for(int i=0; i<nx; i++)
    {
      iptr1 = k*nx + i;     //(i,k)
      x2d_temp[i] = x2d[iptr1];
      z2d_temp[i] = z2d[iptr1];
    }
    // cal arc length
    for(int i=1; i<nx; i++)
    {
      x_len = x2d_temp[i] - x2d_temp[i-1];
      z_len = z2d_temp[i] - z2d_temp[i-1];
      dh_len = sqrt(pow(x_len,2) + pow(z_len,2));
      s[i] = s[i-1] + dh_len;
    }
    // arc length normalized
    for(int i=0; i<nx; i++)
    {
      u[i] = s[i]/s[nx-1];
    }
    for(int i=1; i<nx-1; i++)
    {
      r = arc_len[i];
      for(int m=0; m<nx-1; m++)
      {
        if(r>=u[m] && r<u[m+1]) {
          n=m; 
          break;
        }
      }

      // linear interp
      iptr = k*nx + i;
      x_len = x2d_temp[n+1] - x2d_temp[n];
      z_len = z2d_temp[n+1] - z2d_temp[n];
      ratio = (r - u[n])/(u[n+1]-u[n]);
      x2d[iptr] = x2d_temp[n] + x_len*ratio;
      z2d[iptr] = z2d_temp[n] + z_len*ratio;
    }
  }

  free(x2d_temp);
  free(z2d_temp);
  free(s);
  free(u);

  return 0;
}

// grid sample, linear interpolation  
int 
sample_interp(gd_t *gdcurv_new, gd_t *gdcurv)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int nx_new = gdcurv_new->nx;
  int nz_new = gdcurv_new->nz;

  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *x2d_new = gdcurv_new->x2d;
  float *z2d_new = gdcurv_new->z2d;

  float *x2d_temp;
  float *z2d_temp;
  float *u;
  float *v;
  size_t iptr,iptr1,iptr2;
  float x_len,z_len;
  float r, ratio;
  int n;

  x2d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  z2d_temp = (float *)mem_calloc_1d_float(
              nx, 0.0, "init");
  u = (float *)mem_calloc_1d_float(
              nz, 0.0, "gd_curv_init");
  v = (float *)mem_calloc_1d_float(
              nx, 0.0, "gd_curv_init");

  // first interp zt direction 
  // line by line. i=0 -> i=nx-1 
  for(int i=0; i<nx; i++)
  {
    // point number normalized [0,1]
    for(int k=0; k<nz; k++)
    {
      u[k] = (1.0*k)/(nz-1);
    }

    for(int k_new=0; k_new<nz_new; k_new++)
    {
      r = (1.0*k_new)/(nz_new-1);
      for(int m=0; m<nz-1; m++)
      {
        if(r>=u[m] && r<u[m+1]) {
          n=m; 
          break;
        }
      }

      // linear interp
      iptr = k_new*nx_new + i;
      iptr1 = n*nx + i;
      iptr2 = (n+1)*nx + i;
      x_len = x2d[iptr2] - x2d[iptr1];
      z_len = z2d[iptr2] - z2d[iptr1];
      ratio = (r - u[n])/(u[n+1]-u[n]);
      x2d_new[iptr] = x2d[iptr1] + x_len*ratio;
      z2d_new[iptr] = z2d[iptr1] + z_len*ratio;
    }
  }

  // then interp xi direction 
  // line by line. k=0 -> k=nz_new-1 
  for(int k_new=0; k_new<nz_new; k_new++)
  {
    // copy old coords to temp space
    for(int i=0; i<nx; i++)
    {
      iptr1 = k_new*nx_new + i;     //(i,k)
      x2d_temp[i] = x2d_new[iptr1];
      z2d_temp[i] = z2d_new[iptr1];
    }
    // point number normalized [0,1]
    for(int i=0; i<nx; i++)
    {
      v[i] = (1.0*i)/(nx-1);
    }

    for(int i_new=0; i_new<nx_new; i_new++)
    {
      r = (1.0*i_new)/(nx_new-1);
      for(int m=0; m<nx-1; m++)
      {
        if(r>=v[m] && r<v[m+1]) {
          n=m; 
          break;
        }
      }

      // linear interp
      iptr = k_new*nx_new + i_new;
      x_len = x2d_temp[n+1] - x2d_temp[n];
      z_len = z2d_temp[n+1] - z2d_temp[n];
      ratio = (r - v[n])/(v[n+1]-v[n]);
      x2d_new[iptr] = x2d_temp[n] + x_len*ratio;
      z2d_new[iptr] = z2d_temp[n] + z_len*ratio;
    }
  }

  free(u);
  free(v);
  free(x2d_temp);
  free(z2d_temp);
  // use sample grid gdcurv_new, so free gdcurv space 
  free(gdcurv->v3d);

  return 0;
}

