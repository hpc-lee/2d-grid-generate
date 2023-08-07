#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "elliptic.h"
#include "constants.h"
#include "solver.h"
#include "fdlib_mem.h"

int
higen_gene(gd_t *gdcurv, par_t *par)
{
  float d1 = par->distance[0];
  float d2 = par->distance[1];
  float err = par->i_err;
  int max_iter = par->max_iter;
  float coef = par->coef;

  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;

  float *P = NULL; //source term P
  float *Q = NULL; //source term Q
  P = (float *)fdlib_mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source P");
  Q = (float *)fdlib_mem_calloc_1d_float(gdcurv->siz_icmp, 0.0, "source Q");

  float w; // SOR coef
  // now we default set w=1.0.
  // this is G-S
  w=1.0;

  higen_P_SOR(x2d,z2d,nx,nz,P,Q,d1,d2,coef,err,max_iter,w);

  return 0;
}

int 
higen_P_SOR(float *x2d, float *z2d, int nx, int nz, float *P, float *Q,
            float d1, float d2, float coef, float err, int max_iter, float w)
{
  float *x2d_tmp =NULL;
  float *z2d_tmp =NULL;
  x2d_tmp = (float *)fdlib_mem_calloc_1d_float(nx*nz, 0.0, "x2d_tmp");
  z2d_tmp = (float *)fdlib_mem_calloc_1d_float(nx*nz, 0.0, "z2d_tmp");

  int Niter = 0;
  size_t iptr;
  float resi, resk;
  // copy coord
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++)
    {
      iptr = k*nx + i;
      x2d_tmp[iptr] = x2d[iptr];
      z2d_tmp[iptr] = z2d[iptr];
    }
  }

  int flag_true = 1;
  while(flag_true)
  {
    // update solver
    update_SOR(x2d,z2d,x2d_tmp,z2d_tmp,nx,nz,P,Q,w);
    Niter += 1;

    resi = 0.0;
    resk = 0.0;
    // cal iter error
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr = k*nx + i;
        resi += fabs((x2d_tmp[iptr] - x2d[iptr])/(x2d_tmp[iptr]+pow(10,-6)));
        resk += fabs((z2d_tmp[iptr] - z2d[iptr])/(z2d_tmp[iptr]+pow(10,-6)));
      }
    }

    resi = resi/((nx-2)*(nz-2));
    resk = resk/((nx-2)*(nz-2));

    // copy coord
    for(int k=0; k<nz; k++) {
      for(int i=0; i<nx; i++)
      {
        iptr = k*nx + i;
        x2d[iptr] = x2d_tmp[iptr];
        z2d[iptr] = z2d_tmp[iptr];
      }
    }
    
    set_src_higen(x2d,z2d,P,Q,coef,nx,nz,d1,d2);

    if(Niter>max_iter) {
      flag_true = 0;
    }

    if(resi < err && resk < err) {
      flag_true = 0;
    }

    fprintf(stdout,"number of iter is %d\n", Niter);
    fprintf(stdout,"resi is %f, resk is %f\n", resi, resk);
  }

  free(x2d_tmp);
  free(z2d_tmp);

  return 0;
}

int
set_src_higen(float *x2d, float *z2d, float *P, float *Q,
              float coef, int nx, int nz, float d1, float d2)
{
  float theta0 = PI/2;
  float a = 0.1;
  size_t iptr,iptr1,iptr2,iptr3;
  float dot, len_xi, len_zt, dif_dis;
  float cos_theta, theta, dif_theta;
  float x_xi, z_xi, x_zt, z_zt;
  // bdry z1 zt=0
  for(int i=1; i<nx-1; i++)
  {
    iptr  = i;         //(i,0)
    iptr1 = i+1;       //(i+1,0)
    iptr2 = i-1;       //(i-1,0)
    iptr3 = 1*nx + i;  //(i,1)

    x_xi = 0.5*(x2d[iptr1]-x2d[iptr2]);
    z_xi = 0.5*(z2d[iptr1]-z2d[iptr2]);
    x_zt = x2d[iptr3] - x2d[iptr];
    z_zt = z2d[iptr3] - z2d[iptr];
    
    // cos(theta) = a.b/(|a|*|b|)
    dot = x_xi*x_zt + z_xi*z_zt;
    len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
    len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
    cos_theta = dot/(len_xi*len_zt);
    theta = acos(cos_theta);
    dif_theta = (theta0-theta)/theta0;
    P[iptr] = P[iptr] - a*tanh(dif_theta);

    dif_dis = (d1-len_zt)/d1;
    Q[iptr] = Q[iptr] + a*tanh(dif_dis);
  }

  // NOTE: P and Q update need change sign
  // P is +a and Q is -a
  // bdry z2 zt=0
  for(int i=1; i<nx-1; i++)
  {
    iptr  = (nz-1)*nx + i;      //(i,nz-1)
    iptr1 = (nz-1)*nx + i+1;    //(i+1,nz-1)
    iptr2 = (nz-1)*nx + i-1;    //(i-1,nz-1)
    iptr3 = (nz-2)*nx + i;      //(i,nz-2)

    x_xi = 0.5*(x2d[iptr1]-x2d[iptr2]);
    z_xi = 0.5*(z2d[iptr1]-z2d[iptr2]);
    x_zt = x2d[iptr] - x2d[iptr3];
    z_zt = z2d[iptr] - z2d[iptr3];
    
    // cos(theta) = a.b/(|a|*|b|)
    dot = x_xi*x_zt + z_xi*z_zt;
    len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
    len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
    cos_theta = dot/(len_xi*len_zt);
    theta = acos(cos_theta);
    dif_theta = (theta0-theta)/theta0;
    P[iptr] = P[iptr] + a*tanh(dif_theta);

    dif_dis = (d2-len_zt)/d2;
    Q[iptr] = Q[iptr] - a*tanh(dif_dis);
  }

  // x1 x2 is less important bdry
  // only cal Q to contral bdry orth
  // bdry x1 xi=0
  for(int k=1; k<nz-1; k++)
  {
    iptr  =  k*nx + 0;         //(0,k)
    iptr1 = (k+1)*nx + 0;      //(0,k+1)
    iptr2 = (k-1)*nx + 0;      //(0,k-1)
    iptr3 = k*nx + 1;  //(1,k)

    x_zt = 0.5*(x2d[iptr1]-x2d[iptr2]);
    z_zt = 0.5*(z2d[iptr1]-z2d[iptr2]);
    x_xi = x2d[iptr3] - x2d[iptr];
    z_xi = z2d[iptr3] - z2d[iptr];
    
    // cos(theta) = a.b/(|a|*|b|)
    dot = x_xi*x_zt + z_xi*z_zt;
    len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
    len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
    cos_theta = dot/(len_xi*len_zt);
    theta = acos(cos_theta);
    dif_theta = (theta0-theta)/theta0;
    Q[iptr] = Q[iptr] - a*tanh(dif_theta);
  }
  // bdry x2 xi=1
  for(int k=1; k<nz-1; k++)
  {
    iptr  =  k*nx + (nx-1);         //(nx-1,k)
    iptr1 = (k+1)*nx + (nx-1);      //(nx-1,k+1)
    iptr2 = (k-1)*nx + (nx-1);      //(nx-1,k-1)
    iptr3 = k*nx + (nx-2);  //(nx-2,k)

    x_zt = 0.5*(x2d[iptr1]-x2d[iptr2]);
    z_zt = 0.5*(z2d[iptr1]-z2d[iptr2]);
    x_xi = x2d[iptr] - x2d[iptr3];
    z_xi = z2d[iptr] - z2d[iptr3];
    
    // cos(theta) = a.b/(|a|*|b|)
    dot = x_xi*x_zt + z_xi*z_zt;
    len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
    len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
    cos_theta = dot/(len_xi*len_zt);
    theta = acos(cos_theta);
    dif_theta = (theta0-theta)/theta0;
    Q[iptr] = Q[iptr] + a*tanh(dif_theta);
  }

  // first use bdry z1 z2 to interp inner point
  float xi,zt,c0,c1,r0,r1;
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      zt = (1.0*k)/(nz-1);
      c0 = 1-zt;
      c1 = zt;

      r0 = exp(coef*zt);
      r1 = exp(coef*(1-zt)); 
      
      iptr  = k*nx + i;
      iptr1 = i;
      iptr2 = (nz-1)*nx + i;
      P[iptr] = r0*P[iptr1] + r1*P[iptr2];
      Q[iptr] = c0*Q[iptr1] + c1*Q[iptr2];
    }
  }
  // then use bdry x1 x2 to interp inner point
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      xi = (1.0*i)/(nx-1);

      r0 = exp(coef*xi);
      r1 = exp(coef*(1-xi)); 
      
      iptr  = k*nx + i;
      iptr1 = k*nx + 0;
      iptr2 = k*nx + (nx-1);
      Q[iptr] = Q[iptr] + r0*Q[iptr1] + r1*Q[iptr2];
    }
  }

  return 0;
}
