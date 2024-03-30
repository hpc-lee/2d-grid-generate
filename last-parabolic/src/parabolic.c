#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "parabolic.h"
#include "solver.h"
#include "lib_mem.h"

int 
para_gene(gd_t *gdcurv, par_t *par)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;

  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *step = gdcurv->step;

  float coef = par->coef;
  int o2i = par->o2i;

  if(o2i == 1)
  {
    fprintf(stdout,"parabolic method, march direction from inner bdry(k=0) to outer bdry(nz-1)\n");
    fprintf(stdout,"so if we want from outer bdry(nz-1) to inner bdry(k=0), must be flip\n");
    fprintf(stdout,"flip coord, trans outer bdry(nz-1) to inner bdry(0)\n");
    flip_coord_z(gdcurv);
  }

  // malloc space for thomas method 
  // a,b,c,d_x,d_z,u_x,u_z. 7 vars space
  float *var_th = (float *)mem_calloc_1d_float((nx-2)*7, 0.0, "init");

  // calculate step_length
  float *step_len = (float *)mem_calloc_1d_float(nz, 0.0, "init");
  for(int k=1; k<nz; k++)
  {
    step_len[k] = step_len[k-1] + step[k-1];
  }


  for(int k=1; k<nz-1; k++)
  {
    // predict k+1 layer points
    predict_point(x2d,z2d,nx,nz,k,o2i,coef,step_len);
    // base predict points
    // update k layer points
    update_point(x2d,z2d,var_th,nx,k);
    assign_bdry_coords(x2d,z2d,nx,k);
    fprintf(stdout,"number of layer is %d\n",k+1);
    fflush(stdout);
  }

  if(o2i == 1)
  {
    // flip to return.
    flip_coord_z(gdcurv);
  }

  free(var_th);
  free(step_len);

  return 0;
}

int 
predict_point(float *x2d, float *z2d, int nx, int nz, int k, int o2i, 
              float coef, float *step_len)
{
  // k-1 layer points is know
  // predict points k+1 and k layer
  // first predict further k+1 layer
  // base predict k+1 layer, calculate k layer
  // NOTE: this predict point only used
  // by calculate matrix or coefficient
  // not the final point
  //
  // cal k-1 layer point unit normal vector
  // vt -> vector tangential
  // vn -> vector normal
  
  float zt,cs;
  size_t iptr1,iptr2,iptr3,iptr4;
  size_t iptr5,iptr6;
  int sign1;
  float vt_x,vt_z,vn_x,vn_z,len_vt;
  float R_x,R_z,R;
  float R1, r1, r2;
  float c1,c2;
  float x0,z0,xs,zs;

  // cal switching factor
  zt = (1.0*k)/(nz-1);
  cs = exp(coef*zt);
  // o2i outer bdry to inner bdry
  if(o2i == 1) { 
    sign1 = 1;
  } else {
    sign1 = -1;
  }

  if(k < nz-2)
  {
    for(int i=1; i<nx-1; i++)
    {
      // cal normal vector 
      iptr1 = (k-1)*nx + (i+1); // (i+1,k-1)
      iptr2 = (k-1)*nx + (i-1); // (i-1,k-1)
      vt_x = 0.5*(x2d[iptr1] - x2d[iptr2]); 
      vt_z = 0.5*(z2d[iptr1] - z2d[iptr2]); 
      len_vt = sqrt(pow(vt_x,2)+pow(vt_z,2));
      vn_x =  sign1*vt_z/len_vt;
      vn_z = -sign1*vt_x/len_vt;
      
      // inner point i
      iptr1 = (k-1)*nx + i;  // (i,k-1)
      iptr2 = (nz-1)*nx + i; // (i,nz-1)
      R_x = x2d[iptr1] - x2d[iptr2];
      R_z = z2d[iptr1] - z2d[iptr2];
      R = sqrt(pow(R_x,2)+pow(R_z,2));

      // cal clustering factor
      R1 = step_len[nz-1] - step_len[k-1];
      r1 = step_len[k+1] - step_len[k-1];
      r2 = step_len[k] - step_len[k-1];
      c1 = r1/R1;
      c2 = r2/r1;
      
      iptr1 = (k-1)*nx + i; //(i,k-1)
      // x0,z0 normal point
      x0 = x2d[iptr1] + vn_x*c1*R;
      z0 = z2d[iptr1] + vn_z*c1*R;
      
      // xs,zs linear distance point
      iptr1 = (k-1)*nx + i; //(i,k-1)
      iptr2 = (nz-1)*nx + i;   //(i,nz-1)
      xs = x2d[iptr1] + c1*(x2d[iptr2]-x2d[iptr1]);
      zs = z2d[iptr1] + c1*(z2d[iptr2]-z2d[iptr1]);

      iptr3 = (k+1)*nx + i;     //(i,k+1)
      x2d[iptr3] = cs*x0 + (1-cs)*xs;
      z2d[iptr3] = cs*z0 + (1-cs)*zs;

      iptr4 = k*nx + i;    //(i,k)
      x2d[iptr4] = x2d[iptr1] + c2*(x2d[iptr3]-x2d[iptr1]);
      z2d[iptr4] = z2d[iptr1] + c2*(z2d[iptr3]-z2d[iptr1]);
    }

    // float bdry
    iptr1 = (k-1)*nx+0;   // (0,k-1)
    iptr2 = (k-1)*nx+1;   // (1,k-1)
    iptr3 = k*nx+0;   // (0,k)
    iptr4 = k*nx+1;   // (1,k)
    x2d[iptr3] = x2d[iptr1] + (x2d[iptr4] - x2d[iptr2]);
    z2d[iptr3] = z2d[iptr1] + (z2d[iptr4] - z2d[iptr2]);

    iptr5 = (k+1)*nx+0;   // (0,k+1)
    iptr6 = (k+1)*nx+1;   // (1,k+1)
    x2d[iptr5] = x2d[iptr1] + (x2d[iptr6] - x2d[iptr2]);
    z2d[iptr5] = z2d[iptr1] + (z2d[iptr6] - z2d[iptr2]);

    iptr1 = (k-1)*nx+(nx-1);   // (nx-1,k-1)
    iptr2 = (k-1)*nx+(nx-2);   // (nx-2,k-1)
    iptr3 = k*nx+(nx-1);   // (nx-1,k)
    iptr4 = k*nx+(nx-2);   // (nx-2,k)
    x2d[iptr3] = x2d[iptr1] + (x2d[iptr4] - x2d[iptr2]);
    z2d[iptr3] = z2d[iptr1] + (z2d[iptr4] - z2d[iptr2]);

    iptr5 = (k+1)*nx+(nx-1);   // (nx-1,k+1)
    iptr6 = (k+1)*nx+(nx-2);   // (nx-2,k+1)
    x2d[iptr5] = x2d[iptr1] + (x2d[iptr6] - x2d[iptr2]);
    z2d[iptr5] = z2d[iptr1] + (z2d[iptr6] - z2d[iptr2]);
  }

  if(k == nz-2)
  {
    for(int i=0; i<nx; i++)
    {
      r1 = step_len[k+1] - step_len[k-1];
      r2 = step_len[k] - step_len[k-1];
      c2 = r2/r1;
      iptr1 = (k-1)*nx + i; //(i,nz-3)
      iptr3 = (k+1)*nx + i; //(i,nz-1)
      iptr4 = k*nx + i;     //(i,nz-2)
      x2d[iptr4] = x2d[iptr1] + c2*(x2d[iptr3]-x2d[iptr1]);
      z2d[iptr4] = z2d[iptr1] + c2*(z2d[iptr3]-z2d[iptr1]);
    }
  }

  return 0;
}

int
update_point(float *x2d, float *z2d, float *var_th, int nx, int k)
{
  float *a;
  float *b;
  float *c;
  float *d_x;
  float *d_z;
  float *u_x;
  float *u_z;
  int siz_vec = nx-2;
  a = var_th;
  b = var_th + siz_vec;
  c = var_th + 2*siz_vec;
  d_x = var_th + 3*siz_vec;
  d_z = var_th + 4*siz_vec;
  u_x = var_th + 5*siz_vec;
  u_z = var_th + 6*siz_vec;

  size_t iptr,iptr1,iptr2,iptr3,iptr4;
  float x_xi,z_xi,x_zt,z_zt;
  float temp_x,temp_z;
  float x_xizt,z_xizt;
  float g11,g22,g12;

  for(int i=1; i<nx-1; i++)
  {
    iptr1 =  k*nx + (i-1); // (i-1,k)
    iptr2 =  k*nx + (i+1); // (i+1,k)
    iptr3 = (k-1)*nx + i;  // (i,k-1)
    iptr4 = (k+1)*nx + i;  // (i,k+1)

    x_xi = 0.5*(x2d[iptr2] - x2d[iptr1]);
    z_xi = 0.5*(z2d[iptr2] - z2d[iptr1]);

    x_zt = 0.5*(x2d[iptr4] - x2d[iptr3]);
    z_zt = 0.5*(z2d[iptr4] - z2d[iptr3]);

    temp_x = x2d[iptr4] + x2d[iptr3];
    temp_z = z2d[iptr4] + z2d[iptr3];

    iptr1 = (k-1)*nx + (i-1);   // (i-1,k-1)
    iptr2 = (k-1)*nx + (i+1);   // (i+1,k-1)
    iptr3 = (k+1)*nx + (i-1);   // (i-1,k+1)
    iptr4 = (k+1)*nx + (i+1);   // (i+1,k+1)

    x_xizt = 0.25*(x2d[iptr4] + x2d[iptr1] - x2d[iptr2] - x2d[iptr3]);
    z_xizt = 0.25*(z2d[iptr4] + z2d[iptr1] - z2d[iptr2] - z2d[iptr3]);

    g11 = x_xi*x_xi + z_xi*z_xi;
    g22 = x_zt*x_zt + z_zt*z_zt;
    g12 = x_xi*x_zt + z_xi*z_zt;

    a[i-1] = g22;
    b[i-1] = -2*(g22+g11);
    c[i-1] = g22;

    d_x[i-1] = -g11*temp_x + 2*g12*x_xizt;
    d_z[i-1] = -g11*temp_z + 2*g12*z_xizt;
  }

  // i=1 modify a,d_x,d_z
  iptr = k*nx+0;
  d_x[0] = d_x[0] - a[0]*x2d[iptr];
  d_z[0] = d_z[0] - a[0]*z2d[iptr];
  a[0] = 0;

  // i=nx-2 modify c,d_x,d_z
  iptr = k*nx + (nx-1);
  d_x[nx-3] = d_x[nx-3] - c[nx-3]*x2d[iptr];
  d_z[nx-3] = d_z[nx-3] - c[nx-3]*z2d[iptr];
  c[nx-3] = 0;
  
  // cal coords and update
  thomas(siz_vec,a,b,c,d_x,d_z,u_x,u_z);
  for(int i=1; i<nx-1; i++)
  {
    iptr = k*nx + i;
    x2d[iptr] = u_x[i-1];
    z2d[iptr] = u_z[i-1];
  }

  return 0;
}


int 
assign_bdry_coords(float *x2d, float *z2d, int nx, int k)
{
  size_t iptr1,iptr2,iptr3,iptr4;
  iptr1 = (k-1)*nx+0;   // (0,k-1)
  iptr2 =  k*nx+0;      // (0,k)

  iptr3 = (k-1)*nx+1;   // (1,k-1)
  iptr4 =  k*nx+1;      // (1,k)

  x2d[iptr2] = x2d[iptr1] + (x2d[iptr4]-x2d[iptr3]);
  z2d[iptr2] = z2d[iptr1] + (z2d[iptr4]-z2d[iptr3]);

  iptr1 = (k-1)*nx+(nx-1);   // (nx-1,k-1)
  iptr2 =  k*nx+(nx-1);      // (nx-1,k)

  iptr3 = (k-1)*nx+(nx-2);   // (nx-2,k-1)
  iptr4 =  k*nx+(nx-2);      // (nx-2,k)

  x2d[iptr2] = x2d[iptr1] + (x2d[iptr4]-x2d[iptr3]);
  z2d[iptr2] = z2d[iptr1] + (z2d[iptr4]-z2d[iptr3]);

  return 0;
}

