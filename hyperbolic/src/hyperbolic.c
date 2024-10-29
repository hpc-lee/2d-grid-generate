#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "hyperbolic.h"
#include "solver.h"
#include "lib_mem.h"
#include "lib_math.h"
#include "constants.h"
#include "io_funcs.h"

int 
hyper_gene(gd_t *gdcurv, par_t *par)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int n = nx-2;  // not include bdry 2 points
  float coef = par->coef; 
  int t2b = par->t2b;
  int *bdry_itype = par->bdry_itype;
  float *epsilon = par->epsilon;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *step = gdcurv->step;

  float *bdry1 = NULL;
  float *bdry2 = NULL;
  FILE *fp = NULL;
  char str[500];

  float *coef_e = (float *)mem_calloc_1d_float(nx, 0.0, "init");
  float *area = (float *)mem_calloc_1d_float(nx*2, 0.0, "init");
  // malloc space for thomas_block method 
  double *a = (double *)mem_calloc_1d_double(n*CONST_NDIM*CONST_NDIM, 0.0, "init");
  double *b = (double *)mem_calloc_1d_double(n*CONST_NDIM*CONST_NDIM, 0.0, "init");
  double *c = (double *)mem_calloc_1d_double(n*CONST_NDIM*CONST_NDIM, 0.0, "init");
  double *d = (double *)mem_calloc_1d_double(n*CONST_NDIM, 0.0, "init");
  double *xz = (double *)mem_calloc_1d_double(n*CONST_NDIM, 0.0, "init");
  double *D = (double *)mem_calloc_1d_double(n*CONST_NDIM*CONST_NDIM, 0.0, "init");
  double *y = (double *)mem_calloc_1d_double(n*CONST_NDIM, 0.0, "init");

  // solve first layer  
  int k=1;
  cal_matrix(x2d,z2d,nx,k,step,a,b,c,d,area);
  modify_smooth(x2d,z2d,nx,k,a,b,c,d,coef_e);
  modify_bdry(n,a,b,c,d,epsilon,bdry_itype,bdry1,bdry2,k,nz);
  thomas_block(n,a,b,c,d,xz,D,y);
  assign_coords(xz,x2d,z2d,nx,nz,k,epsilon,bdry_itype,bdry1,bdry2);

  for(int k=1; k<nz; k++)
  {
    cal_smooth_coef(coef,x2d,z2d,nx,nz,k,step,coef_e);
    cal_matrix(x2d,z2d,nx,k,step,a,b,c,d,area);
    modify_smooth(x2d,z2d,nx,k,a,b,c,d,coef_e);
    modify_bdry(n,a,b,c,d,epsilon,bdry_itype,bdry1,bdry2,k,nz);
    thomas_block(n,a,b,c,d,xz,D,y);
    assign_coords(xz,x2d,z2d,nx,nz,k,epsilon,bdry_itype,bdry1,bdry2);

    fprintf(stdout,"number of layers is %d\n",k);
    fflush(stdout);
  }

  if(t2b== 0)
  {
    //fprintf(stdout,"the init bdry is max index, so index must be flip\n");
    flip_coord_z(gdcurv);
  }

  free(coef_e);
  free(area);
  free(a);
  free(b);
  free(c);
  free(d);
  free(xz);
  free(D);
  free(y);
  free(bdry1);
  free(bdry2);

  return 0;
}

int
cal_smooth_coef(float coef, float *x2d, float *z2d, int nx, int nz, int k, float *step, float *coef_e)
{
  float S;
  size_t iptr1, iptr2, iptr3, iptr4;
  float x_xi,z_xi,x_zt,z_zt;
  float xi_len,zt_len,N_xi;
  float x_xi_plus,z_xi_plus,x_xi_minus,z_xi_minus;
  float xi_plus1,xi_minus1,xi_plus2,xi_minus2;
  float d1,d2,delta,delta_mdfy;
  float x_plus,z_plus,x_minus,z_minus;
  float dot,det,theta,alpha;
  float temp;

  S = sqrt((1.0*k)/(nz-1));
  //if(nz<=50) {
  //  S = sqrt((1.0*k)/(nz-1));
  //}

  //if(nz>50) 
  //{
  //  if(k<=50)
  //  {
  //    S = sqrt((1.0*k)/50);
  //  } else {
  //    //S = 1.0;
  //    S = 1.0 + 1.0*(k-50)/(nz-51);
  //  }
  //}

  if(k==1)
  {
    int k1=2;
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (k1-1)*nx + i+1;   // (i+1,k-1)
      iptr2 = (k1-1)*nx + i-1;   // (i-1,k-1)
      iptr3 = (k1-1)*nx + i;     // (i,k-1)
      iptr4 = (k1-2)*nx + i;     // (i,k-2)
      x_xi = 0.5*(x2d[iptr1] - x2d[iptr2]);
      z_xi = 0.5*(z2d[iptr1] - z2d[iptr2]); 
      x_zt = x2d[iptr3] - x2d[iptr4];
      z_zt = z2d[iptr3] - z2d[iptr4];
      xi_len = sqrt(pow(x_xi,2) + pow(z_xi,2));
      zt_len = sqrt(pow(x_zt,2) + pow(z_zt,2));
      N_xi = zt_len/xi_len;

      iptr1 = (k1-2)*nx + i+1;   // (i+1,k-2)
      iptr2 = (k1-2)*nx + i;     // (i,  k-2)
      iptr3 = (k1-2)*nx + i-1;   // (i-1,k-2)
      x_xi_plus = x2d[iptr1] - x2d[iptr2];
      z_xi_plus = z2d[iptr1] - z2d[iptr2];
      xi_plus1 = sqrt(pow(x_xi_plus,2) + pow(z_xi_plus,2));
      x_xi_minus = x2d[iptr3] - x2d[iptr2];
      z_xi_minus = z2d[iptr3] - z2d[iptr2];
      xi_minus1 = sqrt(pow(x_xi_minus,2) + pow(z_xi_minus,2));

      iptr1 = (k1-1)*nx + i+1;   // (i+1,k-1)
      iptr2 = (k1-1)*nx + i;     // (i,  k-1)
      iptr3 = (k1-1)*nx + i-1;   // (i-1,k-1)
      x_xi_plus = x2d[iptr1] - x2d[iptr2];
      z_xi_plus = z2d[iptr1] - z2d[iptr2];
      xi_plus2 = sqrt(pow(x_xi_plus,2) + pow(z_xi_plus,2));
      x_xi_minus = x2d[iptr3] - x2d[iptr2];
      z_xi_minus = z2d[iptr3] - z2d[iptr2];
      xi_minus2 = sqrt(pow(x_xi_minus,2) + pow(z_xi_minus,2));

      d1 = xi_plus1 + xi_minus1;
      d2 = xi_plus2 + xi_minus2;
      delta = d1/d2;
      delta_mdfy = fmax(pow(delta,2/S),0.01);
      // normalization
      x_plus = (x2d[iptr1]-x2d[iptr2])/xi_plus2;
      z_plus = (z2d[iptr1]-z2d[iptr2])/xi_plus2;
      x_minus = (x2d[iptr3]-x2d[iptr2])/xi_minus2;
      z_minus = (z2d[iptr3]-z2d[iptr2])/xi_minus2;

      dot = x_plus*x_minus + z_plus*z_minus;
      det = x_plus*z_minus - z_plus*x_minus;

      // cal two normal vector clockwise angle. 
      // the method from website
      // from right vector to left vector
      // z axis upward, so is -det
      theta = atan2(-det,dot);
      if(theta<0)
      {
        theta = theta + 2*PI;
      }
      if(theta<PI)
      {
        alpha = 1.0/(1-pow(cos(theta/2),2));
      }
      if(theta>=PI)
      {
        //alpha = pow(sin(theta/2),10);
        alpha = 1;
      }
      coef_e[i] = coef*N_xi*S*delta_mdfy*alpha;
    }
  }

  if(k>1)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (k-1)*nx + i+1;   // (i+1,k-1)
      iptr2 = (k-1)*nx + i-1;   // (i-1,k-1)
      iptr3 = (k-1)*nx + i;     // (i,k-1)
      iptr4 = (k-2)*nx + i;     // (i,k-2)
      x_xi = 0.5*(x2d[iptr1] - x2d[iptr2]);
      z_xi = 0.5*(z2d[iptr1] - z2d[iptr2]); 
      x_zt = x2d[iptr3] - x2d[iptr4];
      z_zt = z2d[iptr3] - z2d[iptr4];
      xi_len = sqrt(pow(x_xi,2) + pow(z_xi,2));
      zt_len = sqrt(pow(x_zt,2) + pow(z_zt,2));
      N_xi = zt_len/xi_len;

      iptr1 = (k-2)*nx + i+1;   // (i+1,k-2)
      iptr2 = (k-2)*nx + i;     // (i,  k-2)
      iptr3 = (k-2)*nx + i-1;   // (i-1,k-2)
      x_xi_plus = x2d[iptr1] - x2d[iptr2];
      z_xi_plus = z2d[iptr1] - z2d[iptr2];
      xi_plus1 = sqrt(pow(x_xi_plus,2) + pow(z_xi_plus,2));
      x_xi_minus = x2d[iptr3] - x2d[iptr2];
      z_xi_minus = z2d[iptr3] - z2d[iptr2];
      xi_minus1 = sqrt(pow(x_xi_minus,2) + pow(z_xi_minus,2));

      iptr1 = (k-1)*nx + i+1;   // (i+1,k-1)
      iptr2 = (k-1)*nx + i;     // (i,  k-1)
      iptr3 = (k-1)*nx + i-1;   // (i-1,k-1)
      x_xi_plus = x2d[iptr1] - x2d[iptr2];
      z_xi_plus = z2d[iptr1] - z2d[iptr2];
      xi_plus2 = sqrt(pow(x_xi_plus,2) + pow(z_xi_plus,2));
      x_xi_minus = x2d[iptr3] - x2d[iptr2];
      z_xi_minus = z2d[iptr3] - z2d[iptr2];
      xi_minus2 = sqrt(pow(x_xi_minus,2) + pow(z_xi_minus,2));

      d1 = xi_plus1 + xi_minus1;
      d2 = xi_plus2 + xi_minus2;
      delta = d1/d2;
      delta_mdfy = fmax(pow(delta,2/S),0.01);

      // normalization
      x_plus = (x2d[iptr1]-x2d[iptr2])/xi_plus2;
      z_plus = (z2d[iptr1]-z2d[iptr2])/xi_plus2;
      x_minus = (x2d[iptr3]-x2d[iptr2])/xi_minus2;
      z_minus = (z2d[iptr3]-z2d[iptr2])/xi_minus2;

      dot = x_plus*x_minus + z_plus*z_minus;
      det = x_plus*z_minus - z_plus*x_minus;

      // cal two normal vector clockwise angle. 
      // the method from website
      // from right vector to left vector
      // z axis upward, so is -det
      theta = atan2(-det,dot);
      if(theta<0)
      {
        theta = theta + 2*PI;
      }
      if(theta<PI)
      {
        alpha = 1.0/(1-pow(cos(theta/2),2));
      }
      if(theta>=PI)
      {
        //alpha = pow(sin(theta/2),2);
        alpha = 1;
      }
      coef_e[i] = coef*N_xi*S*delta_mdfy*alpha;
    }
  }

  return 0;
}

int 
cal_matrix(float *x2d, float *z2d, int nx, int k, float *step,
           double *a, double *b, double *c, double *d, float *area)
{
  double A[2][2], B[2][2];
  double mat[2][2], vec[2];
  double mat_b[2][2], vec_d[2];
  size_t iptr1,iptr2,iptr3;
  float x_xi0,z_xi0,x_zt0,z_zt0;
  float x_plus,z_plus,x_minus,z_minus;
  float arc_plus,arc_minus,arc_len,temp;

  if(k==1)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (k-1)*nx + i+1;
      iptr2 = (k-1)*nx + i;
      iptr3 = (k-1)*nx + i-1;
      x_xi0 = 0.5*(x2d[iptr1] - x2d[iptr3]);
      z_xi0 = 0.5*(z2d[iptr1] - z2d[iptr3]);
      x_plus = x2d[iptr1] - x2d[iptr2];
      z_plus = z2d[iptr1] - z2d[iptr2];
      arc_plus = sqrt(pow(x_plus,2) + pow(z_plus,2));
      x_minus = x2d[iptr3] - x2d[iptr2];
      z_minus = z2d[iptr3] - z2d[iptr2];
      arc_minus = sqrt(pow(x_minus,2) + pow(z_minus,2));
      arc_len = 0.5*(arc_plus + arc_minus);
      // arc_length -> area
      // area(i) = A1 
      // area(i+nx) = A0
      area[i] = arc_len * step[k-1];
      area[i+nx] = area[i];
      temp = pow(x_xi0,2) + pow(z_xi0,2);
      x_zt0 = -z_xi0*area[i+nx]/temp;
      z_zt0 =  x_xi0*area[i+nx]/temp;
      // add damping factor, maybe inv(B) singular
      A[0][0] = x_zt0;      A[0][1] = z_zt0;
      A[1][0] = z_zt0;      A[1][1] =-x_zt0; 
      B[0][0] = x_xi0+1e-7; B[0][1] = z_xi0;
      B[1][0] =-z_xi0;      B[1][1] = x_xi0+1e-7; 
      mat_invert2x2(B);
      mat_mul2x2(B,A,mat);
      mat_iden2x2(mat_b);
      vec[0] = 0; vec[1] = area[i];
      mat_mul2x1(B,vec,vec_d);
      iptr1 = (i-1)*CONST_NDIM*CONST_NDIM;
      iptr2 = (i-1)*CONST_NDIM;
      for(int ii=0; ii<2; ii++) {
        for(int jj=0; jj<2; jj++) {
          a[iptr1+2*ii+jj] = -0.5*mat[ii][jj];
          b[iptr1+2*ii+jj] = mat_b[ii][jj];
          c[iptr1+2*ii+jj] = 0.5*mat[ii][jj];
        }
        d[iptr2+ii] = vec_d[ii];
      }
    }
  } else {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (k-1)*nx + i+1;
      iptr2 = (k-1)*nx + i;
      iptr3 = (k-1)*nx + i-1;
      x_xi0 = 0.5*(x2d[iptr1] - x2d[iptr3]);
      z_xi0 = 0.5*(z2d[iptr1] - z2d[iptr3]);
      x_plus = x2d[iptr1] - x2d[iptr2];
      z_plus = z2d[iptr1] - z2d[iptr2];
      arc_plus = sqrt(pow(x_plus,2) + pow(z_plus,2));
      x_minus = x2d[iptr3] - x2d[iptr2];
      z_minus = z2d[iptr3] - z2d[iptr2];
      arc_minus = sqrt(pow(x_minus,2) + pow(z_minus,2));
      arc_len = 0.5*(arc_plus + arc_minus);
      // arc_length -> area
      // area(i) = A1 
      // area(i+nx) = A0
      area[i+nx] = area[i];
      area[i] = arc_len * step[k-1];
      temp = pow(x_xi0,2) + pow(z_xi0,2);
      x_zt0 = -z_xi0*area[i+nx]/temp;
      z_zt0 =  x_xi0*area[i+nx]/temp;
      A[0][0] = x_zt0+1e-7; A[0][1] = z_zt0;
      A[1][0] = z_zt0;      A[1][1] =-x_zt0+1e-7; 
      B[0][0] = x_xi0+1e-7; B[0][1] = z_xi0;
      B[1][0] =-z_xi0;      B[1][1] = x_xi0+1e-7; 
      mat_invert2x2(B);
      mat_mul2x2(B,A,mat);
      mat_iden2x2(mat_b);
      vec[0] = 0; vec[1] = area[i];
      mat_mul2x1(B,vec,vec_d);
      iptr1 = (i-1)*CONST_NDIM*CONST_NDIM;
      iptr2 = (i-1)*CONST_NDIM;
      for(int ii=0; ii<2; ii++) {
        for(int jj=0; jj<2; jj++) {
          a[iptr1+2*ii+jj] = -0.5*mat[ii][jj];
          b[iptr1+2*ii+jj] = mat_b[ii][jj];
          c[iptr1+2*ii+jj] = 0.5*mat[ii][jj];
        }
        d[iptr2+ii] = vec_d[ii];
      }
    }
  }

  return 0;
}

int
modify_smooth(float *x2d, float *z2d, int nx, int k, double *a,
              double *b, double *c, double *d, float *coef_e)
{
  double mat[2][2], vec1[2], vec2[2], vec3[2];
  float coef_i;
  size_t iptr1, iptr2, iptr3, iptr4, iptr5;
  mat_iden2x2(mat);
  for(int i=1; i<nx-1; i++)
  {
    iptr1 = (k-1)*nx + i-1;
    iptr2 = (k-1)*nx + i;
    iptr3 = (k-1)*nx + i+1;
    vec1[0] = x2d[iptr1];
    vec1[1] = z2d[iptr1];
    vec2[0] = x2d[iptr2];
    vec2[1] = z2d[iptr2];
    vec3[0] = x2d[iptr3];
    vec3[1] = z2d[iptr3];

    iptr4 = (i-1)*CONST_NDIM*CONST_NDIM;
    iptr5 = (i-1)*CONST_NDIM;
    
    coef_i = 2*coef_e[i];

    for(int ii=0; ii<2; ii++) {
      for(int jj=0; jj<2; jj++) {
        a[iptr4+2*ii+jj] = a[iptr4+2*ii+jj] - coef_i*mat[ii][jj];
        b[iptr4+2*ii+jj] = b[iptr4+2*ii+jj] + 2*coef_i*mat[ii][jj];
        c[iptr4+2*ii+jj] = c[iptr4+2*ii+jj] - coef_i*mat[ii][jj];
      }
      d[iptr5+ii] = d[iptr5+ii] + coef_e[i]*(vec1[ii]+vec3[ii]-2*vec2[ii]);
    }
  }

  return 0;
}

int
modify_bdry(int n, double *a, double *b, double *c, double *d,
            float *epsilon, int *bdry_itype, 
            float *bdry1, float *bdry2, int k, int nz)
{
  size_t iptr;
  // left bdry
  // float boundry
  if(bdry_itype[0] == 1)
  {
    for(int ii=0; ii<2; ii++) {
      for(int jj=0; jj<2; jj++) {
        // modify i=0
        b[ii*2+jj] = b[ii*2+jj] + (1+epsilon[0])*a[ii*2+jj];
        c[ii*2+jj] = c[ii*2+jj] - epsilon[0]*a[ii*2+jj];
      }
    }
  }
  else if(bdry_itype[0] == 2)  // dx=0 cartesian boundry
  {
    // only modify second column
    int jj = 1;
    for(int ii=0; ii<2; ii++) {
      // modify i=0
      b[ii*2+jj] = b[ii*2+jj] + a[ii*2+jj];
    }
  }
  else if(bdry_itype[0] == 3) // fixed boundry modify d = d-A*delta(r1)
  {
    double A[2][2], vec_1[2], vec[2];
    for(int ii=0; ii<2; ii++) {
      for(int jj=0; jj<2; jj++) {
        A[ii][jj] = a[ii*2+jj];
      }
      vec_1[ii] = bdry1[k+ii*nz]-bdry1[(k-1)+ii*nz];  // coord increment
    }
    mat_mul2x1(A,vec_1,vec);
    for(int ii=0; ii<2; ii++)
    {
        d[ii] = d[ii] - vec[ii];
    }
  }
  // right bdry
  // float boundry
  if(bdry_itype[1] == 1)
  {
    iptr = (n-1)*CONST_NDIM*CONST_NDIM;
    for(int ii=0; ii<2; ii++) {
      for(int jj=0; jj<2; jj++) {
        // modify i=n-1
        b[iptr+ii*2+jj] = b[iptr+ii*2+jj] + (1+epsilon[1])*c[iptr+ii*2+jj];
        a[iptr+ii*2+jj] = a[iptr+ii*2+jj] - epsilon[1]*c[iptr+ii*2+jj];
      }
    }
  }
  else if(bdry_itype[1] == 2)  // dx=0 cartesian boundry
  {
    iptr = (n-1)*CONST_NDIM*CONST_NDIM;
    // only modify second column
    int jj = 1;
    for(int ii=0; ii<2; ii++) {
      // modify i=nx-3
      b[iptr+ii*2+jj] = b[iptr+ii*2+jj] + c[iptr+ii*2+jj];
    }
  }
  else if(bdry_itype[1] == 3) // fixed boundry
  {
    iptr = (n-1)*CONST_NDIM*CONST_NDIM;
    double A[2][2], vec_2[2], vec[2];
    for(int ii=0; ii<2; ii++) {
      for(int jj=0; jj<2; jj++) {
        A[ii][jj] = a[iptr+ii*2+jj];
      }
      vec_2[ii] = bdry2[k+ii*nz]-bdry2[(k-1)+ii*nz];  // coord increment
    }
    mat_mul2x1(A,vec_2,vec);
    iptr = (n-1)*CONST_NDIM;
    for(int ii=0; ii<2; ii++) 
    {
      d[iptr+ii] = d[iptr+ii] - vec[ii];
    }
  }

  return 0;
}

int
assign_coords(double *xz, float *x2d, float *z2d, int nx, int nz, int k,
              float *epsilon, int *bdry_itype, float *bdry1, float *bdry2)
{
  size_t iptr,iptr1,iptr2;
  size_t iptr3,iptr4,iptr5;
  for(int i=1; i<nx-1; i++)
  {
    iptr  =  k*nx + i;
    iptr1 = (k-1)*nx + i;
    iptr2 = (i-1)*CONST_NDIM;
    x2d[iptr] = x2d[iptr1] + xz[iptr2];
    z2d[iptr] = z2d[iptr1] + xz[iptr2+1];
  }

  // left
  // floating boundary
  if(bdry_itype[0] == 1)
  {
    iptr  = k*nx+0;       // (0,k)
    iptr1 = (k-1)*nx+0;   // (0,k-1)
    iptr2 = k*nx+1;       // (1,k)
    iptr3 = (k-1)*nx+1;   // (1,k-1)
    iptr4 = k*nx+2;       // (2,k)
    iptr5 = (k-1)*nx+2;   // (2,k-1)

    x2d[iptr] = x2d[iptr1] + (1+epsilon[0])*(x2d[iptr2]-x2d[iptr3])
               -epsilon[0]*(x2d[iptr4]-x2d[iptr5]);
    z2d[iptr] = z2d[iptr1] + (1+epsilon[0])*(z2d[iptr2]-z2d[iptr3])
               -epsilon[0]*(z2d[iptr4]-z2d[iptr5]);

  }
  else if(bdry_itype[0] == 2) // cartesian boundary
  {
    iptr  = k*nx+0;       // (0,k) 
    iptr1 = k*nx+1;       // (1,k) 
    iptr2 = (k-1)*nx+0;   // (0,k-1)
    iptr3 = (k-1)*nx+1;   // (1,k-1)

    x2d[iptr] = x2d[iptr2];
    z2d[iptr] = z2d[iptr2] + (z2d[iptr1] - z2d[iptr3]);
  }

  if(bdry_itype[1] == 1)
  {
    iptr  = k*nx+(nx-1);       // (nx-1,k)
    iptr1 = (k-1)*nx+(nx-1);   // (nx-1,k-1)
    iptr2 = k*nx+(nx-2);       // (nx-2,k)
    iptr3 = (k-1)*nx+(nx-2);   // (nx-2,k-1)
    iptr4 = k*nx+(nx-3);       // (nx-3,k)
    iptr5 = (k-1)*nx+(nx-3);   // (nx-3,k-1)

    x2d[iptr] = x2d[iptr1] + (1+epsilon[1])*(x2d[iptr2]-x2d[iptr3])
               -epsilon[1]*(x2d[iptr4]-x2d[iptr5]);
    z2d[iptr] = z2d[iptr1] + (1+epsilon[1])*(z2d[iptr2]-z2d[iptr3])
               -epsilon[1]*(z2d[iptr4]-z2d[iptr5]);
  }
  else if(bdry_itype[1] == 2) // cartesian boundary
  {
    iptr  = k*nx+(nx-1);   // (nx-1,k) 
    iptr1 = k*nx+(nx-2);   // (nx-2,k)
    iptr2 = (k-1)*nx+nx-1; // (nx-1,k-1)
    iptr3 = (k-1)*nx+nx-2; // (nx-2,k-1)

    x2d[iptr] = x2d[iptr2];
    z2d[iptr] = z2d[iptr2] + (z2d[iptr1] - z2d[iptr3]);
  }

  return 0;
}

// function: check increment, if too big or to small
// modofy increment
// not used
int
modify_incre(double *xz, float *x2d, float *z2d, int nx, int k)
{
  size_t  iptr1, iptr2, iptr3;
  float x0, x1, z0, z1;
  float x ,z;
  float arc_k0, arc_k1;
  if(k > 1)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (k-1)*nx + i;
      iptr2 = (k-2)*nx + i;
      iptr3 = (i-1)*CONST_NDIM;
      x0 = x2d[iptr1] - x2d[iptr2];
      z0 = z2d[iptr1] - z2d[iptr2];
      x1 = xz[iptr3];
      z1 = xz[iptr3+1];
      arc_k0 = sqrt(pow(x0,2) + pow(z0,2));
      arc_k1 = sqrt(pow(x1,2) + pow(z1,2));

      // normalize get unit vertor (x,z)
      x = x1/arc_k1;
      z = z1/arc_k1;
      if( arc_k1 < 0.5*arc_k0)
      {
        xz[iptr3]   = x * 0.5*arc_k0;
        xz[iptr3+1] = z * 0.5*arc_k0;
      }
      if( arc_k1 > 2.0*arc_k0)
      {
        xz[iptr3]   = x * 2.0*arc_k0;
        xz[iptr3+1] = z * 2.0*arc_k0;
      }
    }
  }
    
  return 0;
}
