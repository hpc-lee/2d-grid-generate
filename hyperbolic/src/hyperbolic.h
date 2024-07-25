#ifndef HYPERBOLIC_H
#define HYPERBOLIC_H

#include "gd_t.h"
#include "par_t.h"
/*************************************************
 * function prototype
 *************************************************/

int 
hyper_gene(gd_t *gdcurv, par_t *par);

int
cal_smooth_coef(float coef, float *x2d, float *z2d,
                int nx, int nz, int k, float *step, float *coef_e);

int 
cal_matrix(float *x2d, float *z2d, int nx, int k, float *step,
           double *a, double *b, double *c, double *d, float *area);

int
modify_smooth(float *x2d, float *z2d, int nx, int k, double *a,
              double *b, double *c, double *d, float *coef_e);

int
modify_bdry(int n, double *a, double *b, double *c, double *d,
            float epsilon, int bdry_itype);

int
assign_coords(double *xz, float *x2d, float *z2d, int nx, int k,
              float epsilon, int bdry_itype);

int
modify_incre(double *xz, float *x2d, float *z2d, int nx, int k);

 #endif
