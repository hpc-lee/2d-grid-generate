#ifndef PARABOLIC_H
#define PARABOLIC_H

#include "gd_t.h"
#include "par_t.h"
/*************************************************
 * function prototype
 *************************************************/

int 
para_gene(gd_t *gdcurv, par_t *par);

int 
predict_point(float *x2d, float *z2d, int nx, int nz, int k, int t2b, 
              float coef, float *step_len, float *x_pre, float *z_pre);

int
update_point(float *x2d, float *z2d, float *var_th, int nx, int k,
             float *x_pre, float *z_pre);

int 
assign_bdry_coords(float *x2d, float *z2d, int nx, int k);

int
flip_step_z(float *step, int nz);

#endif
