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
predict_point(float *x2d, float *z2d, float *step, int nx, int nz, int k,
              float coef, float *x_pre, float *z_pre);

int
update_point(float *x2d, float *z2d, float *thomas, int nx, int k, 
             float *x_pre, float *z_pre);

int 
assign_bdry_coords(float *x2d, float *z2d, int nx, int k);

#endif
