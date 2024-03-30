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
predict_point(float *x2d, float *z2d, int nx, int nz, int k, int o2i, 
              float coef, float *step_len);

int
update_point(float *x2d, float *z2d, float *var_th, int nx, int k);

int 
assign_bdry_coords(float *x2d, float *z2d, int nx, int k);

#endif
