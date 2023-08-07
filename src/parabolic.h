#ifndef PARABOLIC_H
#define PARABOLIC_H

#include "gd_t.h"
/*************************************************
 * function prototype
 *************************************************/

int 
para_gene(gd_t *gdcurv,float coef, int o2i);

int 
predict_point(float *x2d, float *z2d, int nx, int nz, int k, int o2i, float coef);

int
update_point(float *x2d, float *z2d, float *thomas, int nx, int nz, int k);

int
bdry_effct(float *x2d, float *z2d, int nx, int nz, int k);

#endif
