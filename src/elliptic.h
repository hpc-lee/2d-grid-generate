#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include "gd_t.h"
#include "par_t.h"

int
higen_gene(gd_t *gdcurv, par_t *par);

int 
higen_P_SOR(float *x2d, float *z2d, int nx, int nz, float *P, float *Q,
            float d1, float d2, float coef, float err, int max_iter, float w);


int
set_src_higen(float *x2d, float *z2d, float *P, float *Q,
              float coef, int nx, int nz, float d1, float d2);

#endif
