#ifndef SOLVER_H
#define SOLVER_H

/*************************************************
 * function prototype
 *************************************************/

int
update_SOR(float *x2d, float *z2d, float *x2d_tmp, float *z2d_tmp,
           int nx, int nz, float *P, float *Q, float w);

#endif
