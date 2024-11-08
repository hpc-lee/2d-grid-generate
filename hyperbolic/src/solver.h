#ifndef SOLVER_H
#define SOLVER_H

/*************************************************
 * function prototype
 *************************************************/
int
thomas(int n, float *a, float *b, float *c, float *d_x,
       float *d_z, float *u_x, float *u_z);

int
thomas_block(int n, float *a, float *b, float *c, float *d,
             float *xz, float *D, float *y);

#endif
