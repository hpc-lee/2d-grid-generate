#ifndef GD_CURV_H
#define GD_CURV_H

#include "par_t.h"
/*************************************************
 * structure
 *************************************************/

typedef struct {

  int nx;
  int nz;
  int ncmp;
  
  float *v3d; // pointer to var
  float *x2d; 
  float *z2d;
  
  float *step; // for hyperbolic

  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;

} gd_t;

/*************************************************
 * function prototype
 *************************************************/

int
init_gdcurv(gd_t *gdcurv, int nx, int nz);

int 
grid_init_set(gd_t *gdcurv, par_t *par);

int
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, float coef_x, float coef_z);

int
flip_coord_z(gd_t *gdcurv);

#endif
