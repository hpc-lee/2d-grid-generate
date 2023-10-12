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

  size_t siz_iz;
  size_t siz_icmp;


} gd_t;

/*************************************************
 * function prototype
 *************************************************/

int 
init_gdcurv(gd_t *gdcurv, int nx, int nz);

int
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, int coef_x, int coef_z);

int
gd_info_set(gd_t *gdcurv, par_t *par, int iprocx, int iprocz,
           int *global_index, int *count);

#endif
