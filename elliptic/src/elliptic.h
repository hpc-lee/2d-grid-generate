#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include "gd_t.h"
#include "par_t.h"
#include "mympi_t.h"

typedef struct {

  float *P_x1;  
  float *Q_x1;  
  float *P_x2;  
  float *Q_x2;  
  float *P_z1;  
  float *Q_z1;  
  float *P_z2;  
  float *Q_z2;  
  float *P_x1_loc;  
  float *Q_x1_loc;  
  float *P_x2_loc;  
  float *Q_x2_loc;  
  float *P_z1_loc;  
  float *Q_z1_loc;  
  float *P_z2_loc;  
  float *Q_z2_loc;  

  float *P;
  float *Q;

} src_t;

int
init_src(src_t *src, gd_t *gdcurv);

int
diri_gene(gd_t *gdcurv, par_t *par, mympi_t *mympi); 

int
set_src_diri(float *x2d, float *z2d, gd_t *gdcurv,
             src_t *src, float *p_x, float *p_z, int *flag_orth,
             float *g11_x, float *g22_z, mympi_t *mympi);

int
ghost_cal(float *x2d, float *z2d, int nx, int nz,
          float *p_x, float *p_z, float *g11_x,
          float *g22_z, int *neighid);

int
higen_gene(gd_t *gdcurv, par_t *par, mympi_t *mympi);

int
set_src_higen(float *x2d, float *z2d, gd_t *gdcurv, 
              src_t *src, float *dx1, float *dx2,
              float *dz1, float *dz2, mympi_t *mympi);
              
int
dist_cal(gd_t *gdcurv, float *dx1, float *dx2, float *dz1, 
         float *dz2, int *neighid);

int interp_inner_source(src_t *src, gd_t *gdcurv, float coef);

#endif
