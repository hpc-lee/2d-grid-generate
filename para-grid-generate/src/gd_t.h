#ifndef GD_CURV_H
#define GD_CURV_H

#include "par_t.h"
#include "mympi_t.h"
/*************************************************
 * structure
 *************************************************/

typedef struct {

  int ni, nk;
  int nx, nz;
  int ni1, ni2;
  int nk1, nk2;
  // global index
  int gni1, gnk1; // global index, do not accout ghost point
  int gni2, gnk2; // global index

  int total_nx;
  int total_nz;

  int ncmp;
  
  float *v3d; // pointer to var
  float *x2d; 
  float *z2d;

  // temp pointer to var, iter need
  float *v3d_tmp; 
  float *x2d_tmp;   
  float *z2d_tmp;

  size_t siz_iz;
  size_t siz_icmp;

  char fname_part[CONST_MAX_STRLEN];
  char output_dir[CONST_MAX_STRLEN];

} gd_t;

typedef struct {

  float *var;  
  float *x1;
  float *x2;
  float *z1;
  float *z2;
  int number_of_grid_points_x;
  int number_of_grid_points_z;

} bdry_t;
/*************************************************
 * function prototype
 *************************************************/
int
gd_info_set(gd_t *gdcurv, mympi_t *mympi,
            int number_of_grid_points_x,
            int number_of_grid_points_z,
            int verbose);

int
set_output_dir(gd_t *gdcurv,
               mympi_t *mympi,
               char *output_dir,
               const int verbose);

int 
init_gdcurv(gd_t *gdcurv);

int
init_bdry(bdry_t *bdry, par_t *par);

int
read_bdry(int myid, bdry_t *bdry, char *geometry_file);

int
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, float coef_x, float coef_z);

int 
check_bdry(float *x1, float *x2, float *z1, float *z2, int nx, int nz);

int
gd_info_print(gd_t *gdcurv, int myid);

int
gd_curv_coord_exchange(gd_t *gdcurv, int *neighid, MPI_Comm topocomm);

int
grid_mesg_init(mympi_t *mympi, gd_t *gdcurv);

int
grid_pack_mesg(mympi_t *mympi, gd_t *gdcurv, float *x2d, float *z2d);

int
grid_unpack_mesg(mympi_t *mympi, gd_t *gdcurv, float *x2d, float *z2d);

#endif
