#ifndef GD_CURV_H
#define GD_CURV_H

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

  size_t *cmp_pos;
  char  **cmp_name;

} gd_t;

/*************************************************
 * function prototype
 *************************************************/

int 
init_gdcurv(gd_t *gdcurv, int nx, int nz);

int 
grid_init_set(gd_t *gdcurv, char *input_file);

int
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, float coef_x, float coef_z);

int
gd_curv_coord_export(gd_t *gdcurv, char *output_dir);

int 
check_bdry(float *x1, float *x2, float *z1, float *z2, int nx, int nz);

int
io_get_nextline(FILE *fp, char *str, int length);

#endif
