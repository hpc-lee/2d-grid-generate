#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lib_mem.h"
#include "gd_t.h"
#include "constants.h"
#include "io_funcs.h"

int
init_gdcurv(gd_t *gdcurv, int nx, int nz)
{
  gdcurv->nx = nx;
  gdcurv->nz = nz;
  //2 dimension, x and z
  gdcurv->ncmp = CONST_NDIM; 
  gdcurv->siz_iz   = gdcurv->nx;
  gdcurv->siz_icmp = gdcurv->nx * gdcurv->nz;
  // malloc grid space 
  gdcurv->v3d = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp*gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v3d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }

  // position of each v3d
  size_t *cmp_pos = (size_t *) mem_calloc_1d_sizet(gdcurv->ncmp,
                                                         0,
                                                         "gd_curv_init");
  
  // name of each v3d
  char **cmp_name = (char **) mem_malloc_2l_char(gdcurv->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "gd_curv_init");
  
  // set value
  int icmp = 0;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","x");
  gdcurv->x2d = gdcurv->v3d + cmp_pos[icmp];

  icmp += 1;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","z");
  gdcurv->z2d = gdcurv->v3d + cmp_pos[icmp];
  
  // set pointer
  gdcurv->cmp_pos  = cmp_pos;
  gdcurv->cmp_name = cmp_name;

  return 0;
}

int
grid_init_set(gd_t *gdcurv, par_t *par)
{
  FILE *fp = NULL;
  char str[500];
  char *geometry_file = par->geometry_input_file;
  char *step_file = par->step_input_file;
  int dire_itype = par->dire_itype;
  
  int nx;
  int nz;
  int num_step;

  // open step file
  if ((fp = fopen(step_file,"r"))==NULL) {
     fprintf(stderr,"ERROR: fail to open step file=%s\n", step_file);
     fflush(stdout); exit(1);
  }
  // number of step
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&num_step);
  }

  nz = num_step+1; 

  gdcurv->step = (float *)mem_calloc_1d_float(
                          num_step, 0.0, "step length");

  for (int k=0; k<num_step; k++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f",gdcurv->step + k);
    }
  }
  // close step file and free local pointer
  fclose(fp);

  // open geometry file
  if ((fp = fopen(geometry_file,"r"))==NULL) {
     fprintf(stderr,"ERROR: fail to open geometry file=%s\n", geometry_file);
     fflush(stdout); exit(1);
  }
  // nx number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&nx);
  }

  init_gdcurv(gdcurv,nx,nz);
 
  size_t iptr;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  // z1 
  for (int i=0; i<nx; i++)
  {
    iptr = i;  // (i,0)
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",x2d+iptr,z2d+iptr);
    }
  }
  // z2 
  for (int i=0; i<nx; i++)
  {
    iptr = i + (nz-1)*nx;  // (i,nz-1)
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",x2d+iptr,z2d+iptr);
    }
  }
  // close geometry file and free local pointer
  fclose(fp);

  return 0;
}

// 2D array flip z direction.  nz-1->0 0->nz-1 i->(nz-1)-i 
int
flip_coord_z(gd_t *gdcurv)
{
  size_t iptr,iptr1;
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;

  float *tmp_coord_x = NULL;
  tmp_coord_x = (float *) malloc(nx*nz*sizeof(float));
  float *tmp_coord_z = NULL;
  tmp_coord_z = (float *) malloc(nx*nz*sizeof(float));
  // copy data
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++) 
    {
      iptr = k*nx + i;
      tmp_coord_x[iptr] = x2d[iptr];
      tmp_coord_z[iptr] = z2d[iptr];
    }
  }
  // flip coord
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++) 
    {
      iptr = k*nx + i;
      iptr1 = (nz-1-k)*nx + i;
      x2d[iptr] = tmp_coord_x[iptr1];
      z2d[iptr] = tmp_coord_z[iptr1];
    }
  }

  free(tmp_coord_x);
  free(tmp_coord_z);

  return 0;
}
