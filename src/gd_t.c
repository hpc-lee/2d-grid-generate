#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h" 

#include "gd_t.h"
#include "constants.h"
#include "fdlib_mem.h"
#include "algebra.h"

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
  gdcurv->v3d = (float *)fdlib_mem_calloc_1d_float(
                  gdcurv->siz_icmp*gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v3d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }

  // position of each v3d
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(gdcurv->ncmp,
                                                         0,
                                                         "gd_curv_init");
  
  // name of each v3d
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(gdcurv->ncmp,
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
grid_init_set(gd_t *gdcurv, char *input_file)
{
  FILE *fp = NULL;
  char str[500];
  
  int nx;
  int nz;
  // open geometry file
  if ((fp = fopen(input_file,"r"))==NULL) {
     fprintf(stderr,"ERROR: fail to open geometry file=%s\n", input_file);
     fflush(stdout); exit(1);
  }
  // nx number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&nx);
  }
  // nz number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&nz);
  }
  
  init_gdcurv(gdcurv,nx,nz);

  // malloc 4 bdry space, read bdry coods
  float *x1;
  float *x2;
  float *z1;
  float *z2;
  x1 = (float *)fdlib_mem_calloc_1d_float(
            nz*gdcurv->ncmp, 0.0, "bdry_coords");
  x2 = (float *)fdlib_mem_calloc_1d_float(
            nz*gdcurv->ncmp, 0.0, "bdry_coords");
  z1 = (float *)fdlib_mem_calloc_1d_float(
            nx*gdcurv->ncmp, 0.0, "bdry_coords");
  z2 = (float *)fdlib_mem_calloc_1d_float(
            nx*gdcurv->ncmp, 0.0, "bdry_coords");
  // x1 
  for (int k=0; k<nz; k++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",x1+k,(x1+nz)+k);
    }
  }
  // x2 
  for (int k=0; k<nz; k++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",x2+k,(x2+nz)+k);
    }
  }
  // z1 
  for (int i=0; i<nx; i++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",z1+i,(z1+nx)+i);
    }
  }
  // z2 
  for (int i=0; i<nx; i++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",z2+i,(z2+nx)+i);
    }
  }
  // close file and free local pointer
  fclose(fp); 

  check_bdry(x1,x2,z1,z2,nx,nz);

  // read boundry grid coords 
  size_t iptr;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  // x1 i=0
  for (int k=0; k<nz; k++)
  {
    iptr = k*nx;
    x2d[iptr] = x1[k];
    z2d[iptr] = x1[k+nz];
  }
  // x2 i=nx-1
  for (int k=0; k<nz; k++)
  {
    iptr = k*nx + (nx-1);
    x2d[iptr] = x2[k];
    z2d[iptr] = x2[k+nz];
  }
  // z1 k=0
  for (int i=0; i<nx; i++)
  {
    iptr = i;
    x2d[iptr] = z1[i];
    z2d[iptr] = z1[i+nx];
  }
  // z2 k=nz-1
  for (int i=0; i<nx; i++)
  {
    iptr = (nz-1)*nx + i;
    x2d[iptr] = z2[i];
    z2d[iptr] = z2[i+nx];
  }

  free(x1);
  free(x2);
  free(z1);
  free(z2);

  return 0;
}

int
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, float coef_x, float coef_z)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int nx_new = (int) (nx*coef_x);
  int nz_new = (int) (nz*coef_z);

  init_gdcurv(gdcurv_new, nx_new, nz_new);
    
  sample_interp(gdcurv_new, gdcurv); 

  return 0;
}

int
gd_curv_coord_export(gd_t *gdcurv, char *output_dir)
{
  size_t *restrict c3d_pos   = gdcurv->cmp_pos;
  char  **restrict c3d_name  = gdcurv->cmp_name;
  int  nx = gdcurv->nx;
  int  nz = gdcurv->nz;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/coord.nc", output_dir);
  
  // read in nc
  int ncid;
  int varid[gdcurv->ncmp];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    ierr = nc_def_var(ncid, gdcurv->cmp_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
    handle_nc_err(ierr);
  }

  // attribute: index in output snapshot, index w ghost in thread
  int l_count[] = { nx, nz };
  nc_put_att_int(ncid,NC_GLOBAL,"number_of_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    float *ptr = gdcurv->v3d + gdcurv->cmp_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
    handle_nc_err(ierr);
  }
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  return 0;
}

int 
check_bdry(float *x1, float *x2, float *z1, float *z2, int nx, int nz)
{ 
  int ierr = 0;
  float p1_x, p1_z, p2_x, p2_z;
  //  check four corner points
  //  (0,0)
  p1_x = x1[0];  
  p1_z = x1[0+nz];  
  p2_x = z1[0];
  p2_z = z1[0+nx];
  if(p1_x == p2_x && p1_z == p2_z) {
    ierr = 0;
  } else {
    ierr =1;
    fprintf(stdout, "point (0,0) error, please check x1 and z1 boundary");
    exit(1);
  }

  //  (nx-1,0)
  p1_x = x2[0];  
  p1_z = x2[0+nz];  
  p2_x = z1[nx-1];
  p2_z = z1[nx-1+nx];
  if(p1_x == p2_x && p1_z == p2_z) {
    ierr = 0;
  } else {
    ierr =1;
    fprintf(stdout, "point (nx-1,0) error, please check x2 and z1 boundary");
    exit(1);
  }

  //  (nx-1,nz-1)
  p1_x = x2[nz-1];  
  p1_z = x2[nz-1+nz];  
  p2_x = z2[nx-1];
  p2_z = z2[nx-1+nx];
  if(p1_x == p2_x && p1_z == p2_z) {
    ierr = 0;
  } else {
    ierr =1;
    fprintf(stdout, "point (nx-1,nz-1) error, please check x2 and z2 boundary");
    exit(1);
  }

  //  (0,nz-1)
  p1_x = x1[nz-1];  
  p1_z = x1[nz-1+nz];  
  p2_x = z2[0];
  p2_z = z2[0+nx];
  if(p1_x == p2_x && p1_z == p2_z) {
    ierr = 0;
  } else {
    ierr =1;
    fprintf(stdout, "point (0,nz-1) error, please check  x1 and z2  boundary");
    exit(1);
  }

  return 0;
}


int
io_get_nextline(FILE *fp, char *str, int length)
{
  int ierr = 0;

  do
  {
    if (fgets(str, length, fp) == NULL)
    {
      ierr = 1;
      return ierr;
    }
  } while (str[0] == '#' || str[0] == '\n');

  // remove newline char
  int len = strlen(str);

  if (len > 0 && str[len-1] == '\n') {
    str[len-1] = '\0';
  }

  // for debug:
  //fprintf(stdout," --return: %s\n", str);

  return ierr;
}






