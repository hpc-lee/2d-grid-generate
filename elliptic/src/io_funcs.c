#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h" 

#include "lib_mem.h"
#include "constants.h"
#include "io_funcs.h"

int
init_io_quality(io_quality_t *io_quality, gd_t *gdcurv)
{
  io_quality->nx = gdcurv->nx;
  io_quality->nz = gdcurv->nz;
  
  // malloc quality space 
  io_quality->var = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp, 0.0, "quality_init");
  if (io_quality->var == NULL) {
      fprintf(stderr,"Error: failed to alloc quality vars\n");
      fflush(stderr);
  }
   
  return 0;
}

int
gd_curv_coord_export(gd_t *gdcurv, mympi_t *mympi)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int ni = gdcurv->ni;
  int nk = gdcurv->nk;
  int gni1 = gdcurv->gni1;
  int gnk1 = gdcurv->gnk1;
  int gni2 = gdcurv->gni2;
  int gnk2 = gdcurv->gnk2;
  size_t iptr, iptr1;
  size_t siz_iz = gdcurv->siz_iz;

  // first modify ni nk ... , we need export bdry ghost point
  if(mympi->neighid[0] == MPI_PROC_NULL)
  {
    ni = ni+1;
    ni1 = ni1-1;
  }
  if(mympi->neighid[1] == MPI_PROC_NULL)
  {
    ni = ni+1;
    ni2 = ni2+1;
  }
  if(mympi->neighid[2] == MPI_PROC_NULL)
  {
    nk = nk+1;
    nk1 = nk1-1;
  }
  if(mympi->neighid[3] == MPI_PROC_NULL)
  {
    nk = nk+1;
    nk2 = nk2+1;
  }
  if(mympi->topoid[0] != 0)
  {
    gni1 = gni1 + 1; 
  }
  if(mympi->topoid[1] != 0)
  {
    gnk1 = gnk1 + 1; 
  }
  

  float *x2d =  gdcurv->x2d;
  float *z2d =  gdcurv->z2d;
  float *coord_x = (float *) malloc(sizeof(float)*ni*nk);
  float *coord_z = (float *) malloc(sizeof(float)*ni*nk);

  for(int k=nk1; k<=nk2; k++) {
    for(int i=ni1; i<=ni2; i++)
    {
      iptr = k*siz_iz + i;
      iptr1 = (k-nk1)*ni + (i-ni1);
      coord_x[iptr1] = x2d[iptr];
      coord_z[iptr1] = z2d[iptr];
    }
  }

  char *output_dir = gdcurv->output_dir;
  char *fname = gdcurv->fname_part;
  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/coord_%s.nc", output_dir,fname);
  
  // read in nc
  int ncid;
  int xid, zid;
  int dimid[2];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", ni, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nk, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  ierr = nc_def_var(ncid, "x", NC_FLOAT, 2, dimid, &xid);
  handle_nc_err(ierr);
  ierr = nc_def_var(ncid, "z", NC_FLOAT, 2, dimid, &zid);
  handle_nc_err(ierr);

  // attribute: index in output snapshot, index w ghost in thread
  int g_start[] = { gni1, gnk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                   NC_INT,2,g_start);

  int l_count[] = { ni, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,2,l_count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  ierr = nc_put_var_float(ncid, xid, coord_x);  handle_nc_err(ierr);
  ierr = nc_put_var_float(ncid, zid, coord_z);  handle_nc_err(ierr);
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  free(coord_x);
  free(coord_z);

  return 0;
}

int
quality_export(io_quality_t *io_quality, gd_t *gdcurv, mympi_t *mympi, char *var_name)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int ni = gdcurv->ni;
  int nk = gdcurv->nk;
  int gni1 = gdcurv->gni1;
  int gnk1 = gdcurv->gnk1;
  int gni2 = gdcurv->gni2;
  int gnk2 = gdcurv->gnk2;
  size_t iptr, iptr1;
  size_t siz_iz = gdcurv->siz_iz;

  // first modify ni nk ... , we need export bdry ghost point
  if(mympi->neighid[0] == MPI_PROC_NULL)
  {
    ni = ni+1;
    ni1 = ni1-1;
  }
  if(mympi->neighid[1] == MPI_PROC_NULL)
  {
    ni = ni+1;
    ni2 = ni2+1;
  }
  if(mympi->neighid[2] == MPI_PROC_NULL)
  {
    nk = nk+1;
    nk1 = nk1-1;
  }
  if(mympi->neighid[3] == MPI_PROC_NULL)
  {
    nk = nk+1;
    nk2 = nk2+1;
  }
  if(mympi->topoid[0] != 0)
  {
    gni1 = gni1 + 1; 
  }
  if(mympi->topoid[1] != 0)
  {
    gnk1 = gnk1 + 1; 
  }

  float *var =  io_quality->var;
  float *var_out = (float *) malloc(sizeof(float)*ni*nk);

  for(int k=nk1; k<=nk2; k++) {
    for(int i=ni1; i<=ni2; i++)
    {
      iptr = k*siz_iz + i;
      iptr1 = (k-nk1)*ni + (i-ni1);
      var_out[iptr1] = var[iptr];
    }
  }

  char *output_dir = gdcurv->output_dir;
  char *fname = gdcurv->fname_part;
  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/%s_%s.nc", output_dir,var_name,fname);
  
  // read in nc
  int ncid;
  int varid;
  int dimid[2];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", ni, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nk, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  ierr = nc_def_var(ncid, var_name, NC_FLOAT, 2, dimid, &varid);
  handle_nc_err(ierr);

  // attribute: index in output snapshot, index w ghost in thread
  int g_start[] = { gni1, gnk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                   NC_INT,2,g_start);

  int l_count[] = { ni, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,2,l_count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  ierr = nc_put_var_float(ncid, varid, var_out);  handle_nc_err(ierr);
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  free(var_out); 

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
