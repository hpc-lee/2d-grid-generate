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
read_import_coord(gd_t *gdcurv, par_t *par)
{
  int ierr;
  int total_nx = par->number_of_grid_points_x;
  int total_nz = par->number_of_grid_points_z;

  int nprocx_in = par->number_of_mpiprocs_x_in;
  int nprocz_in = par->number_of_mpiprocs_z_in;

  int gni1, gnk1;
  int gni, gnk;
  int ni, nk;
  size_t iptr, iptr1;

  init_gdcurv(gdcurv,total_nx,total_nz);
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;

  float *coord_x = (float *) malloc(sizeof(float)*total_nx*total_nz);
  float *coord_z = (float *) malloc(sizeof(float)*total_nx*total_nz);

  int ncid;
  int xid, zid;
  int global_index[2];
  int count_points[2];
  size_t start[2] = {0, 0};
  size_t count[2];
  char att_global[CONST_MAX_STRLEN] = "global_index_of_first_physical_points";
  char att_count[CONST_MAX_STRLEN] = "count_of_physical_points";

  char fname_coords[CONST_MAX_STRLEN];
  char in_file[CONST_MAX_STRLEN];
  // read coord nc file
  for(int kk=0; kk<nprocz_in; kk++) {
    for(int ii=0; ii<nprocx_in; ii++)
    {
      sprintf(fname_coords,"px%d_pz%d",ii,kk);
      sprintf(in_file,"%s/coord_%s.nc",par->import_dir,fname_coords);

      ierr = nc_open(in_file, NC_NOWRITE, &ncid); handle_nc_err(ierr);

      ierr = nc_get_att_int(ncid,NC_GLOBAL,att_global,global_index);
      ierr = nc_get_att_int(ncid,NC_GLOBAL,att_count,count_points);

      gni1 = global_index[0];
      gnk1 = global_index[1];

      ni = count_points[0];
      nk = count_points[1];

      count[0] = nk;
      count[1] = ni;

      //read vars
      ierr = nc_inq_varid(ncid, "x", &xid);  handle_nc_err(ierr);
      ierr = nc_inq_varid(ncid, "z", &zid);  handle_nc_err(ierr);

      ierr = nc_get_vara_float(ncid, xid, start, count, coord_x);  handle_nc_err(ierr);
      ierr = nc_get_vara_float(ncid, zid, start, count, coord_z);  handle_nc_err(ierr);

      for(int k=0; k<nk; k++) {
        for(int i=0; i<ni; i++)
        {
          gni = gni1 + i;
          gnk = gnk1 + k;
          iptr = gnk*total_nx + gni;

          iptr1 = k*ni + i;

          x2d[iptr] = coord_x[iptr1];
          z2d[iptr] = coord_z[iptr1];
        }
      }

      //close file
      ierr = nc_close(ncid);  handle_nc_err(ierr);
    }
  }

  free(coord_x);
  free(coord_z);

  return 0;
}

int
gd_curv_coord_export(gd_t *gdcurv, par_t *par)
{
  int nprocx_out = par->number_of_mpiprocs_x_out;
  int nprocz_out = par->number_of_mpiprocs_z_out;

  int total_nx = gdcurv->nx;
  int total_nz = gdcurv->nz;

  int gni1, gnk1;
  int gni, gnk;
  int ni, nk;
  size_t iptr, iptr1;

  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *coord_x = (float *) malloc(sizeof(float)*total_nx*total_nz);
  float *coord_z = (float *) malloc(sizeof(float)*total_nx*total_nz);

  char fname_coords[CONST_MAX_STRLEN];
  char ou_file[CONST_MAX_STRLEN];

  // read in nc
  int ierr;
  int ncid;
  int xid, zid;
  int dimid[2];
  int g_start[2];
  int count[2];

  // read coord nc file
  for(int kk=0; kk<nprocz_out; kk++) {
    for(int ii=0; ii<nprocx_out; ii++)
    {
      sprintf(fname_coords,"px%d_pz%d",ii,kk);
      sprintf(ou_file,"%s/coord_%s.nc",par->export_dir,fname_coords);

      gd_info_set(gdcurv, par, ii, kk, g_start, count);
      gni1 = g_start[0];
      gnk1 = g_start[1];
      ni = count[0];
      nk = count[1];

      for(int k=0; k<nk; k++) {
        for(int i=0; i<ni; i++)
        {
          gni = gni1 + i;
          gnk = gnk1 + k;
          iptr = gnk*total_nx + gni;

          iptr1 = k*ni + i;

          coord_x[iptr1] = x2d[iptr];
          coord_z[iptr1] = z2d[iptr];
        }
      }

      ierr = nc_create(ou_file, NC_CLOBBER, &ncid); handle_nc_err(ierr);

      // define dimension
      ierr = nc_def_dim(ncid, "i", ni, &dimid[1]); handle_nc_err(ierr);
      ierr = nc_def_dim(ncid, "k", nk, &dimid[0]); handle_nc_err(ierr);

      // define vars
      ierr = nc_def_var(ncid, "x", NC_FLOAT, 2, dimid, &xid); handle_nc_err(ierr);
      ierr = nc_def_var(ncid, "z", NC_FLOAT, 2, dimid, &zid); handle_nc_err(ierr);

      // attribute: index in output snapshot, index w ghost in thread
      nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
          NC_INT,2,g_start);

      nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
          NC_INT,2,count);

      // end def
      ierr = nc_enddef(ncid);
      handle_nc_err(ierr);

      // add vars
      ierr = nc_put_var_float(ncid, xid, coord_x);  handle_nc_err(ierr);
      ierr = nc_put_var_float(ncid, zid, coord_z);  handle_nc_err(ierr);

      // close file
      ierr = nc_close(ncid);
      handle_nc_err(ierr);
    }
  }

  free(coord_x);
  free(coord_z);

  return 0;
}

int
quality_export(io_quality_t *io_quality, gd_t *gdcurv, par_t *par, char *var_name)
{
  int nprocx_out = par->number_of_mpiprocs_x_out;
  int nprocz_out = par->number_of_mpiprocs_z_out;

  int total_nx = io_quality->nx;
  int total_nz = io_quality->nz;

  int gni1, gnk1;
  int gni, gnk;
  int ni, nk;
  size_t iptr, iptr1;

  float *var = io_quality->var;
  float *var_out = (float *) malloc(sizeof(float)*total_nx*total_nz);

  char fname_coords[CONST_MAX_STRLEN];
  char ou_file[CONST_MAX_STRLEN];

  // read in nc
  int ierr;
  int ncid;
  int varid;
  int dimid[2];
  int g_start[2];
  int count[2];

  // read coord nc file
  for(int kk=0; kk<nprocz_out; kk++) {
    for(int ii=0; ii<nprocx_out; ii++)
    {
      sprintf(fname_coords,"px%d_pz%d",ii,kk);
      sprintf(ou_file,"%s/%s_%s.nc",par->export_dir,var_name,fname_coords);

      gd_info_set(gdcurv, par, ii, kk, g_start, count);
      gni1 = g_start[0];
      gnk1 = g_start[1];
      ni = count[0];
      nk = count[1];

      for(int k=0; k<nk; k++) {
        for(int i=0; i<ni; i++)
        {
          gni = gni1 + i;
          gnk = gnk1 + k;
          iptr = gnk*total_nx + gni;

          iptr1 = k*ni + i;

          var_out[iptr1] = var[iptr];
        }
      }

      ierr = nc_create(ou_file, NC_CLOBBER, &ncid); handle_nc_err(ierr);

      // define dimension
      ierr = nc_def_dim(ncid, "i", ni, &dimid[1]); handle_nc_err(ierr);
      ierr = nc_def_dim(ncid, "k", nk, &dimid[0]); handle_nc_err(ierr);

      // define vars
      ierr = nc_def_var(ncid, var_name, NC_FLOAT, 2, dimid, &varid); handle_nc_err(ierr);

      // attribute: index in output snapshot, index w ghost in thread
      nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
          NC_INT,2,g_start);

      nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
          NC_INT,2,count);

      // end def
      ierr = nc_enddef(ncid);
      handle_nc_err(ierr);

      // add vars
      ierr = nc_put_var_float(ncid, varid, var_out);  handle_nc_err(ierr);

      // close file
      ierr = nc_close(ncid);
      handle_nc_err(ierr);
    }
  }

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
