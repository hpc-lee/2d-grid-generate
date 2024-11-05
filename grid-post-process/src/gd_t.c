#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lib_mem.h"
#include "lib_math.h"
#include "gd_t.h"
#include "algebra.h"
#include "constants.h"
#include "io_funcs.h"

int
init_gdcurv(gd_t *gdcurv, int nx, int nz)
{
  gdcurv->nx = nx;
  gdcurv->nz = nz;
  //2 dimension, x and z
  gdcurv->ncmp = CONST_NDIM; 
  gdcurv->siz_iz   = nx;
  gdcurv->siz_icmp = nx*nz;
  // malloc grid space 
  gdcurv->v3d = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp*gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v3d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }
  
  // set value
  int icmp = 0;
  gdcurv->x2d = gdcurv->v3d + icmp * gdcurv->siz_icmp;

  icmp += 1;
  gdcurv->z2d = gdcurv->v3d + icmp * gdcurv->siz_icmp;

  return 0;
}

int
grid_sample(gd_t *gdcurv_new, gd_t *gdcurv, int coef_x, int coef_z)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int nx_new =  (nx-1)*coef_x + 1;
  int nz_new =  (nz-1)*coef_z + 1;
  if(nx_new < nx || nz_new < nz)
  {
    fprintf(stdout,"only support up sample, \
                    nx_new(ny_new,nz_new) must >= nx(ny,nz)\n");
    exit(1);
  }

  init_gdcurv(gdcurv_new, nx_new, nz_new);
    
  sample_interp(gdcurv_new, gdcurv); 

  return 0;
}

int
gd_info_set(gd_t *gdcurv, par_t *par, int iprocx, int iprocz,
           int *global_index, int *count)
{
  int ierr = 0;

  int total_nx = gdcurv->nx;
  int total_nz = gdcurv->nz;

  int nprocx_out = par->num_of_procs_out[0];
  int nprocz_out = par->num_of_procs_out[1];

  int num_of_pml_x1 = par->num_of_pml_x1;
  int num_of_pml_x2 = par->num_of_pml_x2;
  int num_of_pml_z1 = par->num_of_pml_z1;
  int num_of_pml_z2 = par->num_of_pml_z2;


  int gni1, gnj1, gnk1;
  // determine ni
  int nx_et = total_nx;

  // double cfspml layer, load balance
  nx_et += num_of_pml_x1 + num_of_pml_x2;

  // partition into average plus left at last
  int nx_avg  = nx_et / nprocx_out;
  int nx_left = nx_et % nprocx_out;

  // nx_avg must > pml layers
  if(nx_avg<=num_of_pml_x1 || nx_avg<=num_of_pml_x2)
  {
    fprintf(stdout,"nx must large pml_layers\n");
    fflush(stdout);
    exit(1);
  }

  // default set to average value
  int ni = nx_avg;
  // subtract nlay for pml node
  if (iprocx == 0) {
    ni -= num_of_pml_x1;
  }
  if (iprocx == nprocx_out-1) {
    ni -= num_of_pml_x2;
  }

  // first nx_left node add one more point
  if (iprocx < nx_left) {
    ni++;
  }
  // global index
  if (iprocx==0) {
    gni1 = 0;
  } else {
    gni1 = iprocx * nx_avg - num_of_pml_x1;
  }
  if (nx_left != 0) {
    gni1 += (iprocx < nx_left) ? iprocx : nx_left;
  }

  // determine nk
  int nz_et = total_nz;

  // double cfspml layer, load balance
  nz_et += num_of_pml_z1 + num_of_pml_z2;

  int nz_avg  = nz_et / nprocz_out;
  int nz_left = nz_et % nprocz_out;

  // ny_avg must > pml layers
  if(nz_avg<=num_of_pml_z1 || nz_avg<=num_of_pml_z2)
  {
    fprintf(stdout,"nz must large pml_layers\n");
    fflush(stdout);
    exit(1);
  }

  // default set to average value
  int nk = nz_avg;
  // subtract nlay for pml node
  if (iprocz == 0) {
    nk -= num_of_pml_z1;
  }
  if (iprocz == nprocz_out-1) {
    nk -= num_of_pml_z2;
  }

  // first nz_left node add one more point
  if (iprocz < nz_left) {
    nk++;
  }
  // global index
  if (iprocz==0) {
    gnk1 = 0;
  } else {
    gnk1 = iprocz * nz_avg - num_of_pml_z1;
  }
  if (nz_left != 0) {
    gnk1 += (iprocz < nz_left) ? iprocz : nz_left;
  }

  global_index[0] = gni1;
  global_index[1] = gnk1;

  count[0] = ni;
  count[1] = nk;

  return ierr;
}

int
cal_min_dist(gd_t *gdcurv, int *indx_i,  int *indx_k, float *dL_min)
{
  float dL_min_local = 1e10;
  float dL_min_global = 1e10;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t siz_iz = gdcurv->siz_iz;

  for (int k = 1; k < nz-1; k++)
  {
    for (int i = 1; i < nx-1; i++)
    {
      size_t iptr = i + k * siz_iz;
      float p0[] = { x2d[iptr], z2d[iptr] };

      // min L to 4 adjacent planes
      for (int kk = -1; kk <=1; kk = kk+2)
      {
        for (int ii = -1; ii <= 1; ii = ii+2) 
        {
          float p1[] = { x2d[iptr-ii], z2d[iptr-ii] };
          float p2[] = { x2d[iptr-kk*siz_iz],
                         z2d[iptr-kk*siz_iz] };

          float L = dist_point2line(p0, p1, p2);

          if (dL_min_local > L) dL_min_local = L;
        }
      }

      if (dL_min_global > dL_min_local) 
      {
        dL_min_global = dL_min_local;
        *dL_min = dL_min_global;
        *indx_i = i;
        *indx_k = k;
      }
    }
  }

  return 0;
}
