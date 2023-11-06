#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lib_mem.h"
#include "gd_t.h"
#include "constants.h"
#include "algebra.h"
#include "io_funcs.h"

int
init_gdcurv(gd_t *gdcurv)
{
  //2 dimension, x and z
  gdcurv->ncmp = CONST_NDIM; 
  // malloc grid space 
  gdcurv->v3d = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp*gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v3d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }
  
  // set value
  int icmp;
  icmp = 0;
  gdcurv->x2d = gdcurv->v3d + icmp * gdcurv->siz_icmp;

  icmp = 1;
  gdcurv->z2d = gdcurv->v3d + icmp * gdcurv->siz_icmp;

  // malloc temp grid space 
  gdcurv->v3d_tmp = (float *)mem_calloc_1d_float(
                  gdcurv->siz_icmp*gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v3d_tmp == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }
  
  // set value
  icmp = 0;
  gdcurv->x2d_tmp = gdcurv->v3d_tmp + icmp * gdcurv->siz_icmp;

  icmp = 1;
  gdcurv->z2d_tmp = gdcurv->v3d_tmp + icmp * gdcurv->siz_icmp;

  return 0;
}

int
init_bdry(bdry_t *bdry, par_t *par)
{
  bdry->number_of_grid_points_x = par->number_of_grid_points_x;
  bdry->number_of_grid_points_z = par->number_of_grid_points_z;
  int size = 2*(bdry->number_of_grid_points_x + bdry->number_of_grid_points_z);
  // malloc grid space. x and z coord 
  bdry->var = (float *)mem_calloc_1d_float(size*2, 0.0, "gd_curv_init");
  if (bdry->var == NULL) {
      fprintf(stderr,"Error: failed to alloc bdry vars\n");
      fflush(stderr);
  }

  bdry->x1 = bdry->var;
  bdry->x2 = bdry->x1+2*bdry->number_of_grid_points_z;
  bdry->z1 = bdry->x2+2*bdry->number_of_grid_points_z;
  bdry->z2 = bdry->z1+2*bdry->number_of_grid_points_x;

  return 0;
}

int
read_bdry(int myid, bdry_t *bdry, char *geometry_file)
{
  FILE *fp = NULL;
  char str[500];

  float *x1 = bdry->x1;
  float *x2 = bdry->x2;
  float *z1 = bdry->z1;
  float *z2 = bdry->z2;
  
  int num_point_x;
  int num_point_z;
  // open geometry file
  if ((fp = fopen(geometry_file,"r"))==NULL) {
     fprintf(stderr,"ERROR: fail to open geometry file=%s\n", geometry_file);
     fflush(stdout); exit(1);
  }
  // nx number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&num_point_x);
  }
  // nz number
  if (!io_get_nextline(fp,str,500)) {
    sscanf(str,"%d",&num_point_z);
  }
  // x1 
  for (int k=0; k<num_point_z; k++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",x1+k,(x1+num_point_z)+k);
    }
  }
  // x2 
  for (int k=0; k<num_point_z; k++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",x2+k,(x2+num_point_z)+k);
    }
  }
  // z1 
  for (int i=0; i<num_point_x; i++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",z1+i,(z1+num_point_x)+i);
    }
  }
  // z2 
  for (int i=0; i<num_point_x; i++)
  {
    if (!io_get_nextline(fp,str,500)) {
      sscanf(str,"%f %f",z2+i,(z2+num_point_x)+i);
    }
  }
  // close file and free local pointer
  fclose(fp); 
  if(myid == 0)
  {
    check_bdry(x1,x2,z1,z2,num_point_x,num_point_z);
  }

  return 0;
}

int 
check_bdry(float *x1, float *x2, float *z1, float *z2, int nx, int nz)
{ 
  int ierr = 0;
  float p1_x, p1_z, p2_x, p2_z;
  float dif_x, dif_z, dif;
  //  check four corner points
  //  (0,0)
  p1_x = x1[0];  
  p1_z = x1[0+nz];  
  p2_x = z1[0];
  p2_z = z1[0+nx];
  dif_x = p1_x-p2_x;
  dif_z = p1_z-p2_z;
  dif = fabs(dif_x) + fabs(dif_z);
  if(dif > 1e-6) {
    ierr =1;
    fprintf(stdout, "point (0,0) error, please check x1 and z1 boundary");
    exit(1);
  } else {
    ierr = 0;
  }

  //  (nx-1,0)
  p1_x = x2[0];  
  p1_z = x2[0+nz];  
  p2_x = z1[nx-1];
  p2_z = z1[nx-1+nx];
  dif_x = p1_x-p2_x;
  dif_z = p1_z-p2_z;
  dif = fabs(dif_x) + fabs(dif_z);
  if(dif > 1e-6) {
    ierr =1;
    fprintf(stdout, "point (nx-1,0) error, please check x2 and z1 boundary");
    exit(1);
  } else {
    ierr = 0;
  }

  //  (nx-1,nz-1)
  p1_x = x2[nz-1];  
  p1_z = x2[nz-1+nz];  
  p2_x = z2[nx-1];
  p2_z = z2[nx-1+nx];
  dif_x = p1_x-p2_x;
  dif_z = p1_z-p2_z;
  dif = fabs(dif_x) + fabs(dif_z);
  if(dif > 1e-6) {
    ierr =1;
    fprintf(stdout, "point (nx-1,nz-1) error, please check x2 and z2 boundary");
    exit(1);
  } else {
    ierr = 0;
  }

  //  (0,nz-1)
  p1_x = x1[nz-1];  
  p1_z = x1[nz-1+nz];  
  p2_x = z2[0];
  p2_z = z2[0+nx];
  dif_x = p1_x-p2_x;
  dif_z = p1_z-p2_z;
  dif = fabs(dif_x) + fabs(dif_z);
  if(dif > 1e-6) {
    ierr =1;
    fprintf(stdout, "point (0,nz-1) error, please check  x1 and z2  boundary");
    exit(1);
  } else {
    ierr = 0;
  }

  return 0;
}

int
gd_info_set(gd_t *gdcurv,
            mympi_t *mympi,
            int number_of_grid_points_x,
            int number_of_grid_points_z,
            int verbose)
{
  int ierr = 0;

  gdcurv->total_nx = number_of_grid_points_x;
  gdcurv->total_nz = number_of_grid_points_z;

  // 3 point center difference  
  // ghost number is 1

  // determine ni
  // bbry point to set ghost
  int nx_et = number_of_grid_points_x-2;

  // partition into average plus left at last
  int nx_avg  = nx_et / mympi->nprocx;
  int nx_left = nx_et % mympi->nprocx;

  // default set to average value
  int ni = nx_avg;
  // first nx_left node add one more point
  if (mympi->topoid[0] < nx_left) {
    ni++;
  }
  // global index
  if (mympi->topoid[0]==0) {
    gdcurv->gni1 = 0;
  } else {
    gdcurv->gni1 = mympi->topoid[0] * nx_avg;
  }
  if (nx_left != 0) {
    gdcurv->gni1 += (mympi->topoid[0] < nx_left)? mympi->topoid[0] : nx_left;
  }

  // determine nk
  int nz_et = number_of_grid_points_z-2;

  int nz_avg  = nz_et / mympi->nprocz;
  int nz_left = nz_et % mympi->nprocz;

  int nk = nz_avg;
  // not equal divided points given to first nz_left procs
  if (mympi->topoid[1] < nz_left) {
    nk++;
  }
  // global index
  if (mympi->topoid[1]==0) {
    gdcurv->gnk1 = 0;
  } else {
    gdcurv->gnk1 = mympi->topoid[1] * nz_avg;
  }
  if (nz_left != 0) {
    gdcurv->gnk1 += (mympi->topoid[1] < nz_left)? mympi->topoid[1] : nz_left;
  }

  // add ghost point
  int nx = ni + 2;
  int nz = nk + 2;

  gdcurv->ni = ni;
  gdcurv->nk = nk;

  gdcurv->nx = nx;
  gdcurv->nz = nz;

  gdcurv->ni1 = 1;
  gdcurv->ni2 = gdcurv->ni1 + ni - 1;

  gdcurv->nk1 = 1;
  gdcurv->nk2 = gdcurv->nk1 + nk - 1;

  // global index end
  gdcurv->gni2 = gdcurv->gni1 + gdcurv->ni - 1;
  gdcurv->gnk2 = gdcurv->gnk1 + gdcurv->nk - 1;
  
  // x dimention varies first
  gdcurv->siz_iz = nx; 
  gdcurv->siz_icmp = nx*nz;

  return ierr;
}

int
set_output_dir(gd_t *gdcurv,
               mympi_t *mympi,
               char *output_dir,
               int verbose)
{
  // output file name
  sprintf(gdcurv->fname_part,"px%d_pz%d", mympi->topoid[0],mympi->topoid[1]);

  // output
  sprintf(gdcurv->output_dir, "%s", output_dir);

  return 0;
}

int
gd_info_print(gd_t *gdcurv, int myid)
{    
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> grid info:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "my rank id is %d\n", myid);
  fprintf(stdout, " nx    = %-10d\n", gdcurv->nx);
  fprintf(stdout, " nz    = %-10d\n", gdcurv->nz);
  fprintf(stdout, " ni    = %-10d\n", gdcurv->ni);
  fprintf(stdout, " nk    = %-10d\n", gdcurv->nk);

  fprintf(stdout, " ni1   = %-10d\n", gdcurv->ni1);
  fprintf(stdout, " ni2   = %-10d\n", gdcurv->ni2);
  fprintf(stdout, " nk1   = %-10d\n", gdcurv->nk1);
  fprintf(stdout, " nk2   = %-10d\n", gdcurv->nk2);

  fprintf(stdout, " ni1_to_glob_phys0   = %-10d\n", gdcurv->gni1);
  fprintf(stdout, " ni2_to_glob_phys0   = %-10d\n", gdcurv->gni2);
  fprintf(stdout, " nk1_to_glob_phys0   = %-10d\n", gdcurv->gnk1);
  fprintf(stdout, " nk2_to_glob_phys0   = %-10d\n", gdcurv->gnk2);

  fflush(stdout);

  return(0);
}

int
gd_curv_coord_exchange(gd_t *gdcurv, int *neighid, MPI_Comm topocomm)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  size_t siz_iz = gdcurv->siz_iz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  // extend to ghosts, using mpi exchange
  // NOTE in different myid, nx(or ny) may not equal
  // so send type DTypeXL not equal recv type DTypeXL
  size_t s_iptr;
  size_t r_iptr;

  MPI_Status status;
  MPI_Datatype DTypeXL, DTypeZL;

  MPI_Type_vector(nz,
                  1,
                  nx,
                  MPI_FLOAT,
                  &DTypeXL);
  MPI_Type_vector(1,
                  nx,
                  1,
                  MPI_FLOAT,
                  &DTypeZL);
  MPI_Type_commit(&DTypeXL);
  MPI_Type_commit(&DTypeZL);

  // to X2
  s_iptr = ni1;
  r_iptr = ni2+1;
  MPI_Sendrecv(&x2d[s_iptr],1,DTypeXL,neighid[0],110,
               &x2d[r_iptr],1,DTypeXL,neighid[1],110,
               topocomm,&status);
  MPI_Sendrecv(&z2d[s_iptr],1,DTypeXL,neighid[0],110,
               &z2d[r_iptr],1,DTypeXL,neighid[1],110,
               topocomm,&status);
  // to X1
  s_iptr = ni2;         
  r_iptr = ni1-1;      
  MPI_Sendrecv(&x2d[s_iptr],1,DTypeXL,neighid[1],120,
               &x2d[r_iptr],1,DTypeXL,neighid[0],120,
               topocomm,&status);
  MPI_Sendrecv(&z2d[s_iptr],1,DTypeXL,neighid[1],120,
               &z2d[r_iptr],1,DTypeXL,neighid[0],120,
               topocomm,&status);
  // to Z2
  s_iptr = nk1 * siz_iz;        
  r_iptr = (nk2+1) * siz_iz;    
  MPI_Sendrecv(&x2d[s_iptr],1,DTypeZL,neighid[2],210,
               &x2d[r_iptr],1,DTypeZL,neighid[3],210,
               topocomm,&status);
  MPI_Sendrecv(&z2d[s_iptr],1,DTypeZL,neighid[2],210,
               &z2d[r_iptr],1,DTypeZL,neighid[3],210,
               topocomm,&status);
  // to Z1
  s_iptr = nk2 * siz_iz;   
  r_iptr = (nk1-1) * siz_iz;     
  MPI_Sendrecv(&x2d[s_iptr],1,DTypeZL,neighid[3],220,
               &x2d[r_iptr],1,DTypeZL,neighid[2],220,
               topocomm,&status);
  MPI_Sendrecv(&z2d[s_iptr],1,DTypeZL,neighid[3],220,
               &z2d[r_iptr],1,DTypeZL,neighid[2],220,
               topocomm,&status);


  return 0;
}

int
grid_mesg_init(mympi_t *mympi, gd_t *gdcurv)
{
  int ni = gdcurv->ni;
  int nk = gdcurv->nk;

  mympi->siz_sbuff_x1 = 2*nk;
  mympi->siz_sbuff_x2 = 2*nk;
  mympi->siz_sbuff_z1 = 2*ni;
  mympi->siz_sbuff_z2 = 2*ni;

  mympi->siz_rbuff_x1 = 2*nk;
  mympi->siz_rbuff_x2 = 2*nk;
  mympi->siz_rbuff_z1 = 2*ni;
  mympi->siz_rbuff_z2 = 2*ni;

  mympi->siz_sbuff = mympi->siz_sbuff_x1 + mympi->siz_sbuff_x2
                   + mympi->siz_sbuff_z1 + mympi->siz_sbuff_z2;

  mympi->siz_rbuff = mympi->siz_rbuff_x1 + mympi->siz_rbuff_x2
                   + mympi->siz_rbuff_z1 + mympi->siz_rbuff_z2;

  mympi->sbuff = (float *) malloc(mympi->siz_sbuff*sizeof(MPI_FLOAT));
  mympi->rbuff = (float *) malloc(mympi->siz_rbuff*sizeof(MPI_FLOAT));
  
  int tag[4] = { 11, 12, 21, 22 };

  mympi->s_reqs = (MPI_Request *) malloc(4*sizeof(MPI_Request));
  mympi->r_reqs = (MPI_Request *) malloc(4*sizeof(MPI_Request));
  // send
  float *sbuff_x1 = mympi->sbuff;
  float *sbuff_x2 = sbuff_x1 + mympi->siz_sbuff_x1;
  float *sbuff_z1 = sbuff_x2 + mympi->siz_sbuff_x2;
  float *sbuff_z2 = sbuff_z1 + mympi->siz_sbuff_z1;
  MPI_Send_init(sbuff_x1, mympi->siz_sbuff_x1, MPI_FLOAT, mympi->neighid[0], tag[0],
                mympi->topocomm, &(mympi->s_reqs[0]));
  MPI_Send_init(sbuff_x2, mympi->siz_sbuff_x2, MPI_FLOAT, mympi->neighid[1], tag[1],
                mympi->topocomm, &(mympi->s_reqs[1]));
  MPI_Send_init(sbuff_z1, mympi->siz_sbuff_z1, MPI_FLOAT, mympi->neighid[2], tag[2],
                mympi->topocomm, &(mympi->s_reqs[2]));
  MPI_Send_init(sbuff_z2, mympi->siz_sbuff_z2, MPI_FLOAT, mympi->neighid[3], tag[3],
                mympi->topocomm, &(mympi->s_reqs[3]));

  // recv
  float *rbuff_x1 = mympi->rbuff;
  float *rbuff_x2 = rbuff_x1 + mympi->siz_rbuff_x1;
  float *rbuff_z1 = rbuff_x2 + mympi->siz_rbuff_x2;
  float *rbuff_z2 = rbuff_z1 + mympi->siz_rbuff_z1;
  MPI_Recv_init(rbuff_x1, mympi->siz_rbuff_x1, MPI_FLOAT, mympi->neighid[0], tag[1],
                mympi->topocomm, &(mympi->r_reqs[0]));
  MPI_Recv_init(rbuff_x2, mympi->siz_rbuff_x2, MPI_FLOAT, mympi->neighid[1], tag[0],
                mympi->topocomm, &(mympi->r_reqs[1]));
  MPI_Recv_init(rbuff_z1, mympi->siz_rbuff_z1, MPI_FLOAT, mympi->neighid[2], tag[3],
                mympi->topocomm, &(mympi->r_reqs[2]));
  MPI_Recv_init(rbuff_z2, mympi->siz_rbuff_z2, MPI_FLOAT, mympi->neighid[3], tag[2],
                mympi->topocomm, &(mympi->r_reqs[3]));

  return 0;
}

int
grid_pack_mesg(mympi_t *mympi, gd_t *gdcurv, float *x2d, float *z2d)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int ni = gdcurv->ni;
  int nk = gdcurv->nk;
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t iptr, iptr1;

  int *topoid = mympi->topoid;
  float *sbuff_x1 = mympi->sbuff;
  float *sbuff_x2 = sbuff_x1 + mympi->siz_sbuff_x1;
  float *sbuff_z1 = sbuff_x2 + mympi->siz_sbuff_x2;
  float *sbuff_z2 = sbuff_z1 + mympi->siz_sbuff_z1;
  // x1 
  for(int k=nk1; k<=nk2; k++)
  {
    iptr = k*nx + ni1;
    iptr1 = k-nk1;
    sbuff_x1[iptr1] = x2d[iptr];
    sbuff_x1[iptr1+nk] = z2d[iptr];
  }

  // x2 
  for(int k=nk1; k<=nk2; k++)
  {
    iptr = k*nx + ni2;
    iptr1 = k-nk1;
    sbuff_x2[iptr1] = x2d[iptr];
    sbuff_x2[iptr1+nk] = z2d[iptr];
  }

  // z1 
  for(int i=ni1; i<=ni2; i++)
  {
    iptr = nk1*nx + i;
    iptr1 = i-ni1;
    sbuff_z1[iptr1] = x2d[iptr];
    sbuff_z1[iptr1+ni] = z2d[iptr];
  }

  // z2 
  for(int i=ni1; i<=ni2; i++)
  {
    iptr = nk2*nx + i;
    iptr1 = i-ni1;
    sbuff_z2[iptr1] = x2d[iptr];
    sbuff_z2[iptr1+ni] = z2d[iptr];
  }

  return 0;
}

int
grid_unpack_mesg(mympi_t *mympi, gd_t *gdcurv, float *x2d, float *z2d)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int ni = gdcurv->ni;
  int nk = gdcurv->nk;
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int *neighid = mympi->neighid;
  int *topoid = mympi->topoid;
  size_t iptr, iptr1;

  float *rbuff_x1 = mympi->rbuff;
  float *rbuff_x2 = rbuff_x1 + mympi->siz_rbuff_x1;
  float *rbuff_z1 = rbuff_x2 + mympi->siz_rbuff_x2;
  float *rbuff_z2 = rbuff_z1 + mympi->siz_rbuff_z1;
  // x1 
  if(neighid[0] != MPI_PROC_NULL)
  {
    for(int k=nk1; k<=nk2; k++)
    {
      iptr = k*nx + (ni1-1);
      iptr1 = k-nk1;
      x2d[iptr] = rbuff_x1[iptr1];
      z2d[iptr] = rbuff_x1[iptr1+nk];
    }
  }

  // x2 
  if(neighid[1] != MPI_PROC_NULL)
  {
    for(int k=nk1; k<=nk2; k++)
    {
      iptr = k*nx + (ni2+1);
      iptr1 = k-nk1;
      x2d[iptr] = rbuff_x2[iptr1];
      z2d[iptr] = rbuff_x2[iptr1+nk];
    }
  }

  // z1 
  if(neighid[2] != MPI_PROC_NULL)
  {
    for(int i=ni1; i<=ni2; i++)
    {
      iptr = (nk1-1)*nx + i;
      iptr1 = i-ni1;
      x2d[iptr] = rbuff_z1[iptr1];
      z2d[iptr] = rbuff_z1[iptr1+ni];
    }
  }

  // z2 
  if(neighid[3] != MPI_PROC_NULL)
  {
    for(int i=ni1; i<=ni2; i++)
    {
      iptr = (nk2+1)*nx + i;
      iptr1 = i-ni1;
      x2d[iptr] = rbuff_z2[iptr1];
      z2d[iptr] = rbuff_z2[iptr1+ni];
    }
  }

  return 0;
}
