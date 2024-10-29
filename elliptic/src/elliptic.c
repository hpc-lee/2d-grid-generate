#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "elliptic.h"
#include "constants.h"
#include "solver.h"
#include "gd_t.h"
#include "lib_mem.h"

int
init_src(src_t *src, gd_t *gdcurv)
{
  int ni = gdcurv->ni;
  int nk = gdcurv->nk;
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  
  int total_nx = gdcurv->total_nx;
  int total_nz = gdcurv->total_nz;

  // not include 2 bdry points
  src->P_x1_loc = (float *)mem_calloc_1d_float(total_nz-2, 0.0, "P_x1");
  src->Q_x1_loc = (float *)mem_calloc_1d_float(total_nz-2, 0.0, "P_x1");
  src->P_x2_loc = (float *)mem_calloc_1d_float(total_nz-2, 0.0, "P_x2");
  src->Q_x2_loc = (float *)mem_calloc_1d_float(total_nz-2, 0.0, "P_x2");
  src->P_z1_loc = (float *)mem_calloc_1d_float(total_nx-2, 0.0, "P_z1");
  src->Q_z1_loc = (float *)mem_calloc_1d_float(total_nx-2, 0.0, "P_z1");
  src->P_z2_loc = (float *)mem_calloc_1d_float(total_nx-2, 0.0, "P_z2");
  src->Q_z2_loc = (float *)mem_calloc_1d_float(total_nx-2, 0.0, "P_z2");

  src->P_x1 = (float *)mem_calloc_1d_float(total_nz-2, 0.0, "P_x1");
  src->Q_x1 = (float *)mem_calloc_1d_float(total_nz-2, 0.0, "P_x1");
  src->P_x2 = (float *)mem_calloc_1d_float(total_nz-2, 0.0, "P_x2");
  src->Q_x2 = (float *)mem_calloc_1d_float(total_nz-2, 0.0, "P_x2");
  src->P_z1 = (float *)mem_calloc_1d_float(total_nx-2, 0.0, "P_z1");
  src->Q_z1 = (float *)mem_calloc_1d_float(total_nx-2, 0.0, "P_z1");
  src->P_z2 = (float *)mem_calloc_1d_float(total_nx-2, 0.0, "P_z2");
  src->Q_z2 = (float *)mem_calloc_1d_float(total_nx-2, 0.0, "P_z2");

  // not include 2 ghost points
  src->P = (float *)mem_calloc_1d_float(nx*nz, 0.0, "source P");
  src->Q = (float *)mem_calloc_1d_float(nx*nz, 0.0, "source Q");
}


int
diri_gene(gd_t *gdcurv, par_t *par, mympi_t *mympi)
{
  float err = par->iter_err;
  int max_iter = par->max_iter;
  float *coef = par->coef;

  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t iptr, iptr1, iptr2;

  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *x2d_tmp = gdcurv->x2d_tmp;
  float *z2d_tmp = gdcurv->z2d_tmp;

  MPI_Comm comm = mympi->comm;
  int myid = mympi->myid;
  int *neighid = mympi->neighid;

  int num_of_s_reqs = 4;
  int num_of_r_reqs = 4;
  
  src_t *src = (src_t *) malloc(sizeof(src_t)); 
  init_src(src,gdcurv);
  
  // before update grid. use init grid to
  // calculate ghost points and bdry g11,g22 
  float *p_x; // point_x_bdry
  float *p_z; // point_z_bdry
  float *g11_x;  // g11_x_bdry
  float *g22_z;  // g22_z_bdry
  p_x = (float *)mem_calloc_1d_float(nz*2*2, 0.0, "p_x");
  p_z = (float *)mem_calloc_1d_float(nx*2*2, 0.0, "p_z");
  g11_x = (float *)mem_calloc_1d_float(nz*2, 0.0, "g11_x"); 
  g22_z = (float *)mem_calloc_1d_float(nx*2, 0.0, "g22_z"); 
  ghost_cal(x2d,z2d,nx,nz,p_x,p_z,g11_x,g22_z,neighid);

  set_src_diri(x2d,z2d,gdcurv,src,p_x,p_z,g11_x,g22_z,mympi);
  interp_inner_source(src, gdcurv, coef);

  // copy coord
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++)
    {
      iptr = k*nx + i;
      x2d_tmp[iptr] = x2d[iptr];
      z2d_tmp[iptr] = z2d[iptr];
    }
  }

  // now we default set w=1.0.
  // this is Gauss-Seidel
  float w=1.0; // SOR coef

  int Niter = 0;
  float *ptr_x, *ptr_z;

  float resi, resk;
  float max_resi, max_resk;

  int flag_true = 1;
  float dif1, dif2, dif3, dif_x, dif_z;
  while(flag_true)
  {
    // update solver
    update_SOR(x2d,z2d,x2d_tmp,z2d_tmp,nx,nz,src->P,src->Q,w);
    Niter += 1;

    MPI_Startall(num_of_r_reqs, mympi->r_reqs);

    grid_pack_mesg(mympi,gdcurv,x2d_tmp,z2d_tmp);

    MPI_Startall(num_of_s_reqs, mympi->s_reqs);

    MPI_Waitall(num_of_s_reqs, mympi->s_reqs,MPI_STATUS_IGNORE);
    MPI_Waitall(num_of_r_reqs, mympi->r_reqs,MPI_STATUS_IGNORE);

    grid_unpack_mesg(mympi,gdcurv,x2d_tmp,z2d_tmp);

    max_resi = 0.0;
    max_resk = 0.0;
    // cal iter error
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr  = k*nx + i;
        iptr1 = k*nx + i+1;
        iptr2 = (k+1)*nx + i;
        dif_x = x2d_tmp[iptr]-x2d[iptr];
        dif_z = z2d_tmp[iptr]-z2d[iptr];
        dif1 = sqrt(pow(dif_x,2) + pow(dif_z,2));
        dif_x = x2d[iptr1]-x2d[iptr];
        dif_z = z2d[iptr1]-z2d[iptr];
        dif2 = sqrt(pow(dif_x,2) + pow(dif_z,2));
        dif_x = x2d[iptr2]-x2d[iptr];
        dif_z = z2d[iptr2]-z2d[iptr];
        dif3 = sqrt(pow(dif_x,2) + pow(dif_z,2));
        resi = dif1/dif2;
        resk = dif1/dif3;
        max_resi = fmax(max_resi,resi);
        max_resk = fmax(max_resk,resk);
      }
    }

    float sendbuf = max_resi;
    MPI_Allreduce(&sendbuf, &max_resi, 1, MPI_REAL, MPI_MAX, comm);
    sendbuf = max_resk;
    MPI_Allreduce(&sendbuf, &max_resk, 1, MPI_REAL, MPI_MAX, comm);

    if(myid == 0)
    {
      fprintf(stdout,"number of iter is %d\n", Niter);
      fprintf(stdout,"max_resi is %f, max_resk is %f\n", max_resi, max_resk);
    }
    
    // swap pointer, avoid coping
    ptr_x = x2d;
    ptr_z = z2d;
    x2d = x2d_tmp;
    z2d = z2d_tmp;
    x2d_tmp = ptr_x;
    z2d_tmp = ptr_z;
    
    set_src_diri(x2d,z2d,gdcurv,src,p_x,p_z,g11_x,g22_z,mympi);
    interp_inner_source(src, gdcurv, coef);

    if(Niter>max_iter) {
      flag_true = 0;
    }

    if(max_resi < err && max_resk < err) {
      flag_true = 0;
    }
  }

  gdcurv->v3d = x2d;
  gdcurv->x2d = x2d;
  gdcurv->z2d = z2d;

  gdcurv->v3d_tmp = x2d_tmp;

  free(gdcurv->v3d_tmp); // free temp grid space
  free(p_x);
  free(p_z);
  free(g11_x);
  free(g22_z);

  return 0;
}

int
set_src_diri(float *x2d, float *z2d, gd_t *gdcurv, 
             src_t *src, float *p_x, float *p_z,
             float *g11_x, float *g22_z, mympi_t *mympi)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int ni = gdcurv->ni;
  int nk = gdcurv->nk;
  int gni1 = gdcurv->gni1;
  int gnk1 = gdcurv->gnk1;
  int total_nx = gdcurv->total_nx;
  int total_nz = gdcurv->total_nz;
  int gni,gnk;

  float *p_x1;
  float *p_x2;
  float *p_z1;
  float *p_z2;
  p_x1 = p_x;
  p_x2 = p_x + nz*2;
  p_z1 = p_z;
  p_z2 = p_z + nx*2;
  float *g11_x1, *g11_x2;
  float *g22_z1, *g22_z2;
  g11_x1 = g11_x;
  g11_x2 = g11_x + nz;
  g22_z1 = g22_z;
  g22_z2 = g22_z + nx;

  MPI_Comm topocomm = mympi->topocomm;
  int *neighid = mympi->neighid;

  float *P_x1_loc = src->P_x1_loc;
  float *Q_x1_loc = src->Q_x1_loc;
  float *P_x2_loc = src->P_x2_loc;
  float *Q_x2_loc = src->Q_x2_loc;
  float *P_z1_loc = src->P_z1_loc;
  float *Q_z1_loc = src->Q_z1_loc;
  float *P_z2_loc = src->P_z2_loc;
  float *Q_z2_loc = src->Q_z2_loc;

  float *P_x1 = src->P_x1;
  float *Q_x1 = src->Q_x1;
  float *P_x2 = src->P_x2;
  float *Q_x2 = src->Q_x2;
  float *P_z1 = src->P_z1;
  float *Q_z1 = src->Q_z1;
  float *P_z2 = src->P_z2;
  float *Q_z2 = src->Q_z2;

  size_t iptr, iptr1, iptr2, iptr3;
  float x_xi,z_xi,x_xixi,z_xixi;
  float x_zt,z_zt,x_ztzt,z_ztzt;
  float g11, g22;

  // bdry x1 xi=0
  if(neighid[0] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++)
    {
      iptr  = k*nx;         //(0,k)
      iptr1 = (k+1)*nx;     //(0,k+1)
      iptr2 = (k-1)*nx;     //(0,k-1)
      x_zt = 0.5*(x2d[iptr1] - x2d[iptr2]);
      z_zt = 0.5*(z2d[iptr1] - z2d[iptr2]);
      x_ztzt = x2d[iptr1] + x2d[iptr2] - 2*x2d[iptr];
      z_ztzt = z2d[iptr1] + z2d[iptr2] - 2*z2d[iptr];

      iptr3 = k*nx + 1;  //(1,k)
      x_xi = x2d[iptr3] - x2d[iptr];
      z_xi = z2d[iptr3] - z2d[iptr];
      x_xixi = p_x1[k] + x2d[iptr3] - 2*x2d[iptr];
      z_xixi = p_x1[k+nz] + z2d[iptr3] - 2*z2d[iptr];
      g22 = pow(x_zt,2) + pow(z_zt,2);
      
      gnk = gnk1 + (k-1);
      P_x1_loc[gnk] = -(x_xi*x_xixi + z_xi*z_xixi)/g11_x1[k]  
                      -(x_xi*x_ztzt + z_xi*z_ztzt)/g22;
      Q_x1_loc[gnk] = -(x_zt*x_xixi + z_zt*z_xixi)/g11_x1[k]  
                      -(x_zt*x_ztzt + z_zt*z_ztzt)/g22;
    }
  }

  // bdry x2 xi=1
  if(neighid[1] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++)
    {
      iptr  = k*nx + (nx-1);         //(nx-1,k)
      iptr1 = (k+1)*nx + (nx-1);     //(nx-1,k+1)
      iptr2 = (k-1)*nx + (nx-1);     //(nx-1,k-1)
      x_zt = 0.5*(x2d[iptr1] - x2d[iptr2]);
      z_zt = 0.5*(z2d[iptr1] - z2d[iptr2]);
      x_ztzt = x2d[iptr1] + x2d[iptr2] - 2*x2d[iptr];
      z_ztzt = z2d[iptr1] + z2d[iptr2] - 2*z2d[iptr];

      iptr3 = k*nx + (nx-2);  //(nx-2,k)
      x_xi = x2d[iptr] - x2d[iptr3];
      z_xi = z2d[iptr] - z2d[iptr3];
      x_xixi = p_x2[k] + x2d[iptr3] - 2*x2d[iptr];
      z_xixi = p_x2[k+nz] + z2d[iptr3] - 2*z2d[iptr];
      g22 = pow(x_zt,2) + pow(z_zt,2);

      gnk = gnk1 + (k-1);
      P_x2_loc[gnk] = -(x_xi*x_xixi + z_xi*z_xixi)/g11_x2[k]  
                      -(x_xi*x_ztzt + z_xi*z_ztzt)/g22;
      Q_x2_loc[gnk] = -(x_zt*x_xixi + z_zt*z_xixi)/g11_x2[k]  
                      -(x_zt*x_ztzt + z_zt*z_ztzt)/g22;
    }
  }

  // bdry z1 zt=0
  if(neighid[2] == MPI_PROC_NULL)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = i;         //(i,0)
      iptr1 = i+1;       //(i+1,0)
      iptr2 = i-1;       //(i-1,0)
      x_xi = 0.5*(x2d[iptr1] - x2d[iptr2]);
      z_xi = 0.5*(z2d[iptr1] - z2d[iptr2]);
      x_xixi = x2d[iptr1] + x2d[iptr2] - 2*x2d[iptr];
      z_xixi = z2d[iptr1] + z2d[iptr2] - 2*z2d[iptr];

      iptr3 = 1*nx + i;  //(i,1)
      x_zt = x2d[iptr3] - x2d[iptr];
      z_zt = z2d[iptr3] - z2d[iptr];
      x_ztzt = p_z1[i] + x2d[iptr3] - 2*x2d[iptr];
      z_ztzt = p_z1[i+nx] + z2d[iptr3] - 2*z2d[iptr];
      g11 = pow(x_xi,2) + pow(z_xi,2);
      
      gni = gni1 + (i-1);
      P_z1_loc[gni] = -(x_xi*x_xixi + z_xi*z_xixi)/g11  
                      -(x_xi*x_ztzt + z_xi*z_ztzt)/g22_z1[i];
      Q_z1_loc[gni] = -(x_zt*x_xixi + z_zt*z_xixi)/g11  
                      -(x_zt*x_ztzt + z_zt*z_ztzt)/g22_z1[i];
    }
  }

  // bdry z2 zt=1
  if(neighid[3] == MPI_PROC_NULL)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = (nz-1)*nx + i;      //(i,nz-1)
      iptr1 = (nz-1)*nx + (i+1);  //(i+1,nz-1)
      iptr2 = (nz-1)*nx + (i-1);  //(i-1,nz-1)
      x_xi = 0.5*(x2d[iptr1] - x2d[iptr2]);
      z_xi = 0.5*(z2d[iptr1] - z2d[iptr2]);
      x_xixi = x2d[iptr1] + x2d[iptr2] - 2*x2d[iptr];
      z_xixi = z2d[iptr1] + z2d[iptr2] - 2*z2d[iptr];

      iptr3 = (nz-2)*nx + i;  //(i,nz-2)
      x_zt = x2d[iptr] - x2d[iptr3];
      z_zt = z2d[iptr] - z2d[iptr3];
      x_ztzt = p_z2[i] + x2d[iptr3] - 2*x2d[iptr];
      z_ztzt = p_z2[i+nx] + z2d[iptr3] - 2*z2d[iptr];
      g11 = pow(x_xi,2) + pow(z_xi,2);

      gni = gni1 + (i-1);
      P_z2_loc[gni] = -(x_xi*x_xixi + z_xi*z_xixi)/g11  
                      -(x_xi*x_ztzt + z_xi*z_ztzt)/g22_z2[i];
      Q_z2_loc[gni] = -(x_zt*x_xixi + z_zt*z_xixi)/g11  
                      -(x_zt*x_ztzt + z_zt*z_ztzt)/g22_z2[i];
    }
  }

  MPI_Barrier(topocomm);

  MPI_Allreduce(P_x1_loc, P_x1, total_nz-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_x1_loc, Q_x1, total_nz-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_x2_loc, P_x2, total_nz-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_x2_loc, Q_x2, total_nz-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_z1_loc, P_z1, total_nx-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_z1_loc, Q_z1, total_nx-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_z2_loc, P_z2, total_nx-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_z2_loc, Q_z2, total_nx-2, MPI_FLOAT, MPI_SUM, topocomm);

  return 0;
}

int
ghost_cal(float *x2d, float *z2d, int nx, int nz, float *p_x, float *p_z,
          float *g11_x, float *g22_z, int *neighid)
{
  size_t iptr1,iptr2,iptr3,iptr4;
  float x_xi,z_xi,x_zt,z_zt;
  float x_xi0,z_xi0,x_zt0,z_zt0;
  float vn_xi, vn_zt, vn_xi0, vn_zt0;
  float len_vn, dot;
  float *p_x1;
  float *p_x2;
  float *p_z1;
  float *p_z2;
  p_x1 = p_x;
  p_x2 = p_x + nz*2;
  p_z1 = p_z;
  p_z2 = p_z + nx*2;
  // bdry x1
  if(neighid[0] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++)
    {
      iptr1 = k*nx + 1;    //(1,k)
      iptr2 = k*nx + 0;    //(0,k)
      iptr3 = (k+1)*nx + 0;//(0,k+1)
      iptr4 = (k-1)*nx + 0;//(0,k-1)
      x_xi0 = x2d[iptr1] - x2d[iptr2];
      z_xi0 = z2d[iptr1] - z2d[iptr2];
      x_zt = 0.5*(x2d[iptr3] - x2d[iptr4]);
      z_zt = 0.5*(z2d[iptr3] - z2d[iptr4]);
      // orth vector
      vn_xi =  z_zt;
      vn_zt = -x_zt;
      len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
      // norm
      vn_xi0 = vn_xi/len_vn;
      vn_zt0 = vn_zt/len_vn;
      // projection from r_xi0 to vn
      dot = x_xi0*vn_xi0 + z_xi0*vn_zt0;
      x_xi = dot*vn_xi0;
      z_xi = dot*vn_zt0;

      p_x1[k] = x2d[iptr2] - x_xi;
      p_x1[k+nz] = z2d[iptr2] - z_xi;
      g11_x[k] = pow(x_xi,2) + pow(z_xi,2);
    }
  }
  // bdry x2
  if(neighid[1] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++)
    {
      iptr1 = k*nx + (nx-1);    //(nx-1,k)
      iptr2 = k*nx + (nx-2);    //(nx-2,k)
      iptr3 = (k+1)*nx + (nx-1);//(nx-1,k+1)
      iptr4 = (k-1)*nx + (nx-1);//(nx-1,k-1)
      x_xi0 = x2d[iptr1] - x2d[iptr2];
      z_xi0 = z2d[iptr1] - z2d[iptr2];
      x_zt = 0.5*(x2d[iptr3] - x2d[iptr4]);
      z_zt = 0.5*(z2d[iptr3] - z2d[iptr4]);
      // orth vector
      vn_xi =  z_zt;
      vn_zt = -x_zt;
      len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
      // norm
      vn_xi0 = vn_xi/len_vn;
      vn_zt0 = vn_zt/len_vn;
      // projection from r_xi0 to vn
      dot = x_xi0*vn_xi0 + z_xi0*vn_zt0;
      x_xi = dot*vn_xi0;
      z_xi = dot*vn_zt0;

      p_x2[k] = x2d[iptr1] + x_xi;
      p_x2[k+nz] = z2d[iptr1] + z_xi;
      g11_x[k+nz] = pow(x_xi,2) + pow(z_xi,2);
    }
  }
  // bdry z1
  if(neighid[2] == MPI_PROC_NULL)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = 1*nx + i;    //(i,1)
      iptr2 = 0*nx + i;    //(i,0)
      iptr3 = 0*nx + i+1;  //(i+1,0)
      iptr4 = 0*nx + i-1;  //(i-1,0)
      x_zt0 = x2d[iptr1] - x2d[iptr2];
      z_zt0 = z2d[iptr1] - z2d[iptr2];
      x_xi = 0.5*(x2d[iptr3] - x2d[iptr4]);
      z_xi = 0.5*(z2d[iptr3] - z2d[iptr4]);
      // orth vector
      vn_xi = -z_xi;
      vn_zt =  x_xi;
      len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
      // norm
      vn_xi0 = vn_xi/len_vn;
      vn_zt0 = vn_zt/len_vn;
      // projection from r_zt0 to vn
      dot = x_zt0*vn_xi0 + z_zt0*vn_zt0;
      x_zt = dot*vn_xi0;
      z_zt = dot*vn_zt0;

      p_z1[i] = x2d[iptr2] - x_zt;
      p_z1[i+nx] = z2d[iptr2] - z_zt;
      g22_z[i] = pow(x_zt,2) + pow(z_zt,2);
    }
  }
  // bdry z2
  if(neighid[3] == MPI_PROC_NULL)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (nz-1)*nx + i;    //(i,nz-1)
      iptr2 = (nz-2)*nx + i;    //(i,nz-2)
      iptr3 = (nz-1)*nx + i+1;  //(i+1,nz-1)
      iptr4 = (nz-1)*nx + i-1;  //(i-1,nz-1)
      x_zt0 = x2d[iptr1] - x2d[iptr2];
      z_zt0 = z2d[iptr1] - z2d[iptr2];
      x_xi = 0.5*(x2d[iptr3] - x2d[iptr4]);
      z_xi = 0.5*(z2d[iptr3] - z2d[iptr4]);
      // orth vector
      vn_xi = -z_xi;
      vn_zt =  x_xi;
      len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
      // norm
      vn_xi0 = vn_xi/len_vn;
      vn_zt0 = vn_zt/len_vn;
      // projection from r_zt0 to vn
      dot = x_zt0*vn_xi0 + z_zt0*vn_zt0;
      x_zt = dot*vn_xi0;
      z_zt = dot*vn_zt0;

      p_z2[i] = x2d[iptr1] + x_zt;
      p_z2[i+nx] = z2d[iptr1] + z_zt;
      g22_z[i+nx] = pow(x_zt,2) + pow(z_zt,2);
    }
  }

  return 0;
}

int
higen_gene(gd_t *gdcurv, par_t *par, mympi_t *mympi)
{
  float err = par->iter_err;
  int max_iter = par->max_iter;
  float *coef = par->coef;

  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  size_t iptr, iptr1, iptr2;

  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *x2d_tmp = gdcurv->x2d_tmp;
  float *z2d_tmp = gdcurv->z2d_tmp;

  MPI_Comm comm = mympi->comm;
  int myid = mympi->myid;
  int *neighid = mympi->neighid;

  int num_of_s_reqs = 4;
  int num_of_r_reqs = 4;
  
  src_t *src = (src_t *) malloc(sizeof(src_t)); 
  init_src(src,gdcurv);

  float *dx1 = (float *)mem_calloc_1d_float(nz, 0.0, "dx1"); 
  float *dx2 = (float *)mem_calloc_1d_float(nz, 0.0, "dx2"); 
  float *dz1 = (float *)mem_calloc_1d_float(nx, 0.0, "dz1"); 
  float *dz2 = (float *)mem_calloc_1d_float(nx, 0.0, "dz2"); 
  dist_cal(gdcurv,dx1,dx2,dz1,dz2,neighid);

  set_src_higen(x2d,z2d,gdcurv,src,dx1,dx2,dz1,dz2,mympi);
  interp_inner_source(src, gdcurv, coef);

  // copy coord
  for(int k=0; k<nz; k++) {
    for(int i=0; i<nx; i++)
    {
      iptr = k*nx + i;
      x2d_tmp[iptr] = x2d[iptr];
      z2d_tmp[iptr] = z2d[iptr];
    }
  }

  // now we default set w=1.0.
  // this is Gauss-Seidel
  float w=1.0; // SOR coef

  int Niter = 0;
  float *ptr_x, *ptr_z;

  float resi, resk;
  float max_resi, max_resk;

  int flag_true = 1;
  float dif1, dif2, dif3, dif_x, dif_z;
  while(flag_true)
  {
    // update solver
    update_SOR(x2d,z2d,x2d_tmp,z2d_tmp,nx,nz,src->P,src->Q,w);
    Niter += 1;

    MPI_Startall(num_of_r_reqs, mympi->r_reqs);

    grid_pack_mesg(mympi,gdcurv,x2d_tmp,z2d_tmp);

    MPI_Startall(num_of_s_reqs, mympi->s_reqs);

    MPI_Waitall(num_of_s_reqs, mympi->s_reqs,MPI_STATUS_IGNORE);
    MPI_Waitall(num_of_r_reqs, mympi->r_reqs,MPI_STATUS_IGNORE);

    grid_unpack_mesg(mympi,gdcurv,x2d_tmp,z2d_tmp);

    max_resi = 0.0;
    max_resk = 0.0;
    // cal iter error
    for(int k=1; k<nz-1; k++) {
      for(int i=1; i<nx-1; i++)
      {
        iptr = k*nx + i;
        iptr1 = k*nx + i+1;
        iptr2 = (k+1)*nx + i;
        dif_x = x2d_tmp[iptr]-x2d[iptr];
        dif_z = z2d_tmp[iptr]-z2d[iptr];
        dif1 = sqrt(pow(dif_x,2) + pow(dif_z,2));
        dif_x = x2d[iptr1]-x2d[iptr];
        dif_z = z2d[iptr1]-z2d[iptr];
        dif2 = sqrt(pow(dif_x,2) + pow(dif_z,2));
        dif_x = x2d[iptr2]-x2d[iptr];
        dif_z = z2d[iptr2]-z2d[iptr];
        dif3 = sqrt(pow(dif_x,2) + pow(dif_z,2));
        resi = dif1/dif2;
        resk = dif1/dif3;
        max_resi = fmax(max_resi,resi);
        max_resk = fmax(max_resk,resk);
      }
    }

    float sendbuf = max_resi;
    MPI_Allreduce(&sendbuf, &max_resi, 1, MPI_REAL, MPI_MAX, comm);
    sendbuf = max_resk;
    MPI_Allreduce(&sendbuf, &max_resk, 1, MPI_REAL, MPI_MAX, comm);

    if(myid == 0)
    {
      fprintf(stdout,"number of iter is %d\n", Niter);
      fprintf(stdout,"max_resi is %f, max_resk is %f\n", max_resi, max_resk);
    }

    // swap pointer, avoid coping
    ptr_x = x2d;
    ptr_z = z2d;
    x2d = x2d_tmp;
    z2d = z2d_tmp;
    x2d_tmp = ptr_x;
    z2d_tmp = ptr_z;
    
    set_src_higen(x2d,z2d,gdcurv,src,dx1,dx2,dz1,dz2,mympi);
    interp_inner_source(src, gdcurv, coef);

    fprintf(stdout,"number of iter is %d\n", Niter);
    fprintf(stdout,"max_resi is %f, max_resk is %f\n", max_resi, max_resk);

    if(Niter>max_iter) {
      flag_true = 0;
    }

    if(max_resi < err && max_resk < err) {
      flag_true = 0;
    }

  }

  gdcurv->v3d = x2d;
  gdcurv->x2d = x2d;
  gdcurv->z2d = z2d;

  gdcurv->v3d_tmp = x2d_tmp;

  free(gdcurv->v3d_tmp); // free temp grid space
  free(dx1);
  free(dx2);
  free(dz1);
  free(dz2);

  return 0;
}

int
set_src_higen(float *x2d, float *z2d, gd_t *gdcurv, 
              src_t *src, float *dx1, float *dx2,
              float *dz1, float *dz2, mympi_t *mympi)
             
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int ni = gdcurv->ni;
  int nk = gdcurv->nk;
  int gni1 = gdcurv->gni1;
  int gnk1 = gdcurv->gnk1;
  int total_nx = gdcurv->total_nx;
  int total_nz = gdcurv->total_nz;
  int gni,gnk;

  MPI_Comm topocomm = mympi->topocomm;
  int *neighid = mympi->neighid;

  float *P_x1_loc = src->P_x1_loc;
  float *Q_x1_loc = src->Q_x1_loc;
  float *P_x2_loc = src->P_x2_loc;
  float *Q_x2_loc = src->Q_x2_loc;
  float *P_z1_loc = src->P_z1_loc;
  float *Q_z1_loc = src->Q_z1_loc;
  float *P_z2_loc = src->P_z2_loc;
  float *Q_z2_loc = src->Q_z2_loc;

  float *P_x1 = src->P_x1;
  float *Q_x1 = src->Q_x1;
  float *P_x2 = src->P_x2;
  float *Q_x2 = src->Q_x2;
  float *P_z1 = src->P_z1;
  float *Q_z1 = src->Q_z1;
  float *P_z2 = src->P_z2;
  float *Q_z2 = src->Q_z2;

  float theta0 = PI/2;
  float a = 0.1;
  size_t iptr,iptr1,iptr2,iptr3;
  float dot, len_xi, len_zt, dif_dis;
  float cos_theta, theta, dif_theta;
  float x_xi, z_xi, x_zt, z_zt;

  // bdry x1 xi=0
  if(neighid[0] == MPI_PROC_NULL)
  {
    // bdry x1 xi=0
    for(int k=1; k<nz-1; k++)
    {
      iptr  =  k*nx + 0;         //(0,k)
      iptr1 = (k+1)*nx + 0;      //(0,k+1)
      iptr2 = (k-1)*nx + 0;      //(0,k-1)
      iptr3 = k*nx + 1;  //(1,k)

      x_zt = 0.5*(x2d[iptr1]-x2d[iptr2]);
      z_zt = 0.5*(z2d[iptr1]-z2d[iptr2]);
      x_xi = x2d[iptr3] - x2d[iptr];
      z_xi = z2d[iptr3] - z2d[iptr];
      
      // cos(theta) = a.b/(|a|*|b|)
      dot = x_xi*x_zt + z_xi*z_zt;
      len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
      len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
      cos_theta = dot/(len_xi*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      dif_dis = (dx1[k]-len_xi)/dx1[k];

      gnk = gnk1 + (k-1);
      Q_x1_loc[gnk] = Q_x1_loc[gnk] - a*tanh(dif_theta);
      P_x1_loc[gnk] = P_x1_loc[gnk] + a*tanh(dif_dis);
    }
  }

  // bdry x2 xi=1
  if(neighid[1] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++)
    {
      iptr  =  k*nx + (nx-1);         //(nx-1,k)
      iptr1 = (k+1)*nx + (nx-1);      //(nx-1,k+1)
      iptr2 = (k-1)*nx + (nx-1);      //(nx-1,k-1)
      iptr3 = k*nx + (nx-2);  //(nx-2,k)

      x_zt = 0.5*(x2d[iptr1]-x2d[iptr2]);
      z_zt = 0.5*(z2d[iptr1]-z2d[iptr2]);
      x_xi = x2d[iptr] - x2d[iptr3];
      z_xi = z2d[iptr] - z2d[iptr3];
      
      // cos(theta) = a.b/(|a|*|b|)
      dot = x_xi*x_zt + z_xi*z_zt;
      len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
      len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
      cos_theta = dot/(len_xi*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      dif_dis = (dx2[k]-len_xi)/dx2[k];

      gnk = gnk1 + (k-1);
      Q_x2_loc[gnk] = Q_x2_loc[gnk] + a*tanh(dif_theta);
      P_x2_loc[gnk] = P_x2_loc[gnk] - a*tanh(dif_dis);
    }
  }

  // bdry z1 zt=0
  if(neighid[2] == MPI_PROC_NULL)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = i;         //(i,0)
      iptr1 = i+1;       //(i+1,0)
      iptr2 = i-1;       //(i-1,0)
      iptr3 = 1*nx + i;  //(i,1)

      x_xi = 0.5*(x2d[iptr1]-x2d[iptr2]);
      z_xi = 0.5*(z2d[iptr1]-z2d[iptr2]);
      x_zt = x2d[iptr3] - x2d[iptr];
      z_zt = z2d[iptr3] - z2d[iptr];
      
      // cos(theta) = a.b/(|a|*|b|)
      dot = x_xi*x_zt + z_xi*z_zt;
      len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
      len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
      cos_theta = dot/(len_xi*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      dif_dis = (dz1[i]-len_zt)/dz1[i];
    
      gni = gni1 + (i-1);
      P_z1_loc[gni] = P_z1_loc[gni] - a*tanh(dif_theta);
      Q_z1_loc[gni] = Q_z1_loc[gni] + a*tanh(dif_dis);
    }
  }

  // bdry z2 zt=1
  if(neighid[3] == MPI_PROC_NULL)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr  = (nz-1)*nx + i;      //(i,nz-1)
      iptr1 = (nz-1)*nx + i+1;    //(i+1,nz-1)
      iptr2 = (nz-1)*nx + i-1;    //(i-1,nz-1)
      iptr3 = (nz-2)*nx + i;      //(i,nz-2)

      x_xi = 0.5*(x2d[iptr1]-x2d[iptr2]);
      z_xi = 0.5*(z2d[iptr1]-z2d[iptr2]);
      x_zt = x2d[iptr] - x2d[iptr3];
      z_zt = z2d[iptr] - z2d[iptr3];
      
      // cos(theta) = a.b/(|a|*|b|)
      dot = x_xi*x_zt + z_xi*z_zt;
      len_xi = sqrt(pow(x_xi,2) + pow(z_xi,2));
      len_zt = sqrt(pow(x_zt,2) + pow(z_zt,2));
      cos_theta = dot/(len_xi*len_zt);
      theta = acos(cos_theta);
      dif_theta = (theta0-theta)/theta0;
      dif_dis = (dz2[i]-len_zt)/dz2[i];

      gni = gni1 + (i-1);
      P_z2_loc[gni] = P_z2_loc[gni] + a*tanh(dif_theta);
      Q_z2_loc[gni] = Q_z2_loc[gni] - a*tanh(dif_dis);
    }
  }

  MPI_Barrier(topocomm);

  MPI_Allreduce(P_x1_loc, P_x1, total_nz-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_x1_loc, Q_x1, total_nz-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_x2_loc, P_x2, total_nz-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_x2_loc, Q_x2, total_nz-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_z1_loc, P_z1, total_nx-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_z1_loc, Q_z1, total_nx-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(P_z2_loc, P_z2, total_nx-2, MPI_FLOAT, MPI_SUM, topocomm);
  MPI_Allreduce(Q_z2_loc, Q_z2, total_nx-2, MPI_FLOAT, MPI_SUM, topocomm);

  return 0;
}

int
dist_cal(gd_t *gdcurv, float *dx1, float *dx2, 
         float *dz1, float *dz2, int *neighid)
{
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  size_t iptr1,iptr2,iptr3,iptr4;
  float x_xi,z_xi,x_zt,z_zt;
  float x_xi0,z_xi0,x_zt0,z_zt0;
  float vn_xi, vn_zt, vn_xi0, vn_zt0, len_vn;
  // bdry x1
  if(neighid[0] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++)
    {
      iptr1 = k*nx + 1;    //(1,k)
      iptr2 = k*nx + 0;    //(0,k)
      iptr3 = (k+1)*nx + 0;//(0,k+1)
      iptr4 = (k-1)*nx + 0;//(0,k-1)
      x_xi0 = x2d[iptr1] - x2d[iptr2];
      z_xi0 = z2d[iptr1] - z2d[iptr2];
      x_zt = 0.5*(x2d[iptr3] - x2d[iptr4]);
      z_zt = 0.5*(z2d[iptr3] - z2d[iptr4]);
      // orth vector
      vn_xi =  z_zt;
      vn_zt = -x_zt;
      len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
      // norm
      vn_xi0 = vn_xi/len_vn;
      vn_zt0 = vn_zt/len_vn;
      // projection from r_xi0 to vn
      dx1[k] = x_xi0*vn_xi0 + z_xi0*vn_zt0;
    }
  }
  // bdry x2
  if(neighid[1] == MPI_PROC_NULL)
  {
    for(int k=1; k<nz-1; k++)
    {
      iptr1 = k*nx + (nx-1);    //(nx-1,k)
      iptr2 = k*nx + (nx-2);    //(nx-2,k)
      iptr3 = (k+1)*nx + (nx-1);//(nx-1,k+1)
      iptr4 = (k-1)*nx + (nx-1);//(nx-1,k-1)
      x_xi0 = x2d[iptr1] - x2d[iptr2];
      z_xi0 = z2d[iptr1] - z2d[iptr2];
      x_zt = 0.5*(x2d[iptr3] - x2d[iptr4]);
      z_zt = 0.5*(z2d[iptr3] - z2d[iptr4]);
      // orth vector
      vn_xi =  z_zt;
      vn_zt = -x_zt;
      len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
      // norm
      vn_xi0 = vn_xi/len_vn;
      vn_zt0 = vn_zt/len_vn;
      // projection from r_xi0 to vn
      dx2[k] = x_xi0*vn_xi0 + z_xi0*vn_zt0;
    }
  }
  // bdry z1
  if(neighid[2] == MPI_PROC_NULL)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = 1*nx + i;    //(i,1)
      iptr2 = 0*nx + i;    //(i,0)
      iptr3 = 0*nx + i+1;  //(i+1,0)
      iptr4 = 0*nx + i-1;  //(i-1,0)
      x_zt0 = x2d[iptr1] - x2d[iptr2];
      z_zt0 = z2d[iptr1] - z2d[iptr2];
      x_xi = 0.5*(x2d[iptr3] - x2d[iptr4]);
      z_xi = 0.5*(z2d[iptr3] - z2d[iptr4]);
      // orth vector
      vn_xi = -z_xi;
      vn_zt =  x_xi;
      len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
      // norm
      vn_xi0 = vn_xi/len_vn;
      vn_zt0 = vn_zt/len_vn;
      // projection from r_zt0 to vn
      dz1[i] = x_zt0*vn_xi0 + z_zt0*vn_zt0;
    }
  }
  // bdry z2
  if(neighid[3] == MPI_PROC_NULL)
  {
    for(int i=1; i<nx-1; i++)
    {
      iptr1 = (nz-1)*nx + i;    //(i,nz-1)
      iptr2 = (nz-2)*nx + i;    //(i,nz-2)
      iptr3 = (nz-1)*nx + i+1;  //(i+1,nz-1)
      iptr4 = (nz-1)*nx + i-1;  //(i-1,nz-1)
      x_zt0 = x2d[iptr1] - x2d[iptr2];
      z_zt0 = z2d[iptr1] - z2d[iptr2];
      x_xi = 0.5*(x2d[iptr3] - x2d[iptr4]);
      z_xi = 0.5*(z2d[iptr3] - z2d[iptr4]);
      // orth vector
      vn_xi = -z_xi;
      vn_zt =  x_xi;
      len_vn = sqrt(pow(vn_xi,2)+pow(vn_zt,2));
      // norm
      vn_xi0 = vn_xi/len_vn;
      vn_zt0 = vn_zt/len_vn;
      // projection from r_zt0 to vn
      dz2[i] = x_zt0*vn_xi0 + z_zt0*vn_zt0;
    }
  }

  return 0;
}


int interp_inner_source(src_t *src, gd_t *gdcurv, float *coef)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int nx = gdcurv->nx;
  int nz = gdcurv->nz;
  int gni1 = gdcurv->gni1;
  int gnk1 = gdcurv->gnk1;
  int total_nx = gdcurv->total_nx;
  int total_nz = gdcurv->total_nz;
  
  int gni, gnk;

  float *P_x1 = src->P_x1;
  float *Q_x1 = src->Q_x1;
  float *P_x2 = src->P_x2;
  float *Q_x2 = src->Q_x2;
  float *P_z1 = src->P_z1;
  float *Q_z1 = src->Q_z1;
  float *P_z2 = src->P_z2;
  float *Q_z2 = src->Q_z2;
  float *P = src->P;
  float *Q = src->Q;

  float xi,zt,c0,c1,r0,r1;
  size_t iptr;
  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      gnk = gnk1 + k-1; 
      gni = gni1 + i;
      xi = (1.0*gni)/(total_nx-1);

      c0 = 1-xi;
      c1 = xi;

      r0 = exp(-coef[0]*xi);
      r1 = exp(-coef[1]*(1-xi)); 
      
      iptr  = k*nx + i;
      P[iptr] = c0*P_x1[gnk] + c1*P_x2[gnk];
      Q[iptr] = r0*Q_x1[gnk] + r1*Q_x2[gnk];
    }
  }

  for(int k=1; k<nz-1; k++) {
    for(int i=1; i<nx-1; i++)
    {
      gnk = gnk1 + k; 
      gni = gni1 + i-1;
      zt = (1.0*gnk)/(total_nz-1);
      c0 = 1-zt;
      c1 = zt;

      r0 = exp(-coef[2]*zt);
      r1 = exp(-coef[3]*(1-zt)); 
      
      iptr  = k*nx + i;
      P[iptr] = P[iptr] + r0*P_z1[gni] + r1*P_z2[gni];
      Q[iptr] = Q[iptr] + c0*Q_z1[gni] + c1*Q_z2[gni];
    }
  }

  return 0;
}

