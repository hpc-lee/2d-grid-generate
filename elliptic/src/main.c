#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>
#include <mpi.h>

#include "lib_mem.h"
#include "par_t.h"
#include "constants.h"
#include "gd_t.h"
#include "mympi_t.h"
#include "algebra.h"
#include "io_funcs.h"
#include "quality_check.h"
#include "elliptic.h"

int main(int argc, char** argv)
{
  char *par_fname;
  char err_message[CONST_MAX_STRLEN];

  //-------------------------------------------------------------------------------
  // get commond-line argument
  //-------------------------------------------------------------------------------

  // argc checking
  if (argc < 2) {
    fprintf(stdout,"usage: main_grid_2d <par_file> \n");
    exit(1);
  }

  par_fname = argv[1];

  fprintf(stdout,"par file =  %s\n", par_fname); fflush(stdout);

  // init MPI

  int myid, mpi_size;
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &mpi_size);


  if (myid==0) fprintf(stdout,"comm=%d, size=%d\n", comm, mpi_size); 
  if (myid==0) fprintf(stdout,"par file =  %s\n", par_fname); 

  // read par
  par_t *par = (par_t *) malloc(sizeof(par_t));

  par_mpi_get(par_fname, myid, comm, par);

  if (myid==0) par_print(par);

  gd_t *gdcurv = (gd_t *) malloc(sizeof(gd_t));
  mympi_t *mympi = (mympi_t *) malloc(sizeof(mympi_t));
  // set mpi
  if (myid==0) fprintf(stdout,"set mpi topo ...\n"); 
  mympi_set(mympi,
            par->number_of_mpiprocs_x,
            par->number_of_mpiprocs_z,
            comm,
            myid);

  // set gdinfo
  gd_info_set(gdcurv, mympi,
              par->number_of_grid_points_x,
              par->number_of_grid_points_z);

  gd_info_print(gdcurv,myid);

  init_gdcurv(gdcurv);

  // set str in blk
  set_output_dir(gdcurv, mympi,
                 par->output_dir);
  
  // read bdry and init iter grid 
  bdry_t *bdry = (bdry_t *) malloc(sizeof(bdry_t));
  init_bdry(bdry,par);
  read_bdry(myid,bdry,par->geometry_input_file);

  // linear tfi generate init iter grid
  linear_tfi(gdcurv,bdry,mympi);
  gd_curv_coord_exchange(gdcurv,mympi->neighid,mympi->topocomm);

  grid_mesg_init(mympi, gdcurv);

  time_t t_start = time(NULL);
  switch(par->method_itype)
  {
    case ELLI_DIRI : {
      
      diri_gene(gdcurv,par,mympi);

      break;
    }
    case ELLI_HIGEN : {

      higen_gene(gdcurv,par,mympi);

      break;
    }
  }

  time_t t_end = time(NULL);

  if(myid == 0)
  {
    fprintf(stdout,"\n************************************\n");
    fprintf(stdout,"grid generate running time is :%f s \n", difftime(t_end,t_start));
    fprintf(stdout,"************************************\n \n");
  }

  if(myid == 0) fprintf(stdout,"export coord to file ... \n");
  gd_curv_coord_export(gdcurv,mympi);

  // grid quality check and export quality data
  io_quality_t *io_quality = (io_quality_t *) malloc(sizeof(io_quality_t));
  if(par->grid_check == 1)
  {
    if(myid == 0)
    {
      fprintf(stdout,"****************************************************** \n");
      fprintf(stdout,"***** grid quality check and export quality data ***** \n");
      fprintf(stdout,"****************************************************** \n");
    }
    init_io_quality(io_quality,gdcurv);
    grid_quality_check(io_quality,gdcurv,par,mympi);
  }
  MPI_Finalize();
  return 0;
}

