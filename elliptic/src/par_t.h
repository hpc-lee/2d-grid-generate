#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "cJSON.h"
#include "constants.h"

#define PAR_MAX_STRLEN 1000
#define ELLI_DIRI 1
#define ELLI_HIGEN 2

#define X_DIRE 1
#define Z_DIRE 2

typedef struct{

  //grid size
  int number_of_grid_points_x;
  int number_of_grid_points_z;

  // MPI
  int number_of_mpiprocs_x;
  int number_of_mpiprocs_z;

  int grid_check;
  int check_orth;
  int check_jac;
  int check_ratio;
  int check_step_xi;
  int check_step_zt;
  int check_smooth_xi;
  int check_smooth_zt;

  char geometry_input_file[PAR_MAX_STRLEN];
  char output_dir[PAR_MAX_STRLEN];
  
  // TFI hermite elliptic-dirichlet
  // elliptic-hilgenstock
  // parabolic hyperbolic
  int method_itype;

  float coef;

  float iter_err;   // iteration error
  int max_iter;  // max iterations
} par_t;

int
par_mpi_get(char *par_fname, int myid, MPI_Comm comm, par_t *par, int verbose);

int 
par_read_from_str(const char *str, par_t *par);

int
par_print(par_t *par);

#endif

