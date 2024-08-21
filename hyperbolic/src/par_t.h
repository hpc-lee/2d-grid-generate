#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cJSON.h"
#include "constants.h"
#include "par_t.h"

#define PAR_MAX_STRLEN 1000
#define TFI 1
#define PARABOLIC 2
#define HYPERBOLIC 3

#define X_DIRE 1
#define Z_DIRE 2

typedef struct{

  int number_of_grid_points_x;
  int number_of_grid_points_z;

  int grid_check;
  int check_orth;
  int check_jac;
  int check_ratio;
  int check_step_xi;
  int check_step_zt;
  int check_smooth_xi;
  int check_smooth_zt;

  char geometry_input_file[PAR_MAX_STRLEN];
  char grid_export_dir[PAR_MAX_STRLEN];

  char step_input_file[PAR_MAX_STRLEN];
  char bdry_file1[PAR_MAX_STRLEN];
  char bdry_file2[PAR_MAX_STRLEN];

  int method_itype;
  int dire_itype;
  char direction[PAR_MAX_STRLEN];
  float coef;
  int index_is_min;
  int bdry_itype[2]; 
  float epsilon[2];  
} par_t;

int
par_read_from_file(char *par_fname, par_t *par);

int 
par_read_from_str(const char *str, par_t *par);

int
par_print(par_t *par);

#endif

