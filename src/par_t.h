#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cJSON.h"
#include "par_t.h"

#define PAR_MAX_STRLEN 1000
#define TFI 1
#define HERMITE 2
#define ELLI_DIRI 3
#define ELLI_HIGEN 4
#define PARABOLIC 5
#define HYPERBOLIC 6

typedef struct{

  int grid_check;
  int check_orth;
  int check_jac;
  int check_ratio;
  int check_step_x;
  int check_step_z;
  int check_smooth_x;
  int check_smooth_z;

  int flag_strech_x;
  int flag_strech_z;
  float strech_x_coef;
  float strech_z_coef;

  int flag_sample_x;
  int flag_sample_z;
  float sample_factor_x;
  float sample_factor_z;

  char geometry_input_file[PAR_MAX_STRLEN];
  char step_input_file[PAR_MAX_STRLEN];
  char grid_export_dir[PAR_MAX_STRLEN];
  
  // TFI hermite elliptic-dirichlet
  // elliptic-hilgenstock
  // parabolic hyperbolic
  int  method_itype;

  float coef;
  int o2i;
  int num_layers; 
} par_t;

int
par_read_from_file(char *par_fname, par_t *par, int verbose);

int 
par_read_from_str(const char *str, par_t *par);

int
par_print(par_t *par);

#endif

