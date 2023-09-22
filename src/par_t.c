#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "par_t.h"

/*
 * read from file
 */

int
par_read_from_file(char *par_fname,  par_t *par, int verbose)
{
  //
  // read whole file into str
  //
  FILE *fp = fopen(par_fname,"r");
  if (!fp) {
    fprintf(stderr,"Error: can't open par file: %s\n", par_fname);
    exit(1);
  }

  fseek(fp, 0, SEEK_END);
  long len = ftell(fp);

  fseek(fp, 0, SEEK_SET);
  char *str = (char*)malloc(len+1);
  fread(str, 1, len, fp);
  fclose(fp);

  // read from str
  par_read_from_str(str, par);

  return 0;
}

/*
 * funcs to get par from alread read in str
 */
int 
par_read_from_str(const char *str, par_t *par)
{
  int ierr = 0;

  // convert str to json
  cJSON *root = cJSON_Parse(str);
  if (NULL == root) {
    printf("Error at parsing json!\n");
    exit(1);
  }

  cJSON *item;
  cJSON *subitem, *thirditem;

  // default not check
  par->grid_check = 0;
  par->check_orth  = 0;
  par->check_jac   = 0;
  par->check_ratio = 0;
  par->check_step_xi  = 0;
  par->check_step_zt  = 0;
  par->check_smooth_xi = 0;
  par->check_smooth_zt = 0;
  if (item = cJSON_GetObjectItem(root, "check_orth")) {
    par->check_orth = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_jac")) {
    par->check_jac = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_ratio")) {
    par->check_ratio = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_step_xi")) {
    par->check_step_xi = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_step_zt")) {
    par->check_step_zt = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_smooth_xi")) {
    par->check_smooth_xi = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "check_smooth_zt")) {
    par->check_smooth_zt = item->valueint;
  }
  int check = par->check_orth + par->check_jac /
            + par->check_ratio + par->check_step_xi /
            + par->check_step_zt + par->check_smooth_xi /
            + par->check_smooth_zt; 
  if(check != 0)
  {
    par->grid_check = 1;
  }

  // default not strech
  par->flag_strech_xi = 0;
  par->flag_strech_zt = 0;
  if (item = cJSON_GetObjectItem(root, "flag_strech_xi")) {
    par->flag_strech_xi = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "strech_xi_coef")) {
    par->strech_xi_coef = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "flag_strech_zt")) {
    par->flag_strech_zt = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "strech_zt_coef")) {
    par->strech_zt_coef = item->valuedouble;
  }

  // default not intep
  par->sample_factor_xi = 1.0;
  par->sample_factor_zt = 1.0;
  if (item = cJSON_GetObjectItem(root, "flag_sample_xi")) {
    par->flag_sample_xi = item->valueint;
  }
  if (par->flag_sample_xi == 1) 
  {
    if (item = cJSON_GetObjectItem(root, "sample_factor_xi")) {
      par->sample_factor_xi = item->valuedouble;
      if((par->sample_factor_xi-1) < 0.0)
      {
        fprintf(stdout,"sample_factor_xi must >= 1\n");
        exit(1);
      }
    }
  }
  if (item = cJSON_GetObjectItem(root, "flag_sample_zt")) {
    par->flag_sample_zt = item->valueint;
  }
  if (par->flag_sample_zt == 1) 
  {
    if (item = cJSON_GetObjectItem(root, "sample_factor_zt")) {
      par->sample_factor_zt = item->valuedouble;
      if((par->sample_factor_zt-1) < 0.0)
      {
        fprintf(stdout,"sample_factor_zt must >= 1\n");
        exit(1);
      }
    }
  }

  if (item = cJSON_GetObjectItem(root, "geometry_input_file")) {
    sprintf(par->geometry_input_file, "%s", item->valuestring);
  }
  if (item = cJSON_GetObjectItem(root, "grid_export_dir")) {
    sprintf(par->grid_export_dir, "%s", item->valuestring);
  }

  if (item = cJSON_GetObjectItem(root, "grid_method")) {
    if (subitem = cJSON_GetObjectItem(item, "linear_TFI")) {
      par->method_itype = TFI;
    }
    if (subitem = cJSON_GetObjectItem(item, "hermite")) {
      par->method_itype = HERMITE;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        par->coef = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "direction")) {
        sprintf(par->direction, "%s", thirditem->valuestring);
        if(strcmp(par->direction,"x") == 0)
        {
          par->dire_itype = X_DIRE;
        }
        if(strcmp(par->direction,"z") == 0)
        {
          par->dire_itype = Z_DIRE;
        }
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "elli_diri")) {
      par->method_itype = ELLI_DIRI;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        par->coef = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "iter_err")) {
        par->iter_err = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "max_iter")) {
        par->max_iter = thirditem->valueint;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "first_dire")) {
        sprintf(par->first_dire, "%s", thirditem->valuestring);
        if(strcmp(par->first_dire,"x") == 0)
        {
          par->first_dire_itype = X_DIRE;
        }
        if(strcmp(par->first_dire,"z") == 0)
        {
          par->first_dire_itype = Z_DIRE;
        }
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "elli_higen")) {
      par->method_itype = ELLI_HIGEN;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        par->coef = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "distance")) {
        for (int i = 0; i < 4; i++) {
          par->distance[i] = cJSON_GetArrayItem(thirditem, i)->valuedouble;
        }
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "iter_err")) {
        par->iter_err = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "max_iter")) {
        par->max_iter = thirditem->valueint;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "first_dire")) {
        sprintf(par->first_dire, "%s", thirditem->valuestring);
        if(strcmp(par->first_dire,"x") == 0)
        {
          par->first_dire_itype = X_DIRE;
        }
        if(strcmp(par->first_dire,"z") == 0)
        {
          par->first_dire_itype = Z_DIRE;
        }
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "parabolic")) {
      par->method_itype = PARABOLIC;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        par->coef = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "o2i")) {
        par->o2i = thirditem->valueint;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "direction")) {
        sprintf(par->direction, "%s", thirditem->valuestring);
        if(strcmp(par->direction,"x") == 0)
        {
          par->dire_itype = X_DIRE;
        }
        if(strcmp(par->direction,"z") == 0)
        {
          par->dire_itype = Z_DIRE;
        }
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "hyperbolic")) {
      par->method_itype = HYPERBOLIC;
      if (thirditem = cJSON_GetObjectItem(subitem, "coef")) {
        par->coef = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "epsilon")) {
        par->epsilon = thirditem->valuedouble;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "bdry_type")) {
        par->bdry_itype = thirditem->valueint;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "o2i")) {
        par->o2i = thirditem->valueint;
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "direction")) {
        sprintf(par->direction, "%s", thirditem->valuestring);
        if(strcmp(par->direction,"x") == 0)
        {
          par->dire_itype = X_DIRE;
        }
        if(strcmp(par->direction,"z") == 0)
        {
          par->dire_itype = Z_DIRE;
        }
      }
      if (thirditem = cJSON_GetObjectItem(subitem, "step_input_file")) {
        sprintf(par->step_input_file, "%s", thirditem->valuestring);
      }
    }
  }

  return ierr;
}


int
par_print(par_t *par)
{    
  int ierr = 0;

  fprintf(stdout,"input geometry file is \n %s\n",par->geometry_input_file);
  fprintf(stdout,"export grid dir is \n %s\n",par->grid_export_dir);
  fprintf(stdout, "-------------------------------------------------------\n");
  if (par->grid_check == 1) {
    fprintf(stdout, "------- grid quality check-------\n");
  }
  if (par->check_orth == 1) {
    fprintf(stdout, "------- check grid orthogonality-------\n");
  }
  if (par->check_jac == 1) {
    fprintf(stdout, "------- check grid jacobi-------\n");
  }
  if (par->check_ratio == 1) {
    fprintf(stdout, "------- check grid ratio-------\n");
  }
  if (par->check_step_xi == 1) {
    fprintf(stdout, "------- check grid step xi direction-------\n");
  }
  if (par->check_step_zt == 1) {
    fprintf(stdout, "------- check grid step zt direction-------\n");
  }
  if (par->check_smooth_xi == 1) {
    fprintf(stdout, "------- check grid smooth xi direction-------\n");
  }
  if (par->check_smooth_zt == 1) {
    fprintf(stdout, "------- check grid smooth zt direction-------\n");
  }

  if (par->flag_sample_xi == 1) {
    fprintf(stdout,"sample grid xi direction factor is %f\n",par->sample_factor_xi);
  }
  if (par->flag_sample_zt == 1) {
    fprintf(stdout,"sample grid zt direction factor is %f\n",par->sample_factor_zt);
  }
  if(par->flag_strech_xi == 1) {
    fprintf(stdout, "------- strech xi and strech coef is %f-------\n",par->strech_xi_coef);
  }
  if(par->flag_strech_zt == 1) {
    fprintf(stdout, "------- strech zt and strech coef is %f-------\n",par->strech_zt_coef);
  }

  fprintf(stdout, "------- grid generate method-------\n");
  if(par->method_itype == TFI) {
    fprintf(stdout, "grid generate method is linear TFI\n");
  }
  if(par->method_itype == HERMITE) {
    fprintf(stdout, "grid generate method is unidirection hermite\n");
    fprintf(stdout, "hermite coef is %f\n", par->coef);
    if(par->dire_itype == X_DIRE) {
      fprintf(stdout, "grid generate direction is x\n");
    }
    if(par->dire_itype == Z_DIRE) {
      fprintf(stdout, "grid generate direction is z\n");
    }
  }
  if(par->method_itype == ELLI_DIRI) {
    fprintf(stdout, "grid generate method is elliptic_dirichlet\n");
    fprintf(stdout, "elli_diri coef is %f\n", par->coef);
    fprintf(stdout, "max_iteration is %d\n", par->max_iter);
    fprintf(stdout, "iter_error is %f\n", par->iter_err);
    if(par->first_dire_itype == X_DIRE) {
      fprintf(stdout, "grid generate first direction is x\n");
    }
    if(par->first_dire_itype == Z_DIRE) {
      fprintf(stdout, "grid generate first direction is z\n");
    }
  }
  if(par->method_itype == ELLI_HIGEN) {
    fprintf(stdout, "grid generate method is elliptic_hilgenstock\n");
    fprintf(stdout, "elli_higen coef is %f\n", par->coef);
    fprintf(stdout, "max_iteration is %d\n", par->max_iter);
    fprintf(stdout, "iter_error is %f\n", par->iter_err);
    if(par->first_dire_itype == X_DIRE) {
      fprintf(stdout, "grid generate first direction is x\n");
    }
    if(par->first_dire_itype == Z_DIRE) {
      fprintf(stdout, "grid generate first direction is z\n");
    }
    fprintf(stdout, "expect distance x1 bdry is %f\n",par->distance[0]);
    fprintf(stdout, "expect distance x2 bdry is %f\n",par->distance[1]);
    fprintf(stdout, "expect distance z1 bdry is %f\n",par->distance[2]);
    fprintf(stdout, "expect distance z2 bdry is %f\n",par->distance[3]);
  }
  if(par->method_itype == PARABOLIC) {
    fprintf(stdout, "grid generate method is parabolic\n");
    fprintf(stdout, "parabolic coef is %f\n", par->coef);
    if(par->dire_itype == X_DIRE) {
      fprintf(stdout, "grid generate direction is x\n");
    }
    if(par->dire_itype == Z_DIRE) {
      fprintf(stdout, "grid generate direction is z\n");
    }
    if(par->o2i == 1)
    {
      fprintf(stdout, "outer(bdry_2) to inner(bdry_1)\n");
    } else {
      fprintf(stdout, "inner(bdry_1) to outer(bdry_2)\n");
    }
  }

  if(par->method_itype == HYPERBOLIC) {
    fprintf(stdout, "grid generate method is hyperbolic\n");
    fprintf(stdout, "hyperbolic coef is %f\n", par->coef);
    if(par->bdry_itype == 1) {
      fprintf(stdout, "boundary type is floating boundary\n");
    }
    if(par->bdry_itype == 2) {
      fprintf(stdout, "boundary type is cartesian boundary\n");
    }
    if(par->dire_itype == X_DIRE) {
      fprintf(stdout, "grid generate direction is x\n");
    }
    if(par->dire_itype == Z_DIRE) {
      fprintf(stdout, "grid generate direction is z\n");
    }
    if(par->o2i == 1)
    {
      fprintf(stdout, "outer(bdry_2) to inner(bdry_1)\n");
    } else {
      fprintf(stdout, "inner(bdry_1) to outer(bdry_2)\n");
    }
    fprintf(stdout, "step file is  %s\n",par->step_input_file);
  }

  return ierr;
}
