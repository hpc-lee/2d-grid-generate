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

  if (item = cJSON_GetObjectItem(root, "number_of_grid_points_x")) {
    par->number_of_grid_points_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_grid_points_z")) {
    par->number_of_grid_points_z = item->valueint;
  }

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
  int check = par->check_orth + par->check_jac
            + par->check_ratio + par->check_step_xi
            + par->check_step_zt + par->check_smooth_xi
            + par->check_smooth_zt; 
  if(check != 0)
  {
    par->grid_check = 1;
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
  }

  return ierr;
}


int
par_print(par_t *par)
{    
  int ierr = 0;

  fprintf(stdout,"number of total gird points x is %d\n",par->number_of_grid_points_x);
  fprintf(stdout,"number of total gird points z is %d\n",par->number_of_grid_points_z);

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

  fprintf(stdout, "------- grid generate method-------\n");
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

  return ierr;
}
