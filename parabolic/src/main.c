#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>

#include "par_t.h"
#include "constants.h"
#include "gd_t.h"
#include "io_funcs.h"
#include "quality_check.h"
#include "parabolic.h"

int main(int argc, char** argv)
{
  int verbose;
  char *par_fname;
  char err_message[CONST_MAX_STRLEN];

  //-------------------------------------------------------------------------------
  // get commond-line argument
  //-------------------------------------------------------------------------------

  // argc checking
  if (argc < 3) {
    fprintf(stdout,"usage: main_grid_2d <par_file> <opt: verbose>\n");
    exit(1);
  }

  par_fname = argv[1];

  if (argc >= 3) {
    verbose = atoi(argv[2]); // verbose number
    fprintf(stdout,"verbose=%d\n", verbose); fflush(stdout);
  }

  fprintf(stdout,"par file =  %s\n", par_fname); fflush(stdout);

  // read par
  par_t *par = (par_t *) malloc(sizeof(par_t));

  par_read_from_file(par_fname, par, verbose);

  if (verbose>0) par_print(par);
   
  // generate grid 
  gd_t *gdcurv = (gd_t *) malloc(sizeof(gd_t));

  time_t t_start = time(NULL);
  grid_init_set(gdcurv,par);

  para_gene(gdcurv,par);
  time_t t_end = time(NULL);

  fprintf(stdout,"\n************************************\n");
  fprintf(stdout,"grid generate running time is :%f s \n", difftime(t_end,t_start));
  fprintf(stdout,"************************************\n \n");

  fprintf(stdout,"export coord to file ... \n");
  gd_curv_coord_export(gdcurv, par->grid_export_dir);

  // grid quality check and export quality data
  io_quality_t *io_quality = (io_quality_t *) malloc(sizeof(io_quality_t));
  if(par->grid_check == 1)
  {
    fprintf(stdout,"****************************************************** \n");
    fprintf(stdout,"***** grid quality check and export quality data ***** \n");
    fprintf(stdout,"****************************************************** \n");
    init_io_quality(io_quality,gdcurv);
    grid_quality_check(io_quality,gdcurv,par);
  }
  return 0;
}

