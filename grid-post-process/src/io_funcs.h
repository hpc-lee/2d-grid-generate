#ifndef IO_FUNCS_H
#define IO_FUNCS_H

#include "gd_t.h"
#include "algebra.h"

/*************************************************
 * structure
 *************************************************/
typedef struct
{
  int nx;
  int nz;
  float *var; // pointer to var
} io_quality_t;


/*************************************************
 * function prototype
 *************************************************/

int
init_io_quality(io_quality_t *io_quality, gd_t *gdcurv);

int
read_import_coord(gd_t *gdcurv, par_t *par);

int
gd_curv_coord_export(gd_t *gdcurv, par_t *par);

int
quality_export(io_quality_t *io_quality, gd_t *gdcurv, par_t *par, char *var_name);

int
io_get_nextline(FILE *fp, char *str, int length);

#endif
