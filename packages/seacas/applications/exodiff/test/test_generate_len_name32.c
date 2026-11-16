/*
 * Copyright(C) 2026 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 * Generate two identical databases whose `len_name` dimension is 32 and
 * which contain a 32-character global variable name (stored without a
 * terminating NUL).  Used to verify that exodiff can select that variable
 * by its full name in a command file.
 */
#include "exodusII.h"
#include "netcdf.h"
#include <stdio.h>
#include <string.h>

#define EX_TEST_FILENAME_A "len-name-32-a.exo"
#define EX_TEST_FILENAME_B "len-name-32-b.exo"

static const char *var_name = "LN012345678901234567890123456789";

static int create_file(const char *filename)
{
  int cpu_ws = 8;
  int io_ws  = 8;
  int exoid  = ex_create(filename, EX_CLOBBER, &cpu_ws, &io_ws);
  if (exoid < 0) {
    fprintf(stderr, "ERROR: could not create '%s'\n", filename);
    return 1;
  }

  /* Define `len_name` before exodus does so that it is 32 instead of the default. */
  int dimid = 0;
  nc_redef(exoid);
  if (nc_def_dim(exoid, "len_name", 32, &dimid) != NC_NOERR) {
    fprintf(stderr, "ERROR: could not define len_name in '%s'\n", filename);
    return 1;
  }
  nc_enddef(exoid);

  double x[]    = {0.0, 1.0, 1.0, 0.0};
  double y[]    = {0.0, 0.0, 1.0, 1.0};
  int    conn[] = {1, 2, 3, 4};
  ex_put_init(exoid, "len_name 32 test", 2, 4, 1, 1, 0, 0);
  ex_put_coord(exoid, x, y, NULL);
  ex_put_block(exoid, EX_ELEM_BLOCK, 1, "QUAD4", 1, 4, 0, 0, 0);
  ex_put_conn(exoid, EX_ELEM_BLOCK, 1, conn, NULL, NULL);
  ex_put_variable_param(exoid, EX_GLOBAL, 1);

  /* Write the 32-character name directly; exodus would truncate it to 31. */
  int    varid    = 0;
  size_t start[2] = {0, 0};
  size_t count[2] = {1, strlen(var_name)};
  if (nc_inq_varid(exoid, "name_glo_var", &varid) != NC_NOERR ||
      nc_put_vara_text(exoid, varid, start, count, var_name) != NC_NOERR) {
    fprintf(stderr, "ERROR: could not write global variable name in '%s'\n", filename);
    return 1;
  }

  int max_len = (int)strlen(var_name);
  nc_redef(exoid);
  nc_put_att_int(exoid, NC_GLOBAL, "maximum_name_length", NC_INT, 1, &max_len);
  nc_enddef(exoid);

  double time  = 0.0;
  double value = 1.0;
  ex_put_time(exoid, 1, &time);
  ex_put_var(exoid, 1, EX_GLOBAL, 1, 0, 1, &value);

  ex_close(exoid);
  return 0;
}

int main(void)
{
  ex_opts(EX_VERBOSE);
  int status = create_file(EX_TEST_FILENAME_A);
  status += create_file(EX_TEST_FILENAME_B);
  return status;
}
