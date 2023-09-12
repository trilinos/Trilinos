/*
 * Copyright(C) 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * exggv - ex_get_glob_vars
 *
 * entry conditions -
 *   input parameters:
 *       int     exoid                   exodus file id
 *       int     time_step               time step number
 *       int     num_glob_vars           number of global vars in file
 *
 * exit conditions -
 *       float*  glob_var_vals           array of global variable values
 *
 * revision history -
 *
 *
 *****************************************************************************/

#include "exodusII.h"     // for ex_err, etc
#include "exodusII_int.h" // for ex__comp_ws, EX_FATAL, etc

/*!
 Internal function. Do not use in client code.
 */

int ex__get_glob_vars_multi_time(int exoid, int num_glob_vars, int beg_time_step, int end_time_step,
                                 void *glob_var_vals)
{
  int    varid;
  int    status;
  size_t start[2], count[2];
  char   errmsg[MAX_ERR_LENGTH];

  EX_FUNC_ENTER();
  if (ex__check_valid_file_id(exoid, __func__) == EX_FATAL) {
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* inquire previously defined variable */
  if ((status = nc_inq_varid(exoid, VAR_GLO_VAR, &varid)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "Warning: failed to locate global variables in file id %d",
             exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_WARN);
  }

  /* read values of global variables */
  start[0] = --beg_time_step;
  start[1] = 0;

  count[0] = end_time_step - beg_time_step;
  count[1] = num_glob_vars;

  if (ex__comp_ws(exoid) == 4) {
    status = nc_get_vara_float(exoid, varid, start, count, glob_var_vals);
  }
  else {
    status = nc_get_vara_double(exoid, varid, start, count, glob_var_vals);
  }

  if (status != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get global variable values from file id %d",
             exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }
  EX_FUNC_LEAVE(EX_NOERR);
}
