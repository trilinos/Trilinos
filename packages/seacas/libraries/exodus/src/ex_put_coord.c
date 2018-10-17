/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
/*****************************************************************************
 *
 * expcor - ex_put_coord
 *
 * entry conditions -
 *   input parameters:
 *       int     exoid                   exodus file id
 *       float*  x_coord                 X coord array
 *       float*  y_coord                 y coord array
 *       float*  z_coord                 z coord array
 *
 * exit conditions -
 *
 * revision history -
 *
 *
 *****************************************************************************/

#include "exodusII.h"     // for ex_err, etc
#include "exodusII_int.h" // for EX_FATAL, ex_comp_ws, etc
#include "netcdf.h"       // for NC_NOERR, nc_inq_varid, etc
#include <stddef.h>       // for size_t
#include <stdio.h>

/*!
The function ex_put_coord() writes the nodal coordinates of the nodes
in the model. The function ex_put_init() must be invoked before this
call is made.

Because the coordinates are floating point values, the application
code must declare the arrays passed to be the appropriate type
(\e float or \e double) to match the compute word size passed in
ex_create() or ex_open().

\return In case of an error, ex_put_coord() returns a negative
number; a warning will return a positive number.
Possible causes of errors include:
  -  data file not properly opened with call to ex_create() or ex_open()
  -  data file opened for read only.
  -  data file not initialized properly with call to ex_put_init().

\param[in] exoid   exodus file ID returned from a previous call to ex_create()
or ex_open().
\param[in] x_coor  The X-coordinates of the nodes. If this is NULL, the
                   X-coordinates will not be written.
\param[in] y_coor  The Y-coordinates of the nodes. These are stored only if num_dim
                   > 1; otherwise, pass in NULL. If this is NULL, the
                   Y-coordinates will not be written.
\param[in] z_coor  The Z-coordinates of the nodes. These are stored only if num_dim
                   > 2; otherwise, pass in NULL. If this is NULL, the
                   Z-coordinates will not be written.

The following will write the nodal coordinates to an open
exodus file :

~~~{.c}
int error, exoid;

// \comment{if file opened with compute word size of sizeof(float)}
float x[8], y[8], z[8];

// \comment{write nodal coordinates values to database}
x[0] = 0.0; y[0] = 0.0; z[0] = 0.0;
x[1] = 0.0; y[1] = 0.0; z[1] = 1.0;
x[2] = 1.0; y[2] = 0.0; z[2] = 1.0;
x[3] = 1.0; y[3] = 0.0; z[3] = 0.0;
x[4] = 0.0; y[4] = 1.0; z[4] = 0.0;
x[5] = 0.0; y[5] = 1.0; z[5] = 1.0;
x[6] = 1.0; y[6] = 1.0; z[6] = 1.0;
x[7] = 1.0; y[7] = 1.0; z[7] = 0.0;

error = ex_put_coord(exoid, x, y, z);

\comment{Do the same as the previous call in three separate calls}
error = ex_put_coord(exoid, x,    NULL, NULL);
error = ex_put_coord(exoid, NULL, y,    NULL);
error = ex_put_coord(exoid, NULL, NULL, z);

~~~

 */

int ex_put_coord(int exoid, const void *x_coor, const void *y_coor, const void *z_coor)
{
  int status;
  int coordid;
  int coordidx, coordidy, coordidz;

  int    numnoddim, ndimdim;
  size_t num_nod, num_dim, i;
  char   errmsg[MAX_ERR_LENGTH];

  EX_FUNC_ENTER();
  ex_check_valid_file_id(exoid, __func__);

  /* inquire id's of previously defined dimensions  */

  if ((status = nc_inq_dimid(exoid, DIM_NUM_NODES, &numnoddim)) != NC_NOERR) {
    /* If not found, then this file is storing 0 nodes.
       Return immediately */
    EX_FUNC_LEAVE(EX_NOERR);
  }

  if ((status = nc_inq_dimlen(exoid, numnoddim, &num_nod)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: inquire failed to return number of nodes in file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  if ((status = nc_inq_dimid(exoid, DIM_NUM_DIM, &ndimdim)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to locate number of dimensions in file id %d",
             exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  if ((status = nc_inq_dimlen(exoid, ndimdim, &num_dim)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get number of dimensions in file id %d",
             exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* write out the coordinates  */
  if ((status = nc_inq_varid(exoid, VAR_COORD_X, &coordidx)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to locate x nodal coordinates in file id %d",
             exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  if (num_dim > 1) {
    if ((status = nc_inq_varid(exoid, VAR_COORD_Y, &coordidy)) != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to locate y nodal coordinates in file id %d",
               exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }
  }
  else {
    coordidy = 0;
  }
  if (num_dim > 2) {
    if ((status = nc_inq_varid(exoid, VAR_COORD_Z, &coordidz)) != NC_NOERR) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to locate z nodal coordinates in file id %d",
               exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }
  }
  else {
    coordidz = 0;
  }

  /* write out the coordinates  */
  for (i = 0; i < num_dim; i++) {
    const void *coor  = NULL;
    char *      which = NULL;

    if (i == 0) {
      coor    = x_coor;
      which   = "X";
      coordid = coordidx;
    }
    else if (i == 1) {
      coor    = y_coor;
      which   = "Y";
      coordid = coordidy;
    }
    else if (i == 2) {
      coor    = z_coor;
      which   = "Z";
      coordid = coordidz;
    }

    if (coor != NULL && coordid != 0) {
      if (ex_comp_ws(exoid) == 4) {
        status = nc_put_var_float(exoid, coordid, coor);
      }
      else {
        status = nc_put_var_double(exoid, coordid, coor);
      }

      if (status != NC_NOERR) {
        snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to put %s coord array in file id %d", which,
                 exoid);
        ex_err_fn(exoid, __func__, errmsg, status);
        EX_FUNC_LEAVE(EX_FATAL);
      }
    }
  }
  EX_FUNC_LEAVE(EX_NOERR);
}
