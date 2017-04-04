/*
 * Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
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
 *     * Neither the name of Sandia Corporation nor the names of its
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
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function(s) contained in this file:
 *      ex_put_elem_map()
 *****************************************************************************
 * This function outputs the elemental map.
 *****************************************************************************
 *  Variable Index:
 *      exoid             - The NetCDF ID of an already open NemesisI file.
 *      elem_mapi        - Vector of internal element IDs.
 *      elem_mapb        - Vector of border element IDs.
 *      processor        - The processor ID for which info is to be read.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <exodusII.h>     // for ex_err, exerrval, etc
#include <exodusII_int.h> // for EX_FATAL, DIM_NUM_BOR_ELEMS, etc
#include <netcdf.h>       // for NC_NOERR, nc_inq_varid, etc
#include <stddef.h>       // for size_t
#include <stdio.h>
#include <sys/types.h> // for int64_t

int ex_put_processor_elem_maps(int exoid, void_int *elem_mapi, void_int *elem_mapb, int processor)
{
  const char *func_name = "ex_put_processor_elem_maps";

  char    ftype[2];
  int     status, dimid, varid;
  size_t  start[1], count[1];
  int64_t varidx[2];
  int     nmstat;

  char errmsg[MAX_ERR_LENGTH];
  /*-----------------------------Execution begins-----------------------------*/

  ex_check_valid_file_id(exoid);

  /* Get the file type */
  if (ex_get_file_type(exoid, ftype) != EX_NOERR) {
    exerrval = EX_MSG;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: unable to find file type for file ID %d", exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Get the status of the internal element map */
  if ((status = nc_inq_varid(exoid, VAR_INT_E_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to find variable ID for \"%s\" in file ID %d",
             VAR_INT_E_STAT, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  if (ftype[0] == 's') {
    start[0] = processor;
  }
  else {
    start[0] = 0;
  }

  if ((status = nc_get_var1_int(exoid, varid, start, &nmstat)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get variable \"%s\" from file ID %d",
             VAR_INT_E_STAT, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  if (nmstat == 1) {
    /* get the index */
    if (ex_get_idx(exoid, VAR_ELEM_MAP_INT_IDX, varidx, processor) == -1) {
      exerrval = status;
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to find index variable, \"%s\", in file ID %d", VAR_ELEM_MAP_INT_IDX,
               exoid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    if (varidx[1] == -1) {
      /* Get the size of the internal element map */
      if ((status = nc_inq_dimid(exoid, DIM_NUM_INT_ELEMS, &dimid)) != NC_NOERR) {
        exerrval = status;
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to find dimension ID for \"%s\" in file ID %d", DIM_NUM_INT_ELEMS,
                 exoid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      if ((status = nc_inq_dimlen(exoid, dimid, count)) != NC_NOERR) {
        exerrval = status;
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to find length of dimension \"%s\" in file ID %d",
                 DIM_NUM_INT_ELEMS, exoid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      varidx[1] = count[0];
    }

    if ((status = nc_inq_varid(exoid, VAR_ELEM_MAP_INT, &varid)) != NC_NOERR) {
      exerrval = status;
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to find variable ID for \"%s\" in file ID %d",
               VAR_ELEM_MAP_INT, exoid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* Output the map */
    start[0] = varidx[0];
    count[0] = varidx[1] - varidx[0];
    if (ex_int64_status(exoid) & EX_MAPS_INT64_API) {
      status = nc_put_vara_longlong(exoid, varid, start, count, elem_mapi);
    }
    else {
      status = nc_put_vara_int(exoid, varid, start, count, elem_mapi);
    }
    if (status != NC_NOERR) {
      exerrval = status;
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to output variable \"%s\" in file ID %d",
               VAR_ELEM_MAP_INT, exoid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

  } /* End "if (nmstat == 1)" */

  /* Get the status of the border element map */
  if ((status = nc_inq_varid(exoid, VAR_BOR_E_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to find variable ID for \"%s\" in file ID %d",
             VAR_BOR_E_STAT, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  if (ftype[0] == 's') {
    start[0] = processor;
  }
  else {
    start[0] = 0;
  }

  if ((status = nc_get_var1_int(exoid, varid, start, &nmstat)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get status for \"%s\" from file %d",
             VAR_BOR_E_STAT, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  if (nmstat == 1) {
    /* get the index */
    if (ex_get_idx(exoid, VAR_ELEM_MAP_BOR_IDX, varidx, processor) == -1) {
      exerrval = status;
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to find index variable, \"%s\", in file ID %d", VAR_ELEM_MAP_BOR_IDX,
               exoid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    if (varidx[1] == -1) {
      /* Get the size of the border element map */
      if ((status = nc_inq_dimid(exoid, DIM_NUM_BOR_ELEMS, &dimid)) != NC_NOERR) {
        exerrval = status;
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to find dimension ID for \"%s\" in file ID %d", DIM_NUM_BOR_ELEMS,
                 exoid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      if ((status = nc_inq_dimlen(exoid, dimid, count)) != NC_NOERR) {
        exerrval = status;
        snprintf(errmsg, MAX_ERR_LENGTH,
                 "ERROR: failed to find length of dimension \"%s\" in file ID %d",
                 DIM_NUM_BOR_ELEMS, exoid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      varidx[1] = count[0];
    }

    if ((status = nc_inq_varid(exoid, VAR_ELEM_MAP_BOR, &varid)) != NC_NOERR) {
      exerrval = status;
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to find variable ID for \"%s\" in file ID %d",
               VAR_ELEM_MAP_BOR, exoid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* Output the map */
    start[0] = varidx[0];
    count[0] = varidx[1] - varidx[0];
    if (ex_int64_status(exoid) & EX_MAPS_INT64_API) {
      status = nc_put_vara_longlong(exoid, varid, start, count, elem_mapb);
    }
    else {
      status = nc_put_vara_int(exoid, varid, start, count, elem_mapb);
    }
    if (status != NC_NOERR) {
      exerrval = status;
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to output variable \"%s\" in file ID %d",
               VAR_ELEM_MAP_BOR, exoid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  } /* End "if (nmstat == 1)" */
  return (EX_NOERR);
}
