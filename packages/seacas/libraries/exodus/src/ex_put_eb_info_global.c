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
 *     ex_put_eb_info_global()
 *****************************************************************************
 * This function outputs the global element block IDs of all the element
 * blocks associated with a geometry.
 *****************************************************************************
 *  Variable Index:
 *      exoid            - The NetCDF ID of an already open NemesisI file.
 *      el_blk_ids      - Pointer to vector of global element block IDs.
 *      el_blk_cnts     - Pointer to vector of global element block counts.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include <stdio.h>

#include <netcdf.h>

#include <exodusII.h>
#include <exodusII_int.h>

int ex_put_eb_info_global(int exoid, void_int *el_blk_ids, void_int *el_blk_cnts)
{
  const char *func_name = "ex_put_eb_info_global";

  int  varid, status;
  char errmsg[MAX_ERR_LENGTH];

  ex_check_valid_file_id(exoid);

  /* Find the variable ID for the element block IDs */
  if ((status = nc_inq_varid(exoid, VAR_ELBLK_IDS_GLOBAL, &varid)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to find variable ID for \"%s\" in file ID %d",
             VAR_ELBLK_IDS_GLOBAL, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Output the global element block IDs */
  if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
    status = nc_put_var_longlong(exoid, varid, el_blk_ids);
  }
  else {
    status = nc_put_var_int(exoid, varid, el_blk_ids);
  }
  if (status != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to output variable \"%s\" in file ID %d",
             VAR_ELBLK_IDS_GLOBAL, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Find the variable ID for the element block counts */
  if ((status = nc_inq_varid(exoid, VAR_ELBLK_CNT_GLOBAL, &varid)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to find variable ID for \"%s\" in file ID %d",
             VAR_ELBLK_CNT_GLOBAL, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Output the global element block counts */
  if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
    status = nc_put_var_longlong(exoid, varid, el_blk_cnts);
  }
  else {
    status = nc_put_var_int(exoid, varid, el_blk_cnts);
  }
  if (status != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to output variable \"%s\" in file ID %d",
             VAR_ELBLK_CNT_GLOBAL, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  return (EX_NOERR);
}
