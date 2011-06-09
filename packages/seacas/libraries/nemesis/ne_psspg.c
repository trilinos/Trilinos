/*
 * Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
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
 *     ne_put_ss_param_global()
 *****************************************************************************
 * This function outputs the global side-set parameters.
 *****************************************************************************
 *  Variable Index:
 *      neid            - The NetCDF ID of an already open NemesisI file.
 *      global_ids      - Pointer to a vector of global side-set IDs.
 *      side_cnts       - Pointer to a vector of global side counts in
 *                        each global side set.
 *      df_cnts         - Pointer to a vector of global distribution
 *                        factors in each global side set.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include <stdio.h>

#include <netcdf.h>

#include <exodusII.h>
#include <exodusII_int.h>

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

int ne_put_ss_param_global(int neid,
                           int *global_ids,
                           int *side_cnts,
                           int *df_cnts
                           )
{
  char   *func_name="ne_put_ss_param_global";
  int     varid;

  int     status;
  char    errmsg[MAX_ERR_LENGTH];
/*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* Get the variable ID for the vector of global side set IDs */
  if ((status = nc_inq_varid(neid, VAR_SS_IDS_GLOBAL, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_SS_IDS_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Output the vector of global side set IDs */
  status = nc_put_var_int(neid, varid, global_ids);
  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to output variable \"%s\" to file ID %d",
            VAR_SS_IDS_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Get the variable ID for the vector of global side-set side counts */
  if ((status = nc_inq_varid(neid, VAR_SS_SIDE_CNT_GLOBAL, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_SS_SIDE_CNT_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Output the vector of global side counts in each global side set */
  status = nc_put_var_int(neid, varid, side_cnts);
  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to put variable \"%s\" in file ID %d",
            VAR_SS_SIDE_CNT_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Get the variable ID for the number of dist. factors in each side set */
  if ((status = nc_inq_varid(neid, VAR_SS_DF_CNT_GLOBAL, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_SS_DF_CNT_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Output the vector of dist. factor counts */
  status = nc_put_var_int(neid, varid, df_cnts);
  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to output variable \"%s\" in file ID %d",
            VAR_SS_DF_CNT_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  return (EX_NOERR);
}
