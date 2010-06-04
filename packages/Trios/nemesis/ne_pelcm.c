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
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/* Function(s) contained in this file:
 *     ne_put_elem_cmap()
 ****************************************************************************
 * The function outputs an elemental communication map.
 ****************************************************************************
 * Variable Index:
 *      neid            - The NetCDF ID of an already open NemesisI file.
 *      map_id          - The ID of the elementa communication map to write.
 *      elem_ids        - Pointer to vector of element IDs to output.
 *      side_ids        - Pointer to vector of side IDs for each element
 *                        in "elem_ids".
 *      proc_ids        - Pointer to vector of processor IDs for each
 *                        element in "elem_ids".
 *      processor       - The processor the file being read was written for.
 */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

#include <stdio.h>

#include <netcdf.h>

#include <exodusII.h>
#include <exodusII_int.h>

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

int ne_put_elem_cmap(int  neid,
                     int  map_id,
                     int *elem_ids,
                     int *side_ids,
                     int *proc_ids,
                     int  processor
                     )
{
  char   *func_name="ne_put_elem_cmap";

  int     map_idx, varid, dimid, status;
  size_t  start[1], count[1], ret_val;
  size_t  varidx[2];
  int     value;

  char    errmsg[MAX_ERR_LENGTH];
/*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* get the index for the comm map information variables */
  if (ne_get_idx(neid, VAR_E_COMM_INFO_IDX, varidx, processor) == -1) {
    sprintf(errmsg,
            "Error: failed to find index variable, \"%s\", in file ID %d",
            VAR_E_COMM_INFO_IDX, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Get the index for this map_id */
  if ((map_idx=ne_id_lkup(neid, VAR_E_COMM_IDS, varidx, map_id)) == -1) {
    sprintf(errmsg,
            "Error: failed to find index for variable \"%s\" in file ID %d",
            VAR_E_COMM_IDS, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /*
   * Find out if this is a NULL comm map by checking it's entry in
   * the status vector.
   */
  if ((status = nc_inq_varid(neid, VAR_E_COMM_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_E_COMM_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  start[0] = map_idx;
  if ((status = nc_get_var1_int(neid, varid, start, &value)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get variable \"%s\" from file ID %d",
            VAR_E_COMM_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (value == 0) return(EX_NOERR);   /* NULL set */

  /* now I need to get the comm map data index */
  if (ne_get_idx(neid, VAR_E_COMM_DATA_IDX, varidx, map_idx) == -1) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find index variable, \"%s\", in file ID %d",
            VAR_E_COMM_DATA_IDX, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* check if I need to get the dimension of the cmap data */
  if (varidx[1] == -1) {
    /* Get the size of the comm maps */
    if ((status = nc_inq_dimid(neid, DIM_ECNT_CMAP, &dimid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to get dimension ID for \"%s\" in file ID %d",
              DIM_ECNT_CMAP, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    if ((status = nc_inq_dimlen(neid, dimid, &ret_val)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to get length of dimension \"%s\" in file ID %d",
              DIM_ECNT_CMAP, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    varidx[1] = ret_val;
  } /* "if (varidx[1]==-1)" */

  start[0] = varidx[0];
  count[0] = varidx[1] - varidx[0];

  /* Output the element IDs for this comm map */
  if ((status = nc_inq_varid(neid, VAR_E_COMM_EIDS, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_E_COMM_EIDS, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  status = nc_put_vara_int(neid, varid, start, count, elem_ids);
  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to output vector \"%s\" in file ID %d",
            VAR_E_COMM_EIDS, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Output the processor IDs for this map */
  if ((status = nc_inq_varid(neid, VAR_E_COMM_PROC, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_E_COMM_PROC, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  status = nc_put_vara_int(neid, varid, start, count, proc_ids);
  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to output variable \"%s\" in file ID %d",
            VAR_E_COMM_PROC, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_varid(neid, VAR_E_COMM_SIDS, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_E_COMM_SIDS, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  status = nc_put_vara_int(neid, varid, start, count, side_ids);
  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to ouput variable \"%s\" in file ID %d",
            VAR_E_COMM_SIDS, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}
