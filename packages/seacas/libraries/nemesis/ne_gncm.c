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
/* Function(s) in this file:
 *     ne_get_node_cmap()
 *
 ****************************************************************************
 * Variable Index:
 *      neid            - The NetCDF ID of an already open NemesisI file.
 *      map_id          - The ID of the nodal communication map to retrieve.
 *      node_ids        - Pointer to vector for retrieval of FEM node IDs
 *                        that make up this communication map.
 *      proc_ids        - Pointer to vector for retrieval of the processors
 *                        associated with each of the nodes in this nodal
 *                        communication map.
 *      processor         - The processor the file being read was written for.
 */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <netcdf.h>

#include "exodusII.h"
#include "exodusII_int.h"

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

int ne_get_node_cmap(int  neid,
                     int  map_id,
                     int *node_ids,
                     int *proc_ids,
                     int  processor
                     )
{
  char   *func_name="ne_get_node_cmap";

  int     map_idx, dimid, varid[2], status;
  size_t  start[1], count[1];
  int64_t varidx[2];

  char    errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* get the cmap information variables index */
  if (ne_get_idx(neid, VAR_N_COMM_INFO_IDX, varidx, processor) == -1) {
    sprintf(errmsg,
            "Error: failed to find index variable, \"%s\", in file ID %d",
            VAR_N_COMM_INFO_IDX, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /*
   * no need to check if the second index is -1 that is handled
   * in ne_id_lkup, where the dimension must be looked up anyways
   */

  /* Get the index of the nodal comm map with the given ID */
  if ((map_idx=ne_id_lkup(neid, VAR_N_COMM_IDS, varidx, map_id)) < 0) {
    exerrval = EX_MSG;
    sprintf(errmsg,
            "Error: failed to find nodal comm map with ID %d in file ID %d",
            map_id, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* get the cmap data variables index for this map */
  if (ne_get_idx(neid, VAR_N_COMM_DATA_IDX, varidx, map_idx) == -1) {
    sprintf(errmsg,
            "Error: failed to find index variable, \"%s\", in file ID %d",
            VAR_N_COMM_DATA_IDX, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  if (varidx[1] == -1) {
    /* Get the dimension of this nodal communication map */
    if ((status = nc_inq_dimid(neid, DIM_NCNT_CMAP, &dimid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find dimension ID for \"%s\" in file ID %d",
              DIM_NCNT_CMAP, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    if ((status = nc_inq_dimlen(neid, dimid,count)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find length of dimension \"%s\" in file ID %d",
              DIM_NCNT_CMAP, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
    varidx[1] = count[0];
  }

  /* Get the variable ID for the nodal comm map node IDs */
  if ((status = nc_inq_varid(neid, VAR_N_COMM_NIDS, &varid[0])) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_N_COMM_NIDS, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Get the variable ID for the nodal comm map processor IDs */
  if ((status = nc_inq_varid(neid, VAR_N_COMM_PROC, &varid[1])) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_N_COMM_PROC, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Get the nodal comm map node IDs */
  start[0] = varidx[0];
  count[0] = varidx[1] - varidx[0];
  status = nc_get_vara_int(neid, varid[0], start, count, node_ids);

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get variable \"%s\" from file ID %d",
            VAR_N_COMM_NIDS, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Get the nodal comm map processor IDs */
  status = nc_get_vara_int(neid, varid[1], start, count, proc_ids);

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get variable \"%s\" from file ID %d",
            VAR_N_COMM_PROC, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  return (EX_NOERR);
}
