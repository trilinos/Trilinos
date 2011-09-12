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

/*****************************************************************************
 * This function outputs the nodal map.
 *****************************************************************************
 * Variable Index:
 *      neid            - The NetCDF ID of an already open NemesisI file.
 *      node_mapi       - Pointer to vector containing the internal FEM
 *                        nodal IDs.
 *      node_mapb       - Pointer to vector containing the border FEM
 *                        nodal IDs.
 *      node_mape       - Pointer to vector containing the external FEM
 *                        nodal IDs.
 *      proc_id         - The processor the file is being written for.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>

#include <netcdf.h>

#include <exodusII.h>
#include <exodusII_int.h>

#include "ne_nemesisI.h"
#include "ne_nemesisI_int.h"

int ne_put_node_map(int  neid,
                    int *node_mapi,
                    int *node_mapb,
                    int *node_mape,
                    int  proc_id
                    )
{
  char  *func_name="ne_put_node_map";

  int     status, varid, dimid;
  char    ftype[2];
  size_t  start[1], count[1];
  int64_t varidx[2];
  int  nmstat;

  char   errmsg[MAX_ERR_LENGTH];
  /*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* Get the file type */
  if (ne_get_file_type(neid, ftype) != EX_NOERR) {
    exerrval = EX_MSG;
    sprintf(errmsg,
            "Error: unable to find file type for file ID %d",
            neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Get the status of this node map */
  if ((status = nc_inq_varid(neid, VAR_INT_N_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" from file ID %d",
            VAR_INT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (ftype[0] == 'p')
    start[0] = 0;
  else
    start[0] = proc_id;

  if ((status = nc_get_var1_int(neid, varid, start, &nmstat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get status for \"%s\" from file %d",
            VAR_INT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Write out the internal node-number map */
  if (nmstat == 1) {
    /* get the index */
    if (ne_get_idx(neid, VAR_NODE_MAP_INT_IDX, varidx, proc_id) == -1) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find index variable, \"%s\", in file ID %d",
              VAR_NODE_MAP_INT_IDX, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    /* check if I need to get the dimension */
    if (varidx[1] == -1) {
      if ((status = nc_inq_dimid(neid, DIM_NUM_INT_NODES, &dimid)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to find dimension ID for \"%s\" in file ID %d",
                DIM_NUM_INT_NODES, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      if ((status = nc_inq_dimlen(neid, dimid, count)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
                DIM_NUM_INT_NODES, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      varidx[1] = count[0];
    }

    if ((status = nc_inq_varid(neid, VAR_NODE_MAP_INT, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable ID for \"%s\" in file ID %d",
              VAR_NODE_MAP_INT, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    start[0] = varidx[0];
    count[0] = varidx[1] - varidx[0];
    status = nc_put_vara_int(neid, varid, start, count, node_mapi);
    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output variable \"%s\" in file ID %d",
              VAR_NODE_MAP_INT, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  } /* End "if (nmstat == 1)" */

  /* Get the status of this node map */
  if ((status = nc_inq_varid(neid, VAR_BOR_N_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" from file ID %d",
            VAR_BOR_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (ftype[0] == 'p')
    start[0] = 0;
  else
    start[0] = proc_id;

  if ((status = nc_get_var1_int(neid, varid, start, &nmstat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get status for \"%s\" from file %d",
            VAR_BOR_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (nmstat == 1) {
    /* Write out the border node-number map */
    /* get the index */
    if (ne_get_idx(neid, VAR_NODE_MAP_BOR_IDX, varidx, proc_id) == -1) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find index variable, \"%s\", in file ID %d",
              VAR_NODE_MAP_BOR_IDX, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    /* check if I need to get the dimension */
    if (varidx[1] == -1) {
      if ((status = nc_inq_dimid(neid, DIM_NUM_BOR_NODES, &dimid)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to find dimension ID for \"%s\" in file ID %d",
                DIM_NUM_BOR_NODES, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      if ((status = nc_inq_dimlen(neid, dimid, count)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
                DIM_NUM_BOR_NODES, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      varidx[1] = count[0];
    }

    if ((status = nc_inq_varid(neid, VAR_NODE_MAP_BOR, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable ID for \"%s\" in file ID %d",
              VAR_NODE_MAP_BOR, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    /* Output the map */
    start[0] = varidx[0];
    count[0] = varidx[1] - varidx[0];
    status = nc_put_vara_int(neid, varid, start, count, node_mapb);
    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output variable \"%s\" in file ID %d",
              VAR_NODE_MAP_BOR, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  } /* End "if (nmstat == 1)" */

  /* Get the status of this node map */
  if ((status = nc_inq_varid(neid, VAR_EXT_N_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" from file ID %d",
            VAR_EXT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (ftype[0] == 'p')
    start[0] = 0;
  else
    start[0] = proc_id;

  if ((status = nc_get_var1_int(neid, varid, start, &nmstat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get status for \"%s\" from file %d",
            VAR_EXT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (nmstat == 1) {
    /* Write out the external node-number map */
    if (ne_get_idx(neid, VAR_NODE_MAP_EXT_IDX, varidx, proc_id) == -1) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find index variable, \"%s\", in file ID %d",
              VAR_NODE_MAP_EXT_IDX, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    /* check if I need to get the dimension */
    if (varidx[1] == -1) {
      if ((status = nc_inq_dimid(neid, DIM_NUM_EXT_NODES, &dimid)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to find dimension ID for \"%s\" in file ID %d",
                DIM_NUM_EXT_NODES, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      if ((status = nc_inq_dimlen(neid, dimid, count)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
                DIM_NUM_EXT_NODES, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }
      varidx[1] = count[0];
    }

    if ((status = nc_inq_varid(neid, VAR_NODE_MAP_EXT, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable ID for \"%s\" in file ID %d",
              VAR_NODE_MAP_EXT, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    /* Output the map */
    start[0] = varidx[0];
    count[0] = varidx[1] - varidx[0];
    status = nc_put_vara_int(neid, varid, start, count, node_mape);
    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output variable \"%s\" in file ID %d",
              VAR_NODE_MAP_EXT, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  } /* End "if (nmstat == 1)" */
  return (EX_NOERR);
}
