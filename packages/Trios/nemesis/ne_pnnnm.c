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
*
* ne_pnnnm - ne_put_n_node_num_map
*
* entry conditions - 
*   input parameters:
*	int	neid			exodus file id
*	int	start_ent		first entry in node_map
*	int	num_ents		number of entries in node_map
*       int*    node_map                node numbering map
*
* exit conditions - 
*
* revision history - 
*
*
*****************************************************************************/

#include "exodusII.h"
#include "exodusII_int.h"
#include "ne_nemesisI.h"

/*
 * writes out the node numbering map to the database; allows node numbers
 * to be non-contiguous
 */

int ne_put_n_node_num_map (int  neid,
                           int  start_ent,
                           int  num_ents,
                           int *node_map)
{
  char  *func_name="ne_put_n_node_num_map";

  int numnodedim, dims[1], mapid, status;
  size_t num_nodes, start[1], count[1]; 
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* inquire id's of previously defined dimensions  */

  if ((status = nc_inq_dimid (neid, DIM_NUM_NODES, &numnodedim)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to locate number of nodes in file id %d",
            neid);
    ex_err(func_name,errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_dimlen(neid, numnodedim, &num_nodes)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get number of nodes in file id %d",
            neid);
    ex_err(func_name,errmsg,exerrval);
    return (EX_FATAL);
  }

  /* Check input parameters for a valid range of numbers */
  if (start_ent < 0 || start_ent > num_nodes) {
    exerrval = status;
    sprintf(errmsg,
            "Error: Invalid input to function %s!\n", func_name);
    ex_err(func_name,errmsg,exerrval);
    return (EX_FATAL);
  }

  if (num_ents < 0) {
    exerrval = status;
    sprintf(errmsg, "Error: Invalid number of entries in map!\n");
    ex_err(func_name,errmsg,exerrval);
    return (EX_FATAL);
  }

  /* start_ent now starts at 1, not 0 */
  if ((start_ent + num_ents - 1) > num_nodes) {
    exerrval = status;
    sprintf(errmsg, "Error: request range invalid!\n");
    ex_err(func_name,errmsg,exerrval);
    return (EX_FATAL);
  }

  /* check to see if node number map variable is already defined */
  if ((status = nc_inq_varid (neid, VAR_NODE_NUM_MAP, &mapid)) != NC_NOERR) {
    /* if not, then create it */

    /* put netcdf file into define mode  */
    if ((status = nc_redef (neid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to put file id %d into define mode",
              neid);
      ex_err(func_name,errmsg,exerrval);
      return (EX_FATAL);
    }

    /* create a variable array in which to store the node numbering map  */
    dims[0] = numnodedim;
    if ((status = nc_def_var(neid, VAR_NODE_NUM_MAP, NC_INT, 1, dims, &mapid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to create node numbering map array in file id %d",
              neid);
      ex_err(func_name,errmsg,exerrval);

      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* leave define mode  */
    if ((status = nc_enddef (neid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to complete definition in file id %d",
              neid);
      ex_err(func_name,errmsg,exerrval);
      return (EX_FATAL);
    }
  } 


  /* write out the node numbering map  */
  start[0] = --start_ent;
  count[0] = num_ents;
  status = nc_put_vara_int (neid, mapid, start, count, node_map);

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to store node numbering map in file id %d",
            neid);
    ex_err(func_name,errmsg,exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}
