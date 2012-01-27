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
* ne_gnnnm - ne_get_n_node_num_map
*
* environment - UNIX
*
* entry conditions -
*   input parameters:
*	int	neid			nemesis file id
*
* exit conditions -
*	int*	node_map		node numbering map array
*
* revision history -
*
*
*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "exodusII.h"
#include "exodusII_int.h"

#include "ne_nemesisI.h"

/*
 *  reads the node numbering map from the database
 */

int ne_get_n_node_num_map (int  neid,
                           int  start_ent,
                           int  num_ents,
                           int *node_map)
{
  int     numnodedim, mapid, i, status;
  size_t  num_nodes,  start[1], count[1];
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* inquire id's of previously defined dimensions and variables  */

  if ((status = nc_inq_dimid (neid, DIM_NUM_NODES, &numnodedim)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to locate number of nodes in file id %d",
            neid);
    ex_err("ne_get_n_node_num_map",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_dimlen(neid, numnodedim, &num_nodes)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get number of nodes in file id %d",
            neid);
    ex_err("ne_get_n_node_num_map",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* Check input parameters for a valid range of numbers */
  if (start_ent < 0 || start_ent > num_nodes) {
    fprintf(stderr, "ERROR: Invalid input to function"
                    "  ne_get_n_node_num_map!\n");
    return -1;
  }

  if (num_ents < 0) {
    fprintf(stderr, "ERROR: Invalid number of entries in map!\n");
    return -1;
  }

  /* start_ent now starts at 1, not 0 */
  if ((start_ent + num_ents - 1) > num_nodes) {
    fprintf(stderr, "ERROR: request range invalid!\n");
    return -1;
  }

  if ((status = nc_inq_varid (neid, VAR_NODE_NUM_MAP, &mapid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
  "Warning: node numbering map not stored in file id %d; returning default map",
            neid);
    ex_err("ne_get_n_node_num_map",errmsg,exerrval);

    /* generate default map of 1..n, where n is num_nodes */
    for (i=0; i < num_ents; i++)
      node_map[i] = start_ent+i;

    return (EX_WARN);
  }

  /* read in the node numbering map  */
  start[0] = --start_ent;
  count[0] = num_ents;

  status = nc_get_vara_int(neid, mapid, start, count, node_map);

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get node numbering map in file id %d",
            neid);
    ex_err("ne_get_n_node_num_map",errmsg,exerrval);
    return (EX_FATAL);
  }
  return(EX_NOERR);
}
