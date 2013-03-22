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
/* Function(s) contained in this file:
 *
 *      ex_get_partial_node_set()
 *
 *****************************************************************************
 *
 *  Variable Index:
 *
 *      exoid               - The NetCDF ID of an already open NemesisI file.
 *      node_set_id        - ID of node set to read.
 *      start_node_num     - The starting index of the nodes to be read.
 *      num_nodes          - The number of nodes to read in.
 *      node_set_node_list - List of node IDs in node set.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "exodusII.h"
#include "exodusII_int.h"

/*
 * reads the node list for a single node set
 */

int ex_get_partial_node_set (int   exoid,
                       ex_entity_id node_set_id,
                       int64_t   start_node_num,
                       int64_t   num_nodes,
                       void_int  *node_set_node_list)
{
  int     dimid, node_list_id, node_set_id_ndx, status;
  size_t  num_nodes_in_set, start[1], count[1];
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* first check if any node sets are specified */

  if ((status = nc_inq_dimid (exoid, DIM_NUM_NS, &dimid))  != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Warning: no node sets defined in file id %d",
            exoid);
    ex_err("ex_get_partial_node_set",errmsg,exerrval);
    return (EX_WARN);
  }

  /* Lookup index of node set id in VAR_NS_IDS array */
  if ((node_set_id_ndx = ex_id_lkup(exoid, EX_NODE_SET, node_set_id)) < 0) {
    if (exerrval == EX_NULLENTITY) {
      sprintf(errmsg,
              "Warning: node set %"PRId64" is NULL in file id %d",
              node_set_id,exoid);
      ex_err("ex_get_partial_node_set",errmsg,EX_MSG);
      return (EX_WARN);
    } else {

      sprintf(errmsg,
              "Error: failed to locate node set %"PRId64" in %s in file id %d",
              node_set_id,VAR_NS_IDS,exoid);
      ex_err("ex_get_partial_node_set",errmsg,exerrval);
      return (EX_FATAL);
    }
  }

  /* inquire id's of previously defined dimensions and variables */
  if ((status = nc_inq_dimid (exoid, DIM_NUM_NOD_NS(node_set_id_ndx), &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
         "Error: failed to locate number of nodes in node set %"PRId64" in file id %d",
            node_set_id,exoid);
    ex_err("ex_get_partial_node_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_dimlen(exoid, dimid, &num_nodes_in_set)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get number of nodes in set %"PRId64" in file id %d",
            node_set_id, exoid);
    ex_err("ex_get_partial_node_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* Check input parameters for a valid range of numbers */
  if (start_node_num < 0 || start_node_num > num_nodes_in_set) {
    exerrval = EX_BADPARAM;
    sprintf(errmsg, "Error: Invalid input");
    ex_err("ex_get_partial_node_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  if (num_nodes < 0) {
    exerrval = EX_BADPARAM;
    sprintf(errmsg, "Error: Invalid number of nodes in nodes set!");
    ex_err("ex_get_partial_node_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* start_node_num now starts at 1, not 0 */
  if ((start_node_num + num_nodes - 1) > num_nodes_in_set) {
    exerrval = EX_BADPARAM;
    sprintf(errmsg, "Error: request larger than number of nodes in set!");
    ex_err("ex_get_partial_node_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_varid (exoid, VAR_NODE_NS(node_set_id_ndx), &node_list_id)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to locate node set %"PRId64" node list in file id %d",
            node_set_id,exoid);
    ex_err("ex_get_partial_node_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* read in the node list array */
  start[0] = --start_node_num;
  count[0] = num_nodes;

  if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
    status = nc_get_vara_longlong(exoid, node_list_id, start, count, node_set_node_list);
  } else {
    status = nc_get_vara_int(exoid, node_list_id, start, count, node_set_node_list);
  }

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get node set node list in file id %d",
            exoid);
    ex_err("ex_get_partial_node_set",errmsg,exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}
