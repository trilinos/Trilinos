/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
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
* exgpem - ex_get_partial_elem_map
*
* entry conditions - 
*   input parameters:
*       int     exoid                   exodus file id
*       int     map_id                  element map id
*       int     ent_start               first entry in map
*       int     ent_count               number of entries in map
*
* exit conditions - 
*       int*    elem_map                element map
*
* revision history - 
*
*
*****************************************************************************/

#include <stdlib.h>
#include "exodusII.h"
#include "exodusII_int.h"

/*
 * reads the element map with specified ID
 */

int ex_get_partial_num_map (int  exoid,
			    ex_entity_type map_type,
			    ex_entity_id  map_id,
			    int64_t  ent_start,
			    int64_t  ent_count, 
			    void_int *map)
{
  int dimid, var_id, id_ndx, status;
  size_t num_mobj, start[1], count[1]; 
  char errmsg[MAX_ERR_LENGTH];
  const char* dim_map_size;
  const char* dim_num_maps;

  switch (map_type) {
  case EX_NODE_MAP:
    dim_map_size = DIM_NUM_NODES;
    dim_num_maps = DIM_NUM_NM;
    break;
  case EX_EDGE_MAP:
    dim_map_size = DIM_NUM_EDGE;
    dim_num_maps = DIM_NUM_EDM;
    break;
  case EX_FACE_MAP:
    dim_map_size = DIM_NUM_FACE;
    dim_num_maps = DIM_NUM_FAM;
    break;
  case EX_ELEM_MAP:
    dim_map_size = DIM_NUM_ELEM;
    dim_num_maps = DIM_NUM_EM;
    break;
  default:
    exerrval = EX_BADPARAM;
    sprintf( errmsg, "Bad map type (%d) specified", map_type );
    ex_err( "ex_get_partial_num_map", errmsg, exerrval );
    return (EX_FATAL);
  }

  exerrval = 0; /* clear error code */

  /* See if file contains any elements...*/
  if (nc_inq_dimid (exoid, dim_map_size, &dimid) != NC_NOERR) {
    return (EX_NOERR);
  }

  if ((status = nc_inq_dimlen(exoid, dimid, &num_mobj)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get number of mesh objects in file id %d", exoid);
    ex_err("ex_get_partial_num_map",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* Check input parameters for a valid range of numbers */
  if (ent_start <= 0 || ent_start > num_mobj) {
    exerrval = EX_FATAL;
    sprintf(errmsg,
	    "Error: start count is invalid in file id %d",
	    exoid);
    ex_err("ex_get_partial_num_map",errmsg,exerrval);
    return (EX_FATAL);
  }

  if (ent_count < 0) {
    exerrval = EX_FATAL;
    sprintf(errmsg,
	    "Error: Invalid count value in file id %d",
	    exoid);
    ex_err("ex_get_partial_num_map",errmsg,exerrval);
    return (EX_FATAL);
  }

  if (ent_start+ent_count-1 > num_mobj) {
    exerrval = EX_FATAL;
    sprintf(errmsg,
	    "Error: start+count-1 is larger than element count in file id %d",
	    exoid);
    ex_err("ex_get_partial_num_map",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* first check if any maps have been defined */
  if ((status = nc_inq_dimid (exoid, dim_num_maps, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Warning: no %ss defined in file id %d",
	    ex_name_of_object(map_type), exoid);
    ex_err("ex_get_partial_num_map",errmsg,exerrval);
    return (EX_WARN);
  }

  /* Lookup index of element map id property array */
  id_ndx = ex_id_lkup(exoid,map_type,map_id);
  if (exerrval != 0) {
    sprintf(errmsg,
	    "Error: failed to locate %s id %"PRId64" in id variable in file id %d",
	    ex_name_of_object(map_type),map_id,exoid);
    ex_err("ex_get_partial_num_map",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* inquire id's of previously defined dimensions and variables */
  if ((status = nc_inq_varid(exoid, ex_name_of_map(map_type,id_ndx), &var_id)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to locate %s %"PRId64" in file id %d",
	    ex_name_of_object(map_type),map_id,exoid);
    ex_err("ex_get_partial_num_map",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* read in the map */
  start[0] = ent_start-1;
  count[0] = ent_count;

  if (ex_int64_status(exoid) & EX_MAPS_INT64_API) {
    status = nc_get_vara_longlong(exoid, var_id, start, count, map);
  } else {
    status = nc_get_vara_int(exoid, var_id, start, count, map);
  }

  if (status == -1) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get %s in file id %d",
	    ex_name_of_object(map_type),exoid);
    ex_err("ex_get_partial_num_map",errmsg,exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}
