/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
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
 *     * Neither the name of NTESS nor the names of its
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
 * exgnm - ex_get_map
 *
 * entry conditions -
 *   input parameters:
 *       int     exoid                   exodus file id
 *       int     map_type                type of map (node, edge, face, element)
 *       int     map_id                  map id
 *
 * exit conditions -
 *       int*    map                     map
 *
 * revision history -
 *   20060930 - David Thompson - Adapted from ex_get_node_map
 *
 *****************************************************************************/

#include "exodusII.h"     // for ex_err, etc
#include "exodusII_int.h" // for EX_FATAL, EX_NOERR, etc
#include "netcdf.h"       // for NC_NOERR, nc_inq_dimid, etc
#include <inttypes.h>     // for PRId64
#include <stdio.h>

/*
 * reads the map with specified ID
 */

int ex_get_num_map(int exoid, ex_entity_type map_type, ex_entity_id map_id, void_int *map)
{
  int         dimid, var_id, id_ndx, status;
  char        errmsg[MAX_ERR_LENGTH];
  const char *dim_map_size;
  const char *dim_num_maps;

  EX_FUNC_ENTER();
  ex_check_valid_file_id(exoid, __func__);

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
    snprintf(errmsg, MAX_ERR_LENGTH, "Bad map type (%d) specified", map_type);
    ex_err(__func__, errmsg, EX_BADPARAM);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* See if any entries are stored in this file */
  if (nc_inq_dimid(exoid, dim_map_size, &dimid) != NC_NOERR) {
    EX_FUNC_LEAVE(EX_NOERR);
  }

  /* first check if any maps have been defined */
  if ((status = nc_inq_dimid(exoid, dim_num_maps, &dimid)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "Warning: no %ss defined in file id %d",
             ex_name_of_object(map_type), exoid);
    ex_err(__func__, errmsg, status);
    EX_FUNC_LEAVE(EX_WARN);
  }

  /* Lookup index of map id property array */
  id_ndx = ex_id_lkup(exoid, map_type, map_id);
  if (id_ndx <= 0) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to locate %s id %" PRId64 " in id variable in file id %d",
             ex_name_of_object(map_type), map_id, exoid);
    ex_err(__func__, errmsg, EX_LASTERR);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* inquire id's of previously defined dimensions and variables */
  if ((status = nc_inq_varid(exoid, ex_name_of_map(map_type, id_ndx), &var_id)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to locate %s %" PRId64 " in file id %d",
             ex_name_of_object(map_type), map_id, exoid);
    ex_err(__func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* read in the map */
  if (ex_int64_status(exoid) & EX_MAPS_INT64_API) {
    status = nc_get_var_longlong(exoid, var_id, map);
  }
  else {
    status = nc_get_var_int(exoid, var_id, map);
  }

  if (status != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get %s in file id %d",
             ex_name_of_object(map_type), exoid);
    ex_err(__func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  EX_FUNC_LEAVE(EX_NOERR);
}
