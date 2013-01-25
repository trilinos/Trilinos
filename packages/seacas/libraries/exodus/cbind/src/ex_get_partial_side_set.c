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
 *
 *      ex_get_partial_side_set()
 *
 *****************************************************************************
 *
 *  Variable Index:
 *
 *      exoid               - The NetCDF ID of an already open NemesisI file.
 *      side_set_id        - ID of side set to read.
 *      start_side_num     - The starting index of the sides to be read.
 *      num_sides          - The number of sides to read in.
 *      side_set_elem_list - List of element IDs in side set.
 *      side_set_side_list - List of side IDs in side set.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "exodusII.h"
#include "exodusII_int.h"

/*
 * reads the side set element list and side set side list for a single side set
 */

int ex_get_partial_side_set (int   exoid,
                       ex_entity_id   side_set_id,
                       int64_t   start_side_num,
                       int64_t   num_sides,
                       void_int  *side_set_elem_list,
                       void_int  *side_set_side_list)
{

  int      dimid, elem_list_id, side_list_id, status;
  int      side_set_id_ndx;
  size_t   num_side_in_set, count[1], start[1];
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* first check if any side sets are specified */
  if ((status = nc_inq_dimid (exoid, DIM_NUM_SS, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Warning: no side sets stored in file id %d",
            exoid);
    ex_err("ex_get_partial_side_set",errmsg,exerrval);
    return (EX_WARN);
  }

  /* Lookup index of side set id in VAR_SS_IDS array */
  if ((side_set_id_ndx = ex_id_lkup(exoid, EX_SIDE_SET, side_set_id)) < 0) {
    if (exerrval == EX_NULLENTITY) {
      sprintf(errmsg,
              "Warning: side set %"PRId64" is NULL in file id %d",
              side_set_id,exoid);
      ex_err("ex_get_partial_side_set",errmsg,EX_MSG);
      return (EX_WARN);
    } else {

      sprintf(errmsg,
	      "Error: failed to locate side set id %"PRId64" in %s array in file id %d",
              side_set_id, VAR_SS_IDS, exoid);
      ex_err("ex_get_partial_side_set",errmsg,exerrval);
      return (EX_FATAL);
    }
  }

  /* inquire id's of previously defined dimensions and variables */
  if ((status = nc_inq_dimid (exoid, DIM_NUM_SIDE_SS(side_set_id_ndx), &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to locate number of sides in side set %"PRId64" in file id %d",
            side_set_id,exoid);
    ex_err("ex_get_partial_side_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_dimlen(exoid, dimid, &num_side_in_set)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get number of sides in side set %"PRId64" in file id %d",
            side_set_id,exoid);
    ex_err("ex_get_partial_side_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_varid(exoid, VAR_ELEM_SS(side_set_id_ndx), &elem_list_id)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to locate element list for side set %"PRId64" in file id %d",
            side_set_id,exoid);
    ex_err("ex_get_partial_side_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_varid (exoid, VAR_SIDE_SS(side_set_id_ndx), &side_list_id)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to locate side list for side set %"PRId64" in file id %d",
            side_set_id,exoid);
    ex_err("ex_get_partial_side_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* Check input parameters for a valid range of side numbers */
  if (start_side_num < 0 || start_side_num > num_side_in_set) {
    exerrval = EX_BADPARAM;
    sprintf(errmsg, "Error: Invalid input");
    ex_err("ex_get_partial_side_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  if (num_sides < 0) {
    exerrval = EX_BADPARAM;
    sprintf(errmsg, "Error: Invalid number of elements in side set!");
    ex_err("ex_get_partial_side_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* start_side_num now starts at 1, not 0 */
  if ((start_side_num + num_sides - 1) > num_side_in_set) {
    exerrval = EX_BADPARAM;
    sprintf(errmsg, "Error: request larger than number of elements in set!");
    ex_err("ex_get_partial_side_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* read in the element list and side list arrays */
  start[0] = --start_side_num;
  count[0] = num_sides;

  if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
    status = nc_get_vara_longlong(exoid, elem_list_id, start, count, side_set_elem_list);
  } else {
    status = nc_get_vara_int(exoid, elem_list_id, start, count, side_set_elem_list);
  }

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get element list for side set %"PRId64" in file id %d",
            side_set_id,exoid);
    ex_err("ex_get_partial_side_set",errmsg,exerrval);
    return (EX_FATAL);
  }

  if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
    status = nc_get_vara_longlong(exoid, side_list_id, start, count, side_set_side_list);
  } else {
    status = nc_get_vara_int(exoid, side_list_id, start, count, side_set_side_list);
  }

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get side list for side set %"PRId64" in file id %d",
            side_set_id,exoid);
    ex_err("ex_get_partial_side_set",errmsg,exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}
