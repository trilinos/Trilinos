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
 *      ne_get_n_elem_var()
 *
 *****************************************************************************
 *
 *  Variable Index:
 *
 *      neid               - The NetCDF ID of an already open NemesisI file.
 *      time_step          - The time step to write this data to.
 *      elem_var_index     - The index of this elemental variable.
 *      elem_blk_id        - The ID of the element block being written to.
 *      num_elem_this_blk  - The number of elements in this block.
 *      start_elem_num     - The start point for outputting data.
 *      num_elem           - The number of values to be output.
 *      elem_var_vals      - Pointer to the vector of values to be output.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include "exodusII.h"
#include "exodusII_int.h"

#include "ne_nemesisI.h"
#include "ne_nemesisI_int.h"

/*
 * reads the values of a single element variable for one element block at
 * one time step in the database; assume the first time step and
 * element variable index is 1
 */

int ne_get_n_elem_var (int   neid,
		       int   time_step,
		       int   elem_var_index,
		       int   elem_blk_id,
		       int   num_elem_this_blk,
		       int   start_elem_num,
		       int   num_elem,
		       void *elem_var_vals)
{
  int varid, elem_blk_id_ndx, status;
  size_t start[2], count[2];
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* Determine index of elem_blk_id in VAR_ID_EL_BLK array */
  if ((elem_blk_id_ndx = ex_id_lkup(neid,EX_ELEM_BLOCK,elem_blk_id)) < 0) {
    if (exerrval == EX_NULLENTITY) {
      sprintf(errmsg,
              "Warning: no element variables for NULL block %d in file id %d",
              elem_blk_id,neid);
      ex_err("ne_get_n_elem_var",errmsg,EX_MSG);
      return (EX_WARN);
    } else {
      sprintf(errmsg,
     "Error: failed to locate element block id %d in %s variable in file id %d",
              elem_blk_id, VAR_ID_EL_BLK, neid);
      ex_err("ne_get_n_elem_var",errmsg,exerrval);
      return (EX_FATAL);
    }
  }


  /* inquire previously defined variable */
  if ((status = nc_inq_varid(neid,VAR_ELEM_VAR(elem_var_index,elem_blk_id_ndx), &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
         "Error: failed to locate elem var %d for elem block %d in file id %d",
         elem_var_index,elem_blk_id,neid); /* this msg needs to be improved */
    ex_err("ne_get_n_elem_var",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* read values of element variable */
  start[0] = --time_step;
  start[1] = --start_elem_num;

  count[0] = 1;
  if ((num_elem_this_blk - start_elem_num) > num_elem)
    count[1] = num_elem;
  else
    count[1] = num_elem_this_blk - start_elem_num;

  if (ex_comp_ws(neid) == 4) {
    status = nc_get_vara_float(neid, varid, start, count, elem_var_vals);
  } else {
    status = nc_get_vara_double(neid, varid, start, count, elem_var_vals);
  }
  
  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
        "Error: failed to get elem var %d for block %d in file id %d",
            elem_var_index,elem_blk_id,neid);/*this msg needs to be improved*/
    ex_err("ne_get_n_elem_var",errmsg,exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}
