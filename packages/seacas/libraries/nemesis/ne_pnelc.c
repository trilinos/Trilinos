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
 *      ne_put_n_elem_conn()
 *
 *****************************************************************************
 *
 *  Variable Index:
 *
 *      neid               - The NetCDF ID of an already open NemesisI file.
 *      elem_blk_id        - The element block ID.
 *      start_elem_num     - The starting index of the elements to be
 *                           obtained.
 *      num_elems          - The number of FEM elements to read coords for.
 *      connect            - Pointer to the connectivity vector.
 *
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "exodusII.h"
#include "exodusII_int.h"

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

/*
 * writes the connectivity array for an element block
 */

int ne_put_n_elem_conn (int  neid,
			int  elem_blk_id,
			int  start_elem_num,
			int  num_elems,
			int *connect)
{
  int numelbdim, nelnoddim, connid, elem_blk_id_ndx, status;
  size_t num_elem_this_blk, num_nod_per_elem, start[2], count[2]; 
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* Determine index of elem_blk_id in VAR_ID_EL_BLK array */
  if ((elem_blk_id_ndx = ex_id_lkup(neid, EX_ELEM_BLOCK, elem_blk_id)) == -1)
    {
      if (exerrval == EX_NULLENTITY) {
	sprintf(errmsg,
		"Warning: connectivity array not allowed for NULL element block %d in file id %d",
		elem_blk_id, neid);
	ex_err("ne_put_n_elem_conn",errmsg,EX_MSG);
	return (EX_WARN);
      } else {

	sprintf(errmsg,
		"Error: failed to locate element block id %d in %s array in file id %d",
		elem_blk_id,VAR_ID_EL_BLK, neid);
	ex_err("ne_put_n_elem_conn",errmsg,exerrval);
	return (EX_FATAL);
      }
    }

  /* inquire id's of previously defined dimensions  */

  if ((status = nc_inq_dimid (neid, DIM_NUM_EL_IN_BLK(elem_blk_id_ndx), &numelbdim)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to locate number of elements in block %d in file id %d",
	    elem_blk_id, neid);
    ex_err("ne_put_n_elem_conn",errmsg, exerrval);
    return(EX_FATAL);
  }

  if ((status = nc_inq_dimlen(neid, numelbdim, &num_elem_this_blk)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get number of elements in block %d in file id %d",
            elem_blk_id, neid);
    ex_err("ne_put_n_elem_conn",errmsg,exerrval);
    return(EX_FATAL);
  }

  if ((status = nc_inq_dimid (neid, DIM_NUM_NOD_PER_EL(elem_blk_id_ndx), &nelnoddim)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to locate number of nodes/elem in block %d in file id %d",
            elem_blk_id, neid);
    ex_err("ne_put_n_elem_conn",errmsg,exerrval);
    return(EX_FATAL);
  }

  if ((status = nc_inq_dimlen (neid, nelnoddim, &num_nod_per_elem)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get number of nodes/elem in block %d in file id %d",
            elem_blk_id, neid);
    ex_err("ne_put_n_elem_conn",errmsg,exerrval);
    return(EX_FATAL);
  }


  if ((status = nc_inq_varid (neid, VAR_CONN(elem_blk_id_ndx), &connid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to locate connectivity array for element block %d in file id %d",
            elem_blk_id, neid);
    ex_err("ne_put_n_elem_conn",errmsg, exerrval);
    return(EX_FATAL);
  }

  /* do some error checking */
  if (num_elem_this_blk < (start_elem_num + num_elems - 1)) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: requested connectivity from too many elements in this block, %d",
            elem_blk_id);
    ex_err("ne_put_n_elem_conn",errmsg, exerrval);
    return(EX_FATAL);
  }

  /* write out the connectivity array */
  start[0] = --start_elem_num;
  start[1] = 0;

  count[0] = num_elems;
  count[1] = num_nod_per_elem;

  status = nc_put_vara_int(neid, connid, start, count, connect);

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to write connectivity array for block %d in file id %d",
            elem_blk_id, neid);
    ex_err("ne_put_n_elem_conn",errmsg, exerrval);
    return(EX_FATAL);
  }
  return (EX_NOERR);
}
