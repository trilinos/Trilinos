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
 *      ne_put_elem_var_slab()
 *
 *****************************************************************************
 *
 *  Variable Index:
 *
 *      neid               - The NetCDF ID of an already open NemesisI file.
 *      time_step          - The time step to write this data to.
 *      elem_var_index     - The index of this elemental variable.
 *      elem_blk_id        - The ID of the element block being written to.
 *      start_pos          - The start point for outputting data. The
 *                           first value is 0.
 *      num_vals           - The number of values to be output.
 *      elem_var_vals      - Pointer to the vector of values to be output.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include "exodusII.h"
#include "exodusII_int.h"

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"
#include <stdlib.h>

/*
 * writes the values of a single element variable for one element block,
 * starting at start_pos, at one time step to the database; assume the
 * first time step and element variable index are 1
 */

int ne_put_elem_var_slab (int   neid,
		          int   time_step,
		          int   elem_var_index,
		          int   elem_blk_id,
                          int   start_pos,
		          int   num_vals,
		          void *elem_var_vals)
{
  int status;
  int varid, dimid,time_dim, numelbdim, dims[2], elem_blk_id_ndx;
  size_t num_elem_blk, num_elem_var, start[2], count[2];
  int *elem_var_tab;
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* Determine index of elem_blk_id in VAR_ID_EL_BLK array */
  if ((elem_blk_id_ndx = ex_id_lkup(neid, EX_ELEM_BLOCK, elem_blk_id)) < 0) {
    if (exerrval == EX_NULLENTITY) {
      sprintf(errmsg,
	      "Warning: no variables allowed for NULL block %d in file id %d",
	      elem_blk_id, neid);
      ex_err("ne_put_elem_var_slab", errmsg, EX_MSG);
      return (EX_WARN);
    } else {
      sprintf(errmsg,
	      "Error: failed to locate element block id %d in %s array in file id %d",
	      elem_blk_id, VAR_ID_EL_BLK, neid);
      ex_err("ne_put_elem_var_slab", errmsg, exerrval);
      return (EX_FATAL);
    }
  }

  if ((status = nc_inq_varid (neid,
			      VAR_ELEM_VAR(elem_var_index, elem_blk_id_ndx), &varid)) != NC_NOERR) {
    if (status == NC_ENOTVAR) { /* variable doesn't exist, create it! */

      /*    inquire previously defined dimensions */

      /* check for the existance of an element variable truth table */
      if ((status = nc_inq_varid (neid, VAR_ELEM_TAB, &varid)) == NC_NOERR) {
	/* find out number of element blocks and element variables */
	if ((status = nc_inq_dimid (neid, DIM_NUM_EL_BLK, &dimid)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: failed to locate number of element blocks in file id %d",
		  neid);
	  ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	  return (EX_FATAL);
	}

	if ((status = nc_inq_dimlen(neid, dimid, &num_elem_blk)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: failed to get number of element blocks in file id %d",
		  neid);
	  ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	  return (EX_FATAL);
	}

	if ((status = nc_inq_dimid (neid, DIM_NUM_ELE_VAR, &dimid)) != NC_NOERR) {
	  exerrval = EX_BADPARAM;
	  sprintf(errmsg,
		  "Error: no element variables stored in file id %d",
		  neid);
	  ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	  return (EX_FATAL);
	}

	if ((status = nc_inq_dimlen(neid, dimid, &num_elem_var)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: failed to get number of element variables in file id %d",
		  neid);
	  ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	  return (EX_FATAL);
	}

	if (!(elem_var_tab =
	      (int *)malloc(num_elem_blk*num_elem_var*sizeof(int)))) {
	  exerrval = EX_MEMFAIL;
	  sprintf(errmsg,
		  "Error: failed to allocate memory for element variable truth table in file id %d",
		  neid);
	  ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	  return (EX_FATAL);
	}

	/*   read in the element variable truth table */
	if ((status = nc_get_var_int(neid, varid, elem_var_tab)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: failed to get truth table from file id %d", neid);
	  ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	  return (EX_FATAL);
	}

	if (elem_var_tab[num_elem_var*(elem_blk_id_ndx-1)+elem_var_index-1] == 0L) {
	  free(elem_var_tab);
	  exerrval = EX_BADPARAM;
	  sprintf(errmsg,
		  "Error: Invalid element variable %d, block %d in file id %d",
		  elem_var_index, elem_blk_id, neid);
	  ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	  return (EX_FATAL);
	}
	free(elem_var_tab);
      }

      if ((status = nc_inq_dimid (neid, DIM_TIME, &time_dim)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to locate time dimension in file id %d", neid);
	ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	goto error_ret;		/* exit define mode and return */
      }

      if ((status = nc_inq_dimid(neid, DIM_NUM_EL_IN_BLK(elem_blk_id_ndx), &numelbdim)) != NC_NOERR) {
	if (status == NC_EBADDIM) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: number of elements in element block %d not defined in file id %d",
		  elem_blk_id, neid);
	  ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	} else {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: failed to locate number of elements in element block %d in file id %d",
		  elem_blk_id, neid);
	  ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	}
	goto error_ret;
      }

      /*    variable doesn't exist so put file into define mode  */
      if ((status = nc_redef (neid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to put file id %d into define mode", neid);
	ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	return (EX_FATAL);
      }


      /*    define netCDF variable to store element variable values */
      dims[0] = time_dim;
      dims[1] = numelbdim;
      if ((status = nc_def_var(neid, VAR_ELEM_VAR(elem_var_index, elem_blk_id_ndx),
			       nc_flt_code(neid), 2, dims, &varid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to define element variable %d in file id %d",
		elem_var_index, neid);
	ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	goto error_ret;
      }


      /*    leave define mode  */
      if ((status = nc_enddef(neid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to complete element variable %s definition to file id %d",
		VAR_ELEM_VAR(elem_var_index, elem_blk_id_ndx), neid);
	ex_err("ne_put_elem_var_slab", errmsg, exerrval);
	return (EX_FATAL);
      }
    } else {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to locate element variable %s in file id %d",
	      VAR_ELEM_VAR(elem_var_index, elem_blk_id_ndx),neid);
      ex_err("ne_put_elem_var_slab", errmsg, exerrval);
      return (EX_FATAL);
    }
  }

  /* store element variable values */
  start[0] = --time_step;
  start[1] = --start_pos;

  count[0] = 1;
  count[1] = num_vals;

  if (ex_comp_ws(neid) == 4) {
    status = nc_put_vara_float(neid, varid, start, count, elem_var_vals);
  } else {
    status = nc_put_vara_double(neid, varid, start, count, elem_var_vals);
  }

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to store element variable %d in file id %d", 
	    elem_var_index, neid);
    ex_err("ne_put_elem_var_slab", errmsg, exerrval);
    return (EX_FATAL);
  }

  return (EX_NOERR);

  /* Fatal error: exit definition mode and return */
 error_ret:
  if (nc_enddef (neid) != NC_NOERR)     /* exit define mode */
    {
      sprintf(errmsg,
	      "Error: failed to complete definition for file id %d", neid);
      ex_err("ne_put_elem_var_slab", errmsg, exerrval);
    }
  return (EX_FATAL);
}
