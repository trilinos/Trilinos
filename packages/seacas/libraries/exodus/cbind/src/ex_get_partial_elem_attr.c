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
 *      ex_get_partial_elem_attr()
 *
 *****************************************************************************
 *
 *  Variable Index:
 *
 *      exoid               - The NetCDF ID of an already open NemesisI file.
 *      elem_blk_id        - The element block ID.
 *      start_elem_num     - The starting index of the elements to be
 *                           obtained.
 *      num_elems          - The number of FEM elements to read coords for.
 *      attrib             - Pointer to the attribute vector.
 *
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include <stdio.h>

#include "exodusII.h"
#include "exodusII_int.h"

/*
 * reads the attributes for an element block
 */

int ex_get_partial_elem_attr (int   exoid,
                        ex_entity_id   elem_blk_id,
                        int64_t   start_elem_num,
                        int64_t   num_elems,
                        void *attrib)

{
  int numelbdim, numattrdim, attrid, elem_blk_id_ndx, status;
  size_t num_elem_this_blk, num_attr, start[2], count[2];
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* Determine index of elem_blk_id in VAR_ID_EL_BLK array */
  if ((elem_blk_id_ndx = ex_id_lkup(exoid, EX_ELEM_BLOCK, elem_blk_id)) < 0) {
    if (exerrval == EX_NULLENTITY) {
      sprintf(errmsg,
              "Warning: no attributes found for NULL block %"PRId64" in file id %d",
              elem_blk_id,exoid);
      ex_err("ex_get_partial_elem_attr",errmsg,EX_MSG);
      return (EX_WARN);              /* no attributes for this element block */
    } else {
      sprintf(errmsg,
      "Warning: failed to locate element block %"PRId64" in %s array in file id %d",
              elem_blk_id,VAR_ID_EL_BLK, exoid);
      ex_err("ex_get_partial_elem_attr",errmsg,exerrval);
      return (EX_WARN);
    }
  }


/* inquire id's of previously defined dimensions  */

  if ((status = nc_inq_dimid(exoid, DIM_NUM_EL_IN_BLK(elem_blk_id_ndx), &numelbdim)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
        "Error: failed to locate number of elements for block %"PRId64" in file id %d",
            elem_blk_id, exoid);
    ex_err("ex_get_partial_elem_attr",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_dimlen(exoid, numelbdim, &num_elem_this_blk)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
           "Error: failed to get number of elements for block %"PRId64" in file id %d",
            elem_blk_id,exoid);
    ex_err("ex_get_partial_elem_attr",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_dimid(exoid, DIM_NUM_ATT_IN_BLK(elem_blk_id_ndx), &numattrdim)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Warning: no attributes found for block %"PRId64" in file id %d",
            elem_blk_id,exoid);
    ex_err("ex_get_partial_elem_attr",errmsg,EX_MSG);
    return (EX_WARN);              /* no attributes for this element block */
  }

  if ((status = nc_inq_dimlen(exoid, numattrdim, &num_attr)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
         "Error: failed to get number of attributes for block %"PRId64" in file id %d",
            elem_blk_id,exoid);
    ex_err("ex_get_partial_elem_attr",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_varid(exoid, VAR_ATTRIB(elem_blk_id_ndx), &attrid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to locate attributes for block %"PRId64" in file id %d",
            elem_blk_id,exoid);
    ex_err("ex_get_partial_elem_attr",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* do some error checking */
  if (num_elem_this_blk < (start_elem_num + num_elems - 1)) {
    exerrval = status;
    sprintf(errmsg,
      "Error: requested attributes from too many elements in this block, %"PRId64"",
            elem_blk_id);
    ex_err("ex_get_partial_elem_attr",errmsg, exerrval);
    return(EX_FATAL);
  }

/* read in the attributes */

  start[0] = --start_elem_num;
  start[1] = 0;

  count[0] = num_elems;
  count[1] = num_attr;

  if (ex_comp_ws(exoid) == 4) {
    status = nc_get_vara_float(exoid, attrid, start, count, attrib);
  } else {
    status = nc_get_vara_double(exoid, attrid, start, count, attrib);
  }

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get attributes for block %"PRId64" in file id %d",
            elem_blk_id,exoid);
    ex_err("ex_get_partial_elem_attr",errmsg,exerrval);
    return (EX_FATAL);
  }
  return(EX_NOERR);
}
