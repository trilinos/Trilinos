/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
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
* expecpp - ex_put_entity_count_per_polyhedra
*
*****************************************************************************/

#include "exodusII.h"     // for ex_err, exerrval, etc
#include "exodusII_int.h" // for EX_FATAL, ex_id_lkup, etc
#include "netcdf.h"       // for nc_inq_varid, NC_NOERR, etc
#include <inttypes.h>     // for PRId64
#include <stdio.h>

/*!
 * writes out the number of entities (nodes/faces) per polyhedra
 * (nsided/nfaced) in this element block to the database.
 * \param  exoid                exodus file id
 * \param  blk_type             type of block (face, or element)
 * \param  blk_id               block identifer
 * \param  entity_counts        entity_per_polyhedra count array
 */

int ex_put_entity_count_per_polyhedra(int exoid, ex_entity_type blk_type, ex_entity_id blk_id,
                                      const int *entity_counts)
{
  int  npeid = -1, blk_id_ndx, status;
  char errmsg[MAX_ERR_LENGTH];

  ex_check_valid_file_id(exoid);

  exerrval = 0; /* clear error code */

  blk_id_ndx = ex_id_lkup(exoid, blk_type, blk_id);
  if (exerrval != 0) {
    if (exerrval == EX_NULLENTITY) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "Warning: entity_counts array not allowed for NULL %s block %" PRId64
               " in file id %d",
               ex_name_of_object(blk_type), blk_id, exoid);
      ex_err("ex_put_entity_count_per_polyhedra", errmsg, EX_NULLENTITY);
      return (EX_WARN);
    }

    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to locate %s block id %" PRId64 " in id array in file id %d",
             ex_name_of_object(blk_type), blk_id, exoid);
    ex_err("ex_put_entity_count_per_polyhedra", errmsg, exerrval);
    return (EX_FATAL);
  }

  /* inquire id's of previously defined dimensions  */
  switch (blk_type) {
  case EX_ELEM_BLOCK: status = nc_inq_varid(exoid, VAR_EBEPEC(blk_id_ndx), &npeid); break;
  case EX_FACE_BLOCK: status = nc_inq_varid(exoid, VAR_FBEPEC(blk_id_ndx), &npeid); break;
  default:
    exerrval = 1005;
    snprintf(errmsg, MAX_ERR_LENGTH,
             "Internal ERROR: unrecognized block type in switch: %d in file id %d", blk_type,
             exoid);
    ex_err("ex_put_entity_count_per_polyhedra", errmsg, EX_MSG);
    return (EX_FATAL);
  }
  if (status != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to locate entity_counts array for %s block %" PRId64 " in file id %d",
             ex_name_of_object(blk_type), blk_id, exoid);
    ex_err("ex_put_entity_count_per_polyhedra", errmsg, exerrval);
    return (EX_FATAL);
  }

  status = nc_put_var_int(exoid, npeid, entity_counts);
  if (status != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to write node counts array for %s block %" PRId64 " in file id %d",
             ex_name_of_object(blk_type), blk_id, exoid);
    ex_err("ex_put_entity_count_per_polyhedra", errmsg, exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}
