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
 *      ex_put_n_side_set_df()
 *
 *****************************************************************************
 *
 *  Variable Index:
 *
 *      exoid               - The NetCDF ID of an already open NemesisI file.
 *      side_set_id        - ID of side set to written.
 *      start_num          - The starting index of the df's to be written.
 *      num_df_to_get      - The number of sides to write.
 *      side_set_dist_fact - List of distribution factors for the side set.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "exodusII.h"
#include "exodusII_int.h"

/*
 * writes the distribution factors for a single side set
 */

int ex_put_n_side_set_df (int   exoid,
                          ex_entity_id   side_set_id,
                          int64_t   start_num,
                          int64_t   num_df_to_get,
                          void *side_set_dist_fact)
{
  int status;
  int dimid, side_set_id_ndx;
  int dist_id;
  size_t num_df_in_set,  start[1], count[1];
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* first check if any side sets are specified */

  if ((status = nc_inq_dimid (exoid, DIM_NUM_SS, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: no side sets specified in file id %d",
            exoid);
    ex_err("ex_put_n_side_set_df",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* Lookup index of side set id in VAR_SS_IDS array */

  if ((side_set_id_ndx = ex_id_lkup(exoid, EX_SIDE_SET, side_set_id)) < 0)
  {
    if (exerrval == EX_NULLENTITY) {
      sprintf(errmsg,
              "Warning: no data allowed for NULL side set %"PRId64" in file id %d",
              side_set_id, exoid);
      ex_err("ex_put_side_set_fact",errmsg,EX_MSG);
      return (EX_WARN);
    } else {
      sprintf(errmsg,
     "Error: failed to locate side set id %"PRId64" in VAR_SS_IDS array in file id %d",
              side_set_id, exoid);
      ex_err("ex_put_n_side_set_df",errmsg,exerrval);
      return (EX_FATAL);
    }
  }

  /* inquire id's of previously defined dimension and variable */

  if ((status = nc_inq_dimid (exoid, DIM_NUM_DF_SS(side_set_id_ndx), &dimid)) != NC_NOERR) {
    if (status == NC_EBADDIM) {
      exerrval = EX_BADPARAM;
      sprintf(errmsg,
              "Warning: no dist factors defined for side set %"PRId64" in file id %d",
              side_set_id, exoid);
      ex_err("ex_put_n_side_set_df",errmsg,exerrval);
      return (EX_WARN);

    } else {
      exerrval = status;
      sprintf(errmsg,
  "Error: failed to locate number of dist factors in side set %"PRId64" in file id %d",
              side_set_id, exoid);
      ex_err("ex_put_n_side_set_df",errmsg,exerrval);
      return (EX_FATAL);
    }
  }

  if ((status = nc_inq_dimlen(exoid, dimid, &num_df_in_set)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
     "Error: failed to get number of dist factors in side set %"PRId64" in file id %d",
            side_set_id, exoid);
    ex_err("ex_put_n_side_set_df",errmsg,exerrval);
    return (EX_FATAL);
  }

  /* Check input parameters for a valid range of numbers */
  if (start_num < 0 || (num_df_to_get > 0 && start_num > num_df_in_set)) {
    exerrval = EX_BADPARAM;
    sprintf(errmsg, "Error: Invalid input");
    ex_err("ex_put_n_side_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  if (num_df_to_get < 0) {
    exerrval = EX_BADPARAM;
    sprintf(errmsg, "Error: Invalid number of df's to put!");
    ex_err("ex_put_n_side_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  /* start_num now starts at 1, not 0 */
  if ((start_num + num_df_to_get) > num_df_in_set+1) {
    exerrval = EX_BADPARAM;
    sprintf(errmsg, "Error: request larger than number of df's in set!");
    ex_err("ex_put_n_side_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_varid (exoid, VAR_FACT_SS(side_set_id_ndx), &dist_id)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
      "Error: failed to locate dist factors list for side set %"PRId64" in file id %d",
            side_set_id, exoid);
    ex_err("ex_put_n_side_set_df",errmsg,exerrval);
    return (EX_FATAL);
  }


  /* write out the distribution factors array */
  start[0] = --start_num;
  count[0] = num_df_to_get;

  if (ex_comp_ws(exoid) == 4) {
    status = nc_put_vara_float(exoid, dist_id, start, count, side_set_dist_fact);
  } else {
    status = nc_put_vara_double(exoid, dist_id, start, count, side_set_dist_fact);
  }

  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to store dist factors for side set %"PRId64" in file id %d",
            side_set_id, exoid);
    ex_err("ex_put_n_side_set_df",errmsg,exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}
