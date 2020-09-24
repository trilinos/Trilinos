/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * exgssd - ex_get_set_dist_fact
 *
 * entry conditions -
 *   input parameters:
 *       int     exoid                   exodus file id
 *       int     set_type                type of set
 *       int     set_id                  set id
 *
 * exit conditions -
 *       float*  set_dist_fact           array of dist factors for set
 *
 * revision history -
 *
 *
 *****************************************************************************/

#include "exodusII.h"     // for ex_err, ex_name_of_object, etc
#include "exodusII_int.h" // for EX_FATAL, EX_WARN, etc

/*!
 * reads the distribution factors for a single set
 */

int ex_get_set_dist_fact(int exoid, ex_entity_type set_type, ex_entity_id set_id,
                         void *set_dist_fact)
{

  int   dimid, dist_id, set_id_ndx;
  int   status;
  char  errmsg[MAX_ERR_LENGTH];
  char *factptr = NULL;

  EX_FUNC_ENTER();
  if (ex__check_valid_file_id(exoid, __func__) == EX_FATAL) {
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* first check if any sets are specified */
  if ((status = nc_inq_dimid(exoid, ex__dim_num_objects(set_type), &dimid)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "Warning: no %s sets stored in file id %d",
             ex_name_of_object(set_type), exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_WARN);
  }

  /* Lookup index of set id in VAR_*S_IDS array */
  set_id_ndx = ex__id_lkup(exoid, set_type, set_id);
  if (set_id_ndx <= 0) {
    ex_get_err(NULL, NULL, &status);

    if (status != 0) {
      if (status == EX_NULLENTITY) {
        snprintf(errmsg, MAX_ERR_LENGTH, "Warning: %s set %" PRId64 " is NULL in file id %d",
                 ex_name_of_object(set_type), set_id, exoid);
        ex_err_fn(exoid, __func__, errmsg, EX_NULLENTITY);
        EX_FUNC_LEAVE(EX_WARN);
      }
      snprintf(errmsg, MAX_ERR_LENGTH,
               "ERROR: failed to locate %s set %" PRId64 " in VAR_*S_IDS array in file id %d",
               ex_name_of_object(set_type), set_id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_FATAL);
    }
  }

  /* setup more pointers based on set_type */
  if (set_type == EX_NODE_SET) {
    factptr = VAR_FACT_NS(set_id_ndx);
  }
  else if (set_type == EX_EDGE_SET) {
    factptr = VAR_FACT_ES(set_id_ndx);
  }
  else if (set_type == EX_FACE_SET) {
    factptr = VAR_FACT_FS(set_id_ndx);
  }
  else if (set_type == EX_SIDE_SET) {
    factptr = VAR_FACT_SS(set_id_ndx);
  }
  if (set_type == EX_ELEM_SET) {
    factptr = VAR_FACT_ELS(set_id_ndx);
  }

  /* inquire id's of previously defined dimensions and variables */
  if ((status = nc_inq_varid(exoid, factptr, &dist_id)) != NC_NOERR) {
    /* not an error for node sets because this is how we check that df's exist
     */
    if (set_type == EX_NODE_SET) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "Warning: dist factors not stored for %s set %" PRId64 " in file id %d",
               ex_name_of_object(set_type), set_id, exoid);
      ex_err_fn(exoid, __func__, errmsg, status);
      EX_FUNC_LEAVE(EX_WARN); /* complain - but not too loud */
    }
    /* is an error for other sets */

    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to locate dist factors list for %s set %" PRId64 " in file id %d",
             ex_name_of_object(set_type), set_id, exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* read in the distribution factors array */
  if (ex__comp_ws(exoid) == 4) {
    status = nc_get_var_float(exoid, dist_id, set_dist_fact);
  }
  else {
    status = nc_get_var_double(exoid, dist_id, set_dist_fact);
  }

  if (status != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to get dist factors list for %s set %" PRId64 " in file id %d",
             ex_name_of_object(set_type), set_id, exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_FATAL);
  }
  EX_FUNC_LEAVE(EX_NOERR);
}
