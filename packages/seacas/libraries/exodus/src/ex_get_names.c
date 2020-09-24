/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * exgnam - ex_get_names
 *
 * entry conditions -
 *   input parameters:
 *       int     exoid          exodus file id
 *       int    obj_type,
 *
 * exit conditions -
 *       char*   names[]           ptr array of names
 *
 * revision history -
 *
 *
 *****************************************************************************/

#include "exodusII.h"     // for ex_err, etc
#include "exodusII_int.h" // for ex__get_dimension, EX_NOERR, etc

/*
 * reads the entity names from the database
 */

int ex_get_names(int exoid, ex_entity_type obj_type, char **names)
{
  int    status;
  int    varid, temp;
  size_t num_entity, i;
  char   errmsg[MAX_ERR_LENGTH];

  EX_FUNC_ENTER();
  if (ex__check_valid_file_id(exoid, __func__) == EX_FATAL) {
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* inquire previously defined dimensions and variables  */

  switch (obj_type) {
    /* ======== ASSEMBLY ========= */
  case EX_ASSEMBLY:
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: Assembly names are read using `ex_get_assembly()` function");
    ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
    EX_FUNC_LEAVE(EX_FATAL);
    break;
  /*  ======== BLOCKS ========= */
  case EX_EDGE_BLOCK:
    ex__get_dimension(exoid, DIM_NUM_ED_BLK, "edge block", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_ED_BLK, &varid);
    break;
  case EX_FACE_BLOCK:
    ex__get_dimension(exoid, DIM_NUM_FA_BLK, "face block", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_FA_BLK, &varid);
    break;
  case EX_ELEM_BLOCK:
    ex__get_dimension(exoid, DIM_NUM_EL_BLK, "element block", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_EL_BLK, &varid);
    break;

  /*  ======== SETS ========= */
  case EX_NODE_SET:
    ex__get_dimension(exoid, DIM_NUM_NS, "nodeset", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_NS, &varid);
    break;
  case EX_EDGE_SET:
    ex__get_dimension(exoid, DIM_NUM_ES, "edgeset", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_ES, &varid);
    break;
  case EX_FACE_SET:
    ex__get_dimension(exoid, DIM_NUM_FS, "faceset", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_FS, &varid);
    break;
  case EX_SIDE_SET:
    ex__get_dimension(exoid, DIM_NUM_SS, "sideset", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_SS, &varid);
    break;
  case EX_ELEM_SET:
    ex__get_dimension(exoid, DIM_NUM_ELS, "elemset", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_ELS, &varid);
    break;

  /*  ======== MAPS ========= */
  case EX_NODE_MAP:
    ex__get_dimension(exoid, DIM_NUM_NM, "node map", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_NM, &varid);
    break;
  case EX_EDGE_MAP:
    ex__get_dimension(exoid, DIM_NUM_EDM, "edge map", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_EDM, &varid);
    break;
  case EX_FACE_MAP:
    ex__get_dimension(exoid, DIM_NUM_FAM, "face map", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_FAM, &varid);
    break;
  case EX_ELEM_MAP:
    ex__get_dimension(exoid, DIM_NUM_EM, "element map", &num_entity, &temp, __func__);
    status = nc_inq_varid(exoid, VAR_NAME_EM, &varid);
    break;

  /* invalid variable type */
  default:
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: Invalid type specified in file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  if (status == NC_NOERR) {
    if ((status = ex__get_names(exoid, varid, num_entity, names, obj_type, "ex_get_names")) !=
        EX_NOERR) {
      EX_FUNC_LEAVE(status);
    }
  }
  else {
    /* Names variable does not exist on the database; probably since this is an
     * older version of the database.  Return an empty array...
     */
    for (i = 0; i < num_entity; i++) {
      names[i][0] = '\0';
    }
  }
  EX_FUNC_LEAVE(EX_NOERR);
}
