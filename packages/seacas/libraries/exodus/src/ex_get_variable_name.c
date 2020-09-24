/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * exgvnm - ex_get_variable_name
 *
 * entry conditions -
 *   input parameters:
 *       int   exoid                   exodus file id
 *       int   obj_type                variable type
 *       int   var_num                 variable index to read 1..num_var
 *
 * exit conditions -
 *       char*   var_name                ptr to variable name
 *
 * revision history -
 *
 *
 *****************************************************************************/

#include "exodusII.h"     // for ex_err, etc
#include "exodusII_int.h" // for EX_FATAL, etc

/*!
 * \ingroup ResultsData
 *
 * reads the name of a particular results variable from the database
 */

int ex_get_variable_name(int exoid, ex_entity_type obj_type, int var_num, char *var_name)
{
  int         status;
  int         varid;
  char        errmsg[MAX_ERR_LENGTH];
  const char *vname = NULL;

  EX_FUNC_ENTER();
  if (ex__check_valid_file_id(exoid, __func__) == EX_FATAL) {
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /* inquire previously defined variables  */

  switch (obj_type) {
  case EX_GLOBAL: vname = VAR_NAME_GLO_VAR; break;
  case EX_NODAL: vname = VAR_NAME_NOD_VAR; break;
  case EX_ASSEMBLY: vname = VAR_NAME_ASSEMBLY_VAR; break;
  case EX_BLOB: vname = VAR_NAME_BLOB_VAR; break;
  case EX_EDGE_BLOCK: vname = VAR_NAME_EDG_VAR; break;
  case EX_FACE_BLOCK: vname = VAR_NAME_FAC_VAR; break;
  case EX_ELEM_BLOCK: vname = VAR_NAME_ELE_VAR; break;
  case EX_NODE_SET: vname = VAR_NAME_NSET_VAR; break;
  case EX_EDGE_SET: vname = VAR_NAME_ESET_VAR; break;
  case EX_FACE_SET: vname = VAR_NAME_FSET_VAR; break;
  case EX_SIDE_SET: vname = VAR_NAME_SSET_VAR; break;
  case EX_ELEM_SET: vname = VAR_NAME_ELSET_VAR; break;
  default:
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: Invalid variable type (%d) given for file id %d",
             obj_type, exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_BADPARAM);
    EX_FUNC_LEAVE(EX_FATAL);
  }

  if ((status = nc_inq_varid(exoid, vname, &varid)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "Warning: no %s variable names stored in file id %d",
             ex_name_of_object(obj_type), exoid);
    ex_err_fn(exoid, __func__, errmsg, status);
    EX_FUNC_LEAVE(EX_WARN);
  }

  /* read the variable name */
  {
    int db_name_size  = ex_inquire_int(exoid, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);
    int api_name_size = ex_inquire_int(exoid, EX_INQ_MAX_READ_NAME_LENGTH);
    int name_size     = db_name_size < api_name_size ? db_name_size : api_name_size;

    status = ex__get_name(exoid, varid, var_num - 1, var_name, name_size, obj_type, __func__);
    if (status != NC_NOERR) {
      EX_FUNC_LEAVE(EX_FATAL);
    }
  }
  EX_FUNC_LEAVE(EX_NOERR);
}
