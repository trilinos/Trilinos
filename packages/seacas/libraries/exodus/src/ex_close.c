/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * exclos - ex_close
 *
 * entry conditions -
 *   input parameters:
 *       int     exoid                   exodus file id
 *
 * exit conditions -
 *
 * revision history -
 *
 *
 *****************************************************************************/

#include "exodusII.h"     // for ex_err, etc
#include "exodusII_int.h" // for ex__get_counter_list, etc

/*!
\ingroup Utilities

The function ex_close() updates and then closes an open exodus file.

\return In case of an error, ex_close() returns a negative number; a
        warning will return a positive number. Possible causes of errors
        include:
 - data file not properly opened with call to ex_create() or ex_open()

 \param exoid      exodus file ID returned from a previous call to ex_create()
or ex_open().

The following code segment closes an open exodus file:

~~~{.c}
int error,exoid;
error = ex_close (exoid);
~~~

 */
int ex_close(int exoid)
{
  char errmsg[MAX_ERR_LENGTH];
  int  status;
  int  status1;
  int  status2;
#if NC_HAS_HDF5
  int parent_id = 0;
#endif

  EX_FUNC_ENTER();

  if (ex__check_valid_file_id(exoid, __func__) == EX_FATAL) {
    EX_FUNC_LEAVE(EX_FATAL);
  }

  /*
   * NOTE: If using netcdf-4, exoid must refer to the root group.
   * Need to determine whether there are any groups and if so,
   * call ex__rm_file_item and ex__rm_stat_ptr on each group.
   */

#if NC_HAS_HDF5
  /* nc_inq_grp_parent() will return NC_ENOGRP error if exoid
   * refers to the root group (which is what we want)
   */
  if ((status = nc_inq_grp_parent(exoid, &parent_id)) != NC_ENOGRP) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: file id %d does not refer to root group.", exoid);
    ex_err_fn(exoid, __func__, errmsg, EX_NOTROOTID);
    EX_FUNC_LEAVE(EX_FATAL);
  }
#endif

  if ((status1 = nc_sync(exoid)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to update file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, status1);
  }

  if ((status2 = nc_close(exoid)) != NC_NOERR) {
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to close file id %d", exoid);
    ex_err_fn(exoid, __func__, errmsg, status2);
  }

  /* Even if we have failures above due to nc_sync() or nc_close(), we still need to clean up our
   * internal datastructures.
   */

  ex__rm_file_item(exoid, ex__get_counter_list(EX_ELEM_BLOCK));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_FACE_BLOCK));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_EDGE_BLOCK));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_NODE_SET));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_EDGE_SET));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_FACE_SET));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_SIDE_SET));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_ELEM_SET));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_NODE_MAP));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_EDGE_MAP));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_FACE_MAP));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_ELEM_MAP));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_ASSEMBLY));
  ex__rm_file_item(exoid, ex__get_counter_list(EX_BLOB));

  ex__rm_stat_ptr(exoid, &exoII_ed);
  ex__rm_stat_ptr(exoid, &exoII_fa);
  ex__rm_stat_ptr(exoid, &exoII_eb);
  ex__rm_stat_ptr(exoid, &exoII_ns);
  ex__rm_stat_ptr(exoid, &exoII_es);
  ex__rm_stat_ptr(exoid, &exoII_fs);
  ex__rm_stat_ptr(exoid, &exoII_ss);
  ex__rm_stat_ptr(exoid, &exoII_els);
  ex__rm_stat_ptr(exoid, &exoII_nm);
  ex__rm_stat_ptr(exoid, &exoII_edm);
  ex__rm_stat_ptr(exoid, &exoII_fam);
  ex__rm_stat_ptr(exoid, &exoII_em);

  ex__conv_exit(exoid);

  status = EX_NOERR;
  if (status1 != NC_NOERR || status2 != NC_NOERR) {
    status = EX_FATAL;
  }
  EX_FUNC_LEAVE(status);
}
