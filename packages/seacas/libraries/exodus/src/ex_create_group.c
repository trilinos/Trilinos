/*
 * Copyright (c) 2013 Sandia Corporation. Under the terms of Contract
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

#include "exodusII.h"     // for exerrval, ex_err, etc
#include "exodusII_int.h" // for EX_FATAL
#include "netcdf.h"       // for NC_NOERR, nc_def_grp, etc
#include <stdio.h>

int ex_create_group(int parent_id, const char *group_name)
{
  char errmsg[MAX_ERR_LENGTH];
#if NC_HAS_HDF5
  int exoid = -1;
  int status;

  ex_check_valid_file_id(parent_id);

  if ((status = nc_redef(parent_id)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to put file id %d into define mode", parent_id);
    ex_err("ex_create_group", errmsg, exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_def_grp(parent_id, group_name, &exoid)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: group create failed for %s in file id %d", group_name,
             parent_id);
    ex_err("ex_create_group", errmsg, exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_enddef(parent_id)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to complete definition for file id %d", exoid);
    ex_err("ex_create", errmsg, exerrval);
    return (EX_FATAL);
  }
  return (exoid);
#else
  exerrval = NC_ENOTNC4;
  snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: Group capabilities are not available in this netcdf "
                                   "version--not netcdf4");
  ex_err("ex_create_group", errmsg, exerrval);
  return (EX_FATAL);
#endif
}
