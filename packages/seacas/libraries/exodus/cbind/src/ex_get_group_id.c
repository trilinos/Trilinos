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

#include <stdio.h>                      // for sprintf
#include <stdlib.h>                     // for NULL
#include <string.h>                     // for strchr
#include "exodusII.h"                   // for exerrval, ex_err, etc
#include "exodusII_int.h"               // for EX_FATAL, EX_FILE_ID_MASK, etc
#include "netcdf.h"                     // for NC_NOERR, etc

/**
 * Given an exoid and group name (NULL gets root group), return id of that group.
 * If the name is NULL, return the root group.
 * If the name contains "/", then the name is assumed to be a full path name
 * and all groups in the file are searched.
 * Otherwise, the name is assumed to be the name of a child group of exoid
 */
int ex_get_group_id (int parent_id, const char *group_name, int *group_id)
{
  int status;
  char errmsg[MAX_ERR_LENGTH];
   
  exerrval = 0; /* clear error code */

#if defined(NOT_NETCDF4)
  exerrval = NC_ENOTNC4;
  sprintf(errmsg,
	  "Error: Group capabilities are not available in this netcdf version--not netcdf4");
  ex_err("ex_get_group_id",errmsg,exerrval);
  return (EX_FATAL);
#else
  /* See if name contains "/" indicating it is a full path name... */
  if (group_name == NULL) {
    /* Return root */
    *group_id = (unsigned)parent_id & EX_FILE_ID_MASK;
  }
  else if (strchr(group_name, '/') == NULL) {
    /* Local child */
    status = nc_inq_grp_ncid(parent_id, group_name, group_id);
    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: Failed to locate group with name %s as child group in file id %d",
	      group_name, parent_id);
      ex_err("ex_get_group_id",errmsg,exerrval);
      return (EX_FATAL);
    }
  }
  else {
    /* Full path name */
    status = nc_inq_grp_full_ncid(parent_id, group_name, group_id);
    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: Failed to locate group with full path name %s in file id %d",
	      group_name, parent_id);
      ex_err("ex_get_group_id",errmsg,exerrval);
      return (EX_FATAL);
    }
  }
  return (EX_NOERR);
#endif
}
