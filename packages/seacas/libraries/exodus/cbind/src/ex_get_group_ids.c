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
#include "exodusII.h"                   // for exerrval, ex_err, etc
#include "exodusII_int.h"               // for EX_FATAL, EX_NOERR
#include "netcdf.h"                     // for nc_inq_grps, NC_NOERR

/**
 * Given a file or group 'parent' id, return the number of child groups and the ids
 * of the child groups below the parent.  If num_groups is NULL, do not return
 * count; if group_ids is NULL, do not return ids.
 */
int ex_get_group_ids (int parent_id, int *num_groups, int *group_ids)
{
  int status;
  char errmsg[MAX_ERR_LENGTH];
   
  exerrval = 0; /* clear error code */

#if defined(NOT_NETCDF4)
  exerrval = NC_ENOTNC4;
  sprintf(errmsg,
	  "Error: Group capabilities are not available in this netcdf version--not netcdf4");
  ex_err("ex_get_group_ids",errmsg,exerrval);
  return (EX_FATAL);
#else
  status = nc_inq_grps(parent_id, num_groups, group_ids);
  if (status != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: Failed to get child group ids in file id %d",
	    parent_id);
    ex_err("ex_get_group_ids",errmsg,exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
#endif
}
