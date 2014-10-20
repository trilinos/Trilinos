/*
 * Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
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
/*****************************************************************************/
/*****************************************************************************/
/* Function(s) contained in this file:
 * 	ex_get_init_info()
 *****************************************************************************
 * This function reads information about the processors for which the
 * decomposition was performed.
 *****************************************************************************
 * Variable Index:
 *	exoid		  - The NetCDF ID of an already open NemesisI file.
 *	num_proc	  - The number of processors in the decomposition.
 *	num_proc_in_f	  - The number of processors the file contains
 *			    information for.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <exodusII.h>                   // for exerrval, ex_err, EX_MSG, etc
#include <exodusII_int.h>               // for EX_FATAL, DIM_NUM_PROCS, etc
#include <netcdf.h>                     // for NC_NOERR, nc_inq_dimid, etc
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for sprintf



int ex_get_init_info(int   exoid,
                     int  *num_proc,
                     int  *num_proc_in_f,
                     char *ftype
		     )
{
  const char *func_name="ex_get_init_info";

  int   dimid, status;
  size_t ltempsv;

  char   errmsg[MAX_ERR_LENGTH];
  /*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* Get the file type */
  if (ex_get_file_type(exoid, ftype) != EX_NOERR) {
    exerrval = EX_MSG;
    sprintf(errmsg,
            "Error: failed to get file type for file ID %d",
            exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  if ((status = nc_inq_dimid(exoid, DIM_NUM_PROCS, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find dimension ID for \"%s\" in file ID %d",
            DIM_NUM_PROCS, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Get the value of the number of processors */
  if ((status = nc_inq_dimlen(exoid, dimid, &ltempsv)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find length of dimension \"%s\" in file ID %d",
            DIM_NUM_PROCS, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }
  *num_proc = ltempsv;

  /* Get the dimension ID of processors that have info in this file */
  if ((status = nc_inq_dimid(exoid, DIM_NUM_PROCS_F, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find dimension ID for \"%s\" in file ID %d",
            DIM_NUM_PROCS_F, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Get the value of the number of processors that have info in this file */
  if ((status = nc_inq_dimlen(exoid, dimid, &ltempsv)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find length of dimension \"%s\" in file ID %d",
            DIM_NUM_PROCS_F, exoid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }
  *num_proc_in_f = ltempsv;

  return (EX_NOERR);
}
