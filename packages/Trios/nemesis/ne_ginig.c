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
/* Function(s) contained in this file:
 *     ne_get_init_global()
 *****************************************************************************
 * This function reads the global initial information.
 *****************************************************************************
 *  Variable Index:
 *      neid            - The NetCDF ID of an already open NemesisI file.
 *      num_nodes_g     - The number of global FEM nodes. This is output as
 *                        a NetCDF variable.
 *      num_elems_g     - The number of global FEM elements. This is output
 *                        as a NetCDF variable.
 *      num_elem_blks_g - The number of global element blocks. This is output
 *                        as a NetCDF dimension.
 *      num_node_sets_g - The number of global node sets. This is output as
 *                        a NetCDF dimension.
 *      num_side_sets_g - The number of global side sets. This is output as
 *                        a NetCDF dimension.
 */    
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>

#include <netcdf.h>

#include "exodusII.h"
#include "exodusII_int.h"

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

int ne_get_init_global(int   neid,
                       int  *num_nodes_g,
                       int  *num_elems_g,
                       int  *num_elem_blks_g,
                       int  *num_node_sets_g,
                       int  *num_side_sets_g
                       )
{
  char  *func_name="ne_get_init_global";
  int    dimid, status;
  size_t ltempsv;

  char    errmsg[MAX_ERR_LENGTH];
  /*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error flag */

  /* Check the file version information */
  if ((dimid=ne_check_file_version(neid)) != EX_NOERR) return (dimid);

  /* Get the dimension ID for the number of global FEM nodes */
  if ((status = nc_inq_dimid(neid, DIM_NUM_NODES_GLOBAL, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find dimension ID for \"%s\" in file ID %d",
            DIM_NUM_NODES_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    
    return (EX_FATAL);
  }

  /* Get the value of the number of global FEM nodes */
  if ((status = nc_inq_dimlen(neid, dimid, &ltempsv)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find length of dimension \"%s\" in file ID %d",
            DIM_NUM_NODES_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    
    return (EX_FATAL);
  }
  *num_nodes_g = ltempsv;

  /* Get the dimension ID for the number of global FEM elements */
  if ((status = nc_inq_dimid(neid, DIM_NUM_ELEMS_GLOBAL, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find dimension ID for \"%s\" in file ID %d",
            DIM_NUM_ELEMS_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    
    return (EX_FATAL);
  }

  /* Get the value of the number of global FEM elements */
  if ((status = nc_inq_dimlen(neid, dimid, &ltempsv)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find length of dimension \"%s\" in file ID %d",
            DIM_NUM_ELEMS_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    
    return (EX_FATAL);
  }
  *num_elems_g = ltempsv;

  /* Get the dimension ID for the number of global element blocks */
  if ((status = nc_inq_dimid(neid, DIM_NUM_ELBLK_GLOBAL, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find dimension ID for \"%s\" in file ID %d",
            DIM_NUM_ELBLK_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    
    return (EX_FATAL);
  }

  /* Get the value of the number of global element blocks */
  if ((status = nc_inq_dimlen(neid, dimid, &ltempsv)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to find length of dimension \"%s\" in file ID %d",
	    DIM_NUM_ELBLK_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    
    return (EX_FATAL);
  }
  *num_elem_blks_g = ltempsv;

  /* Get the dimension ID for the number of global node sets */
  if ((status = nc_inq_dimid(neid, DIM_NUM_NS_GLOBAL, &dimid)) != NC_NOERR) {
    *num_node_sets_g = (int) 0;
  }
  else {
    /* Get the value of the number of global node sets */
    if ((status = nc_inq_dimlen(neid, dimid, &ltempsv)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find length of dimension \"%s\" in file ID %d",
	      DIM_NUM_NS_GLOBAL, neid);
      ex_err(func_name, errmsg, exerrval);
    
      return (EX_FATAL);
    }
    *num_node_sets_g = ltempsv;
  }

  /* Get the dimension ID for the number of global side sets */
  if ((status = nc_inq_dimid(neid, DIM_NUM_SS_GLOBAL, &dimid)) != NC_NOERR) {
    *num_side_sets_g = (int) 0;
  }
  else {
    /* Get the value of the number of global side sets */
    if ((status = nc_inq_dimlen(neid, dimid, &ltempsv)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find length of dimension \"%s\" in file ID %d",
	      DIM_NUM_SS_GLOBAL, neid);
      ex_err(func_name, errmsg, exerrval);
      
      return (EX_FATAL);
    }
    *num_side_sets_g = ltempsv;
  }

  return (EX_NOERR);
}
