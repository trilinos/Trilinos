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
 *      ne_get_n_coord()
 *****************************************************************************
 * This function retrieves "n" coordinate values.
 *****************************************************************************
 *  Variable Index:
 *      neid               - The NetCDF ID of an already open NemesisI file.
 *      start_node_num     - The starting index of the coords to be obtained.
 *      num_nodes          - The number of FEM nodes to read coords for.
 *      x_coor             - Pointer to location of X coordinate vector.
 *      y_coor             - Pointer to location of Y coordinate vector.
 *      z_coor             - Pointer to location of Z coordinate vector.
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

/*
 * reads the coordinates of the nodes
 */

int ne_get_n_coord (int   neid,
                    int   start_node_num,
                    int   num_nodes,
                    void *x_coor,
                    void *y_coor,
                    void *z_coor)
{
  int coordidx, coordidy, coordidz;
  int coordid;
  int status;
  char *which = NULL;

  int  ndimdim, i;
  size_t num_dim, start[2], count[2];
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0;

  if ((status = nc_inq_dimid(neid, DIM_NUM_DIM, &ndimdim)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find dimension ID for \"%s\" in file id %d",
	    DIM_NUM_DIM, neid);
    ex_err("ne_get_n_coord",errmsg,exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_dimlen(neid, ndimdim, &num_dim)) != NC_NOERR) {
    sprintf(errmsg,
            "Error: failed to get number of dimensions in file id %d",
	    neid);
    ex_err("ne_get_n_coord",errmsg,exerrval);
    return (EX_FATAL);
  }


  /* read in the coordinates  */
  if (ex_large_model(neid) == 0) {
    if ((status = nc_inq_varid (neid, VAR_COORD, &coordid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to locate nodal coordinates in file id %d", neid);
      ex_err("ne_get_n_coord",errmsg,exerrval);
      return (EX_FATAL);
    }

    start_node_num--; /* do this here so it is only decremented once */
    for (i=0; i<num_dim; i++) {
      start[0] = i;
      start[1] = start_node_num;

      count[0] = 1;
      count[1] = num_nodes;

      if (i == 0 && x_coor != NULL) {
	which = "X";
	if (ex_comp_ws(neid) == 4) {
	  status = nc_get_vara_float(neid, coordid, start, count, x_coor);
	} else {
	  status = nc_get_vara_double(neid, coordid, start, count, x_coor);
	}
      }
      else if (i == 1 && y_coor != NULL) {
	which = "Y";
	if (ex_comp_ws(neid) == 4) {
	  status = nc_get_vara_float(neid, coordid, start, count, y_coor);
	} else {
	  status = nc_get_vara_double(neid, coordid, start, count, y_coor);
	}
      }
      else if (i == 2 && z_coor != NULL) {
	which = "Z";
	if (ex_comp_ws(neid) == 4) {
	  status = nc_get_vara_float(neid, coordid, start, count, z_coor);
	} else {
	  status = nc_get_vara_double(neid, coordid, start, count, z_coor);
	}
      }
      if (status != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to get %s coord array in file id %d", which, neid);
	ex_err("ne_get_n_coord",errmsg,exerrval);
	return (EX_FATAL);
      }
    }
  } else {

    start_node_num--; /* do this here so it is only decremented once */
    for (i=0; i<num_dim; i++) {
      start[0] = start_node_num;
      count[0] = num_nodes;

      if (i == 0 && x_coor != NULL) {
	which = "X";
	if ((status = nc_inq_varid (neid, VAR_COORD_X, &coordidx)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Warning: failed to locate %s-nodal coordinates in file id %d", which , neid);
	  ex_err("ne_get_n_coord",errmsg,exerrval);
	  return (EX_WARN);
	}
	if (ex_comp_ws(neid) == 4) {
	  status = nc_get_vara_float(neid, coordidx, start, count, x_coor);
	} else {
	  status = nc_get_vara_double(neid, coordidx, start, count, x_coor);
	}
      }
      else if (i == 1 && y_coor != NULL) {
	which = "Y";
	if ((status = nc_inq_varid (neid, VAR_COORD_Y, &coordidy)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Warning: failed to locate %s-nodal coordinates in file id %d", which , neid);
	  ex_err("ne_get_n_coord",errmsg,exerrval);
	  return (EX_WARN);
	}
	if (ex_comp_ws(neid) == 4) {
	  status = nc_get_vara_float(neid, coordidy, start, count, y_coor);
	} else {
	  status = nc_get_vara_double(neid, coordidy, start, count, y_coor);
	}
      }
      else if (i == 2 && z_coor != NULL) {
	which = "Z";
	if ((status = nc_inq_varid (neid, VAR_COORD_Z, &coordidz)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Warning: failed to locate %s-nodal coordinates in file id %d", which , neid);
	  ex_err("ne_get_n_coord",errmsg,exerrval);
	  return (EX_WARN);
	}
	if (ex_comp_ws(neid) == 4) {
	  status = nc_get_vara_float(neid, coordidz, start, count, z_coor);
	} else {
	  status = nc_get_vara_double(neid, coordidz, start, count, z_coor);
	}
      }
      if (status != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to get %s coord array in file id %d", which, neid);
	ex_err("ne_get_n_coord",errmsg,exerrval);
	return (EX_FATAL);
      }
    }
  }
  return (EX_NOERR);
}
