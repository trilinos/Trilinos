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
 *     ne_get_elem_map()
 *
 *****************************************************************************
 * Variable Index:
 *      neid            - The NetCDF ID of an already open NemesisI file.
 *      elem_mapi       - Pointer to vector for retrieval of internal
 *                        FEM element IDs.
 *      elem_mapb       - Pointer to vector for retrieval of border
 *                        FEM element IDs.
 *      processor       - The processor the file being read was written for.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <netcdf.h>

#include "exodusII.h"
#include "exodusII_int.h"

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

int ne_get_elem_map(int  neid,
                    int *elem_mapi,
                    int *elem_mapb,
                    int  processor
                    )
{
  char  *func_name="ne_get_elem_map";

  char    ftype[2];
  int     dimid, varid, status;
  size_t  start[1], count[1];
  size_t  varidx[2];
  int  emstat;

  char   errmsg[MAX_ERR_LENGTH];

  /* Get the file type */
  if (ne_get_file_type(neid, ftype) != EX_NOERR) {
    exerrval = EX_MSG;
    sprintf(errmsg,
            "Error: unable to find file type for file ID %d",
            neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Check the status of the internal element map */
  if ((status = nc_inq_varid(neid, VAR_INT_E_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_INT_E_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (ftype[0] == 'p')
    start[0] = 0;
  else
    start[0] = processor;

  if ((status = nc_get_var1_int(neid, varid, start, &emstat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get status for \"%s\" from file ID %d",
            VAR_INT_E_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (emstat == 1) {
    /* get the index */
    if (ne_get_idx(neid, VAR_ELEM_MAP_INT_IDX, varidx, processor) == -1) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find index variable, \"%s\", in file ID %d",
              VAR_ELEM_MAP_INT_IDX, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    if (varidx[1] == -1) {
      /* Get the size of the internal element map */
      if ((status = nc_inq_dimid(neid, DIM_NUM_INT_ELEMS, &dimid)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to find dimension ID for \"%s\" in file ID %d",
                DIM_NUM_INT_ELEMS, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      if ((status = nc_inq_dimlen(neid, dimid, count)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
                DIM_NUM_INT_ELEMS, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      varidx[1] = count[0];
    }

    /* Get the map */
    if ((status = nc_inq_varid(neid, VAR_ELEM_MAP_INT, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable ID for \"%s\" in file ID %d",
              VAR_ELEM_MAP_INT, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    start[0] = varidx[0];
    count[0] = varidx[1] - varidx[0];
    status = nc_get_vara_int(neid, varid, start, count, elem_mapi);

    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to get variable \"%s\" from file ID %d",
              VAR_ELEM_MAP_INT, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

  } /* End "if (emstat == 1)" */

  /* Check the status of the internal element map */
  if ((status = nc_inq_varid(neid, VAR_BOR_E_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_BOR_E_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (ftype[0] == 'p')
    start[0] = 0;
  else
    start[0] = processor;

  if ((status = nc_get_var1_int(neid, varid, start, &emstat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get status for \"%s\" from file ID %d",
	    VAR_INT_E_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (emstat == 1) {
    /* get the index */
    if (ne_get_idx(neid, VAR_ELEM_MAP_BOR_IDX, varidx, processor) == -1) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find index variable, \"%s\", in file ID %d",
	      VAR_ELEM_MAP_BOR_IDX, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    if (varidx[1] == -1) {
      /* Get the size of the border element map */
      if ((status = nc_inq_dimid(neid, DIM_NUM_BOR_ELEMS, &dimid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find dimension ID for \"%s\" in file ID %d",
		DIM_NUM_BOR_ELEMS, neid);
	ex_err(func_name, errmsg, exerrval);
	return (EX_FATAL);
      }

      if ((status = nc_inq_dimlen(neid, dimid, count)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
		DIM_NUM_BOR_ELEMS, neid);
	ex_err(func_name, errmsg, exerrval);
	return (EX_FATAL);
      }

      varidx[1] = count[0];
    }

    /* Get the map */
    if ((status = nc_inq_varid(neid, VAR_ELEM_MAP_BOR, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find variable ID for \"%s\" in file ID %d",
	      VAR_ELEM_MAP_BOR, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    start[0] = varidx[0];
    count[0] = varidx[1] - varidx[0];
    status = nc_get_vara_int(neid, varid, start, count, elem_mapb);

    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to get variable \"%s\" from file ID %d",
	      VAR_ELEM_MAP_BOR, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  } /* End "if (emstat == 1)" */

  return (EX_NOERR);
}
