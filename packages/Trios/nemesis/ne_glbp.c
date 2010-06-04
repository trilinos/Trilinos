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
 *      ne_get_loadbal_param()
 *****************************************************************************
 * This function retrieves the load balance parameters.
 *****************************************************************************
 *  Variable Index:
 *      neid             - The NetCDF ID of an already open NemesisI file.
 *      num_int_nodes    - The number of internal FEM nodes.
 *      num_bor_nodes    - The number of border FEM nodes.
 *      num_ext_nodes    - The number of external FEM nodes.
 *      num_int_elems    - The number of internal FEM elements.
 *      num_bor_elems    - The number of border FEM elements.
 *      num_node_cmaps   - The number of nodal communication maps.
 *      num_elem_cmaps   - The number of elemental communication maps.
 *      processor        - The processor the file being read was written for.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>

#include <netcdf.h>

#include "exodusII.h"
#include "exodusII_int.h"

#include "ne_nemesisI.h"
#include "ne_nemesisI_int.h"

int ne_get_loadbal_param(int   neid,
                         int  *num_int_nodes,
                         int  *num_bor_nodes,
                         int  *num_ext_nodes,
                         int  *num_int_elems,
                         int  *num_bor_elems,
                         int  *num_node_cmaps,
                         int  *num_elem_cmaps,
                         int   processor
                         )
{
  char  *func_name="ne_get_loadbal_param";

  int    dimid, varid, status;
  size_t start[1], ltempsv, ltempsv2;
  size_t varidx[2];
  char   ftype[2];
  int nmstat;

  char   errmsg[MAX_ERR_LENGTH];
  /*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* Check the file version information */
  if ((dimid=ne_check_file_version(neid)) != EX_NOERR) return (dimid);

  /* Get the file type */
  if (ne_get_file_type(neid, ftype) != EX_NOERR) {
    exerrval = EX_MSG;
    sprintf(errmsg,
            "Error: unable to find file type for file ID %d",
            neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Get the status for this node map */
  if ((status = nc_inq_varid(neid, VAR_INT_N_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" from file ID %d",
            VAR_INT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  if (ftype[0] == 's')
    start[0] = processor;
  else
    start[0] = 0;

  if ((status = nc_get_var1_int(neid, varid, start, &nmstat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get status for \"%s\" from file ID %d",
	    VAR_INT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (nmstat == 1) {
    if (ne_get_idx(neid, VAR_NODE_MAP_INT_IDX, varidx, processor) == -1) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find index variable, \"%s\", in file ID %d",
	      VAR_NODE_MAP_INT_IDX, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    /* check if I need to get the dimension of the internal node map */
    if (varidx[1] == -1) {
      /* Get the dimension ID for the number of internal nodes */
      if ((status = nc_inq_dimid(neid, DIM_NUM_INT_NODES, &dimid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find dimension ID for \"%s\" in file ID %d",
		DIM_NUM_INT_NODES, neid);
	ex_err(func_name, errmsg, exerrval);

	return (EX_FATAL);
      }

      /*
       * Get the value of the dimension representing the total number of
       * internal FEM nodes.
       */
      if ((status = nc_inq_dimlen(neid, dimid, &ltempsv)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
		DIM_NUM_INT_NODES, neid);
	ex_err(func_name, errmsg, exerrval);

	return (EX_FATAL);
      }

      /* set the end value for the node map */
      varidx[1] = ltempsv;
    }  /* End "if (varidx[1] = -1)" */

    /* now get the number of nodes */
    *num_int_nodes = varidx[1] - varidx[0];
  }
  else { /* if (nmstat != 1) */
    *num_int_nodes = 0;

  } /* End "if (nmstat == 1)" */

  /* Get the status for this node map */
  if ((status = nc_inq_varid(neid, VAR_BOR_N_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to find variable ID for \"%s\" from file ID %d",
	    VAR_BOR_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  if (ftype[0] == 's')
    start[0] = processor;
  else
    start[0] = 0;

  if ((status = nc_get_var1_int(neid, varid, start, &nmstat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get status for \"%s\" from file ID %d",
	    VAR_BOR_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (nmstat == 1) {
    if (ne_get_idx(neid, VAR_NODE_MAP_BOR_IDX, varidx, processor) == -1) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find index variable, \"%s\", in file ID %d",
	      VAR_NODE_MAP_BOR_IDX, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    /* check if I need to get the dimension of the border node map */
    if (varidx[1] == -1) {
      /* Get the dimension ID for the number of border nodes */
      if ((status = nc_inq_dimid(neid, DIM_NUM_BOR_NODES, &dimid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find dimension ID for \"%s\" in file ID %d",
		DIM_NUM_BOR_NODES, neid);
	ex_err(func_name, errmsg, exerrval);

	return (EX_FATAL);
      }

      /*
       * Get the value of the dimension representing the number of border
       * FEM nodes.
       */
      if ((status = nc_inq_dimlen(neid, dimid, &ltempsv)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
		DIM_NUM_BOR_NODES, neid);
	ex_err(func_name, errmsg, exerrval);

	return (EX_FATAL);
      }

      /* set the end value for the node map */
      varidx[1] = ltempsv;
    }  /* End "if (varidx[1] == -1)" */

    /* now calculate the number of nodes */
    *num_bor_nodes = varidx[1] - varidx[0];
  }
  else { /* if (nmstat != 1) */
    *num_bor_nodes = 0;
  }

  /* Get the status for this node map */
  if ((status = nc_inq_varid(neid, VAR_EXT_N_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to find variable ID for \"%s\" from file ID %d",
	    VAR_EXT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  if (ftype[0] == 's')
    start[0] = processor;
  else
    start[0] = 0;

  if ((status = nc_get_var1_int(neid, varid, start, &nmstat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get status for \"%s\" from file ID %d",
	    VAR_EXT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (nmstat == 1) {
    if (ne_get_idx(neid, VAR_NODE_MAP_EXT_IDX, varidx, processor) == -1) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find index variable, \"%s\", in file ID %d",
	      VAR_NODE_MAP_EXT_IDX, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    /* check if I need to get the dimension of the external node map */
    if (varidx[1] == -1) {
      /* Get the dimension ID for the number of external nodes */
      if ((status = nc_inq_dimid(neid, DIM_NUM_EXT_NODES, &dimid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find dimension ID for \"%s\" in file ID %d",
		DIM_NUM_EXT_NODES, neid);
	ex_err(func_name, errmsg, exerrval);

	return (EX_FATAL);
      }

      /*
       * Get the value of the dimension representing the number of external
       * FEM nodes.
       */
      if ((status = nc_inq_dimlen(neid, dimid, &ltempsv)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
		DIM_NUM_EXT_NODES, neid);
	ex_err(func_name, errmsg, exerrval);

	return (EX_FATAL);
      }
      /* set the end value for the node map */
      varidx[1] = ltempsv;
    }  /* End "if (varidx[1] == -1)" */

    /* now get the number of nodes */
    *num_ext_nodes = varidx[1] - varidx[0];
  }
  else { /* if (nmstat != 1) */
    *num_ext_nodes = 0;

  } /* End "if (nmstat == 1)" */

  /* Get the status for this element map */
  if ((status = nc_inq_varid(neid, VAR_INT_E_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to find variable ID for \"%s\" from file ID %d",
	    VAR_INT_E_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (ftype[0] == 's')
    start[0] = processor;
  else
    start[0] = 0;

  if ((status = nc_get_var1_int(neid, varid, start, &nmstat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get status for \"%s\" from file ID %d",
	    VAR_INT_E_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (nmstat == 1) {
    if (ne_get_idx(neid, VAR_ELEM_MAP_INT_IDX, varidx, processor) == -1) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find index variable, \"%s\", in file ID %d",
	      VAR_ELEM_MAP_INT_IDX, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    /* check if I need to get the dimension of the internal element map */
    if (varidx[1] == -1) {
      /* Get the dimension ID for the number of internal elements */
      if ((status = nc_inq_dimid(neid, DIM_NUM_INT_ELEMS, &dimid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find dimension ID for \"%s\" from file ID %d",
		DIM_NUM_INT_ELEMS, neid);
	ex_err(func_name, errmsg, exerrval);

	return (EX_FATAL);
      }

      /*
       * Get the value of the dimension representing the number of internal
       * FEM elements.
       */
      if ((status = nc_inq_dimlen(neid, dimid, &ltempsv)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find length of dimesion \"%s\" in file ID %d",
		DIM_NUM_INT_ELEMS, neid);
	ex_err(func_name, errmsg, exerrval);
	return (EX_FATAL);
      }

      /* set the end value for the node map */
      varidx[1] = ltempsv;
    }  /* End "if (varidx[1] == -1)" */

    /* now get the number of elements */
    *num_int_elems = varidx[1] - varidx[0];
  }
  else { /* if (nmstat != 1) */
    *num_int_elems = 0;

  } /* End "if (nmstat == 1)" */

  /* Get the status for this element map */
  if ((status = nc_inq_varid(neid, VAR_BOR_E_STAT, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to find variable ID for \"%s\" from file ID %d",
	    VAR_BOR_E_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (ftype[0] == 's')
    start[0] = processor;
  else
    start[0] = 0;

  if ((status = nc_get_var1_int(neid, varid, start, &nmstat)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to get status for \"%s\" from file ID %d",
	    VAR_BOR_E_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if (nmstat == 1) {
    if (ne_get_idx(neid, VAR_ELEM_MAP_BOR_IDX, varidx, processor) == -1) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find index variable, \"%s\", in file ID %d",
	      VAR_ELEM_MAP_BOR_IDX, neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

    /* check if I need to get the dimension of the border element map */
    if (varidx[1] == -1) {
      /* Get the dimension ID for the number of border elements */
      if ((status = nc_inq_dimid(neid, DIM_NUM_BOR_ELEMS, &dimid)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find dimension ID for \"%s\" from file ID %d",
		DIM_NUM_BOR_ELEMS, neid);
	ex_err(func_name, errmsg, exerrval);

	return (EX_FATAL);
      }

      /*
       * Get the value of the dimension representing the number of internal
       * FEM elements.
       */
      if ((status = nc_inq_dimlen(neid, dimid, &ltempsv)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find length of dimesion \"%s\" in file ID %d",
		DIM_NUM_BOR_ELEMS, neid);
	ex_err(func_name, errmsg, exerrval);
	return (EX_FATAL);
      }

      /* set the end value for the node map */
      varidx[1] = ltempsv;
    }  /* End "if (varidx[1] == -1)" */

    /* now get the number of nodes */
    *num_bor_elems = varidx[1] - varidx[0];
  }
  else { /* if (nmstat != 1) */
    *num_bor_elems = 0;

  } /* End "if (nmstat == 1)" */

  if (ne_get_idx(neid, VAR_N_COMM_INFO_IDX, varidx, processor) == -1) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to find index variable, \"%s\", in file ID %d",
	    VAR_N_COMM_INFO_IDX, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* check if I need to get the dimension of the nodal comm map */
  if (varidx[1] == -1) {
    /* Get the nodal comm map information */
    if ((status = nc_inq_dimid(neid, DIM_NUM_N_CMAPS, &dimid)) != NC_NOERR)
      varidx[1] = 0;
    else {
      if ((status = nc_inq_dimlen(neid, dimid, &ltempsv2)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
		DIM_NUM_N_CMAPS, neid);
	ex_err(func_name, errmsg, exerrval);
	return (EX_FATAL);
      }
      /* set the end value for the node map */
      varidx[1] = ltempsv2;
    }
  }  /* End "if (varidx[1] == -1)" */

  *num_node_cmaps = varidx[1] - varidx[0];

  if (ne_get_idx(neid, VAR_E_COMM_INFO_IDX, varidx, processor) == -1) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to find index variable, \"%s\", in file ID %d",
	    VAR_E_COMM_INFO_IDX, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* check if I need to get the dimension of the elemental comm map */
  if (varidx[1] == -1) {
    /* Get the elemental comm map information */
    if ((status = nc_inq_dimid(neid, DIM_NUM_E_CMAPS, &dimid)) != NC_NOERR)
      varidx[1] = 0;
    else {
      if ((status = nc_inq_dimlen(neid, dimid, &ltempsv2)) != NC_NOERR) {
	exerrval = status;
	sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
		DIM_NUM_E_CMAPS, neid);
	ex_err(func_name, errmsg, exerrval);
	return (EX_FATAL);
      }

      /* set the end value for the node map */
      varidx[1] = ltempsv2;
    }
  }  /* End "if (varidx[1] == -1)" */

  *num_elem_cmaps = varidx[1] - varidx[0];

  return (EX_NOERR);
}
