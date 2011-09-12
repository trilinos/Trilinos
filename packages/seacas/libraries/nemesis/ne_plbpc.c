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
/*****************************************************************************/
/* Function(s) contained in this file:
 *      ne_put_loadbal_param_cc()
 *****************************************************************************
 * This function outputs the concatenated list of load-balance parameters.
 *****************************************************************************
 *  Variable Index:
 *      neid             - The NetCDF ID of an already open NemesisI file.
 *      num_int_nodes    - Vector of number of internal FEM nodes for
 *			   "num_proc_in_f" processors.
 *      num_bor_nodes    - Vector of number of border FEM nodes for
 *			   "num_proc_in_f" processors.
 *      num_ext_nodes    - Vector of number of external FEM nodes for
 *			   "num_proc_in_f" processors.
 *      num_int_elems    - Vector of number of internal FEM elems for
 *			   "num_proc_in_f" processors.
 *      num_bor_elems    - Vector of number of border FEM elems for
 *			   "num_proc_in_f" processors.
 *      num_node_cmaps   - Vector of number of nodal communication maps
 *                         for "num_proc_in_f" processors.
 *      num_elem_cmaps   - Vector of number of elemental communication maps
 *                         for "num_proc_in_f" processors.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>

#include <exodusII.h>
#include <exodusII_int.h>
#include <netcdf.h>

#include "ne_nemesisI.h"
#include "ne_nemesisI_int.h"

int ne_put_loadbal_param_cc(int   neid,
                            int  *num_int_nodes,
                            int  *num_bor_nodes,
                            int  *num_ext_nodes,
                            int  *num_int_elems,
                            int  *num_bor_elems,
                            int  *num_node_cmaps,
                            int  *num_elem_cmaps
                            )
{
  char  *func_name="ne_put_loadbal_param_cc";

  int     status;
  int     iproc, varid, dimid_npf, dimid[3];
  int     num_proc, num_proc_in_f;
  int     varid_nm[3], varid_em[2];
  int     varid_idx[7] = {0, 0, 0, 0, 0, 0, 0};
  size_t  start[1], ltempsv;
  char    ftype[2];
  int  oldfill;
  int  num_int_elem = 0, num_int_node = 0, num_bor_elem = 0;
  int  num_bor_node = 0, num_ext_node = 0;
  int  num_n_cmaps = 0, num_e_cmaps = 0;
  int  nmstat;

  char   errmsg[MAX_ERR_LENGTH];
  /*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* Get the processor information from the file */
  if (ne_get_init_info(neid, &num_proc, &num_proc_in_f, ftype) != EX_NOERR) {
    exerrval = EX_MSG;
    sprintf(errmsg,
            "Error: Unable to get processor info from file ID %d",
            neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /*
   * Get the dimension ID for the number of processors storing
   * information in this file.
   */
  if ((status = nc_inq_dimid(neid, DIM_NUM_PROCS_F, &dimid_npf)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find dimension ID for \"%s\" in file ID %d",
            DIM_NUM_PROCS_F, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Put NetCDF file into define mode */
  if ((status = nc_redef(neid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to put file id %d into define mode", neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Set the fill mode */
  if ((status = nc_set_fill(neid, NC_NOFILL, &oldfill)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to put file ID %d into no-fill mode",
            neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Output the file version */
  if ((status=ne_put_version(neid)) < 0) return (status);

  /* Output the file type */
  if (nc_inq_varid(neid, VAR_FILE_TYPE, &varid) != NC_NOERR) {
    if ((status = nc_def_var(neid, VAR_FILE_TYPE, NC_INT, 0, NULL, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define file type in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);

      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }
  }

  /* Define the status variables for the nodal vectors */
  if (nc_inq_varid(neid, VAR_INT_N_STAT, &varid) != NC_NOERR) {
    if ((status = nc_def_var(neid, VAR_INT_N_STAT, NC_INT, 1, &dimid_npf, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_INT_N_STAT, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }
  }

  /* Set the dimension for status vectors */
  if (nc_inq_varid(neid, VAR_BOR_N_STAT, &varid) != NC_NOERR) {
    if ((status = nc_def_var(neid, VAR_BOR_N_STAT, NC_INT, 1, &dimid_npf, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_BOR_N_STAT, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }
  }

  if (nc_inq_varid(neid, VAR_EXT_N_STAT, &varid) != NC_NOERR) {
    if ((status = nc_def_var(neid, VAR_EXT_N_STAT, NC_INT, 1, &dimid_npf, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_EXT_N_STAT, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }
  }

  /* Define the variable IDs for the elemental status vectors */
  if (nc_inq_varid(neid, VAR_INT_E_STAT, &varid) != NC_NOERR) {
    if ((status = nc_def_var(neid, VAR_INT_E_STAT, NC_INT, 1, &dimid_npf, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_INT_E_STAT, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }
  }

  if (nc_inq_varid(neid, VAR_BOR_E_STAT, &varid) != NC_NOERR) {
    if ((status = nc_def_var(neid, VAR_BOR_E_STAT, NC_INT, 1, &dimid_npf, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: Failed to define variable \"%s\" in file ID %d",
              VAR_BOR_E_STAT, neid);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }
  }

  /* Get the variable ID for the nodal status vectors */
  if ((status = nc_inq_varid(neid, VAR_INT_N_STAT, &varid_nm[0])) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to find variable ID for \"%s\" in file ID %d",
	    VAR_INT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

    return (EX_FATAL);
  }

  if ((status = nc_inq_varid(neid, VAR_BOR_N_STAT, &varid_nm[1])) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_BOR_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

    return (EX_FATAL);
  }

  if ((status = nc_inq_varid(neid, VAR_EXT_N_STAT, &varid_nm[2])) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_EXT_N_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

    return (EX_FATAL);
  }

  /* Get the variable IDs for the elemental status vectors */
  if ((status = nc_inq_varid(neid, VAR_INT_E_STAT, &varid_em[0])) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_INT_E_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

    return (EX_FATAL);
  }

  if ((status = nc_inq_varid(neid, VAR_BOR_E_STAT, &varid_em[1])) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find variable ID for \"%s\" in file ID %d",
            VAR_BOR_E_STAT, neid);
    ex_err(func_name, errmsg, exerrval);
    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

    return (EX_FATAL);
  }

  /* mms: NEW STUFF HERE */
  /*
  ** first need to loop through the processors in this
  ** file and get the counts of the element and cmap lists
  */
  for(iproc=0; iproc < num_proc_in_f; iproc++) {
    num_int_elem += num_int_elems[iproc];
    num_int_node += num_int_nodes[iproc];
    num_bor_elem += num_bor_elems[iproc];
    num_bor_node += num_bor_nodes[iproc];
    num_ext_node += num_ext_nodes[iproc];
    num_e_cmaps += num_elem_cmaps[iproc];
    num_n_cmaps += num_node_cmaps[iproc];
  }

  /* Define variable for the internal element information */
  if (num_int_elem > 0) {
    ltempsv = num_int_elem;
    if ((status = nc_def_dim(neid, DIM_NUM_INT_ELEMS, ltempsv, &dimid[0])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to dimension \"%s\" in file id %d",
              DIM_NUM_INT_ELEMS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_ELEM_MAP_INT, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_ELEM_MAP_INT, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* and the index variable */
    if ((status = nc_def_var(neid, VAR_ELEM_MAP_INT_IDX, NC_INT64, 1,
			     &dimid_npf, &varid_idx[0])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_ELEM_MAP_INT_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }
  } /* End "if (num_int_elem > 0)" */

  /* Define variable for the border element information */
  if (num_bor_elem > 0) {
    ltempsv = num_bor_elem;
    if ((status = nc_def_dim(neid, DIM_NUM_BOR_ELEMS, ltempsv, &dimid[0])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to dimension \"%s\" in file id %d",
              DIM_NUM_BOR_ELEMS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_ELEM_MAP_BOR, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_ELEM_MAP_BOR, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* and the index variable */
    if ((status = nc_def_var(neid, VAR_ELEM_MAP_BOR_IDX, NC_INT64, 1,
			     &dimid_npf, &varid_idx[1])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_ELEM_MAP_BOR_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

  } /* End "if (num_bor_elem > 0)" */

  if (num_int_node > 0) {
    /* Define variable for vector of internal FEM node IDs */
    ltempsv = num_int_node;
    if ((status = nc_def_dim(neid, DIM_NUM_INT_NODES, ltempsv, &dimid[0])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to dimension \"%s\" in file id %d",
              DIM_NUM_INT_NODES, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_NODE_MAP_INT, NC_INT, 1, &dimid[0], &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_NODE_MAP_INT, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* and the index variable */
    if ((status = nc_def_var(neid, VAR_NODE_MAP_INT_IDX, NC_INT64, 1,
			     &dimid_npf, &varid_idx[2])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_NODE_MAP_INT_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

  } /* End "if (num_int_node > 0)" */

  if (num_bor_node > 0) {
    /* Define variable for vector of border FEM node IDs */
    ltempsv = num_bor_node;
    if ((status = nc_def_dim(neid, DIM_NUM_BOR_NODES, ltempsv, &dimid[1])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to dimension \"%s\" in file id %d",
              DIM_NUM_BOR_NODES, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_NODE_MAP_BOR, NC_INT, 1, &dimid[1], &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_NODE_MAP_BOR, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* and the index variable */
    if ((status = nc_def_var(neid, VAR_NODE_MAP_BOR_IDX, NC_INT64, 1,
			     &dimid_npf, &varid_idx[3])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_NODE_MAP_BOR_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

  } /* End "if (num_bor_node > 0)" */

  if (num_ext_node > 0) {
    /* Define dimension for vector of external FEM node IDs */
    ltempsv = num_ext_node;
    if ((status = nc_def_dim(neid, DIM_NUM_EXT_NODES, ltempsv, &dimid[2])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to dimension \"%s\" in file id %d",
              DIM_NUM_EXT_NODES, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_NODE_MAP_EXT, NC_INT, 1, &dimid[2], &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_NODE_MAP_EXT, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* and the index variable */
    if ((status = nc_def_var(neid, VAR_NODE_MAP_EXT_IDX, NC_INT64, 1,
			     &dimid_npf, &varid_idx[4])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_NODE_MAP_EXT_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

  } /* End "if (num_ext_node > 0)" */

  /* Output the communication map dimensions */
  if (num_n_cmaps > 0) {
    ltempsv = num_n_cmaps;
    if ((status = nc_def_dim(neid, DIM_NUM_N_CMAPS, ltempsv, &dimid[0])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to add dimension \"%s\" in file ID %d",
              DIM_NUM_N_CMAPS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* Add variables for communication maps */
    if ((status = nc_def_var(neid, VAR_N_COMM_IDS, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_N_COMM_IDS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_N_COMM_STAT, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_N_COMM_STAT, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* and the index variable */
    if ((status = nc_def_var(neid, VAR_N_COMM_INFO_IDX, NC_INT64, 1,
			     &dimid_npf, &varid_idx[5])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_N_COMM_INFO_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

  } /* End "if (num_n_cmaps > 0)" */

  if (num_e_cmaps > 0) {
    ltempsv = num_e_cmaps;
    if ((status = nc_def_dim(neid, DIM_NUM_E_CMAPS, ltempsv, &dimid[0])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to add dimension \"%s\" in file ID %d",
              DIM_NUM_E_CMAPS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* Add variables for elemental communication maps */
    if ((status = nc_def_var(neid, VAR_E_COMM_IDS, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_E_COMM_IDS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_E_COMM_STAT, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_E_COMM_STAT, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* and the index variable */
    if ((status = nc_def_var(neid, VAR_E_COMM_INFO_IDX, NC_INT64, 1,
			     &dimid_npf, &varid_idx[6])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to define variable \"%s\" in file ID %d",
              VAR_E_COMM_INFO_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

  } /* End "if (num_e_cmaps > 0)" */

  /* Leave define mode */
  if (ne_leavedef(neid, func_name) != EX_NOERR)
    return (EX_FATAL);

  /* need to reset these counters */
  num_int_elem = 0;
  num_int_node = 0;
  num_bor_elem = 0;
  num_bor_node = 0;
  num_ext_node = 0;
  num_n_cmaps = 0;
  num_e_cmaps = 0;

  /* Update the status vectors */
  for(iproc=0; iproc < num_proc_in_f; iproc++) {
    start[0] = iproc;

    if (num_int_nodes[iproc] > 0)
      nmstat = 1;
    else
      nmstat = 0;

    if ((status = nc_put_var1_int(neid, varid_nm[0], start, &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status int node map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    if (num_bor_nodes[iproc] > 0)
      nmstat = 1;
    else
      nmstat = 0;

    if ((status = nc_put_var1_int(neid, varid_nm[1], start, &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status bor node map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    if (num_ext_nodes[iproc] > 0)
      nmstat = 1;
    else
      nmstat = 0;

    if ((status = nc_put_var1_int(neid, varid_nm[2], start, &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status ext node map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    if (num_int_elems[iproc] > 0)
      nmstat = 1;
    else
      nmstat = 0;

    if ((status = nc_put_var1_int(neid, varid_em[0], start, &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status int elem map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    if (num_bor_elems[iproc] > 0)
      nmstat = 1;
    else
      nmstat = 0;

    if ((status = nc_put_var1_int(neid, varid_em[1], start, &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status bor elem map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* now fill the index variables */
    if (varid_idx[0] > 0) {
      /* increment to the next starting position */
      num_int_elem += num_int_elems[iproc];
      if ((status = nc_put_var1_int(neid, varid_idx[0], start, &num_int_elem)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to output int elem map index in file ID %d",
                neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }
    }

    if (varid_idx[1] > 0) {
      /* increment to the next starting position */
      num_bor_elem += num_bor_elems[iproc];
      if ((status = nc_put_var1_int(neid, varid_idx[1], start, &num_bor_elem)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to output bor elem map index in file ID %d",
                neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }
    }

    if (varid_idx[2] > 0) {
      /* increment to the next starting position */
      num_int_node += num_int_nodes[iproc];
      if ((status = nc_put_var1_int(neid, varid_idx[2], start, &num_int_node)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to output int node map index in file ID %d",
                neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }
    }

    if (varid_idx[3] > 0) {
      /* increment to the next starting position */
      num_bor_node += num_bor_nodes[iproc];
      if ((status = nc_put_var1_int(neid, varid_idx[3], start, &num_bor_node)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to output bor node map index in file ID %d",
                neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }
    }

    if (varid_idx[4] > 0) {
      /* increment to the next starting position */
      num_ext_node += num_ext_nodes[iproc];
      if ((status = nc_put_var1_int(neid, varid_idx[4], start, &num_ext_node)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to output ext node map index in file ID %d",
                neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }
    }

    if (varid_idx[5] > 0) {
      /* increment to the next starting position */
      num_n_cmaps += num_node_cmaps[iproc];
      if ((status = nc_put_var1_int(neid, varid_idx[5], start, &num_n_cmaps)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to output node comm map index in file ID %d",
                neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }
    }

    if (varid_idx[6] > 0) {
      /* increment to the next starting position */
      num_e_cmaps += num_elem_cmaps[iproc];
      if ((status = nc_put_var1_int(neid, varid_idx[6], start, &num_e_cmaps)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to output elem cmap index in file ID %d",
                neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }
    }

  } /* End "for(iproc=0; iproc < num_proc_in_f; iproc++)" */
  return (EX_NOERR);
}
