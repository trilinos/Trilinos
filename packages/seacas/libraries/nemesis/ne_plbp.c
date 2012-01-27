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
 *      ne_put_loadbal_param()
 *****************************************************************************
 * This function outputs the load balance parameters.
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

#include <exodusII.h>
#include <exodusII_int.h>

#include "ne_nemesisI.h"
#include "ne_nemesisI_int.h"

int ne_put_loadbal_param(int   neid,
                         int   num_int_nodes,
                         int   num_bor_nodes,
                         int   num_ext_nodes,
                         int   num_int_elems,
                         int   num_bor_elems,
                         int   num_node_cmaps,
                         int   num_elem_cmaps,
                         int   processor
                         )
{
  char  *func_name="ne_put_loadbal_param";

  int    status, varid;
  int    dimid_npf, dimid[3];
  int    varid_nm[3], varid_em[2];
  char   ftype[2];

  int nmstat, ltempsv;

  char   errmsg[MAX_ERR_LENGTH];
  /*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* Get the file type */
  if (ne_get_file_type(neid, ftype) != EX_NOERR) {
    exerrval = EX_MSG;
    sprintf(errmsg,
            "Error: failed to get file type from file ID %d\n",
            neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* make sure that this is a parallel file */
  if (ftype[0] != 'p') {
    exerrval = EX_MSG;
    sprintf(errmsg,
            "Error: function for use with parallel files only, file ID %d\n",
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

    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

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

  /* Output the file version */
  if ((status=ne_put_version(neid)) < 0) return (status);

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

  /* Define variable for the internal element information */
  if (num_int_elems > 0) {
    ltempsv = num_int_elems;
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

  } /* End "if (num_int_elems > 0)" */

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

  /* Define variable for the border element information */
  if (num_bor_elems > 0) {
    ltempsv = num_bor_elems;
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

  } /* End "if (num_bor_elems > 0)" */

  if (num_int_nodes > 0) {
    /* Define variable for vector of internal FEM node IDs */
    ltempsv = num_int_nodes;
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

  } /* End "if (num_int_nodes > 0)" */

  if (num_bor_nodes > 0) {
    /* Define variable for vector of border FEM node IDs */
    ltempsv = num_bor_nodes;
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

  } /* End "if (num_bor_nodes > 0)" */

  if (num_ext_nodes > 0) {
    /* Define dimension for vector of external FEM node IDs */
    ltempsv = num_ext_nodes;
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

  } /* End "if (num_ext_nodes > 0)" */

  /* Add the nodal communication map count */
  if (num_node_cmaps > 0) {
    ltempsv = num_node_cmaps;
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

    /* Add the ID vector */
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

    /* Add the status vector */
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

  } /* End "if (num_node_cmaps > 0)" */

  /* Add the nodal communication map count */
  if (num_elem_cmaps > 0) {
    ltempsv = num_elem_cmaps;
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

  } /* End "if (num_elem_cmaps > 0)" */

  /* Leave define mode */
  if (ne_leavedef(neid, func_name) != EX_NOERR)
    return (EX_FATAL);

  /*
  ** Set up status vector for internal node map
  ** NOTE(9/26/96): this function is no longer valid
  ** for scaler files, so no need to check for file type
  */
  if (num_int_nodes == 0) {
    /* NULL set for internal nodes */
    nmstat = 0;
    if ((status = nc_put_var_int(neid, varid_nm[0], &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status for int node map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  }
  else {
    nmstat = 1;
    if ((status = nc_put_var_int(neid, varid_nm[0], &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status for int node map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

  } /* End "if (num_int_nodes == 0)" */

  /* Set up status vector for border node map */
  if (num_bor_nodes == 0) {
    /* NULL set for border nodes */
    nmstat = 0;
    if ((status = nc_put_var_int(neid, varid_nm[1], &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status for bor node map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  }
  else {
    /* Set the status indicating non-NULL size */
    nmstat = 1;
    if ((status = nc_put_var_int(neid, varid_nm[1], &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status for bor node map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

  } /* End "if (num_bor_nodes == 0)" */

  /* Set up status vector for external node map */
  if (num_ext_nodes == 0) {
    /* NULL set for external nodes */
    nmstat = 0;
    if ((status = nc_put_var_int(neid, varid_nm[2], &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status for ext node map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  }
  else {
    /* Set the status indicating non-NULL size */
    nmstat = 1;
    if ((status = nc_put_var_int(neid, varid_nm[2], &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status for ext node map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

  } /* End "if (num_ext_nodes == 0)" */

  /* Set up status vector for internal element map */
  if (num_int_elems == 0) {
    /* NULL set for internal elements */
    nmstat = 0;
    if ((status = nc_put_var_int(neid, varid_em[0], &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status for int elem map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  }
  else {
    /* Set the status indicating non-NULL size */
    nmstat = 1;
    if ((status = nc_put_var_int(neid, varid_em[0], &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status for int elem map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

  } /* End "if (num_int_elems == 0)" */

  /* Set up status vector for border element map */
  if (num_bor_elems == 0) {
    /* NULL set for internal elements */
    nmstat = 0;
    if ((status = nc_put_var_int(neid, varid_em[1], &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status for bor elem map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  }
  else {
    /* Set the status indicating non-NULL size */
    nmstat = 1;
    if ((status = nc_put_var_int(neid, varid_em[1], &nmstat)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output status for bor elem map in file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

  } /* End "if (num_bor_elems == 0)" */

  return (EX_NOERR);
}
