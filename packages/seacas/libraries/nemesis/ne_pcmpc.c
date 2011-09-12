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
 *     ne_put_cmap_params_cc()
 *****************************************************************************
 * This function outputs the concantenated list of communication map
 * parameters.
 *****************************************************************************
 *  Variable Index:
 *      neid                - The NetCDF ID of an already open NemesisI file.
 *      node_cmap_ids       - Pointer to vector of nodal communication
 *                            set IDs.
 *      node_cmap_node_cnts - Pointer to a vector which contains a count of
 *                            the number of FEM nodes for each nodal
 *                            communication map.
 *      elem_cmap_ids       - Pointer to vector for retrieval of elemental
 *                            communication set IDs.
 *      elem_cmap_elem_cnts - Pointer to a vector which contains a count of
 *                            the number of FEM elements for each elemental
 *                            communication map.
 *      proc_ids            - The processor the file being read was written
 *                            for.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <netcdf.h>

#include <exodusII.h>
#include <exodusII_int.h>

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

int ne_put_cmap_params_cc(int  neid,
                          int *node_cmap_ids,
                          int *node_cmap_node_cnts,
                          int *node_proc_ptrs,
                          int *elem_cmap_ids,
                          int *elem_cmap_elem_cnts,
                          int *elem_proc_ptrs
                          )
{
  char   *func_name="ne_put_cmap_params_cc";

  size_t  num_n_comm_maps, num_e_comm_maps, num_procs_in_file;
  int     status, icm, n_varid[2], e_varid[2], iproc;
  int     varid, n_dimid[1], e_dimid[1];
  int     n_varid_idx, e_varid_idx;
  int     num_icm;
  size_t start[1], count[1];
  size_t ecnt_cmap, ncnt_cmap;
  int64_t  nl_ecnt_cmap, nl_ncnt_cmap;
  int64_t *n_var_idx = NULL;
  int64_t *e_var_idx = NULL;
  int  nmstat;

  char    errmsg[MAX_ERR_LENGTH];
  /*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* Get the number of processors in the file */
  if ((status = nc_inq_dimid(neid, DIM_NUM_PROCS_F, &n_dimid[0])) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to get dimension ID for \"%s\" in file ID %d",
            DIM_NUM_PROCS_F, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_dimlen(neid, n_dimid[0], &num_procs_in_file)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find length of dimension \"%s\" in file ID %d",
            DIM_NUM_PROCS_F, neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /*
   * since I cannot get variables while in define mode, I need to
   * get the cmap information index variables before I go into
   * define mode
   */

  /* Check to see if there are nodal communications maps in the file */
  if (nc_inq_dimid(neid, DIM_NUM_N_CMAPS, &n_dimid[0]) != NC_NOERR) {
    num_n_comm_maps = 0;
  }
  else {
    if ((status = nc_inq_dimlen(neid, n_dimid[0], &num_n_comm_maps)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find length of dimension \"%s\" in \
file ID %d",
              DIM_NUM_N_CMAPS, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  }

  if (num_n_comm_maps > 0) {
    /* Get the variable ID for the comm map index vector */
    if ((status = nc_inq_varid(neid, VAR_N_COMM_INFO_IDX, &n_varid_idx)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable ID for \"%s\" in file ID %d",
              VAR_N_COMM_INFO_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* allocate space for the index variable */
    n_var_idx = (int64_t*) malloc((num_procs_in_file + 1) * sizeof(int64_t));
    if (!n_var_idx) {
      exerrval = EX_MSG;
      sprintf(errmsg,
              "Error: insufficient memory to read index variable from file ID %d",
              neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* and set the last value of the index */
    n_var_idx[0] = 0;

    /* get the communication map info index */
    if ((status = nc_get_var_longlong(neid, n_varid_idx, &(n_var_idx[1]))) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to get variable \"%s\" from file ID %d",
	      VAR_N_COMM_INFO_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  } /* "if (num_n_comm_maps > 0)" */

    /* Check to see if there are elemental communications maps in the file */
  if ((status = nc_inq_dimid(neid, DIM_NUM_E_CMAPS, &e_dimid[0])) != NC_NOERR) {
    num_e_comm_maps = 0;
  }
  else {
    if ((status = nc_inq_dimlen(neid, e_dimid[0], &num_e_comm_maps)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find length of dimension \"%s\" in \
file ID %d",
	      DIM_NUM_E_CMAPS, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  }

  if (num_e_comm_maps > 0) {
    /* Get the variable ID for the comm map index vector */
    if ((status = nc_inq_varid(neid, VAR_E_COMM_INFO_IDX, &e_varid_idx)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find variable ID for \"%s\" in file ID %d",
	      VAR_E_COMM_INFO_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* allocate space for the index variable */
    e_var_idx = (int64_t*) malloc((num_procs_in_file + 1) * sizeof(int64_t));
    if (!e_var_idx) {
      exerrval = EX_MSG;
      sprintf(errmsg,
	      "Error: insufficient memory to read index variable from file ID %d",
	      neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* and set the first value of the index */
    e_var_idx[0] = 0;

    /* get the communication map info index */
    if ((status = nc_get_var_longlong(neid, e_varid_idx, &(e_var_idx[1]))) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to get variable \"%s\" from file ID %d",
	      VAR_E_COMM_INFO_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  } /* "if (num_e_comm_maps >0)" */

    /* Put NetCDF file into define mode */
  if ((status = nc_redef(neid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
	    "Error: failed to put file ID %d into define mode", neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /*
   * Add dimensions for the size of the number of nodal
   * communication maps.
   */
  if (num_n_comm_maps > 0) {
    /* add the communications data index variable */
    if ((status = nc_def_var(neid, VAR_N_COMM_DATA_IDX, NC_INT64, 1,
			     n_dimid, &n_varid_idx)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to add variable \"%s\" in file ID %d",
	      VAR_N_COMM_DATA_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* now add up all of the nodal communications maps */
    ncnt_cmap = 0;
    for(iproc=0; iproc < num_procs_in_file; iproc++) {
      num_icm = n_var_idx[iproc+1] - n_var_idx[iproc];
      for(icm=0; icm < num_icm; icm++)
	ncnt_cmap += node_cmap_node_cnts[node_proc_ptrs[iproc]+icm];
    }

    if ((status = nc_def_dim(neid, DIM_NCNT_CMAP, ncnt_cmap, &n_dimid[0])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to add dimension for \"%s\" in file ID %d",
	      DIM_NCNT_CMAP, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* Define variables for the nodal IDS and processor vectors */
    if ((status = nc_def_var(neid, VAR_N_COMM_NIDS, NC_INT, 1, n_dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to add variable \"%s\" in file ID %d",
	      VAR_N_COMM_NIDS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_N_COMM_PROC, NC_INT, 1, n_dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to add variable \"%s\" in file ID %d",
	      VAR_N_COMM_PROC, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

  } /* End "if (num_n_comm_maps > 0)" */

    /*
     * Add dimensions for the size of the number of elemental
     * communication maps.
     */
  if (num_e_comm_maps > 0) {
    /* add the communications data index variable */
    if ((status = nc_def_var(neid, VAR_E_COMM_DATA_IDX, NC_INT64, 1,
			     e_dimid, &e_varid_idx)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to add variable \"%s\" in file ID %d",
	      VAR_E_COMM_DATA_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* now add up all of the nodal communications maps */
    ecnt_cmap = 0;
    for(iproc=0; iproc < num_procs_in_file; iproc++) {
      num_icm = e_var_idx[iproc+1] - e_var_idx[iproc];
      for(icm=0; icm < num_icm; icm++)
	ecnt_cmap += elem_cmap_elem_cnts[elem_proc_ptrs[iproc]+icm];
    }

    /* Add dimensions for elemental communications maps */
    if ((status = nc_def_dim(neid, DIM_ECNT_CMAP, ecnt_cmap, &e_dimid[0])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to add dimension for \"%s\" in file ID %d",
	      DIM_ECNT_CMAP, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* Define variables for the nodal IDS and processor vectors */
    if ((status = nc_def_var(neid, VAR_E_COMM_EIDS, NC_INT, 1, e_dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to add variable \"%s\" in file ID %d",
	      VAR_E_COMM_EIDS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_E_COMM_PROC, NC_INT, 1, e_dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to add variable \"%s\" in file ID %d",
	      VAR_E_COMM_PROC, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_E_COMM_SIDS, NC_INT, 1, e_dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to add variable \"%s\" in file ID %d",
	      VAR_E_COMM_SIDS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

  } /* End "if (num_e_comm_maps > 0)" */


    /* Exit define mode */
  ne_leavedef(neid, func_name);

  /* Set the status of the nodal communication maps */
  if (num_n_comm_maps > 0) {

    /* need to get the two "n_comm_*" variable ids */

    if ((status = nc_inq_varid(neid, VAR_N_COMM_STAT, &n_varid[0])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find variable ID for \"%s\" in file ID %d",
	      VAR_N_COMM_STAT, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* Get the variable ID for the comm map IDs vector */
    if ((status = nc_inq_varid(neid, VAR_N_COMM_IDS, &n_varid[1])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find variable ID for \"%s\" in file ID %d",
	      VAR_N_COMM_IDS, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* reset the index variable */
    nl_ncnt_cmap = 0;

    for(iproc=0; iproc < num_procs_in_file; iproc++) {
      num_icm = n_var_idx[iproc+1] - n_var_idx[iproc];
      for(icm=0; icm < num_icm; icm++) {

	start[0] = n_var_idx[iproc] + icm;
	if (node_cmap_node_cnts[node_proc_ptrs[iproc]+icm] > 0)
	  nmstat = 1;
	else
	  nmstat = 0;

	if ((status = nc_put_var1_int(neid, n_varid[0], start, &nmstat)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: unable to output variable in file ID %d", neid);
	  ex_err(func_name, errmsg, exerrval);
	  return (EX_FATAL);
	}

	/* increment to the next starting position */
	nl_ncnt_cmap += node_cmap_node_cnts[node_proc_ptrs[iproc]+icm];

	/* fill the data index variable */
	if ((status = nc_put_var1_longlong(neid, n_varid_idx, start, &nl_ncnt_cmap)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: failed to output int elem map index in file ID %d",
		  neid);
	  ex_err(func_name, errmsg, exerrval);
	  return (EX_FATAL);
	}
      } /* End "for(icm=0; icm < num_icm; icm++)" */

      if (num_icm > 0) {
	/* Output the nodal comm map IDs */
	start[0] = n_var_idx[iproc];
	count[0] =  num_icm;
	status = nc_put_vara_int(neid, n_varid[1], start, count,
			 &node_cmap_ids[node_proc_ptrs[iproc]]);
	if (status != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: failed to output variable in file ID %d", neid);
	  ex_err(func_name, errmsg, exerrval);

	  return (EX_FATAL);
	}
      }
    } /* End "for(iproc=0; iproc < num_procs_in_file; iproc++)" */

      /* free up memory for index */
    free(n_var_idx);

  } /* End "if (num_n_comm_maps > 0)" */


    /* Set the status of the elemental communication maps */
  if (num_e_comm_maps > 0) {

    /* need to get the two "e_comm_*" variables" */

    /* Get variable ID for elemental status vector */
    if ((status = nc_inq_varid(neid, VAR_E_COMM_STAT, &e_varid[0])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find variable ID for \"%s\" in file ID %d",
	      VAR_E_COMM_STAT, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* Get the variable ID for the elemental comm map IDs vector */
    if ((status = nc_inq_varid(neid, VAR_E_COMM_IDS, &e_varid[1])) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
	      "Error: failed to find variable ID for \"%s\" in file ID %d",
	      VAR_E_COMM_IDS, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* reset the index variable */
    nl_ecnt_cmap = 0;

    for(iproc=0; iproc < num_procs_in_file; iproc++) {
      num_icm = e_var_idx[iproc+1] - e_var_idx[iproc];
      for(icm=0; icm < num_icm; icm++) {

	start[0] = e_var_idx[iproc] + icm;
	if (elem_cmap_elem_cnts[elem_proc_ptrs[iproc]+icm] > 0)
	  nmstat = 1;
	else
	  nmstat = 0;

	if ((status = nc_put_var1_int(neid, e_varid[0], start, &nmstat)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: unable to output variable in file ID %d", neid);
	  ex_err(func_name, errmsg, exerrval);
	  return (EX_FATAL);
	}

	/* increment to the next starting position */
	nl_ecnt_cmap += elem_cmap_elem_cnts[elem_proc_ptrs[iproc]+icm];

	/* fill the data index variable */
	if ((status = nc_put_var1_longlong(neid, e_varid_idx, start, &nl_ecnt_cmap)) != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: failed to output int elem map index in file ID %d",
		  neid);
	  ex_err(func_name, errmsg, exerrval);
	  return (EX_FATAL);
	}
      } /* End "for(icm=0; icm < num_icm; icm++)" */

      if (num_icm > 0) {
	/* Output the elemental comm map IDs */
	start[0] = e_var_idx[iproc];
	count[0] = num_icm;
	status = nc_put_vara_int(neid, e_varid[1], start, count,
			 &elem_cmap_ids[elem_proc_ptrs[iproc]]);
	if (status != NC_NOERR) {
	  exerrval = status;
	  sprintf(errmsg,
		  "Error: failed to output variable in file ID %d", neid);
	  ex_err(func_name, errmsg, exerrval);
	  return (EX_FATAL);
	}
      }
    } /* End "for(iproc=0; iproc < num_procs_in_file; iproc++)" */

    free(e_var_idx);

  } /* End "if (num_e_comm_maps > 0)" */

  return (EX_NOERR);
}
