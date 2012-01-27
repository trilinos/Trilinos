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
 *     ne_put_cmap_params()
 *****************************************************************************
 * This function outputs the communication map parameters.
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
 *      processor           - The processor the file being read was written
 *                            for.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>
#include <stdint.h>

#include <netcdf.h>

#include <exodusII.h>
#include <exodusII_int.h>

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

int ne_put_cmap_params(int  neid,
                       int *node_cmap_ids,
                       int *node_cmap_node_cnts,
                       int *elem_cmap_ids,
                       int *elem_cmap_elem_cnts,
                       int  processor
                       )
{
  char   *func_name="ne_put_cmap_params";

  size_t  num_n_comm_maps, num_e_comm_maps;
  size_t  ncnt_cmap, ecnt_cmap;
  int     icm, varid, dimid[1], n_varid, e_varid, status;
  int     n_varid_idx, e_varid_idx;
  size_t  start[1];
  char    ftype[2];
  int  nl_ncnt_cmap, nl_ecnt_cmap;
  int  nmstat;

  char    errmsg[MAX_ERR_LENGTH];
/*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /*
  ** with the new database format, this function sould only
  ** be used for writing a parallel file
  */
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

  /* Put NetCDF file into define mode */
  if ((status = nc_redef(neid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to file ID %d into define mode", neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Check to see if there are nodal communications maps in the file */
  if ((status = nc_inq_dimid(neid, DIM_NUM_N_CMAPS, &dimid[0])) != NC_NOERR) {
    num_n_comm_maps = 0;
  }
  else {
    if ((status = nc_inq_dimlen(neid, dimid[0], &num_n_comm_maps)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find length of dimension \"%s\" in file ID %d",
              DIM_NUM_N_CMAPS, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }
  }

  /*
   * Add dimensions for the size of the number of nodal
   * communication maps.
   */
  if (num_n_comm_maps > 0) {

    /* add the communications data index variable */
    if ((status = nc_def_var(neid, VAR_N_COMM_DATA_IDX, NC_INT, 1,
			     dimid, &n_varid_idx)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to add variable \"%s\" in file ID %d",
              VAR_N_COMM_DATA_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* Add dimensions for all of the nodal communication maps */
    ncnt_cmap = 0;
    for(icm=0; icm < num_n_comm_maps; icm++) {
      ncnt_cmap += node_cmap_node_cnts[icm];
    }

    if ((status = nc_def_dim(neid, DIM_NCNT_CMAP, ncnt_cmap, &dimid[0])) != NC_NOERR) {
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
    if ((status = nc_def_var(neid, VAR_N_COMM_NIDS, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to add variable \"%s\" in file ID %d",
              VAR_N_COMM_NIDS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_N_COMM_PROC, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
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

  /* Check to see if there are elemental communications maps in the file */
  if ((status = nc_inq_dimid(neid, DIM_NUM_E_CMAPS, &dimid[0])) != NC_NOERR) {
    num_e_comm_maps = 0;
  }
  else{
    if ((status = nc_inq_dimlen(neid, dimid[0], &num_e_comm_maps)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find length of dimension \"%s\" in file ID %d",
              DIM_NUM_E_CMAPS, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

  }

  /*
   * Add dimensions for the size of the number of elemental
   * communication maps.
   */
  if (num_e_comm_maps > 0) {

    /* add the communications data index variable */
    if ((status = nc_def_var(neid, VAR_E_COMM_DATA_IDX, NC_INT, 1,
			     dimid, &e_varid_idx)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to add variable \"%s\" in file ID %d",
              VAR_E_COMM_DATA_IDX, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* Add dimensions for each of the nodal communication maps */
    ecnt_cmap = 0;
    for(icm=0; icm < num_e_comm_maps; icm++)
      ecnt_cmap += elem_cmap_elem_cnts[icm];

    if ((status = nc_def_dim(neid, DIM_ECNT_CMAP, ecnt_cmap, &dimid[0])) != NC_NOERR) {
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
    if ((status = nc_def_var(neid, VAR_E_COMM_EIDS, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to add variable \"%s\" in file ID %d",
              VAR_E_COMM_EIDS, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_E_COMM_PROC, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to add variable \"%s\" in file ID %d",
              VAR_E_COMM_PROC, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    if ((status = nc_def_var(neid, VAR_E_COMM_SIDS, NC_INT, 1, dimid, &varid)) != NC_NOERR) {
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

    if ((status = nc_inq_varid(neid, VAR_N_COMM_STAT, &n_varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable ID for \"%s\" in file ID %d",
              VAR_N_COMM_STAT, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    nl_ncnt_cmap = 0; /* reset this for index */
    for(icm=0; icm < num_n_comm_maps; icm++) {

      start[0] = icm;
      if (node_cmap_node_cnts[icm] > 0)
        nmstat = 1;
      else
        nmstat = 0;

      if ((status = nc_put_var1_int(neid, n_varid, start, &nmstat)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: unable to output variable in file ID %d", neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      /* increment to the next starting position */
      nl_ncnt_cmap += node_cmap_node_cnts[icm];

      /* fill the cmap data index */
      if ((status = nc_put_var1_int(neid, n_varid_idx, start, &nl_ncnt_cmap)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to output int elem map index in file ID %d",
                neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }
    } /* End "for(icm=0; icm < num_n_comm_maps; icm++)" */

    /* Get the variable ID for the comm map IDs vector */
    if ((status = nc_inq_varid(neid, VAR_N_COMM_IDS, &n_varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable ID for \"%s\" in file ID %d",
              VAR_N_COMM_IDS, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* Output the nodal comm map IDs */
    status  = nc_put_var_int(neid, n_varid, node_cmap_ids);
    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output variable in file ID %d", neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

  } /* End "if (num_n_comm_maps > 0)" */

  /* Set the status of the elemental communication maps */
  if (num_e_comm_maps > 0) {

    /* Get variable ID for elemental status vector */
    if ((status = nc_inq_varid(neid, VAR_E_COMM_STAT, &e_varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable ID for \"%s\" in file ID %d",
              VAR_E_COMM_STAT, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    nl_ecnt_cmap = 0; /* reset this for index */
    for(icm=0; icm < num_e_comm_maps; icm++) {

      start[0] = icm;
      if (elem_cmap_elem_cnts[icm] > 0)
        nmstat = 1;
      else
        nmstat = 0;

      if ((status = nc_put_var1_int(neid, e_varid, start, &nmstat)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: unable to output variable in file ID %d", neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      /* increment to the next starting position */
      nl_ecnt_cmap += elem_cmap_elem_cnts[icm];

      /* fill the cmap data index */
      if ((status = nc_put_var1_int(neid, e_varid_idx, start, &nl_ecnt_cmap)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to output int elem map index in file ID %d",
                neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }
     } /* End "for(icm=0; icm < num_e_comm_maps; icm++)" */

    /* Get the variable ID for the elemental comm map IDs vector */
    if ((status = nc_inq_varid(neid, VAR_E_COMM_IDS, &e_varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to find variable ID for \"%s\" in file ID %d",
              VAR_E_COMM_IDS, neid);
      ex_err(func_name, errmsg, exerrval);
      return (EX_FATAL);
    }

    /* Output the elemental comm map IDs */
    status = nc_put_var_int(neid, e_varid, elem_cmap_ids);
    if (status != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to output variable in file ID %d", neid);
      ex_err(func_name, errmsg, exerrval);

      return (EX_FATAL);
    }

  } /* End "if (num_e_comm_maps > 0)" */

  return (EX_NOERR);
}
