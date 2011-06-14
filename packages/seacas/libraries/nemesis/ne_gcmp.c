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
 *     ne_get_cmap_params()
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
#include <stdlib.h>

#include <netcdf.h>

#include <exodusII.h>
#include <exodusII_int.h>

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

int ne_get_cmap_params(int  neid,
                       int *node_cmap_ids,
                       int *node_cmap_node_cnts,
                       int *elem_cmap_ids,
                       int *elem_cmap_elem_cnts,
                       int  processor
		       )
{
  char   *func_name="ne_get_cmap_params";

  size_t  cnt, num_n_comm_maps, num_e_comm_maps, start[1], count[1];
  size_t  cmap_info_idx[2], cmap_data_idx[2];
  int     nmstat;
  int     status, map_idx, varid, dimid;

  char    errmsg[MAX_ERR_LENGTH];
  /*-----------------------------Execution begins-----------------------------*/

  exerrval = 0;	/* clear error code */

  /*****************************************************************************/
  /*****************************************************************************/
  /*                    Nodal communication map(s)                             */
  /*****************************************************************************/
  /*****************************************************************************/

  /* get the cmap information variables index */
  if (ne_get_idx(neid, VAR_N_COMM_INFO_IDX, cmap_info_idx, processor) == -1) {
    sprintf(errmsg,
            "Error: failed to find index variable, \"%s\", in file ID %d",
            VAR_N_COMM_INFO_IDX, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Get the number of nodal communications maps in the file */
  if ((status = nc_inq_dimid(neid, DIM_NUM_N_CMAPS, &dimid)) == NC_NOERR) {
    /* check if I need to get the dimension of the nodal comm map */
    if (cmap_info_idx[1] == -1) {
      if ((status = nc_inq_dimlen(neid, dimid, count)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
                DIM_NUM_N_CMAPS, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      /* set the end value for the node map */
      cmap_info_idx[1] = count[0];
    }  /* End "if (cmap_info_idx[1] == -1) */

    num_n_comm_maps = cmap_info_idx[1] - cmap_info_idx[0];

    if (num_n_comm_maps > 0) {
      count[0] = num_n_comm_maps;

      /* Get the variable ID for the vector of nodal comm map IDs */
      if ((status = nc_inq_varid(neid, VAR_N_COMM_IDS, &varid)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to get variable ID for \"%s\" in file ID %d",
                VAR_N_COMM_IDS, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      /* Get the vector of nodal communication map IDs */
      if (node_cmap_ids != NULL) {
        start[0] = cmap_info_idx[0];
	status = nc_get_vara_int(neid, varid, start, count, node_cmap_ids);

        if (status != NC_NOERR) {
          exerrval = status;
          sprintf(errmsg,
                  "Error: failed to get variable \"%s\" from file ID %d",
                  VAR_N_COMM_IDS, neid);
          ex_err(func_name, errmsg, exerrval);
          return (EX_FATAL);
        }

        if ((status = nc_inq_varid(neid, VAR_N_COMM_STAT, &varid)) != NC_NOERR) {
          exerrval = status;
          sprintf(errmsg,
		  "Error: failed to find variable ID for \"%s\" from file ID %d",
                  VAR_N_COMM_STAT, neid);
          ex_err(func_name, errmsg, exerrval);
          return (EX_FATAL);
        }

        if (node_cmap_node_cnts != NULL) {

          /* Get the node counts in each of the nodal communication maps */
          for(cnt=0; cnt < num_n_comm_maps; cnt++) {

            if ((map_idx=ne_id_lkup(neid, VAR_N_COMM_IDS, cmap_info_idx,
				    node_cmap_ids[cnt])) < 0) {
              exerrval = EX_MSG;
              sprintf(errmsg,
		      "Error: failed to find nodal comm map with ID %d in file ID %d",
                      node_cmap_ids[cnt], neid);
              ex_err(func_name, errmsg, exerrval);
              return (EX_FATAL);
            }

            /* Check the status of the node map */
            start[0] = map_idx;
            if ((status = nc_get_var1_int(neid, varid, start, &nmstat)) != NC_NOERR) {
              exerrval = status;
              sprintf(errmsg,
		      "Error: failed to get status for \"%s\" from file ID %d",
                      VAR_N_COMM_STAT, neid);
              ex_err(func_name, errmsg, exerrval);
              return (EX_FATAL);
            }

            if (nmstat == 1) {

              /* get the cmap information variables index */
              if (ne_get_idx(neid, VAR_N_COMM_DATA_IDX, cmap_data_idx,
                             map_idx) == -1) {
                exerrval = status;
                sprintf(errmsg,
			"Error: failed to find index variable, \"%s\", in file ID %d",
                        VAR_N_COMM_DATA_IDX, neid);
                ex_err(func_name, errmsg, exerrval);

                return (EX_FATAL);
              }

              if (cmap_data_idx[1] == -1) {
                /*
                 * Find the dimension ID of the variable containing the
                 * node count
                 */
                if ((status = nc_inq_dimid(neid, DIM_NCNT_CMAP, &dimid)) != NC_NOERR) {
                  exerrval = status;
                  sprintf(errmsg,
			  "Error: failed to find dimension ID for \"%s\" in file ID %d",
                          DIM_NCNT_CMAP, neid);
                  ex_err(func_name, errmsg, exerrval);
                  return (EX_FATAL);
                }

		/* Find the value of the number of nodes in this nodal comm map */
                if ((status = nc_inq_dimlen(neid, dimid, count)) != NC_NOERR) {
                  exerrval = status;
                  sprintf(errmsg,
			  "Error: failed to find length of dimension \"%s\" in file ID %d",
                          DIM_NCNT_CMAP, neid);
                  ex_err(func_name, errmsg, exerrval);
                  return (EX_FATAL);
                }

                cmap_data_idx[1] = count[0];
              }

              node_cmap_node_cnts[cnt] = cmap_data_idx[1] - cmap_data_idx[0];
            }
            else
              node_cmap_node_cnts[cnt] = 0;
          }  /* "for(cnt=0; cnt < num_n_comm_maps; cnt++)" */
        }  /* "if (node_cmap_node_cnts != NULL)" */
      }  /* "if (node_cmap_ids != NULL)" */
    }  /* "if (num_n_comm_maps > 0)" */
  } /* End "if ((dimid = nc_inq_dimid(neid, DIM_NUM_N_CMAPS)) != -1)" */

  /*****************************************************************************/
  /*****************************************************************************/
  /*                Elemental communication map(s)                             */
  /*****************************************************************************/
  /*****************************************************************************/

  /* get the cmap information variables index */
  if (ne_get_idx(neid, VAR_E_COMM_INFO_IDX, cmap_info_idx, processor) == -1) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to find index variable, \"%s\", in file ID %d",
            VAR_E_COMM_INFO_IDX, neid);
    ex_err(func_name, errmsg, exerrval);

    return (EX_FATAL);
  }

  /* Get the number of elemental communications maps in the file */
  if ((status = nc_inq_dimid(neid, DIM_NUM_E_CMAPS, &dimid)) == NC_NOERR) {
    /* check if I need to get the dimension of the nodal comm map */
    if (cmap_info_idx[1] == -1) {
      if ((status = nc_inq_dimlen(neid, dimid, count)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
		"Error: failed to find length of dimension \"%s\" in file ID %d",
                DIM_NUM_E_CMAPS, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      /* set the end value for the node map */
      cmap_info_idx[1] = count[0];
    }  /* End "if (cmap_info_idx[1] == -1) */

    num_e_comm_maps = cmap_info_idx[1] - cmap_info_idx[0];

    if (num_e_comm_maps > 0) {
      count[0] = num_e_comm_maps;

      /* Get the variable ID for the vector of nodal comm map IDs */
      if ((status = nc_inq_varid(neid, VAR_E_COMM_IDS, &varid)) != NC_NOERR) {
        exerrval = status;
        sprintf(errmsg,
                "Error: failed to get variable ID for \"%s\" in file ID %d",
                VAR_E_COMM_IDS, neid);
        ex_err(func_name, errmsg, exerrval);
        return (EX_FATAL);
      }

      /* Get the vector of elemental communication map IDs */
      if (elem_cmap_ids != NULL) {
        start[0] = cmap_info_idx[0];
	status = nc_get_vara_int(neid, varid, start, count, elem_cmap_ids);
        if (status != NC_NOERR) {
          exerrval = status;
          sprintf(errmsg,
                  "Error: failed to get variable \"%s\" from file ID %d",
                  VAR_E_COMM_IDS, neid);
          ex_err(func_name, errmsg, exerrval);
          return (EX_FATAL);
        }

        if ((status = nc_inq_varid(neid, VAR_E_COMM_STAT, &varid)) != NC_NOERR) {
          exerrval = status;
          sprintf(errmsg,
		  "Error: failed to find variable ID for \"%s\" from file ID %d",
                  VAR_E_COMM_STAT, neid);
          ex_err(func_name, errmsg, exerrval);
          return (EX_FATAL);
        }

        if (elem_cmap_elem_cnts != NULL) {
          /*
           * Get the element counts in each of the elemental
           * communication maps
           */
          for(cnt=0; cnt < num_e_comm_maps; cnt++) {

            if ((map_idx=ne_id_lkup(neid, VAR_E_COMM_IDS, cmap_info_idx,
				    elem_cmap_ids[cnt])) < 0) {
              exerrval = EX_MSG;
              sprintf(errmsg,
		      "Error: failed to find elemental comm map with ID %d in file ID %d",
                      elem_cmap_ids[cnt], neid);
              ex_err(func_name, errmsg, exerrval);
              return (EX_FATAL);
            }

            /* Check the status of the requested elemental map */
            start[0] = map_idx;
            if ((status = nc_get_var1_int(neid, varid, start, &nmstat)) != NC_NOERR) {
              exerrval = status;
              sprintf(errmsg,
                      "Error: failed to get status for \"%s\" from file ID %d",
                      VAR_E_COMM_STAT, neid);
              ex_err(func_name, errmsg, exerrval);
              return (EX_FATAL);
            }

            if (nmstat == 1) {

              /* get the cmap information variables index */
              if (ne_get_idx(neid, VAR_E_COMM_DATA_IDX, cmap_data_idx,
                             map_idx) == -1) {
                exerrval = status;
                sprintf(errmsg,
			"Error: failed to find index variable, \"%s\", in file ID %d",
			VAR_E_COMM_DATA_IDX, neid);
                ex_err(func_name, errmsg, exerrval);

                return (EX_FATAL);
              }

              if (cmap_data_idx[1] == -1) {
                /*
                 * Find the dimension ID of the variable containing the
                 * element count
                 */
                if ((status = nc_inq_dimid(neid, DIM_ECNT_CMAP, &dimid)) != NC_NOERR) {
                  exerrval = status;
                  sprintf(errmsg,
			  "Error: failed to find dimension ID for \"%s\" in file ID %d",
                          DIM_ECNT_CMAP, neid);
                  ex_err(func_name, errmsg, exerrval);
                  return (EX_FATAL);
                }

                /*
                 * Find the value of the number of elements in this elemental
                 * comm map
                 */
                if ((status = nc_inq_dimlen(neid, dimid, count)) != NC_NOERR) {
                  exerrval = status;
                  sprintf(errmsg,
			  "Error: failed to find length of dimension \"%s\" in file ID %d",
                          DIM_ECNT_CMAP, neid);
                  ex_err(func_name, errmsg, exerrval);
                  return (EX_FATAL);
                }
                cmap_data_idx[1] = count[0];
              }
              elem_cmap_elem_cnts[cnt] = cmap_data_idx[1] - cmap_data_idx[0];
            }
            else
              elem_cmap_elem_cnts[cnt] = 0;
          }  /* "for(cnt=0; cnt < num_e_comm_maps; cnt++)" */
        }  /* "if (elem_cmap_elem_cnts != NULL)" */
      }  /* "if (elem_cmap_ids != NULL)" */
    }  /* "if (num_e_comm_maps > 0)" */
  } /* End "if ((dimid = nc_inq_dimid(neid, DIM_NUM_E_CMAPS(processor))) != -1)" */
  return (EX_NOERR);
}
