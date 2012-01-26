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
 *     ne_put_init_global()
 *****************************************************************************
 * This function outputs the initial global information.
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

#include <exodusII.h>
#include <exodusII_int.h>

#include "ne_nemesisI_int.h"
#include "ne_nemesisI.h"

int ne_put_init_global(int   neid,
                       int   num_nodes_g,
                       int   num_elems_g,
                       int   num_elem_blks_g,
                       int   num_node_sets_g,
                       int   num_side_sets_g
                       )
{
  char   *func_name="ne_put_init_global";

  int     varid, dimid, status;
  int     ltempsv;

  char    errmsg[MAX_ERR_LENGTH];
/*-----------------------------Execution begins-----------------------------*/

  exerrval = 0; /* clear error code */

  /* Put NetCDF file into define mode */
  if ((status = nc_redef(neid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to put file ID %d into define mode", neid);
    ex_err(func_name, errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Output the file version */
  if ((status=ne_put_version(neid)) < 0) return (status);

  /* Define dimension for number of global nodes */
  ltempsv = num_nodes_g;
  if ((status = nc_def_dim(neid, DIM_NUM_NODES_GLOBAL, ltempsv, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to dimension \"%s\" in file ID %d",
            DIM_NUM_NODES_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

    return (EX_FATAL);
  }

  /* Define dimension for number of global elements */
  ltempsv = num_elems_g;
  if ((status = nc_def_dim(neid, DIM_NUM_ELEMS_GLOBAL, ltempsv, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to dimension \"%s\" in file ID %d",
            DIM_NUM_ELEMS_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

    return (EX_FATAL);
  }

  /*
   * Output the number of global element blocks. This is output as a
   * dimension since the vector of global element block IDs is sized
   * by this quantity.
   */
  ltempsv = num_elem_blks_g;
  if ((status = nc_def_dim(neid, DIM_NUM_ELBLK_GLOBAL, ltempsv, &dimid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to dimension \"%s\" in file ID %d",
            DIM_NUM_ELBLK_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

    return (EX_FATAL);
  }

  /* Define the element block IDs variable. */
  if ((status = nc_def_var(neid, VAR_ELBLK_IDS_GLOBAL, NC_INT, 1, &dimid, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to put variable definition for \"%s\" into \
file ID %d",
            VAR_ELBLK_IDS_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

    return (EX_FATAL);
  }

  /* Define the element block counts variable. */
  if ((status = nc_def_var(neid, VAR_ELBLK_CNT_GLOBAL, NC_INT, 1, &dimid, &varid)) != NC_NOERR) {
    exerrval = status;
    sprintf(errmsg,
            "Error: failed to put variable definition for \"%s\" into \
file ID %d",
            VAR_ELBLK_CNT_GLOBAL, neid);
    ex_err(func_name, errmsg, exerrval);
    /* Leave define mode before returning */
    ne_leavedef(neid, func_name);

    return (EX_FATAL);
  }

  /*
   * Output the number of global node sets. This is output as a
   * dimension since the vector of global element block IDs is sized
   * by this quantity.
   */
  if (num_node_sets_g > 0) {
    ltempsv = num_node_sets_g;
    if ((status = nc_def_dim(neid, DIM_NUM_NS_GLOBAL, ltempsv, &dimid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to dimension \"%s\" in file ID %d",
              DIM_NUM_NS_GLOBAL, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* Define the variable for output of global node set IDs */
    if ((status = nc_def_var(neid, VAR_NS_IDS_GLOBAL, NC_INT, 1, &dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to put variable definition for \"%s\" into \
file ID %d",
              VAR_NS_IDS_GLOBAL, neid);
      ex_err(func_name, errmsg, exerrval);

      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* Define variable for global node counts in each global node set */
    if ((status = nc_def_var(neid, VAR_NS_NODE_CNT_GLOBAL, NC_INT, 1, &dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to put variable definition for \"%s\" into \
file ID %d",
              VAR_NS_NODE_CNT_GLOBAL, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returing */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /*
     * Define variable for global dist. factor count in each global
     * node set
     */
    if ((status = nc_def_var(neid, VAR_NS_DF_CNT_GLOBAL, NC_INT, 1, &dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to put variable definition for \"%s\" into \
file ID %d",
              VAR_NS_DF_CNT_GLOBAL, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returing */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

  } /* End "if (num_node_sets_g > 0)" */

  /*
   * Output the number of global side sets. This is output as a
   * dimension since the vector of global element block IDs is sized
   * by this quantity.
   */
  if (num_side_sets_g > 0) {
    ltempsv = num_side_sets_g;
    if ((status = nc_def_dim(neid, DIM_NUM_SS_GLOBAL, ltempsv, &dimid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to dimension \"%s\" in file id %d",
              DIM_NUM_SS_GLOBAL, neid);
      ex_err(func_name, errmsg, exerrval);
      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /* Define the variable for output of global side set IDs */
    if ((status = nc_def_var(neid, VAR_SS_IDS_GLOBAL, NC_INT, 1, &dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to put variable definition for \"%s\" into \
file id %d",
              VAR_SS_IDS_GLOBAL, neid);
      ex_err(func_name, errmsg, exerrval);

      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /*
     * Define the variable for count of global number of sides in the
     * global side sets.
     */
    if ((status = nc_def_var(neid, VAR_SS_SIDE_CNT_GLOBAL, NC_INT, 1, &dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to put variable definition for \"%s\" into \
file id %d",
              VAR_SS_SIDE_CNT_GLOBAL, neid);
      ex_err(func_name, errmsg, exerrval);

      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

    /*
     * Define the variable for count of global dist. factors in the
     * global side sets.
     */
    if ((status = nc_def_var(neid, VAR_SS_DF_CNT_GLOBAL, NC_INT, 1, &dimid, &varid)) != NC_NOERR) {
      exerrval = status;
      sprintf(errmsg,
              "Error: failed to put variable definition for \"%s\" into \
file id %d",
              VAR_SS_DF_CNT_GLOBAL, neid);
      ex_err(func_name, errmsg, exerrval);

      /* Leave define mode before returning */
      ne_leavedef(neid, func_name);

      return (EX_FATAL);
    }

  } /* End "if (num_side_sets_g > 0)" */

  /* End define mode */
  if (ne_leavedef(neid, func_name) != EX_NOERR)
    return (EX_FATAL);

  return (EX_NOERR);
}
