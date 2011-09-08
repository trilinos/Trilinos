/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "exodusII.h"
#include "ne_nemesisI.h"

#include "rf_salsa.h"
#include "rf_comm.h"
#include "rf_message.h"

#include "pe_common.h"
#include "el_geom_const.h"
#include "ps_pario_const.h"
#include "rf_io_const.h"
#include "rf_mp_const.h"
#include "rf_allo.h"
#include "rf_util.h"

/************ R O U T I N E S   I N   T H I S   F I L E ***********************
*
*  Name_of_Routine              type                 Called by
*  ---------------------------------------------------------------
*
*  read_restart_params ()                         main:pe_exoII_2par_exoII.c
*  read_restart_data ()                           main:pe_exoII_2par_exoII.c
*  read_var_param ()             int              read_restart_params
*  broadcast_var_param ()        int              read_restart_params
*  read_vars()                   int              read_restart_data
*  read_elem_vars ()             int              read_vars
*  read_nodal_vars ()            int              read_vars
*  compare_mesh_param ()         int              read_restart_params
*  find_gnode_inter ()           int              read_nodal_vars
*
******************************************************************************/

/******** P R O T O T Y P E S   O F    F U N C T I O N S *********************/
static int get_free_descriptor_count(void);
extern int  check_monot       (int vector[], int length);
extern void check_exodus_error(
                               int  error,           /* error code         */
                               char function_name[]
                                    /* EXODUS function returning the error   */
                               );

static int read_var_param (int exoid, int max_name_length);
static int broadcast_var_param (RESTART_PTR restart, int max_name_length);
static int read_vars(int exoid, int index, int blk_cnt, int *eb_ids,
                     int *eb_cnts, int ***eb_map_ptr, int **eb_cnts_local,
		     int *ss_ids, int *ss_cnts, int *ns_ids, int *ns_cnts,
                     int io_ws);
static int read_elem_vars(int exoid, int index, int blk_cnt, int *eb_ids,
                          int *eb_cnts, int ***eb_map_ptr,
                          int **eb_cnts_local, int io_ws);
static int read_elem_vars_1(int exoid, int index, int blk_cnt, int *eb_ids,
			    int *eb_cnts, int ***eb_map_ptr,
			    int **eb_cnts_local, int io_ws, int iblk,
			    int eb_offset, int *local_offset);
static int read_elem_vars_n(int exoid, int index, int blk_cnt, int *eb_ids,
			    int *eb_cnts, int ***eb_map_ptr,
			    int **eb_cnts_local, int io_ws, int iblk,
			    int eb_offset, int *local_offset,
			    int max_elem_per_mesg, int num_mesgs, int *mesg_start);

static int read_sset_vars(int exoid, int index, int blk_cnt, int *ss_ids,
                          int *ss_cnts, int io_ws);
static int read_sset_vars_1(int exoid, int index, int blk_cnt, int *ss_ids,
			    int *ss_cnts, int io_ws, int iblk);
static int read_sset_vars_n(int exoid, int index, int blk_cnt, int *ss_ids,
			    int *ss_cnts, int io_ws, int iblk,
			    int max_sset_per_mesg, int num_mesgs, int *mesg_start);

static int read_nset_vars(int exoid, int index, int blk_cnt, int *ns_ids,
                          int *ns_cnts, int io_ws);
static int read_nset_vars_1(int exoid, int index, int blk_cnt, int *ns_ids,
			    int *ns_cnts, int io_ws, int iblk);
static int read_nset_vars_n(int exoid, int index, int blk_cnt, int *ns_ids,
			    int *ns_cnts, int io_ws, int iblk,
			    int max_elem_per_mesg, int num_mesgs,
			    int *mesg_start);

static int read_nodal_vars (int exoid, int index, int blk_cnt, int io_ws);
static int compare_mesh_param (int exoid);
static int find_gnode_inter(int *intersect, int num_glob, int *glob_vec,
			    int num_int, int num_bor, int num_ext,
			    int *loc_vec);


/*****************************************************************************/
/*****************************************************************************/

void read_restart_params(int io_ws)

/* Function which reads the restart variable parameters for the EXODUS II
 * database which contains the results information. Allocate necessary
 * memory on each processor.
 *
 *----------------------------------------------------------------------------
 *
 * Functions called:
 *
 * compare_mesh_param -- function which checks that parameters in
 *                       the restart EXODUS II file are the same as in
 *                       the mesh EXODUS II file
 * read_var_param -- function which reads the time indicies, number
 *                   of variables, and their names from the restart file
 * broadcast_var_param -- function which broadcasts parameters to all of
 *                        the processors
 *
 *----------------------------------------------------------------------------
 */

{
  char  *yo="read_restart_params";

  int    exoid, cpu_ws=0;
  float  vers;
  int    max_name_length = 0;
  
  if (Proc == 0) {
    /* Open the ExodusII file */
    if ((exoid=ex_open(Exo_Res_File, EX_READ, &cpu_ws, &io_ws, &vers)) < 0) {
      fprintf(stderr, "%s: Could not open file %s for restart info\n",
              yo, Exo_Res_File);
      exit(1);
    }

    max_name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
    ex_set_max_name_length(exoid, max_name_length);
    
    /*
     * Just do a rudimentary check to figure out if the mesh parameters
     * in the results file are the same as the mesh parameters in the
     * mesh file.
     */
    if (strcmp(ExoFile, Exo_Res_File) != 0)
      if (!compare_mesh_param(exoid)) {
        fprintf(stderr, "%s: Mesh parameters in mesh and result files"
                " differ\n", yo);
        exit(1);
      }

    /* get the time, and the variable names */
    if (read_var_param(exoid, max_name_length) < 0) {
      fprintf(stderr, "%s: Error occured while reading variable parameters\n",
              yo);
      exit(1);
    }

    /* Close the ExodusII file */
    ex_close(exoid);

  }

  /* now broadcast the variable parameters to all of the processors */
  if (broadcast_var_param(&Restart_Info, max_name_length) < 0) {
    fprintf(stderr, "%s: Error occured while broadcasting variable params\n",
            yo);
    exit(1);
  }

  return;
}

void read_restart_data (int io_ws)

/* Function which reads the restart variable data from the EXODUS II
 * database which contains the results information. Then distribute
 * it to the processors, and write it to the parallel exodus files.
 *
 *----------------------------------------------------------------------------
 *
 * Functions called:
 *
 * read_vars -- function which reads the variable values from the restart
 *              file, and then distributes them to the processors
 *
 * write_var_timestep -- function which writes out the variables for a
 *                       to a parallel ExodusII file.
 *
 *----------------------------------------------------------------------------
 */

{
  char  *yo="read_restart_data";

  int    cnt, ilocal, ifound, offset, iproc, cpu_ws;
  int   *eb_ids_global = NULL, *eb_cnts_global = NULL;
  int   *ss_ids_global = NULL, *ss_cnts_global = NULL;
  int   *ns_ids_global = NULL, *ns_cnts_global = NULL;
  int ***eb_map_ptr = NULL, **eb_cnts_local = NULL;
  int    num_blocks, times_in_blk, iblk, time_idx;
  int    array_size;
  int    dum1, dum2;
  int    exoid=0, *par_exoid = NULL;
  
  int    open_file_count;
  double start_t, end_t;
  float  vers;

  void  *ptr;

  char   cTemp[512];

  /* computing precision should be the same as the database precision
   *
   * EXCEPTION: if the io_ws is smaller than the machine precision,
   * ie - database with io_ws == 4 on a Cray (sizeof(float) == 8),
   * then the cpu_ws must be the machine precision.
   */
  if (io_ws < sizeof(float)) cpu_ws = sizeof(float);
  else                       cpu_ws = io_ws;

  if (Proc == 0) {
    /* Open the ExodusII file */
    if ((exoid=ex_open(Exo_Res_File, EX_READ, &cpu_ws, &io_ws, &vers)) < 0) {
      fprintf(stderr, "%s: Could not open file %s for restart info\n",
              yo, Exo_Res_File);
      exit(1);
    }
  }

  /* allocate memory for the time values */
  ptr = array_alloc (__FILE__, __LINE__, 1, Restart_Info.Block_Size, cpu_ws);
  if (cpu_ws == sizeof(float))
    Restart_Info.Time_sp = (float *) ptr;
  else
    Restart_Info.Time_dp = (double *) ptr;

  /* allocate space for the global variables */
  if (Restart_Info.NVar_Glob > 0) {
    ptr = array_alloc (__FILE__, __LINE__, 2, Restart_Info.Block_Size,
		       Restart_Info.NVar_Glob, cpu_ws);

    if (cpu_ws == sizeof(float))
      Restart_Info.Glob_Vals_sp = (float **) ptr;
    else
      Restart_Info.Glob_Vals_dp = (double **) ptr;
  }

  if (Restart_Info.NVar_Elem > 0 ) {

    /* allocate storage space */
    if (cpu_ws == sizeof(float))
      Restart_Info.Elem_Vals_sp = (float ***) array_alloc (__FILE__, __LINE__,
                                                           2, Proc_Info[2],
							   Restart_Info.Block_Size,
                                                           sizeof(float *));
    else
      Restart_Info.Elem_Vals_dp = (double ***) array_alloc (__FILE__, __LINE__,
                                                            2, Proc_Info[2],
							    Restart_Info.Block_Size,
                                                            sizeof(double *));

    /* now allocate storage for the values */
    for (iproc = 0; iproc < Proc_Info[2]; iproc++) {
      array_size = Restart_Info.NVar_Elem *
	(Num_Internal_Elems[iproc] + Num_Border_Elems[iproc]);
      ptr = array_alloc (__FILE__, __LINE__, 1,
                         (Restart_Info.Block_Size * array_size) , cpu_ws);

      if (cpu_ws == sizeof(float))
        Restart_Info.Elem_Vals_sp[iproc][0] = (float *) ptr;
      else
        Restart_Info.Elem_Vals_dp[iproc][0] = (double *) ptr;

      for (cnt = 1; cnt < Restart_Info.Block_Size; cnt++) {
        if (cpu_ws == sizeof(float))
          Restart_Info.Elem_Vals_sp[iproc][cnt] =
            Restart_Info.Elem_Vals_sp[iproc][cnt-1] + array_size;
        else
          Restart_Info.Elem_Vals_dp[iproc][cnt] =
            Restart_Info.Elem_Vals_dp[iproc][cnt-1] + array_size;
      }
    }

    /*
     * at this point, I need to broadcast the global element block ids
     * and counts to the processors. I know that this is redundant data
     * since they will all receive this information in read_mesh, but
     * the variables which contain that information are static in
     * el_exoII_io.c, and cannot be used here. So, take a second and
     * broadcast all of this out.
     *
     * I want to do this here so that it is done only once no matter
     * how many time steps are retrieved
     */

    /* need to get the element block ids and counts */
    eb_ids_global  = (int *) array_alloc (__FILE__, __LINE__, 1,
                                          (2 * Num_Elem_Blk), sizeof(int));
    if (!eb_ids_global) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n", Proc, yo);
      exit(1);
    }
    eb_cnts_global = eb_ids_global + Num_Elem_Blk;

    if (Proc == 0) {
      /* Get the Element Block IDs from the input file */
      if (ex_get_elem_blk_ids (exoid, eb_ids_global) < 0)
	{
	  fprintf(stderr, "%s: unable to get element block IDs", yo);
	  exit(1);
	}

      /* Get the count of elements in each element block */
      for (cnt = 0; cnt < Num_Elem_Blk; cnt++) {
	if (ex_get_elem_block(exoid, eb_ids_global[cnt], cTemp,
			      &(eb_cnts_global[cnt]), &dum1, &dum2) < 0) {
	  fprintf(stderr, "%s: unable to get element count for block id %d",
		  yo, eb_ids_global[cnt]);
	  exit(1);
	}
      }
    }

    brdcst(Proc, Num_Proc, (char *) eb_ids_global, (2*Num_Elem_Blk*sizeof(int)), 0);

    /*
     * in order to speed up finding matches in the global element
     * number map, set up an array of pointers to the start of
     * each element block's global element number map. That way
     * only entries for the current element block have to be searched
     */
    eb_map_ptr = (int ***) array_alloc (__FILE__, __LINE__, 2, Proc_Info[2],
                                        Num_Elem_Blk, sizeof(int *));
    if (!eb_map_ptr) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n", Proc, yo);
      exit(1);
    }
    eb_cnts_local = (int **) array_alloc (__FILE__, __LINE__, 2, Proc_Info[2],
                                          Num_Elem_Blk, sizeof(int));
    if (!eb_cnts_local) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n", Proc, yo);
      exit(1);
    }

    /*
     * for now, assume that element blocks have been
     * stored in the same order as the global blocks
     */
    for (iproc = 0; iproc < Proc_Info[2]; iproc++) {
      ifound = 0;
      offset = 0;
      for (cnt = 0; cnt < Num_Elem_Blk; cnt++) {
        for (ilocal = ifound; ilocal < Proc_Num_Elem_Blk[iproc]; ilocal++) {
          if (Proc_Elem_Blk_Ids[iproc][ilocal] == eb_ids_global[cnt])
            break;
        }

        if (ilocal < Proc_Num_Elem_Blk[iproc]) {
          eb_map_ptr[iproc][cnt] = &GElems[iproc][offset];
          eb_cnts_local[iproc][cnt] = Proc_Num_Elem_In_Blk[iproc][ilocal];
          offset += Proc_Num_Elem_In_Blk[iproc][ilocal];
          ifound = ilocal; /* don't search the same part of the list over */
        }
        else {
          eb_map_ptr[iproc][cnt] = NULL;
          eb_cnts_local[iproc][cnt] = 0;
        }
      }
    }

  } /* End: "if (Restart_Info.NVar_Elem > 0 )" */

  if (Restart_Info.NVar_Node > 0 ) {
    /* allocate storage space */
    if (cpu_ws == sizeof(float))
      Restart_Info.Node_Vals_sp = (float ***) array_alloc (__FILE__, __LINE__,
                                                           2, Proc_Info[2],
							   Restart_Info.Block_Size,
                                                           sizeof(float *));
    else
      Restart_Info.Node_Vals_dp = (double ***) array_alloc (__FILE__, __LINE__,
                                                            2, Proc_Info[2],
							    Restart_Info.Block_Size,
                                                            sizeof(double *));

    /* now allocate storage for the values */
    for (iproc = 0; iproc < Proc_Info[2]; iproc++) {
      array_size = Restart_Info.NVar_Node * (Num_Internal_Nodes[iproc] +
					     Num_Border_Nodes[iproc] + Num_External_Nodes[iproc]);
      ptr = array_alloc (__FILE__, __LINE__, 1,
                         (Restart_Info.Block_Size * array_size), cpu_ws);

      if (cpu_ws == sizeof(float))
        Restart_Info.Node_Vals_sp[iproc][0] = (float *) ptr;
      else
        Restart_Info.Node_Vals_dp[iproc][0] = (double *) ptr;

      for (cnt = 1; cnt < Restart_Info.Block_Size; cnt++) {
        if (cpu_ws == sizeof(float))
          Restart_Info.Node_Vals_sp[iproc][cnt] =
            Restart_Info.Node_Vals_sp[iproc][cnt-1] + array_size;
        else
          Restart_Info.Node_Vals_dp[iproc][cnt] =
            Restart_Info.Node_Vals_dp[iproc][cnt-1] + array_size;
      }
    }
  }

  if (Restart_Info.NVar_Sset > 0 ) {

    /* allocate storage space */
    if (cpu_ws == sizeof(float))
      Restart_Info.Sset_Vals_sp = (float ***) array_alloc (__FILE__, __LINE__,
                                                           2, Proc_Info[2],
							   Restart_Info.Block_Size,
                                                           sizeof(float *));
    else
      Restart_Info.Sset_Vals_dp = (double ***) array_alloc (__FILE__, __LINE__,
                                                            2, Proc_Info[2],
							    Restart_Info.Block_Size,
                                                            sizeof(double *));

    /* now allocate storage for the values */
    for (iproc = 0; iproc < Proc_Info[2]; iproc++) {
      array_size = Restart_Info.NVar_Sset * Proc_SS_Elem_List_Length[iproc];
                   
      ptr = array_alloc (__FILE__, __LINE__, 1,
                         (Restart_Info.Block_Size * array_size) , cpu_ws);

      if (cpu_ws == sizeof(float))
        Restart_Info.Sset_Vals_sp[iproc][0] = (float *) ptr;
      else
        Restart_Info.Sset_Vals_dp[iproc][0] = (double *) ptr;

      for (cnt = 1; cnt < Restart_Info.Block_Size; cnt++) {
        if (cpu_ws == sizeof(float))
          Restart_Info.Sset_Vals_sp[iproc][cnt] =
            Restart_Info.Sset_Vals_sp[iproc][cnt-1] + array_size;
        else
          Restart_Info.Sset_Vals_dp[iproc][cnt] =
            Restart_Info.Sset_Vals_dp[iproc][cnt-1] + array_size;
      }
    }

    /*
     * at this point, I need to broadcast the ids and counts to the
     * processors. I know that this is redundant data since they will
     * all receive this information in read_mesh, but the variables
     * which contain that information are static in el_exoII_io.c, and
     * cannot be used here. So, take a second and broadcast all of
     * this out.
     *
     * I want to do this here so that it is done only once no matter
     * how many time steps are retrieved
     */

    /* need to get the sideset ids and counts */
    ss_ids_global  = (int *) array_alloc (__FILE__, __LINE__, 1,
                                          (2 * Num_Side_Set), sizeof(int));
    if (!ss_ids_global) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n", Proc, yo);
      exit(1);
    }
    ss_cnts_global = ss_ids_global + Num_Side_Set;

    if (Proc == 0) {
      /* Get the Sideset IDs from the input file */
      if (ex_get_ids (exoid, EX_SIDE_SET, ss_ids_global) < 0) {
	fprintf(stderr, "%s: unable to get sideset IDs", yo);
	exit(1);
      }

      /* Get the count of elements in each sideset */
      for (cnt = 0; cnt < Num_Side_Set; cnt++) {
	if (ex_get_side_set_param(exoid, ss_ids_global[cnt],
				  &(ss_cnts_global[cnt]), &dum1) < 0) {
	  fprintf(stderr, "%s: unable to get element count for sideset id %d",
		  yo, ss_ids_global[cnt]);
	  exit(1);
	}
      }
    }

    brdcst(Proc, Num_Proc, (char *) ss_ids_global, (2*Num_Side_Set*sizeof(int)), 0);

  } /* End: "if (Restart_Info.NVar_Sset > 0 )" */


  if (Restart_Info.NVar_Nset > 0 ) {

    /* allocate storage space */
    if (cpu_ws == sizeof(float))
      Restart_Info.Nset_Vals_sp = (float ***) array_alloc (__FILE__, __LINE__,
                                                           2, Proc_Info[2],
							   Restart_Info.Block_Size,
                                                           sizeof(float *));
    else
      Restart_Info.Nset_Vals_dp = (double ***) array_alloc (__FILE__, __LINE__,
                                                            2, Proc_Info[2],
							    Restart_Info.Block_Size,
                                                            sizeof(double *));

    /* now allocate storage for the values */
    for (iproc = 0; iproc < Proc_Info[2]; iproc++) {
      array_size = Restart_Info.NVar_Nset * Proc_NS_List_Length[iproc];
                   
      ptr = array_alloc (__FILE__, __LINE__, 1,
                         (Restart_Info.Block_Size * array_size) , cpu_ws);

      if (cpu_ws == sizeof(float))
        Restart_Info.Nset_Vals_sp[iproc][0] = (float *) ptr;
      else
        Restart_Info.Nset_Vals_dp[iproc][0] = (double *) ptr;

      for (cnt = 1; cnt < Restart_Info.Block_Size; cnt++) {
        if (cpu_ws == sizeof(float))
          Restart_Info.Nset_Vals_sp[iproc][cnt] =
            Restart_Info.Nset_Vals_sp[iproc][cnt-1] + array_size;
        else
          Restart_Info.Nset_Vals_dp[iproc][cnt] =
            Restart_Info.Nset_Vals_dp[iproc][cnt-1] + array_size;
      }
    }

    /*
     * at this point, I need to broadcast the ids and counts to the
     * processors. I know that this is redundant data since they will
     * all receive this information in read_mesh, but the variables
     * which contain that information are static in el_exoII_io.c, and
     * cannot be used here. So, take a second and broadcast all of
     * this out.
     *
     * I want to do this here so that it is done only once no matter
     * how many time steps are retrieved
     */

    /* need to get the nodeset ids and counts */
    ns_ids_global  = (int *) array_alloc (__FILE__, __LINE__, 1,
                                          (2 * Num_Node_Set), sizeof(int));
    if (!ns_ids_global) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n", Proc, yo);
      exit(1);
    }
    ns_cnts_global = ns_ids_global + Num_Node_Set;

    if (Proc == 0) {
      /* Get the Nodeset IDs from the input file */
      if (ex_get_ids (exoid, EX_NODE_SET, ns_ids_global) < 0) {
	fprintf(stderr, "%s: unable to get nodeset IDs", yo);
	exit(1);
      }

      /* Get the count of elements in each nodeset */
      for (cnt = 0; cnt < Num_Node_Set; cnt++) {
	if (ex_get_node_set_param(exoid, ns_ids_global[cnt],
				  &(ns_cnts_global[cnt]), &dum1) < 0) {
	  fprintf(stderr, "%s: unable to get element count for nodeset id %d",
		  yo, ns_ids_global[cnt]);
	  exit(1);
	}
      }
    }

    brdcst(Proc, Num_Proc, (char *) ns_ids_global, (2*Num_Node_Set*sizeof(int)), 0);

  } /* End: "if (Restart_Info.NVar_Nset > 0 )" */


  /*
   * NOTE: A possible place to speed this up would be to
   * get the global node and element lists here, and broadcast
   * them out only once.
   */

  /* figure out how many blocks need to be looped over */
  num_blocks = Restart_Info.Num_Times / Restart_Info.Block_Size;
  if ((Restart_Info.Num_Times % Restart_Info.Block_Size) != 0)
    num_blocks++;
  
  par_exoid = malloc(Proc_Info[2] * sizeof(int));
  if(!par_exoid) {
    fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n",
	    Proc, yo);
    exit(1);
  }

  /* See if any '/' in the name.  IF present, isolate the basename of the file */
  if (strrchr(PIO_Info.Scalar_LB_File_Name, '/') != NULL) {
    /* There is a path separator.  Get the portion after the
     * separator
     */
    strcpy(cTemp, strrchr(PIO_Info.Scalar_LB_File_Name, '/')+1);
  } else {
    /* No separator; this is already just the basename... */
    strcpy(cTemp, PIO_Info.Scalar_LB_File_Name);
  }    
  
  if (strlen(PIO_Info.Exo_Extension) == 0)
    add_fname_ext(cTemp, ".par");
  else
    add_fname_ext(cTemp, PIO_Info.Exo_Extension);
  
  open_file_count = get_free_descriptor_count();
  if (open_file_count > Proc_Info[2]) {
    fprintf(stderr, "All output files opened simultaneously.\n");
    for (iproc=0; iproc < Proc_Info[2]; iproc++) {
      
      gen_par_filename(cTemp, Par_Nem_File_Name, Proc_Ids[iproc],
		       Proc_Info[0]);
      
      /* Open the parallel Exodus II file for writing */
      if ((par_exoid[iproc]=ex_open(Par_Nem_File_Name, EX_WRITE, &cpu_ws,
				    &io_ws, &vers)) < 0) {
	fprintf(stderr,"[%d] %s Could not open parallel Exodus II file: %s\n",
		iproc, yo, Par_Nem_File_Name);
	exit(1);
      }
    }
  } else {
    fprintf(stderr, "All output files opened one-at-a-time.\n");
  }

  /* Now loop over the number of time steps */
  for (iblk = 0; iblk < num_blocks; iblk++) {

    /* now figure out how many times are in this block */
    if (((iblk + 1) * Restart_Info.Block_Size) > Restart_Info.Num_Times)
      times_in_blk = Restart_Info.Num_Times % Restart_Info.Block_Size;
    else
      times_in_blk = Restart_Info.Block_Size;

    for (cnt = 0; cnt < times_in_blk; cnt++) {
      time_idx = iblk * Restart_Info.Block_Size + cnt;
      start_t = second ();

      /* read and distribute the variables for this time step */
      if (read_vars(exoid, Restart_Info.Time_Idx[iblk], 0,
		    eb_ids_global, eb_cnts_global, eb_map_ptr,
		    eb_cnts_local,
		    ss_ids_global, ss_cnts_global,
		    ns_ids_global, ns_cnts_global, io_ws) < 0) {
	fprintf(stderr, "%s: Error occured while reading variables\n",
		yo);
	exit(1);
      }
      end_t   = second () - start_t;
      if (Proc == 0) printf ("\tTime to read  vars for timestep %d: %f (sec.)\n", (time_idx+1), end_t);

      start_t = second ();
      for (iproc=0; iproc < Proc_Info[2]; iproc++) {

	if (open_file_count < Proc_Info[2]) {
	  gen_par_filename(cTemp, Par_Nem_File_Name, Proc_Ids[iproc],
			   Proc_Info[0]);
	  
	  /* Open the parallel Exodus II file for writing */
	  if ((par_exoid[iproc]=ex_open(Par_Nem_File_Name, EX_WRITE, &cpu_ws,
					&io_ws, &vers)) < 0) {
	    fprintf(stderr,"[%d] %s Could not open parallel Exodus II file: %s\n",
		    iproc, yo, Par_Nem_File_Name);
	    exit(1);
	  }
	}

	/*
	 * Write out the variable data for the time steps in this
	 * block to each parallel file.
	 */
	for (cnt = 0; cnt < times_in_blk; cnt++) {
	  time_idx = iblk * Restart_Info.Block_Size + cnt;
	  write_var_timestep(par_exoid[iproc], iproc, (time_idx+1), cnt,
			     eb_ids_global, ss_ids_global, ns_ids_global,
			     io_ws);
	}

	if (Proc == 0) {
	  if (iproc%10 == 0 || iproc == Proc_Info[2]-1)
	    fprintf(stderr, "%d", iproc);
	  else
	    fprintf(stderr, ".");
	}

	if (open_file_count < Proc_Info[2]) {
	  if (ex_close(par_exoid[iproc]) == -1) {
	    fprintf(stderr, "[%d] %s Could not close the parallel Exodus II file.\n",
		    iproc, yo);
	    exit(1);
	  }
	}
      } /* End "for (iproc=0; iproc < Proc_Info[2]; iproc++)" */

      end_t   = second () - start_t;
      if (Proc == 0) printf ("\n\tTime to write vars for timestep %d: %f (sec.)\n", (time_idx+1), end_t);

    } /* End: "for (iblk = 0; iblk < num_blocks; iblk++)" */
  }
  if (Restart_Info.NVar_Elem > 0 ) {
    safe_free((void **) &eb_ids_global);
    safe_free((void **) &eb_map_ptr);
    safe_free((void **) &eb_cnts_local);
  }

  if (Restart_Info.NVar_Sset > 0 ) {
    safe_free((void **) &ss_ids_global);
  }
  if (Restart_Info.NVar_Nset > 0 ) {
    safe_free((void **) &ns_ids_global);
  }

  /* Close the restart exodus II file */
  if (Proc == 0) {
    if (ex_close(exoid) == -1) {
      fprintf(stderr, "%sCould not close the restart Exodus II file\n",
              yo);
      exit(1);
    }
  }

  for (iproc=0; iproc < Proc_Info[2]; iproc++) {
    /* Close the parallel exodus II file */
    if (ex_close(par_exoid[iproc]) == -1) {
      fprintf(stderr, "[%d] %s Could not close the parallel Exodus II file.\n",
	      iproc, yo);
      exit(1);
    }
  }
  if (par_exoid != NULL) {
    free(par_exoid);
    par_exoid = NULL;
  }
}

static int read_var_param (int exoid, int max_name_length)
{
  char  *yo="read_var_param";

  int    ret_int, cnt;

  /* Get the number of time indices contained in the file */
  ret_int = ex_inquire_int(exoid, EX_INQ_TIME);

  /* see if the user want to get all of the time indices */
  if (Restart_Info.Num_Times == -1) {

    Restart_Info.Num_Times = ret_int;

    if (ret_int > 0) {
      /* allocate array space */
      Restart_Info.Time_Idx = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                  ret_int, sizeof(int));

      for (cnt = 0; cnt < Restart_Info.Num_Times; cnt++)
        Restart_Info.Time_Idx[cnt] = cnt + 1;
    }
  }
  else {
    /* Check to see if the requested indeces are valid */
    for (cnt = 0; cnt < Restart_Info.Num_Times; cnt++) {

      /* if the user wants the last time, then set it */
      if (Restart_Info.Time_Idx[cnt] == 0)
        Restart_Info.Time_Idx[cnt] = ret_int;

      if (Restart_Info.Time_Idx[cnt] > ret_int) {
        fprintf(stderr, "%s: Requested time index, %d, out of range.\n",
                yo, Restart_Info.Time_Idx[cnt]);
        fprintf(stderr, "%s: Valid time indices in %s are from 1 to %d.\n",
                yo, Exo_Res_File, ret_int);
        return -1;
      }

    }
  }

  /* if there are not any time steps, then return here without an error */
  if (Restart_Info.Num_Times == 0) {
    Restart_Info.Flag = 0;
    Restart_Info.NVar_Glob = 0;
    Restart_Info.NVar_Node = 0;
    Restart_Info.NVar_Elem = 0;
    return 0;
  }


  /***************** Global Variables ********************/
  if (ex_get_var_param(exoid, "g", &(Restart_Info.NVar_Glob)) < 0) {
    fprintf(stderr, "%s: Could not get global variable parameter from file\n",
            yo);
    return -1;
  }

  /* allocate space for the global variable names */
  if (Restart_Info.NVar_Glob > 0) {
    Restart_Info.GV_Name = (char **) array_alloc(__FILE__, __LINE__, 2,
                                                 Restart_Info.NVar_Glob,
                                                 max_name_length+1,
                                                 sizeof(char));

    /* get the global variable names */
    if (ex_get_var_names(exoid, "g", Restart_Info.NVar_Glob,
                         Restart_Info.GV_Name) < 0) {
      fprintf(stderr, "%s: Could not get global variable names from file\n",
              yo);
      return -1;
    }
  }

  /***************** Elemental Variables ********************/
  if (ex_get_var_param(exoid, "e", &(Restart_Info.NVar_Elem)) < 0) {
    fprintf(stderr, "%s: Could not get elemental variable param from file\n",
            yo);
    return -1;
  }

  /* allocate space for the elemental variable names */
  if (Restart_Info.NVar_Elem > 0) {
    Restart_Info.EV_Name = (char **) array_alloc(__FILE__, __LINE__, 2,
                                                 Restart_Info.NVar_Elem,
                                                 max_name_length+1,
                                                 sizeof(char));

    /* get the elemental variable names */
    if (ex_get_var_names(exoid, "e", Restart_Info.NVar_Elem,
                         Restart_Info.EV_Name) < 0) {
      fprintf(stderr, "%s: Could not get elemental variable names from file\n",
              yo);
      return -1;
    }

    /* and get the truth table */
    Restart_Info.GElem_TT = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                (Num_Elem_Blk
                                                * Restart_Info.NVar_Elem),
                                                sizeof(int));

    check_exodus_error(ex_get_var_tab(exoid, "e", 
				      Num_Elem_Blk,
				      Restart_Info.NVar_Elem,
				      Restart_Info.GElem_TT),
                       "ex_get_elem_var_tab");
  }

  /******************* Nodal Variables **********************/
  if (ex_get_var_param(exoid, "n", &(Restart_Info.NVar_Node)) < 0) {
    fprintf(stderr, "%s: Could not get nodal variable param from file\n",
            yo);
    return -1;
  }

  /* allocate space for the nodal variable names */
  if (Restart_Info.NVar_Node > 0) {
    Restart_Info.NV_Name = (char **) array_alloc(__FILE__, __LINE__, 2,
                                                 Restart_Info.NVar_Node,
                                                 max_name_length+1,
                                                 sizeof(char));

    /* get the nodal variable names */
    if (ex_get_var_names(exoid, "n", Restart_Info.NVar_Node,
                         Restart_Info.NV_Name) < 0) {
      fprintf(stderr, "%s: Could not get nodal variable names from file\n",
              yo);
      return -1;
    }
  }

  /******************* Sideset Variables **********************/
  if (ex_get_var_param(exoid, "s", &(Restart_Info.NVar_Sset)) < 0) {
    fprintf(stderr, "%s: Could not get sideset variable param from file\n",
            yo);
    return -1;
  }

  /* allocate space for the variable names */
  if (Restart_Info.NVar_Sset > 0) {
    Restart_Info.SSV_Name = (char **) array_alloc(__FILE__, __LINE__, 2,
						  Restart_Info.NVar_Sset,
						  max_name_length+1,
						  sizeof(char));

    /* get the variable names */
    if (ex_get_var_names(exoid, "s", Restart_Info.NVar_Sset,
                         Restart_Info.SSV_Name) < 0) {
      fprintf(stderr, "%s: Could not get sideset variable names from file\n",
              yo);
      return -1;
    }

    /* and get the truth table */
    Restart_Info.GSset_TT = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                (Num_Side_Set
                                                * Restart_Info.NVar_Sset),
                                                sizeof(int));

    check_exodus_error(ex_get_var_tab(exoid, "s", Num_Side_Set,
				      Restart_Info.NVar_Sset,
				      Restart_Info.GSset_TT),
                       "ex_get_var_tab");
  }

  /******************* Nodeset Variables **********************/
  if (ex_get_var_param(exoid, "m", &(Restart_Info.NVar_Nset)) < 0) {
    fprintf(stderr, "%s: Could not get nodeset variable param from file\n",
            yo);
    return -1;
  }

  /* allocate space for the variable names */
  if (Restart_Info.NVar_Nset > 0) {
    Restart_Info.NSV_Name = (char **) array_alloc(__FILE__, __LINE__, 2,
						  Restart_Info.NVar_Nset,
						  max_name_length+1,
						  sizeof(char));

    /* get the variable names */
    if (ex_get_var_names(exoid, "m", Restart_Info.NVar_Nset,
                         Restart_Info.NSV_Name) < 0) {
      fprintf(stderr, "%s: Could not get nodeset variable names from file\n",
              yo);
      return -1;
    }

    /* and get the truth table */
    Restart_Info.GNset_TT = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                (Num_Node_Set
                                                * Restart_Info.NVar_Nset),
                                                sizeof(int));

    check_exodus_error(ex_get_var_tab(exoid, "m", Num_Node_Set,
				      Restart_Info.NVar_Nset,
				      Restart_Info.GNset_TT),
                       "ex_get_var_tab");
  }


#ifdef DEBUG
  if (Debug_Flag >= 2) {
    printf("\n\nRestart Parameters:\n");
    printf("\tNumber of time indices: %d\n", Restart_Info.Num_Times);
    for (cnt = 0; cnt < Restart_Info.Num_Times; cnt++)
      printf("\t\tTime index: %d\n", Restart_Info.Time_Idx[cnt]);
    printf("\tNumber of global variables: %d\n", Restart_Info.NVar_Glob);
    for (cnt = 0; cnt < Restart_Info.NVar_Glob; cnt++)
      printf("\t\tGlobal variable %d: %s\n", (cnt+1),
             Restart_Info.GV_Name[cnt]);
    printf("\tNumber of elental variables: %d\n", Restart_Info.NVar_Elem);
    for (cnt = 0; cnt < Restart_Info.NVar_Elem; cnt++)
      printf("\t\tElemental variable %d: %s\n", (cnt+1),
             Restart_Info.EV_Name[cnt]);
    printf("\tNumber of nodal variables: %d\n", Restart_Info.NVar_Node);
    for (cnt = 0; cnt < Restart_Info.NVar_Node; cnt++)
      printf("\t\tNodal variable %d: %s\n", (cnt+1), Restart_Info.NV_Name[cnt]);
  }
#endif

  return 0;

}

static int broadcast_var_param(RESTART_PTR restart, int max_name_length)
{
  int    iproc, cnt1, cnt2;

  /* broadcast the restart structure so that all processors have the lengths */
  brdcst(Proc, Num_Proc, (char *) restart, sizeof(RESTART), 0);

  /* now allocate memory on all of the other processors */
  if (Proc != 0) {
    /* allocate the space for the time indices, and the times */
    Restart_Info.Time_Idx = (int *) array_alloc (__FILE__, __LINE__, 1,
                                                 restart->Num_Times,
                                                 sizeof(int));

    /* and the space for each of the variable names */
    if (Restart_Info.NVar_Glob > 0)
      restart->GV_Name = (char **) array_alloc(__FILE__, __LINE__, 2,
                                               restart->NVar_Glob,
                                               max_name_length+1,
                                               sizeof(char));

    if (Restart_Info.NVar_Elem > 0) {
      restart->EV_Name = (char **) array_alloc(__FILE__, __LINE__, 2,
                                               restart->NVar_Elem,
                                               max_name_length+1,
                                               sizeof(char));

      Restart_Info.GElem_TT = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                  (Num_Elem_Blk
                                                  * Restart_Info.NVar_Elem),
                                                  sizeof(int));
    }


    if (Restart_Info.NVar_Node > 0)
      restart->NV_Name = (char **) array_alloc(__FILE__, __LINE__, 2,
                                               restart->NVar_Node,
                                               max_name_length+1,
                                               sizeof(char));
  }  /* End "if (Proc != 0)" */

  /* broadcast the time indices */
  brdcst(Proc, Num_Proc, (char *) restart->Time_Idx,
         (restart->Num_Times*sizeof(int)), 0);

  /* now broadcast the variable names
   *
   * (Note: arrays with dimension greater than 1, that were allocated with the
   * array_alloc function must be communicated with either a series of
   * broadcasts due to potential memory positioning conflicts or a broadcast
   * starting at an address after the pointer information.
   */
  if (restart->NVar_Glob > 0)
    brdcst(Proc, Num_Proc, restart->GV_Name[0],
           (restart->NVar_Glob*(max_name_length+1)*sizeof(char)), 0);
  if (restart->NVar_Elem > 0) {
    brdcst(Proc, Num_Proc, restart->EV_Name[0],
           (restart->NVar_Elem*(max_name_length+1)*sizeof(char)), 0);

    /* with the elemental variables, broadcast the truth table as well */
    brdcst(Proc, Num_Proc, (char *) Restart_Info.GElem_TT,
           (Num_Elem_Blk*Restart_Info.NVar_Elem*sizeof(int)), 0);

    Restart_Info.Elem_TT = (int **) array_alloc(__FILE__, __LINE__, 2,
                                                Proc_Info[2],
                                                 (Num_Elem_Blk
                                                 * Restart_Info.NVar_Elem),
                                                 sizeof(int));

    /* and copy it for each proc that this processor is responsible for */
    for (iproc = 0; iproc < Proc_Info[2]; iproc++)
      for (cnt1 = 0; cnt1 < Num_Elem_Blk; cnt1++)
        for (cnt2 = 0; cnt2 < Restart_Info.NVar_Elem; cnt2++)
          Restart_Info.Elem_TT[iproc][cnt1*Restart_Info.NVar_Elem+cnt2] =
            Restart_Info.GElem_TT[cnt1*Restart_Info.NVar_Elem+cnt2];
  }

  if (restart->NVar_Node > 0)
    brdcst(Proc, Num_Proc, restart->NV_Name[0],
           (restart->NVar_Node*(max_name_length+1)*sizeof(char)), 0);

  if (restart->NVar_Nset > 0) {
    brdcst(Proc, Num_Proc, restart->NSV_Name[0],
           (restart->NVar_Nset*(max_name_length+1)*sizeof(char)), 0);

    /* with the nodeset variables, broadcast the truth table as well */
    brdcst(Proc, Num_Proc, (char *) Restart_Info.GNset_TT,
           (Num_Node_Set*Restart_Info.NVar_Nset*sizeof(int)), 0);

    Restart_Info.Nset_TT = (int **) array_alloc(__FILE__, __LINE__, 2,
                                                Proc_Info[2],
						(Num_Node_Set * Restart_Info.NVar_Nset),
                                                 sizeof(int));

    /* and copy it for each proc that this processor is responsible for */
    for (iproc = 0; iproc < Proc_Info[2]; iproc++)
      for (cnt1 = 0; cnt1 < Num_Node_Set; cnt1++)
        for (cnt2 = 0; cnt2 < Restart_Info.NVar_Nset; cnt2++)
          Restart_Info.Nset_TT[iproc][cnt1*Restart_Info.NVar_Nset+cnt2] =
            Restart_Info.GNset_TT[cnt1*Restart_Info.NVar_Nset+cnt2];
  }

  if (restart->NVar_Sset > 0) {
    brdcst(Proc, Num_Proc, restart->SSV_Name[0],
           (restart->NVar_Sset*(max_name_length+1)*sizeof(char)), 0);

    /* with the sideset variables, broadcast the truth table as well */
    brdcst(Proc, Num_Proc, (char *) Restart_Info.GSset_TT,
           (Num_Side_Set*Restart_Info.NVar_Sset*sizeof(int)), 0);

    Restart_Info.Sset_TT = (int **) array_alloc(__FILE__, __LINE__, 2,
                                                Proc_Info[2],
						(Num_Side_Set * Restart_Info.NVar_Sset),
                                                 sizeof(int));

    /* and copy it for each proc that this processor is responsible for */
    for (iproc = 0; iproc < Proc_Info[2]; iproc++)
      for (cnt1 = 0; cnt1 < Num_Side_Set; cnt1++)
        for (cnt2 = 0; cnt2 < Restart_Info.NVar_Sset; cnt2++)
          Restart_Info.Sset_TT[iproc][cnt1*Restart_Info.NVar_Sset+cnt2] =
            Restart_Info.GSset_TT[cnt1*Restart_Info.NVar_Sset+cnt2];
  }
  return 0;

}


static int read_vars(int exoid, int index, int blk_cnt, int *eb_ids,
                     int *eb_cnts, int ***eb_map_ptr, int **eb_cnts_local,
		     int *ss_ids, int *ss_cnts, int *ns_ids, int *ns_cnts,
                     int io_ws)
{
  char  *yo="read_vars";

  void  *ptr;

  /* check to see if the io_ws is smaller than the machine precision */
  if (io_ws < sizeof(float)) io_ws = sizeof(float);

  if (io_ws == sizeof(float)) ptr = (void *) &(Restart_Info.Time_sp[blk_cnt]);
  else                        ptr = (void *) &(Restart_Info.Time_dp[blk_cnt]);
  /* first read the time */
  if (Proc == 0) {
    if (ex_get_time(exoid, index, ptr) < 0) {
      fprintf(stderr, "%s: ERROR, unable to get time for restart index %d!\n",
              yo, index);
      return -1;
    }
  }

  brdcst(Proc, Num_Proc, (char *) ptr, io_ws, 0);

  /***************** Global Variables ********************/
  /* allocate space for the global variables */
  if (Restart_Info.NVar_Glob > 0) {
    if (io_ws == sizeof(float))
      ptr = (void *) Restart_Info.Glob_Vals_sp[blk_cnt];
    else
      ptr = (void *) Restart_Info.Glob_Vals_dp[blk_cnt];

    if (Proc == 0) {
      /* get the global variables */
      if (ex_get_glob_vars(exoid, index, Restart_Info.NVar_Glob, ptr) < 0) {
        fprintf(stderr, "%s: Could not get global variables from file\n", yo);
        return -1;
      }
    }

    brdcst(Proc, Num_Proc, (char *) ptr, (Restart_Info.NVar_Glob*io_ws), 0);
  }

#ifdef DEBUG
  if (Debug_Flag >= 7) {
    int i, j;
    print_sync_start(Proc, Num_Proc, TRUE);
    printf("\n\nRestart global variables for Processor %d\n", Proc);
    printf("Time Index: %d\n", index);
    for (i = 0; i < Proc_Info[2]; i++) {
      printf("Values for processor %d\n", (i*Num_Proc+Proc));
      for (j = 0; j < Restart_Info.NVar_Glob; j++) {
        if (io_ws == sizeof(float))
          printf("Value for variable %d: %f\n", (j+1),
                 Restart_Info.Glob_Vals_sp[blk_cnt][j]);
        else
          printf("Value for variable %d: %lf\n", (j+1),
                 Restart_Info.Glob_Vals_dp[blk_cnt][j]);
      }
    }
    print_sync_end(Proc, Num_Proc, TRUE);
  }
#endif

  if (Restart_Info.NVar_Elem > 0 ) {
    printf("Reading %d element variables...\n", Restart_Info.NVar_Elem);
    if (read_elem_vars(exoid, index, blk_cnt, eb_ids, eb_cnts, eb_map_ptr,
                       eb_cnts_local, io_ws) < 0) {
      fprintf(stderr, "%s: Error distributing elemental variables.\n", yo);
      return -1;
    }
  }

  if (Restart_Info.NVar_Node > 0 ) {
    printf("Reading %d nodal variables...\n", Restart_Info.NVar_Node);
    if (read_nodal_vars(exoid, index, blk_cnt, io_ws) < 0) {
      fprintf(stderr, "%s: Error distributing nodal variables.\n", yo);
      return -1;
    }
  }

  if (Restart_Info.NVar_Sset > 0 ) {
    printf("Reading %d sideset variables...\n", Restart_Info.NVar_Sset);
    if (read_sset_vars(exoid, index, blk_cnt, ss_ids, ss_cnts, io_ws) < 0) {
      fprintf(stderr, "%s: Error distributing sideset variables.\n", yo);
      return -1;
    }
  }


  if (Restart_Info.NVar_Nset > 0 ) {
    printf("Reading %d nodeset variables...\n", Restart_Info.NVar_Nset);
    if (read_nset_vars(exoid, index, blk_cnt, ns_ids, ns_cnts, io_ws) < 0) {
      fprintf(stderr, "%s: Error distributing nodeset variables.\n", yo);
      return -1;
    }
  }

  return 0;

}

static int read_elem_vars(int exoid, int index, int blk_cnt, int *eb_ids,
                          int *eb_cnts, int ***eb_map_ptr,
                          int **eb_cnts_local, int io_ws)
{
  int     iblk, iproc;
  int     num_mesgs, *mesg_start, max_elem_per_mesg;
  int     eb_offset=0, *local_offset;

  /* check to see if the io_ws is smaller than the machine precision */
  if (io_ws < sizeof(float)) io_ws = sizeof(float);

  /* to speed up searches, keep track of element blocks offset on each proc */
  local_offset = (int *) array_alloc (__FILE__, __LINE__, 1, Proc_Info[2],
                                      sizeof(int));
  for (iproc = 0; iproc < Proc_Info[2]; iproc++)
    local_offset[iproc] = 0;

  /* loop over the number of element blocks */
  for (iblk = 0; iblk < Num_Elem_Blk; iblk++) {

    /* calculate the message length for this element block */
    num_mesgs = break_message_up(io_ws, eb_cnts[iblk], MAX_CHUNK_SIZE/2,
                                 &mesg_start);
    max_elem_per_mesg = mesg_start[1] - mesg_start[0];

#ifdef DEBUG
    if (Debug_Flag >= 2 && Proc == 0) {
      printf("\n\nMessage summary for elemental variables:\n");
      printf("\tTime Index: %d\n", index);
      printf("\tElemental Block Id: %d\n", eb_ids[iblk]);
      printf("\tNumber of elemental variables: %d\n", Restart_Info.NVar_Elem);
      printf("\tNumber of messages per elemental variable: %d\n", num_mesgs);
      printf("\tMax message size per elemental variable: %d\n",
             max_elem_per_mesg);
      printf("\tMin message size per elemental variable: %d\n",
             mesg_start[num_mesgs]-mesg_start[num_mesgs-1]);
    }
#endif

    if (Num_Proc == 1) {
      read_elem_vars_1(exoid, index, blk_cnt, eb_ids,
		       eb_cnts, eb_map_ptr, eb_cnts_local, io_ws,
		       iblk, eb_offset, local_offset);
    } else {
      read_elem_vars_n(exoid, index, blk_cnt, eb_ids,
		       eb_cnts, eb_map_ptr, eb_cnts_local, io_ws,
		       iblk, eb_offset, local_offset,
		       max_elem_per_mesg, num_mesgs, mesg_start);
    }

    /* need to keep track of this for the element number map */
    eb_offset += eb_cnts[iblk];

    /* need to set up local offsets for next block */
    for (iproc = 0; iproc < Proc_Info[2]; iproc++)
      local_offset[iproc] += eb_cnts_local[iproc][iblk];

    /* free up memory for the next element block */
    safe_free((void **) &mesg_start);

  } /* End "for (iblk = 0; iblk < Num_Elem_Blk; iblk++)" */

  safe_free((void **) &local_offset);


#ifdef DEBUG
  if (Debug_Flag >= 7) {
    int ielem, num_elem, ivar, var_offset, elem_loc;
    print_sync_start(Proc, Num_Proc, TRUE);
    printf("\n\nRestart elemental variables for Processor %d\n", Proc);
    printf("Time Index: %d\n", index);
    for (iproc = 0; iproc < Proc_Info[2]; iproc++) {
      printf("Values for processor %d\n", (iproc*Num_Proc+Proc));
      num_elem = Num_Internal_Elems[iproc] + Num_Border_Elems[iproc];
      for (ivar = 0; ivar < Restart_Info.NVar_Elem; ivar++) {
        var_offset = ivar * num_elem;
        printf("Values for variable %d\n", (ivar+1));
        for (ielem = 0; ielem < num_elem; ielem++) {
          elem_loc = var_offset + ielem;
          if (io_ws == sizeof(float))
            printf("%f  ",
                   Restart_Info.Elem_Vals_sp[iproc][blk_cnt][elem_loc]);
          else
            printf("%lf  ",
                   Restart_Info.Elem_Vals_dp[iproc][blk_cnt][elem_loc]);
        }
        printf("\n");
      }
    }
    print_sync_end(Proc, Num_Proc, TRUE);
  }
#endif

  return 0;
}

static int read_elem_vars_1(int exoid, int index, int blk_cnt, int *eb_ids,
			    int *eb_cnts, int ***eb_map_ptr,
			    int **eb_cnts_local, int io_ws, int iblk,
			    int eb_offset, int *local_offset)
{
  int     i1, ivar, var_offset, iproc;
  int     elem_loc;
  int    *elem_map, num_elem;
  float  *sp_vals;
  double *dp_vals;
  void   *ptr;

  
  /* Allocate memory for temporary storage */
  ptr  = array_alloc(__FILE__, __LINE__, 1, eb_cnts[iblk], io_ws);
  /* now sort out if this is single or double precision */
  if (io_ws == sizeof(float)) sp_vals = (float *) ptr;
  else                        dp_vals = (double *) ptr;
  
  /* now loop over each variable */
  for (ivar = 0; ivar < Restart_Info.NVar_Elem; ivar++) {

    /* check if this variable exists for this element block */
    if (Restart_Info.GElem_TT[iblk*Restart_Info.NVar_Elem+ivar]) {

      /*
       * Read in the specified element variable values and their associated
       * global FEM element numbers.
       */

      check_exodus_error(ex_get_elem_var(exoid,
					 index,
					 (ivar+1),
					 eb_ids[iblk],
					 eb_cnts[iblk],
					 ptr),
			 "ex_get_elem_var");
      
      /*
       * Find out which FEM elements belong on this processor and copy
       * them to the restart vector.
       */
      for (iproc = 0; iproc < Proc_Info[2]; iproc++) {
	  
	/* check to see if this element block needs this variable */
	if (Restart_Info.Elem_TT[iproc][iblk*Restart_Info.NVar_Elem+ivar]) {
	    
	  /* calculate the offset for this variable */
	  var_offset = ivar * (Num_Internal_Elems[iproc] + Num_Border_Elems[iproc]);
	    
	  elem_map = eb_map_ptr[iproc][iblk];
	  num_elem = eb_cnts_local[iproc][iblk];
	    
	  for (i1 = 0; i1 < num_elem; i1++) {
	    elem_loc = var_offset + i1 + local_offset[iproc];
		
	    if (io_ws == sizeof(float))
	      Restart_Info.Elem_Vals_sp[iproc][blk_cnt][elem_loc]
		= sp_vals[elem_map[i1]-eb_offset];
	    else
	      Restart_Info.Elem_Vals_dp[iproc][blk_cnt][elem_loc]
		= dp_vals[elem_map[i1]-eb_offset];
	  }
	}
      }
    } /* End "if (Restart_Info.GElem_TT[...])" */
  } 
  safe_free((void **) &ptr);
  return 0;
}
    
static int read_elem_vars_n(int exoid, int index, int blk_cnt, int *eb_ids,
			    int *eb_cnts, int ***eb_map_ptr,
			    int **eb_cnts_local, int io_ws, int iblk,
			    int eb_offset,int *local_offset,
			    int max_elem_per_mesg, int num_mesgs,
			    int *mesg_start)
{
  int     i1, ivar, var_offset, imsg, iproc;
  int     elem_loc;
  int    *glob_elem, *elem_map, *proc_elem_inter, num_elem;
  int     istart_elem, iend_elem, num_ev_in_mesg;
  float  *sp_vals;
  double *dp_vals;
  void   *ptr;

  /* Allocate memory for temporary storage */
  ptr  = array_alloc(__FILE__, __LINE__, 1, max_elem_per_mesg, io_ws);
  glob_elem = (int *)   array_alloc(__FILE__, __LINE__, 1,
				    2*max_elem_per_mesg, sizeof(int));
  proc_elem_inter = glob_elem + max_elem_per_mesg;

  /* now loop over each variable */
  for (ivar = 0; ivar < Restart_Info.NVar_Elem; ivar++) {


    /* check if this variable exists for this element block */
    if (Restart_Info.GElem_TT[iblk*Restart_Info.NVar_Elem+ivar]) {

      /* loop over the number of messages */
      for (imsg = 0; imsg < num_mesgs; imsg++) {

	istart_elem = mesg_start[imsg];
	iend_elem   = mesg_start[imsg+1];
	num_ev_in_mesg = iend_elem - istart_elem;

        /*
         * Read in the specified element variable values and their associated
         * global FEM element numbers.
         */

        if (Proc == 0) {
	  for (i1 = 0; i1 < num_ev_in_mesg; i1++)
	    glob_elem[i1] = istart_elem+eb_offset+i1;

	  check_exodus_error(ne_get_n_elem_var(exoid,
					       index,
					       (ivar+1),
					       eb_ids[iblk],
					       eb_cnts[iblk],
					       istart_elem+1,
					       num_ev_in_mesg,
					       ptr),
			     "ne_get_n_nodal_var");

	}

	brdcst_maxlen(Proc, Num_Proc, (char *)glob_elem,
		      num_ev_in_mesg * sizeof(int), 0);

	brdcst_maxlen(Proc, Num_Proc, (char *)ptr,
		      num_ev_in_mesg * io_ws, 0);

	psync(Proc, Num_Proc);

	/* now sort out if this is single or double precision */
	if (io_ws == sizeof(float)) sp_vals = (float *) ptr;
	else                        dp_vals = (double *) ptr;

	/* Check for monotonic listing of elements in glob_elem */
	if (!check_monot(glob_elem, num_ev_in_mesg)) {

	  if (io_ws == sizeof(float))
	    sortN_int_float(num_ev_in_mesg, glob_elem, 1, sp_vals);
	  else
	    sortN_int_double(num_ev_in_mesg, glob_elem, 1, dp_vals);
	}


	/*
	 * Find out which FEM elements belong on this processor and copy
	 * them to the restart vector.
	 */
	for (iproc = 0; iproc < Proc_Info[2]; iproc++) {

	  /* check to see if this element block needs this variable */
	  if (Restart_Info.Elem_TT[iproc][
					  iblk*Restart_Info.NVar_Elem+ivar]) {

	    /* calculate the offset for this variable */
	    var_offset = ivar * (Num_Internal_Elems[iproc]
				 + Num_Border_Elems[iproc]);

	    /*
	     * Get the part of the global map for this element block
	     * on this processor. This portion of the map is sorted
	     * when it is read in, and so it is monotonic.
	     */
	    elem_map = eb_map_ptr[iproc][iblk];
	    num_elem = eb_cnts_local[iproc][iblk];

	    find_inter_pos(proc_elem_inter, num_ev_in_mesg,
			   glob_elem, num_elem, elem_map, 2);

	    for (i1 = 0; i1 < num_ev_in_mesg; i1++) {
	      if (proc_elem_inter[i1] >= 0) {
		elem_loc = var_offset + proc_elem_inter[i1] +
		  local_offset[iproc];

		if (io_ws == sizeof(float))
		  Restart_Info.Elem_Vals_sp[iproc][blk_cnt][elem_loc]
		    = sp_vals[i1];
		else
		  Restart_Info.Elem_Vals_dp[iproc][blk_cnt][elem_loc]
		    = dp_vals[i1];
	      }
	    }
	  }
	}

      } /* End "for (imsg = 0; imsg < num_mesgs; imsg++)" */
    } /* End "if (Restart_Info.GElem_TT[...])" */
  } /* End "for (ivar = 0; ivar < Restart_Info.NVar_Elem; ivar++)" */
  safe_free((void **) &ptr);
  safe_free((void **) &glob_elem);
  return 0;
}

static int read_sset_vars(int exoid, int index, int blk_cnt, int *ss_ids,
                          int *ss_cnts, int io_ws)
{
  int     iset, num_mesgs, *mesg_start, max_elem_per_mesg;

  /* loop over the number of side sets */
  for (iset = 0; iset < Num_Side_Set; iset++) {

    /* calculate the message length for this sideset */
    num_mesgs = break_message_up(io_ws, ss_cnts[iset], MAX_CHUNK_SIZE/2,
                                 &mesg_start);
    max_elem_per_mesg = mesg_start[1] - mesg_start[0];

    if (Num_Proc == 1) {
      read_sset_vars_1(exoid, index, blk_cnt, ss_ids, ss_cnts, io_ws, iset);
    } else {
      printf("Read_sset_vars currently only supports single message \n");
      abort();
      read_sset_vars_n(exoid, index, blk_cnt, ss_ids, ss_cnts, io_ws,
		       iset, max_elem_per_mesg, num_mesgs, mesg_start);
    }
    /* free up memory... */
    safe_free((void **) &mesg_start);
  }
  return 0;
}

static int read_sset_vars_1(int exoid, int index, int blk_cnt, int *ss_ids,
			    int *ss_cnts, int io_ws, int iset)
{
  int     i, i1, ivar, var_offset, iproc;
  float  *sp_vals;
  double *dp_vals;
  void   *ptr;

  
  /* Allocate memory for temporary storage */
  ptr  = array_alloc(__FILE__, __LINE__, 1, ss_cnts[iset], io_ws);
  /* now sort out if this is single or double precision */
  if (io_ws == sizeof(float)) sp_vals = (float *) ptr;
  else                        dp_vals = (double *) ptr;
  
  /* now loop over each variable */
  var_offset = 0;
  for (ivar = 0; ivar < Restart_Info.NVar_Sset; ivar++) {

    /* check if this variable exists for this set */
    if (Restart_Info.GSset_TT[iset*Restart_Info.NVar_Sset+ivar]) {

      /* Read in the specified variable values */
      check_exodus_error(ex_get_sset_var(exoid, index, (ivar+1), ss_ids[iset], ss_cnts[iset],
					 ptr), "ex_get_sset_var");
      
      for (iproc = 0; iproc < Proc_Info[2]; iproc++) {
	int ss_offset = 0;
	var_offset = ivar * Proc_SS_Elem_List_Length[iproc];
	for (i = 0; i < Proc_Num_Side_Sets[iproc]; i++) {
	  if (Proc_SS_Ids[iproc][i] == ss_ids[iset]) {
	
	    int num_elem = Proc_SS_Elem_Count[iproc][i];
	    for (i1 = 0; i1 < num_elem; i1++) {
	      int gelem_loc = Proc_SS_GEMap_List[iproc][i1+ss_offset];
	      assert(gelem_loc < ss_cnts[iset]);
	      if (io_ws == sizeof(float))
		Restart_Info.Sset_Vals_sp[iproc][blk_cnt][i1+ss_offset+var_offset]
		  = sp_vals[gelem_loc];
	      else
		Restart_Info.Sset_Vals_dp[iproc][blk_cnt][i1+ss_offset+var_offset]
		  = dp_vals[gelem_loc];
	    }
	    break;
	  }
	  ss_offset += Proc_SS_Elem_Count[iproc][i];
	}
      } 
    } 
  }
  safe_free((void **) &ptr);
  return 0;
}
    
static int read_sset_vars_n(int exoid, int index, int blk_cnt, int *ss_ids,
			    int *ss_cnts, int io_ws, int iset,
			    int max_elem_per_mesg, int num_mesgs,
			    int *mesg_start)
{
  return -1;
}

static int read_nset_vars(int exoid, int index, int blk_cnt, int *ns_ids,
                          int *ns_cnts, int io_ws)
{
  int     iset, num_mesgs, *mesg_start, max_elem_per_mesg;

  /* loop over the number of node sets */
  for (iset = 0; iset < Num_Node_Set; iset++) {

    /* calculate the message length for this nodeset */
    num_mesgs = break_message_up(io_ws, ns_cnts[iset], MAX_CHUNK_SIZE/2,
                                 &mesg_start);
    max_elem_per_mesg = mesg_start[1] - mesg_start[0];

    if (Num_Proc == 1) {
      read_nset_vars_1(exoid, index, blk_cnt, ns_ids, ns_cnts, io_ws, iset);
    } else {
      printf("Read_nset_vars currently only supports single message \n");
      abort();
      read_nset_vars_n(exoid, index, blk_cnt, ns_ids, ns_cnts, io_ws,
		       iset, max_elem_per_mesg, num_mesgs, mesg_start);
    }
    /* free up memory... */
    safe_free((void **) &mesg_start);
  }
  return 0;
}

static int read_nset_vars_1(int exoid, int index, int blk_cnt, int *ns_ids,
			    int *ns_cnts, int io_ws, int iset)
{
  int     i, i1, ivar, var_offset, iproc;
  float  *sp_vals;
  double *dp_vals;
  void   *ptr;

  
  /* Allocate memory for temporary storage */
  ptr  = array_alloc(__FILE__, __LINE__, 1, ns_cnts[iset], io_ws);
  /* now sort out if this is single or double precision */
  if (io_ws == sizeof(float)) sp_vals = (float *) ptr;
  else                        dp_vals = (double *) ptr;
  
  /* now loop over each variable */
  var_offset = 0;
  for (ivar = 0; ivar < Restart_Info.NVar_Nset; ivar++) {

    /* check if this variable exists for this set */
    if (Restart_Info.GNset_TT[iset*Restart_Info.NVar_Nset+ivar]) {

      /* Read in the specified variable values */
      check_exodus_error(ex_get_nset_var(exoid, index, (ivar+1), ns_ids[iset], ns_cnts[iset],
					 ptr), "ex_get_nset_var");
      
      for (iproc = 0; iproc < Proc_Info[2]; iproc++) {
	int ns_offset = 0;
	var_offset = ivar * Proc_NS_List_Length[iproc];
	for (i = 0; i < Proc_Num_Node_Sets[iproc]; i++) {
	  if (Proc_NS_Ids[iproc][i] == ns_ids[iset]) {
	
	    int num_elem = Proc_NS_Count[iproc][i];
	    for (i1 = 0; i1 < num_elem; i1++) {
	      int gelem_loc = Proc_NS_GNMap_List[iproc][i1+ns_offset];
	      assert(gelem_loc < ns_cnts[iset]);
	      if (io_ws == sizeof(float))
		Restart_Info.Nset_Vals_sp[iproc][blk_cnt][i1+ns_offset+var_offset]
		  = sp_vals[gelem_loc];
	      else
		Restart_Info.Nset_Vals_dp[iproc][blk_cnt][i1+ns_offset+var_offset]
		  = dp_vals[gelem_loc];
	    }
	    break;
	  }
	  ns_offset += Proc_NS_Count[iproc][i];
	}
      } 
    } 
  }
  safe_free((void **) &ptr);
  return 0;
}
    
static int read_nset_vars_n(int exoid, int index, int blk_cnt, int *ns_ids,
			    int *ns_cnts, int io_ws, int iset,
			    int max_elem_per_mesg, int num_mesgs,
			    int *mesg_start)
{
  return -1;
}

static int read_nodal_vars_1 (int exoid, int index, int blk_cnt, int io_ws)
{
  /* Same routine as read_nodal_vars except there is a single message.
     This is a first step at removing some complexity to get some better
     speed.
  */
     
  int    var_num, node_loc, loc_count;
  int     iproc, i2;
  int     var_offset;
  float  *sp_vals;
  double *dp_vals;
  void   *ptr;
  /*---------------------------Begin Execution--------------------------------*/

  /* Allocate memory for temporary storage */
  ptr  = array_alloc(__FILE__, __LINE__, 1, Num_Node, io_ws);

  /* figure out which pointer to use */
  if (io_ws == sizeof(float)) sp_vals = ptr;
  else                       dp_vals = ptr;
  
  /* Loop over each auxiliary variable */
  for (var_num = 0; var_num < Restart_Info.NVar_Node; var_num++) {
    /*
     * Read in the specified nodal variable values and their associated
     * global FEM node numbers.
     */
    check_exodus_error(ne_get_n_nodal_var(exoid, index, (var_num+1),
					  1, Num_Node, ptr),
		       "ne_get_n_nodal_var");
    
    /*
     * Find out which FEM nodes belong on this processor and copy
     * them to the restart vector.
     */
    for (iproc = 0; iproc < Proc_Info[2]; iproc++) {

      /* calculate the offset for this variable */
      loc_count = Num_Internal_Nodes[iproc]+Num_Border_Nodes[iproc]+
	Num_External_Nodes[iproc];

      var_offset = var_num * loc_count;
	
      if (io_ws == sizeof(float)) {
	for (i2 = 0; i2 < loc_count; i2++) {
	  node_loc = var_offset + i2;
	  Restart_Info.Node_Vals_sp[iproc][blk_cnt][node_loc] =
	    sp_vals[GNodes[iproc][i2]-1];
	}
      } else {
	for (i2 = 0; i2 < loc_count; i2++) {
	  node_loc = var_offset + i2;
	  Restart_Info.Node_Vals_dp[iproc][blk_cnt][node_loc] =
	    dp_vals[GNodes[iproc][i2]-1];
	}
      }
    }

  } /* End "for (var_num = 0; var_num < Restart_Info.NVar_Node; var_num++)" */

  safe_free((void **) &ptr);
  return 0;
}

static int read_nodal_vars (int exoid, int index, int blk_cnt, int io_ws)
{

  int    *mesg_start, num_mesgs, max_nv_per_mesg, var_num, node_loc;
  int     iproc, i1, i2, istart_node, iend_node;
  int    *glob_node, *proc_node_inter, num_nv_in_mesg;
  int     var_offset;
  float  *sp_vals;
  double *dp_vals;
  void   *ptr;

  /*---------------------------Begin Execution--------------------------------*/

  /* check to see if the io_ws is smaller than the machine precision */
  if (io_ws < sizeof(float)) io_ws = sizeof(float);

  /*
   * Find out how many pieces to break the read of the nodal variables up
   * into. Note that we're sending both the values AND their global node IDs,
   * thus the MAX_CHUNK_SIZE/2.
   */

  num_mesgs = break_message_up(io_ws, Num_Node, MAX_CHUNK_SIZE/2, &mesg_start);

  if (Num_Proc == 1) {
    safe_free((void **) &mesg_start);
    return read_nodal_vars_1(exoid, index, blk_cnt, io_ws);
  }
  /* Below here we are in parallel or have more than 1 message... */

  max_nv_per_mesg = mesg_start[1] - mesg_start[0];

#ifdef DEBUG
  if (Debug_Flag >= 2 && Proc == 0) {
    printf("\n\nMessage summary for nodal variables:\n");
    printf("\tTime Index: %d\n", index);
    printf("\tNumber of nodal variables: %d\n", Restart_Info.NVar_Node);
    printf("\tNumber of messages per nodal variable: %d\n", num_mesgs);
    printf("\tMax message size per nodal variable: %d\n",
	   max_nv_per_mesg);
    printf("\tMin message size per nodal variable: %d\n",
	   mesg_start[num_mesgs]-mesg_start[num_mesgs-1]);
  }
#endif

  /* Allocate memory for temporary storage */

  ptr  = array_alloc(__FILE__, __LINE__, 1, max_nv_per_mesg, io_ws);
  glob_node = (int *)   array_alloc(__FILE__, __LINE__, 1,
				    2*max_nv_per_mesg, sizeof(int));
  proc_node_inter = glob_node + max_nv_per_mesg;

  /* Loop over each auxiliary variable */

  for (var_num = 0; var_num < Restart_Info.NVar_Node; var_num++) {

    /* Loop over each message */

    for (i1 = 0; i1 < num_mesgs; i1++) {

      istart_node = mesg_start[i1];
      iend_node   = mesg_start[i1+1];
      num_nv_in_mesg = iend_node - istart_node;

      /*
       * Read in the specified nodal variable values and their associated
       * global FEM node numbers.
       */

      if (Proc == 0) {
	for (i2 = 0; i2 < num_nv_in_mesg; i2++) {
	  glob_node[i2] = istart_node+1+i2;
	}

	check_exodus_error(ne_get_n_nodal_var(exoid,
					      index,
					      (var_num+1),
					      istart_node+1,
					      num_nv_in_mesg,
					      ptr),
			   "ne_get_n_nodal_var");

      }

      /* Broadcast the information to the other processors */

      brdcst_maxlen(Proc, Num_Proc, (char *)glob_node,
		    num_nv_in_mesg * sizeof(int), 0);

      psync(Proc, Num_Proc);

      brdcst_maxlen(Proc, Num_Proc, (char *)ptr, num_nv_in_mesg * io_ws, 0);

      psync(Proc, Num_Proc);

      /* figure out which pointer to use */
      if (io_ws == sizeof(float)) sp_vals = ptr;
      else                       dp_vals = ptr;


      /*
       * Find out which FEM nodes belong on this processor and copy
       * them to the restart vector.
       */
      for (iproc = 0; iproc < Proc_Info[2]; iproc++) {

	/* calculate the offset for this variable */
	var_offset = var_num * (Num_Internal_Nodes[iproc]
				+ Num_Border_Nodes[iproc] + Num_External_Nodes[iproc]);


	/*
	 * Find the intersection of the nodes contained on this processor
	 * and the nodes in this message.
	 */

	find_gnode_inter(proc_node_inter, num_nv_in_mesg,
			 glob_node, Num_Internal_Nodes[iproc],
			 Num_Border_Nodes[iproc],
			 Num_External_Nodes[iproc],
			 GNodes[iproc]);

	for (i2=0; i2 < num_nv_in_mesg; i2++) {
	  if (proc_node_inter[i2] >= 0) {
	    node_loc = var_offset + proc_node_inter[i2];

	    if (io_ws == sizeof(float)) {
	      Restart_Info.Node_Vals_sp[iproc][blk_cnt][node_loc] =
		sp_vals[i2];
	    }
	    else {
	      Restart_Info.Node_Vals_dp[iproc][blk_cnt][node_loc] =
		dp_vals[i2];
	    }
	  }
	}
      }
    } /* End "for (i1 = 0; i1 < num_mesgs; i1++)" */
  } /* End "for (var_num = 0; var_num < Restart_Info.NVar_Node; var_num++)" */

#ifdef DEBUG
  if (Debug_Flag >= 7) {
    int num_node;
    int i3;
    print_sync_start(Proc, Num_Proc, TRUE);
    printf("\n\nRestart nodal variables for Processor %d\n", Proc);
    printf("Time Index: %d\n", index);
    for (i1 = 0; i1 < Proc_Info[2]; i1++) {
      printf("Values for processor %d\n", (i1*Num_Proc+Proc));
      num_node = Num_Internal_Nodes[i1] + Num_Border_Nodes[i1]
	+ Num_External_Nodes[i1];
      for (i2 = 0; i2 < Restart_Info.NVar_Node; i2++) {
	var_offset = i2 * num_node;
	printf("Values for variable %d\n", (i2+1));
	for (i3 = 0; i3 < num_node; i3++) {
	  node_loc = var_offset + i3;
	  if (io_ws == sizeof(float))
	    printf("%f  ", Restart_Info.Node_Vals_sp[i1][blk_cnt][node_loc]);
	  else
	    printf("%f  ", Restart_Info.Node_Vals_dp[i1][blk_cnt][node_loc]);
	}
	printf("\n");
      }
    }
    print_sync_end(Proc, Num_Proc, TRUE);
  }
#endif

  safe_free((void **) &mesg_start);
  safe_free((void **) &ptr);
  safe_free((void **) &glob_node);

  return 0;
}


static int compare_mesh_param(int exoid)
{
  char    Title[MAX_LINE_LENGTH+1];

  int     ndim, nnode, nelem, nelem_blk, nns, nss;
  int     ret = 1, error;

  memset(Title, '\0', MAX_LINE_LENGTH*sizeof(char));
  error = ex_get_init (exoid, Title, &ndim, &nnode, &nelem, &nelem_blk,
		       &nns, &nss);
  check_exodus_error (error, "ex_get_init");

  /* now check that the parameters match those retrieved from the mesh file */
  if (ndim           != Num_Dim)      ret = 0;
  else if (nnode     != Num_Node)     ret = 0;
  else if (nelem     != Num_Elem)     ret = 0;
  else if (nelem_blk != Num_Elem_Blk) ret = 0;
  else if (nns       != Num_Node_Set) ret = 0;
  else if (nss       != Num_Side_Set) ret = 0;

  return (ret);
}

  /*****************************************************************************/
int find_gnode_inter(int *intersect, int num_g_nodes, int *glob_vec,
		     int num_int_nodes, int num_bor_nodes,
		     int num_ext_nodes, int *loc_vec)

/*
 * This function assumes that glob_vec is monotonic and that loc_vec is
 * monotonic for each of the internal, border and external node IDs it
 * contains.
 */
{
  int i1, i2, min_set1, max_set1, min_set2, max_set2, offset;
  int count=0;

  /* Initialize the intersect vector */
  for (i1=0; i1 < num_g_nodes; i1++)
    intersect[i1] = -1;

  /* Check for the possibility of an intersection */
  min_set1 = glob_vec[0];
  max_set1 = glob_vec[num_g_nodes-1];

  /* Search through the internal nodes */
  if (num_int_nodes > 0) {
    min_set2 = loc_vec[0];
    max_set2 = loc_vec[num_int_nodes-1];

    if ( (max_set2 >= min_set1) && (min_set2 <= max_set1) ) {
      for (i1=0, i2=0; i1 < num_g_nodes; i1++) {
	while ( (i2 < (num_int_nodes-1)) && (glob_vec[i1] > loc_vec[i2]) ) i2++;
	if (glob_vec[i1] == loc_vec[i2]) {
	  intersect[i1] = i2;
	  count++;
	}
      }
    }
  }

  /* Search through the border nodes */
  if (num_bor_nodes > 0) {
    min_set2 = loc_vec[num_int_nodes];
    max_set2 = loc_vec[num_int_nodes+num_bor_nodes-1];

    offset = num_int_nodes;

    if ( (max_set2 >= min_set1) && (min_set2 <= max_set1) ) {
      for (i1=0, i2=0; i1 < num_g_nodes; i1++) {
	while ( (i2 < (num_bor_nodes-1)) &&
		(glob_vec[i1] > loc_vec[offset+i2]) ) i2++;

	if (glob_vec[i1] == loc_vec[offset+i2]) {
	  intersect[i1] = offset+i2;
	  count++;
	}
      }
    }
  }

  /* Search through the external nodes */
  if (num_ext_nodes > 0) {
    min_set2 = loc_vec[num_int_nodes+num_bor_nodes];
    max_set2 = loc_vec[num_int_nodes+num_bor_nodes+num_ext_nodes-1];

    offset = num_int_nodes + num_bor_nodes;

    if ( (max_set2 >= min_set1) && (min_set2 <= max_set1) ) {
      for (i1=0, i2=0; i1 < num_g_nodes; i1++) {
	while ( (i2 < (num_ext_nodes-1)) &&
		(glob_vec[i1] > loc_vec[offset+i2]) ) i2++;

	if (glob_vec[i1] == loc_vec[offset+i2]) {
	  intersect[i1] = offset+i2;
	  count++;
	}
      }
    }
  }

#ifdef DEBUG
  assert(count < num_g_nodes &&
	 "find_gnode_inter ERROR: Catastrophic error in node structure observed\n");
#endif

  return count;
}

#if defined(__PUMAGON__)
#include <stdio.h>
#else
#include <unistd.h>
#endif
#include <limits.h>

int get_free_descriptor_count(void)
{
  /* Returns maximum number of files that one process can have open
   * at one time. (POSIX)
   */
#if defined(__PUMAGON__)
  int fdmax = FOPEN_MAX;
#else
  int fdmax = sysconf(_SC_OPEN_MAX);
  if (fdmax == -1) {
    /* POSIX indication that there is no limit on open files... */
    fdmax = INT_MAX;
  }
#endif
  /* File descriptors are assigned in order (0,1,2,3,...) on a per-process
   * basis.
   *
   * Assume that we have stdin, stdout, stderr, and input exodus
   * file (4 total).
   */

  return fdmax - 4;

  /* Could iterate from 0..fdmax and check for the first EBADF (bad
   * file descriptor) error return from fcntl, but that takes too long
   * and may cause other problems.  There is a isastream(filedes) call
   * on Solaris that can be used for system-dependent code.
   *
   * Another possibility is to do an open and see which descriptor is
   * returned -- take that as 1 more than the current count of open files.
   */
}
