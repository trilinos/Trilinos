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
#include <assert.h>                     // for assert
#include <limits.h>                     // for INT_MAX
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for fprintf, stderr, NULL, etc
#include <stdlib.h>                     // for exit, free, malloc
#include <string.h>                     // for strcpy, strrchr, memset, etc
#include <unistd.h>                     // for sysconf, _SC_OPEN_MAX
#include <vector>                       // for vector
#include "exodusII.h"                   // for ex_close, etc
#include "nem_spread.h"                 // for NemSpread, etc
#include "pe_common.h"                  // for MAX_CHUNK_SIZE
#include "ps_pario_const.h"             // for Par_Nem_File_Name, PIO_Info, etc
#include "rf_allo.h"                    // for array_alloc, safe_free
#include "rf_io_const.h"                // for Exo_Res_File, ExoFile, etc
#include "rf_util.h"                    // for break_message_up

namespace {
  int get_free_descriptor_count(void);

  template <typename INT>
  size_t find_gnode_inter(INT *intersect, size_t num_g_nodes, INT *glob_vec,
			size_t num_int_nodes, size_t num_bor_nodes,
			size_t num_ext_nodes, INT *loc_vec);
}

#define TOPTR(x) (x.empty() ? NULL : &x[0])

/*****************************************************************************/
/*****************************************************************************/
template void NemSpread<double,int>::read_restart_params();
template void NemSpread<float, int>::read_restart_params();

template void NemSpread<double,int64_t>::read_restart_params();
template void NemSpread<float, int64_t>::read_restart_params();

template <typename T, typename INT>
void NemSpread<T,INT>::read_restart_params()

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
 *
 *----------------------------------------------------------------------------
 */

{
  const char  *yo="read_restart_params";

  int    exoid, cpu_ws=0;
  float  vers;
  int    max_name_length = 0;
  
  /* Open the ExodusII file */
  cpu_ws = io_ws;
  int mode = EX_READ | int64api;
  if ((exoid=ex_open(Exo_Res_File, mode, &cpu_ws, &io_ws, &vers)) < 0) {
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

  if (broadcast_var_param(&Restart_Info, max_name_length) < 0) {
    fprintf(stderr, "%s: Error occured while broadcasting variable params\n",
            yo);
    exit(1);
  }

  return;
}

template void NemSpread<double,int>::read_restart_data();
template void NemSpread<float, int>::read_restart_data();

template void NemSpread<double,int64_t>::read_restart_data();
template void NemSpread<float, int64_t>::read_restart_data();

template <typename T, typename INT>
void NemSpread<T,INT>::read_restart_data ()

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
  const char  *yo="read_restart_data";

  /* need to get the element block ids and counts */
  std::vector<INT> eb_ids_global(globals.Num_Elem_Blk);
  std::vector<INT> eb_cnts_global(globals.Num_Elem_Blk);
  std::vector<INT> ss_ids_global(globals.Num_Side_Set);
  std::vector<INT> ss_cnts_global(globals.Num_Side_Set);
  std::vector<INT> ns_ids_global(globals.Num_Node_Set);
  std::vector<INT> ns_cnts_global(globals.Num_Node_Set);

  INT ***eb_map_ptr = NULL, **eb_cnts_local = NULL;
  int    num_blocks, times_in_blk, iblk, time_idx;
  int    exoid=0, *par_exoid = NULL;

  float  vers;
  char   cTemp[512];

  /* computing precision should be the same as the database precision
   *
   * EXCEPTION: if the io_ws is smaller than the machine precision,
   * ie - database with io_ws == 4 on a Cray (sizeof(float) == 8),
   * then the cpu_ws must be the machine precision.
   */
  int cpu_ws;
  if (io_ws < (int)sizeof(float)) cpu_ws = sizeof(float);
  else                            cpu_ws = io_ws;

  /* Open the ExodusII file */
  {
    cpu_ws = io_ws;
    int mode = EX_READ | int64api;
    if ((exoid=ex_open(Exo_Res_File, mode, &cpu_ws, &io_ws, &vers)) < 0) {
      fprintf(stderr, "%s: Could not open file %s for restart info\n",
	      yo, Exo_Res_File);
      exit(1);
    }
  }

  /* allocate space for the global variables */
  Restart_Info.Glob_Vals.resize(Restart_Info.NVar_Glob);

  if (Restart_Info.NVar_Elem > 0 ) {

    /* allocate storage space */
    Restart_Info.Elem_Vals = (T **) array_alloc (__FILE__, __LINE__,
						 1, Proc_Info[2],
						 sizeof(T *));

    /* now allocate storage for the values */
    for (int iproc = 0; iproc <Proc_Info[2]; iproc++) {
      size_t array_size = Restart_Info.NVar_Elem *
	(globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc]);
      Restart_Info.Elem_Vals[iproc] = (T*)array_alloc (__FILE__, __LINE__, 1,
			     array_size , cpu_ws);
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
    
    /* Get the Element Block IDs from the input file */
    if (ex_get_ids (exoid, EX_ELEM_BLOCK, TOPTR(eb_ids_global)) < 0)
      {
	fprintf(stderr, "%s: unable to get element block IDs", yo);
	exit(1);
      }

    /* Get the count of elements in each element block */
    for (int cnt = 0; cnt < globals.Num_Elem_Blk; cnt++) {
      if (ex_get_block(exoid, EX_ELEM_BLOCK, eb_ids_global[cnt], cTemp,
		       &(eb_cnts_global[cnt]), NULL, NULL, NULL, NULL) < 0) {
	fprintf(stderr, "%s: unable to get element count for block id %lu",
		yo, (size_t)eb_ids_global[cnt]);
	exit(1);
      }
    }

    /*
     * in order to speed up finding matches in the global element
     * number map, set up an array of pointers to the start of
     * each element block's global element number map. That way
     * only entries for the current element block have to be searched
     */
    eb_map_ptr = (INT ***) array_alloc (__FILE__, __LINE__, 2,Proc_Info[2],
                                        globals.Num_Elem_Blk, sizeof(INT *));
    if (!eb_map_ptr) {
      fprintf(stderr, "[%s]: ERROR, insufficient memory!\n", yo);
      exit(1);
    }
    eb_cnts_local = (INT **) array_alloc (__FILE__, __LINE__, 2,Proc_Info[2],
                                          globals.Num_Elem_Blk, sizeof(INT));
    if (!eb_cnts_local) {
      fprintf(stderr, "[%s]: ERROR, insufficient memory!\n", yo);
      exit(1);
    }

    /*
     * for now, assume that element blocks have been
     * stored in the same order as the global blocks
     */
    for (int iproc = 0; iproc <Proc_Info[2]; iproc++) {
      int    ifound = 0;
      size_t offset = 0;
      int    ilocal;
      for (int cnt = 0; cnt < globals.Num_Elem_Blk; cnt++) {
        for (ilocal = ifound; ilocal < globals.Proc_Num_Elem_Blk[iproc]; ilocal++) {
          if (globals.Proc_Elem_Blk_Ids[iproc][ilocal] == eb_ids_global[cnt])
            break;
        }

        if (ilocal < globals.Proc_Num_Elem_Blk[iproc]) {
          eb_map_ptr[iproc][cnt] = &globals.GElems[iproc][offset];
          eb_cnts_local[iproc][cnt] = globals.Proc_Num_Elem_In_Blk[iproc][ilocal];
          offset += globals.Proc_Num_Elem_In_Blk[iproc][ilocal];
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
    Restart_Info.Node_Vals = (T **) array_alloc (__FILE__, __LINE__,
                                                 1,Proc_Info[2],
                                                 sizeof(T *));

    /* now allocate storage for the values */
    for (int iproc = 0; iproc <Proc_Info[2]; iproc++) {
      size_t array_size = Restart_Info.NVar_Node * (globals.Num_Internal_Nodes[iproc] +
					     globals.Num_Border_Nodes[iproc] + globals.Num_External_Nodes[iproc]);
      Restart_Info.Node_Vals[iproc] = (T*)array_alloc (__FILE__, __LINE__, 1,
			     array_size, cpu_ws);

    }
  }

  if (Restart_Info.NVar_Sset > 0 ) {

    /* allocate storage space */
    Restart_Info.Sset_Vals = (T **) array_alloc (__FILE__, __LINE__,
						  1,Proc_Info[2],
						  sizeof(T *));

    /* now allocate storage for the values */
    for (int iproc = 0; iproc <Proc_Info[2]; iproc++) {
      size_t array_size = Restart_Info.NVar_Sset * globals.Proc_SS_Elem_List_Length[iproc];
                   
      Restart_Info.Sset_Vals[iproc] = (T*)array_alloc (__FILE__, __LINE__, 1,
			     array_size , cpu_ws);
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

    /* Get the Sideset IDs from the input file */
    if (ex_get_ids (exoid, EX_SIDE_SET, TOPTR(ss_ids_global)) < 0) {
      fprintf(stderr, "%s: unable to get sideset IDs", yo);
      exit(1);
    }

    /* Get the count of elements in each sideset */
    for (int cnt = 0; cnt < globals.Num_Side_Set; cnt++) {
      if (ex_get_set_param(exoid, EX_SIDE_SET,
			   ss_ids_global[cnt],
			   &(ss_cnts_global[cnt]), NULL) < 0) {
	fprintf(stderr, "%s: unable to get element count for sideset id %lu",
		yo, (size_t)ss_ids_global[cnt]);
	exit(1);
      }
    }
  } /* End: "if (Restart_Info.NVar_Sset > 0 )" */


  if (Restart_Info.NVar_Nset > 0 ) {

    /* allocate storage space */
    Restart_Info.Nset_Vals = (T **) array_alloc (__FILE__, __LINE__,
						  1,Proc_Info[2],
						  sizeof(T *));

    /* now allocate storage for the values */
    for (int iproc = 0; iproc <Proc_Info[2]; iproc++) {
      size_t array_size = Restart_Info.NVar_Nset * globals.Proc_NS_List_Length[iproc];
                   
      Restart_Info.Nset_Vals[iproc] = (T*)array_alloc (__FILE__, __LINE__, 1,
			     array_size , cpu_ws);

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

    /* Get the Nodeset IDs from the input file */
    if (ex_get_ids (exoid, EX_NODE_SET, TOPTR(ns_ids_global)) < 0) {
      fprintf(stderr, "%s: unable to get nodeset IDs", yo);
      exit(1);
    }

    /* Get the count of elements in each nodeset */
    for (int cnt = 0; cnt < globals.Num_Node_Set; cnt++) {
      if (ex_get_set_param(exoid, EX_NODE_SET,
			   ns_ids_global[cnt],
			   &(ns_cnts_global[cnt]), NULL) < 0) {
	fprintf(stderr, "%s: unable to get element count for nodeset id %lu",
		yo, (size_t)ns_ids_global[cnt]);
	exit(1);
      }
    }
  } /* End: "if (Restart_Info.NVar_Nset > 0 )" */


  /*
   * NOTE: A possible place to speed this up would be to
   * get the global node and element lists here, and broadcast
   * them out only once.
   */

  /* figure out how many blocks need to be looped over */
  num_blocks = Restart_Info.Num_Times;
  
  par_exoid = (int*)malloc(Proc_Info[2] * sizeof(int));
  if(!par_exoid) {
    fprintf(stderr, "[%s]: ERROR, insufficient memory!\n",
	    yo);
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
  
  int open_file_count = get_free_descriptor_count();
  if (open_file_count >Proc_Info[5]) {
    fprintf(stderr, "All output files opened simultaneously.\n");
    for (int iproc=Proc_Info[4]; iproc <Proc_Info[4]+Proc_Info[5]; iproc++) {
     
      gen_par_filename(cTemp, Par_Nem_File_Name, Proc_Ids[iproc],
		       Proc_Info[0]);
      
      /* Open the parallel Exodus II file for writing */
      cpu_ws = io_ws;
      int mode = EX_WRITE | int64api | int64db;
      if ((par_exoid[iproc]=ex_open(Par_Nem_File_Name, mode, &cpu_ws,
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
    times_in_blk = 1;

    for (int cnt = 0; cnt < times_in_blk; cnt++) {
      time_idx = iblk; 
      double start_t = second ();

      /* read and distribute the variables for this time step */
      if (read_vars(exoid, Restart_Info.Time_Idx[iblk],
		    TOPTR(eb_ids_global), TOPTR(eb_cnts_global), eb_map_ptr,
		    eb_cnts_local,
		    TOPTR(ss_ids_global), TOPTR(ss_cnts_global),
		    TOPTR(ns_ids_global), TOPTR(ns_cnts_global)) < 0) {
	fprintf(stderr, "%s: Error occured while reading variables\n",
		yo);
	exit(1);
      }
      double end_t   = second () - start_t;
      printf ("\tTime to read  vars for timestep %d: %f (sec.)\n", (time_idx+1), end_t);

      start_t = second ();
      for (int iproc=Proc_Info[4]; iproc <Proc_Info[4]+Proc_Info[5]; iproc++) {

	if (open_file_count <Proc_Info[5]) {
	  gen_par_filename(cTemp, Par_Nem_File_Name, Proc_Ids[iproc],
			   Proc_Info[0]);
	  
	  /* Open the parallel Exodus II file for writing */
	  cpu_ws = io_ws;
	  int mode = EX_WRITE | int64api | int64db;
	  if ((par_exoid[iproc]=ex_open(Par_Nem_File_Name, mode, &cpu_ws,
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
	time_idx = iblk;
	write_var_timestep(par_exoid[iproc], iproc, (time_idx+1),
			   TOPTR(eb_ids_global), TOPTR(ss_ids_global), TOPTR(ns_ids_global));

	if (iproc%10 == 0 || iproc ==Proc_Info[2]-1)
	  fprintf(stderr, "%d", iproc);
	else
	  fprintf(stderr, ".");

	if (open_file_count <Proc_Info[5]) {
	  if (ex_close(par_exoid[iproc]) == -1) {
	    fprintf(stderr, "[%d] %s Could not close the parallel Exodus II file.\n",
		    iproc, yo);
	    exit(1);
	  }
	}
      } /* End "for (iproc=0; iproc <Proc_Info[2]; iproc++)" */

      end_t   = second () - start_t;
      printf ("\n\tTime to write vars for timestep %d: %f (sec.)\n", (time_idx+1), end_t);

    } /* End: "for (iblk = 0; iblk < num_blocks; iblk++)" */
  }
  if (Restart_Info.NVar_Elem > 0 ) {
    safe_free((void **) &eb_map_ptr);
    safe_free((void **) &eb_cnts_local);
  }

  /* Close the restart exodus II file */
  if (ex_close(exoid) == -1) {
    fprintf(stderr, "%sCould not close the restart Exodus II file\n",
	    yo);
    exit(1);
  }

  if (open_file_count >Proc_Info[5]) {
    for (int iproc=Proc_Info[4]; iproc <Proc_Info[4]+Proc_Info[5]; iproc++) {
      /* Close the parallel exodus II file */
      if (ex_close(par_exoid[iproc]) == -1) {
	fprintf(stderr, "[%d] %s Could not close the parallel Exodus II file.\n",
		iproc, yo);
	exit(1);
      }
    }
  }
  if (par_exoid != NULL) {
    free(par_exoid);
    par_exoid = NULL;
  }
}

template <typename T, typename INT>
int NemSpread<T,INT>::read_var_param (int exoid, int max_name_length)
{
  const char  *yo="read_var_param";

  /* Get the number of time indices contained in the file */
  int ret_int = ex_inquire_int(exoid, EX_INQ_TIME);

  /* see if the user want to get all of the time indices */
  if (Restart_Info.Num_Times == -1) {

    Restart_Info.Num_Times = ret_int;

    if (ret_int > 0) {
      /* allocate array space */
      Restart_Info.Time_Idx = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                  ret_int, sizeof(int));

      for (int cnt = 0; cnt < Restart_Info.Num_Times; cnt++)
        Restart_Info.Time_Idx[cnt] = cnt + 1;
    }
  }
  else {
    /* Check to see if the requested indeces are valid */
    for (int cnt = 0; cnt < Restart_Info.Num_Times; cnt++) {

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
  if (ex_get_variable_param(exoid, EX_GLOBAL, &(Restart_Info.NVar_Glob)) < 0) {
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
    if (ex_get_variable_names(exoid, EX_GLOBAL, Restart_Info.NVar_Glob,
                         Restart_Info.GV_Name) < 0) {
      fprintf(stderr, "%s: Could not get global variable names from file\n",
              yo);
      return -1;
    }
  }

  /***************** Elemental Variables ********************/
  if (ex_get_variable_param(exoid, EX_ELEM_BLOCK, &(Restart_Info.NVar_Elem)) < 0) {
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
    if (ex_get_variable_names(exoid, EX_ELEM_BLOCK, Restart_Info.NVar_Elem,
                         Restart_Info.EV_Name) < 0) {
      fprintf(stderr, "%s: Could not get elemental variable names from file\n",
              yo);
      return -1;
    }

    /* and get the truth table */
    Restart_Info.GElem_TT = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                (globals.Num_Elem_Blk
						 * Restart_Info.NVar_Elem),
                                                sizeof(int));

    check_exodus_error(ex_get_truth_table(exoid, EX_ELEM_BLOCK, 
					  globals.Num_Elem_Blk,
					  Restart_Info.NVar_Elem,
					  Restart_Info.GElem_TT),
                       "ex_get_truth_table");
  }

  /******************* Nodal Variables **********************/
  if (ex_get_variable_param(exoid, EX_NODAL, &(Restart_Info.NVar_Node)) < 0) {
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
    if (ex_get_variable_names(exoid, EX_NODAL, Restart_Info.NVar_Node,
                         Restart_Info.NV_Name) < 0) {
      fprintf(stderr, "%s: Could not get nodal variable names from file\n",
              yo);
      return -1;
    }
  }

  /******************* Sideset Variables **********************/
  if (ex_get_variable_param(exoid, EX_SIDE_SET, &(Restart_Info.NVar_Sset)) < 0) {
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
    if (ex_get_variable_names(exoid, EX_SIDE_SET, Restart_Info.NVar_Sset,
                         Restart_Info.SSV_Name) < 0) {
      fprintf(stderr, "%s: Could not get sideset variable names from file\n",
              yo);
      return -1;
    }

    /* and get the truth table */
    Restart_Info.GSset_TT = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                (globals.Num_Side_Set
						 * Restart_Info.NVar_Sset),
                                                sizeof(int));

    check_exodus_error(ex_get_truth_table(exoid, EX_SIDE_SET,
					  globals.Num_Side_Set,
					  Restart_Info.NVar_Sset,
					  Restart_Info.GSset_TT),
                       "ex_get_truth_table");
  }

  /******************* Nodeset Variables **********************/
  if (ex_get_variable_param(exoid, EX_NODE_SET, &(Restart_Info.NVar_Nset)) < 0) {
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
    if (ex_get_variable_names(exoid, EX_NODE_SET, Restart_Info.NVar_Nset,
                         Restart_Info.NSV_Name) < 0) {
      fprintf(stderr, "%s: Could not get nodeset variable names from file\n",
              yo);
      return -1;
    }

    /* and get the truth table */
    Restart_Info.GNset_TT = (int *) array_alloc(__FILE__, __LINE__, 1,
                                                (globals.Num_Node_Set
						 * Restart_Info.NVar_Nset),
                                                sizeof(int));

    check_exodus_error(ex_get_truth_table(exoid, EX_NODE_SET,
					  globals.Num_Node_Set,
					  Restart_Info.NVar_Nset,
					  Restart_Info.GNset_TT),
                       "ex_get_var_tab");
  }


#ifdef DEBUG
  if (Debug_Flag >= 2) {
    printf("\n\nRestart Parameters:\n");
    printf("\tNumber of time indices: %d\n", Restart_Info.Num_Times);
    for (int cnt = 0; cnt < Restart_Info.Num_Times; cnt++)
      printf("\t\tTime index: %d\n", Restart_Info.Time_Idx[cnt]);
    printf("\tNumber of global variables: %d\n", Restart_Info.NVar_Glob);
    for (int cnt = 0; cnt < Restart_Info.NVar_Glob; cnt++)
      printf("\t\tGlobal variable %d: %s\n", (cnt+1),
             Restart_Info.GV_Name[cnt]);
    printf("\tNumber of elental variables: %d\n", Restart_Info.NVar_Elem);
    for (int cnt = 0; cnt < Restart_Info.NVar_Elem; cnt++)
      printf("\t\tElemental variable %d: %s\n", (cnt+1),
             Restart_Info.EV_Name[cnt]);
    printf("\tNumber of nodal variables: %d\n", Restart_Info.NVar_Node);
    for (int cnt = 0; cnt < Restart_Info.NVar_Node; cnt++)
      printf("\t\tNodal variable %d: %s\n", (cnt+1), Restart_Info.NV_Name[cnt]);
  }
#endif

  return 0;

}

template <typename T, typename INT>
int NemSpread<T,INT>::broadcast_var_param(Restart_Description<T> *restart, int max_name_length)
{
  if (restart->NVar_Elem > 0) {
    Restart_Info.Elem_TT = (int **) array_alloc(__FILE__, __LINE__, 2,
                                               Proc_Info[2],
                                                 (globals.Num_Elem_Blk
                                                 * Restart_Info.NVar_Elem),
                                                 sizeof(int));

    /* and copy it for each proc that this processor is responsible for */
    for (int iproc = 0; iproc <Proc_Info[2]; iproc++)
      for (int cnt1 = 0; cnt1 < globals.Num_Elem_Blk; cnt1++)
        for (int cnt2 = 0; cnt2 < Restart_Info.NVar_Elem; cnt2++)
          Restart_Info.Elem_TT[iproc][cnt1*Restart_Info.NVar_Elem+cnt2] =
            Restart_Info.GElem_TT[cnt1*Restart_Info.NVar_Elem+cnt2];
  }

  if (restart->NVar_Nset > 0) {
    Restart_Info.Nset_TT = (int **) array_alloc(__FILE__, __LINE__, 2,
                                               Proc_Info[2],
						(globals.Num_Node_Set * Restart_Info.NVar_Nset),
                                                 sizeof(int));

    /* and copy it for each proc that this processor is responsible for */
    for (int iproc = 0; iproc <Proc_Info[2]; iproc++)
      for (int cnt1 = 0; cnt1 < globals.Num_Node_Set; cnt1++)
        for (int cnt2 = 0; cnt2 < Restart_Info.NVar_Nset; cnt2++)
          Restart_Info.Nset_TT[iproc][cnt1*Restart_Info.NVar_Nset+cnt2] =
            Restart_Info.GNset_TT[cnt1*Restart_Info.NVar_Nset+cnt2];
  }

  if (restart->NVar_Sset > 0) {
    Restart_Info.Sset_TT = (int **) array_alloc(__FILE__, __LINE__, 2,
                                               Proc_Info[2],
						(globals.Num_Side_Set * Restart_Info.NVar_Sset),
                                                 sizeof(int));

    /* and copy it for each proc that this processor is responsible for */
    for (int iproc = 0; iproc <Proc_Info[2]; iproc++)
      for (int cnt1 = 0; cnt1 < globals.Num_Side_Set; cnt1++)
        for (int cnt2 = 0; cnt2 < Restart_Info.NVar_Sset; cnt2++)
          Restart_Info.Sset_TT[iproc][cnt1*Restart_Info.NVar_Sset+cnt2] =
            Restart_Info.GSset_TT[cnt1*Restart_Info.NVar_Sset+cnt2];
  }
  return 0;

}


template <typename T, typename INT>
int NemSpread<T,INT>::read_vars(int exoid, int index, INT *eb_ids,
				INT *eb_cnts, INT ***eb_map_ptr, INT **eb_cnts_local,
				INT *ss_ids, INT *ss_cnts, INT *ns_ids, INT *ns_cnts)
{
  const char  *yo="read_vars";

  T *ptr = &Restart_Info.Time;
  /* first read the time */
  if (ex_get_time(exoid, index, ptr) < 0) {
    fprintf(stderr, "%s: ERROR, unable to get time for restart index %d!\n",
	    yo, index);
    return -1;
  }

  /***************** Global Variables ********************/
  /* allocate space for the global variables */
  if (Restart_Info.NVar_Glob > 0) {
    ptr = &Restart_Info.Glob_Vals[0];

    /* get the global variables */
    if (ex_get_glob_vars(exoid, index, Restart_Info.NVar_Glob, ptr) < 0) {
      fprintf(stderr, "%s: Could not get global variables from file\n", yo);
      return -1;
    }
  }

  if (Restart_Info.NVar_Elem > 0 ) {
    printf("Reading %d element variables...\n", Restart_Info.NVar_Elem);
    if (read_elem_vars(exoid, index, eb_ids, eb_cnts, eb_map_ptr,
                       eb_cnts_local) < 0) {
      fprintf(stderr, "%s: Error distributing elemental variables.\n", yo);
      return -1;
    }
  }

  if (Restart_Info.NVar_Node > 0 ) {
    printf("Reading %d nodal variables...\n", Restart_Info.NVar_Node);
    if (read_nodal_vars(exoid, index) < 0) {
      fprintf(stderr, "%s: Error distributing nodal variables.\n", yo);
      return -1;
    }
  }

  if (Restart_Info.NVar_Sset > 0 ) {
    printf("Reading %d sideset variables...\n", Restart_Info.NVar_Sset);
    if (read_sset_vars(exoid, index, ss_ids, ss_cnts) < 0) {
      fprintf(stderr, "%s: Error distributing sideset variables.\n", yo);
      return -1;
    }
  }


  if (Restart_Info.NVar_Nset > 0 ) {
    printf("Reading %d nodeset variables...\n", Restart_Info.NVar_Nset);
    if (read_nset_vars(exoid, index, ns_ids, ns_cnts) < 0) {
      fprintf(stderr, "%s: Error distributing nodeset variables.\n", yo);
      return -1;
    }
  }

  return 0;

}

template <typename T, typename INT>
int NemSpread<T,INT>::read_elem_vars(int exoid, int index, INT *eb_ids,
				     INT *eb_cnts, INT ***eb_map_ptr,
				     INT **eb_cnts_local)
{

  /* to speed up searches, keep track of element blocks offset on each proc */
  INT *local_offset = (INT *) array_alloc (__FILE__, __LINE__, 1,Proc_Info[2],
                                      sizeof(INT));
  for (int iproc = 0; iproc <Proc_Info[2]; iproc++)
    local_offset[iproc] = 0;

  /* loop over the number of element blocks */
  INT     eb_offset=0;
  for (int iblk = 0; iblk < globals.Num_Elem_Blk; iblk++) {
    read_elem_vars_1(exoid, index, eb_ids,
		     eb_cnts, eb_map_ptr, eb_cnts_local,
		     iblk, eb_offset, local_offset);

    /* need to keep track of this for the element number map */
    eb_offset += eb_cnts[iblk];

    /* need to set up local offsets for next block */
    for (int iproc = 0; iproc <Proc_Info[2]; iproc++)
      local_offset[iproc] += eb_cnts_local[iproc][iblk];

  } /* End "for (iblk = 0; iblk < globals.Num_Elem_Blk; iblk++)" */

  safe_free((void **) &local_offset);
  return 0;
}

template <typename T, typename INT>
int NemSpread<T,INT>::read_elem_vars_1(int exoid, int index, INT *eb_ids,
				       INT *eb_cnts, INT ***eb_map_ptr,
				       INT **eb_cnts_local, int iblk,
				       int eb_offset, INT *local_offset)
{
  /* Allocate memory for temporary storage */
  T *vals  = (T*)array_alloc(__FILE__, __LINE__, 1, eb_cnts[iblk], sizeof(T));
  
  /* now loop over each variable */
  for (int ivar = 0; ivar < Restart_Info.NVar_Elem; ivar++) {

    /* check if this variable exists for this element block */
    if (Restart_Info.GElem_TT[iblk*Restart_Info.NVar_Elem+ivar]) {

      /*
       * Read in the specified element variable values and their associated
       * global FEM element numbers.
       */

      check_exodus_error(ex_get_var(exoid,
				    index,
				    EX_ELEM_BLOCK,
				    (ivar+1),
				    eb_ids[iblk],
				    eb_cnts[iblk],
				    vals),
			 "ex_get_var");
      
      /*
       * Find out which FEM elements belong on this processor and copy
       * them to the restart vector.
       */
      for (int iproc = 0; iproc <Proc_Info[2]; iproc++) {
	  
	/* check to see if this element block needs this variable */
	if (Restart_Info.Elem_TT[iproc][iblk*Restart_Info.NVar_Elem+ivar]) {
	    
	  /* calculate the offset for this variable */
	  size_t var_offset = ivar * (globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc]);
	    
	  INT *elem_map = eb_map_ptr[iproc][iblk];
	  size_t num_elem = eb_cnts_local[iproc][iblk];
	    
	  for (size_t i1 = 0; i1 < num_elem; i1++) {
	    size_t elem_loc = var_offset + i1 + local_offset[iproc];
		
	    Restart_Info.Elem_Vals[iproc][elem_loc] = vals[elem_map[i1]-eb_offset];
	  }
	}
      }
    } /* End "if (Restart_Info.GElem_TT[...])" */
  } 
  safe_free((void **) &vals);
  return 0;
}
    
template <typename T, typename INT>
int NemSpread<T,INT>::read_sset_vars(int exoid, int index, INT *ss_ids, INT *ss_cnts)
{
  /* loop over the number of side sets */
  for (int iset = 0; iset < globals.Num_Side_Set; iset++) {
    read_sset_vars_1(exoid, index, ss_ids, ss_cnts, iset);
  }
  return 0;
}

template <typename T, typename INT>
int NemSpread<T,INT>::read_sset_vars_1(int exoid, int index, INT *ss_ids,
				       INT *ss_cnts, int iset)
{
  /* Allocate memory for temporary storage */
  T *vals  = (T*)array_alloc(__FILE__, __LINE__, 1, ss_cnts[iset], sizeof(T));
  
  /* now loop over each variable */
  for (int ivar = 0; ivar < Restart_Info.NVar_Sset; ivar++) {

    /* check if this variable exists for this set */
    if (Restart_Info.GSset_TT[iset*Restart_Info.NVar_Sset+ivar]) {

      /* Read in the specified variable values */
      check_exodus_error(ex_get_var(exoid, index, EX_SIDE_SET,
				    (ivar+1), ss_ids[iset], ss_cnts[iset],
				    vals), "ex_get_var");
      
      for (int iproc = 0; iproc <Proc_Info[2]; iproc++) {
	size_t ss_offset = 0;
	size_t var_offset = ivar * globals.Proc_SS_Elem_List_Length[iproc];
	for (int i = 0; i < globals.Proc_Num_Side_Sets[iproc]; i++) {
	  if (globals.Proc_SS_Ids[iproc][i] == ss_ids[iset]) {
	
	    size_t num_elem = globals.Proc_SS_Elem_Count[iproc][i];
	    for (size_t i1 = 0; i1 < num_elem; i1++) {
	      INT gelem_loc = globals.Proc_SS_GEMap_List[iproc][i1+ss_offset];
	      assert(gelem_loc < ss_cnts[iset]);
	      Restart_Info.Sset_Vals[iproc][i1+ss_offset+var_offset] = vals[gelem_loc];
	    }
	    break;
	  }
	  ss_offset += globals.Proc_SS_Elem_Count[iproc][i];
	}
      } 
    } 
  }
  safe_free((void **) &vals);
  return 0;
}
    
template <typename T, typename INT>
int NemSpread<T,INT>::read_nset_vars(int exoid, int index, INT *ns_ids, INT *ns_cnts)
{
  /* loop over the number of node sets */
  for (int iset = 0; iset < globals.Num_Node_Set; iset++) {
    read_nset_vars_1(exoid, index, ns_ids, ns_cnts, iset);
  }
  return 0;
}

template <typename T, typename INT>
int NemSpread<T,INT>::read_nset_vars_1(int exoid, int index, INT *ns_ids,
			    INT *ns_cnts, int iset)
{
  /* Allocate memory for temporary storage */
  T *vals  = (T*)array_alloc(__FILE__, __LINE__, 1, ns_cnts[iset], sizeof(T));
  
  /* now loop over each variable */
  for (int ivar = 0; ivar < Restart_Info.NVar_Nset; ivar++) {

    /* check if this variable exists for this set */
    if (Restart_Info.GNset_TT[iset*Restart_Info.NVar_Nset+ivar]) {

      /* Read in the specified variable values */
      check_exodus_error(ex_get_var(exoid, index, EX_NODE_SET, (ivar+1), ns_ids[iset], ns_cnts[iset],
					 vals), "ex_get_nset_var");
      
      for (int iproc = 0; iproc <Proc_Info[2]; iproc++) {
	size_t ns_offset = 0;
	size_t var_offset = ivar * globals.Proc_NS_List_Length[iproc];
	for (int i = 0; i < globals.Proc_Num_Node_Sets[iproc]; i++) {
	  if (globals.Proc_NS_Ids[iproc][i] == ns_ids[iset]) {
	
	    size_t num_elem = globals.Proc_NS_Count[iproc][i];
	    for (size_t i1 = 0; i1 < num_elem; i1++) {
	      INT gelem_loc = globals.Proc_NS_GNMap_List[iproc][i1+ns_offset];
	      assert(gelem_loc < ns_cnts[iset]);
	      Restart_Info.Nset_Vals[iproc][i1+ns_offset+var_offset] = vals[gelem_loc];
	    }
	    break;
	  }
	  ns_offset += globals.Proc_NS_Count[iproc][i];
	}
      } 
    } 
  }
  safe_free((void **) &vals);
  return 0;
}
    
template <typename T, typename INT>
int NemSpread<T,INT>::read_nodal_vars(int exoid, int index)
{
  /* Allocate memory for temporary storage */
  T *vals = (T*)array_alloc(__FILE__, __LINE__, 1, globals.Num_Node, sizeof(T));

  /* Loop over each auxiliary variable */
  for (int var_num = 0; var_num < Restart_Info.NVar_Node; var_num++) {
    /*
     * Read in the specified nodal variable values and their associated
     * global FEM node numbers.
     */
    check_exodus_error(ex_get_var(exoid, index, EX_NODAL, (var_num+1),
				  1, globals.Num_Node, vals),
		       "ex_get_var");
    
    /*
     * Find out which FEM nodes belong on this processor and copy
     * them to the restart vector.
     */
    for (int iproc = 0; iproc <Proc_Info[2]; iproc++) {

      /* calculate the offset for this variable */
      size_t loc_count = globals.Num_Internal_Nodes[iproc]+globals.Num_Border_Nodes[iproc]+
	globals.Num_External_Nodes[iproc];

      size_t var_offset = var_num * loc_count;
	
      for (size_t i2 = 0; i2 < loc_count; i2++) {
	size_t node_loc = var_offset + i2;
	Restart_Info.Node_Vals[iproc][node_loc] =
	  vals[globals.GNodes[iproc][i2]-1];
      }
    }

  } /* End "for (var_num = 0; var_num < Restart_Info.NVar_Node; var_num++)" */

  safe_free((void **) &vals);
  return 0;
}

template <typename T, typename INT>
int NemSpread<T,INT>::compare_mesh_param(int exoid)
{
  int     ret = 1;
  
  ex_init_params info;
  info.title[0] = '\0';
  int error = ex_get_init_ext (exoid, &info);
  check_exodus_error (error, "ex_get_init");

  /* now check that the parameters match those retrieved from the mesh file */
  if (info.num_dim                    != globals.Num_Dim)      ret = 0;
  else if ((size_t)info.num_nodes     != globals.Num_Node)     ret = 0;
  else if ((size_t)info.num_elem      != globals.Num_Elem)     ret = 0;
  else if (info.num_elem_blk  != globals.Num_Elem_Blk) ret = 0;
  else if (info.num_node_sets != globals.Num_Node_Set) ret = 0;
  else if (info.num_side_sets != globals.Num_Side_Set) ret = 0;

  return (ret);
}

  /*****************************************************************************/
namespace {
  template <typename INT>
  size_t find_gnode_inter(INT *intersect, size_t num_g_nodes, INT *glob_vec,
			  size_t num_int_nodes, size_t num_bor_nodes,
			  size_t num_ext_nodes, INT *loc_vec)

  /*
   * This function assumes that glob_vec is monotonic and that loc_vec is
   * monotonic for each of the internal, border and external node IDs it
   * contains.
   */
  {
    size_t count=0;

    /* Initialize the intersect vector */
    for (size_t i1=0; i1 < num_g_nodes; i1++)
      intersect[i1] = -1;

    /* Check for the possibility of an intersection */
    size_t min_set1 = glob_vec[0];
    size_t max_set1 = glob_vec[num_g_nodes-1];

    /* Search through the internal nodes */
    if (num_int_nodes > 0) {
      size_t min_set2 = loc_vec[0];
      size_t max_set2 = loc_vec[num_int_nodes-1];

      if ( (max_set2 >= min_set1) && (min_set2 <= max_set1) ) {
	for (size_t i1=0, i2=0; i1 < num_g_nodes; i1++) {
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
      size_t min_set2 = loc_vec[num_int_nodes];
      size_t max_set2 = loc_vec[num_int_nodes+num_bor_nodes-1];

      size_t offset = num_int_nodes;

      if ( (max_set2 >= min_set1) && (min_set2 <= max_set1) ) {
	for (size_t i1=0, i2=0; i1 < num_g_nodes; i1++) {
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
      size_t min_set2 = loc_vec[num_int_nodes+num_bor_nodes];
      size_t max_set2 = loc_vec[num_int_nodes+num_bor_nodes+num_ext_nodes-1];

      size_t offset = num_int_nodes + num_bor_nodes;

      if ( (max_set2 >= min_set1) && (min_set2 <= max_set1) ) {
	for (size_t i1=0, i2=0; i1 < num_g_nodes; i1++) {
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
}
