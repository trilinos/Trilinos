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
#include <assert.h>
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for printf, fprintf, stderr, etc
#include <stdlib.h>                     // for exit, malloc

#include "exodusII.h"                   // for ex_inquire, ex_opts, etc
#include "globals.h"                    // for ELEM_COMM_MAP, etc
#include "nem_spread.h"                 // for NemSpread, etc
#include "pe_common.h"                  // for PEX_MAX
#include "rf_allo.h"                    // for array_alloc, safe_free
#include "rf_io_const.h"                // for Debug_Flag, Exo_LB_File
#include "rf_util.h"                    // for print_line
#include "sort_utils.h"                 // for gds_qsort

char **qa_record_ptr, **inf_record_ptr;
int num_inf_rec=0, num_qa_rec=0, length_qa=0;

/****************************************************************************
 ****************************************************************************
 ****************************************************************************
 * This function reads load-balance information for salsa on a parallel
 * computer.
 *
 * Author: Gary L. Hennigan (1421)
 *---------------------------------------------------------------------------
 * Revision History.
 *
 *   Gary Hennigan:
 *      08 November 1993 - Modified scalar version to write information out
 *                         to the parallel disk array once processed.
 ****************************************************************************
 ****************************************************************************
 ****************************************************************************/
template void NemSpread<float,int>::load_lb_info(void);
template void NemSpread<double,int>::load_lb_info(void);
template void NemSpread<float,int64_t>::load_lb_info(void);
template void NemSpread<double,int64_t>::load_lb_info(void);

template <typename T, typename INT>
void NemSpread<T,INT>::load_lb_info(void)

/*
 *       load_lb_info:
 *
 *            This function reads and distributes the load balance information.
 *    This information is defined as the following scalars, which are specific
 *    to each processor:
 *
 *        Num_Internal_Nodes
 *        Num_Border_Nodes
 *        Num_External_Nodes
 *        globals.Num_Internal_Elems
 *        globals.Num_Border_Elems
 *
 *    and the following vectors, also specific to each processor:
 *
 *        globals.GNodes [Num_Internal_Nodes+Num_Border_Nodes+Num_External_Nodes]
 *        globals.GElems [globals.Num_Internal_Elems+globals.Num_Border_Elems]
 *
 *    For the 1 processor case, this routine fills in appropriate values
 *    for these quantities.
 */
{

  int    lb_exoid=0;
  INT    cmap_max_size=0, *comm_vec;
  const char  *yo = "load_lb_info";
  char   Title[MAX_LINE_LENGTH+1];
  float  version;

  int    cpu_ws=0;

  /******************************** START EXECUTION ****************************/
  if(Debug_Flag)
    printf ("\nStart to read in and distribute the load balance info\n");

  /* Open the Load Balance exoII file for reading */

  printf ("EXODUS II load-balance file: %s\n", Exo_LB_File);
  cpu_ws = io_ws;
  int mode = EX_READ | int64api;
  int iio_ws = 0; // Don't interfere with exodus files; this is the nemesis file.
  if((lb_exoid = ex_open(Exo_LB_File, mode, &cpu_ws,
			 &iio_ws, &version)) == -1) {
    fprintf(stderr, "%sERROR: Couldn\'t open lb file, %s\n", yo,
	    Exo_LB_File);
    exit(1);
  }

  /* Read information about the processor configuration */
  read_proc_init(lb_exoid,Proc_Info, &Proc_Ids);

  /* Allocate space for the counts */
  globals.Num_Internal_Nodes = (INT *)array_alloc(__FILE__, __LINE__, 1,
						  7*Proc_Info[2], sizeof(INT));
  globals.Num_Border_Nodes   = globals.Num_Internal_Nodes +Proc_Info[2];
  globals.Num_External_Nodes = globals.Num_Border_Nodes   +Proc_Info[2];
  globals.Num_Internal_Elems = globals.Num_External_Nodes +Proc_Info[2];
  globals.Num_Border_Elems   = globals.Num_Internal_Elems +Proc_Info[2];
  globals.Num_N_Comm_Maps    = globals.Num_Border_Elems   +Proc_Info[2];
  globals.Num_E_Comm_Maps    = globals.Num_N_Comm_Maps    +Proc_Info[2];

  /* Allocate space for each processor entity */
  globals.GNodes   = (INT **)array_alloc(__FILE__, __LINE__, 1, 3*Proc_Info[2],
					 sizeof(INT *));
  globals.GElems   = globals.GNodes +Proc_Info[2];
  globals.Elem_Map = globals.GElems +Proc_Info[2];

  /* Allocate contiguous space for the pointer vectors on all processors */
  INT *Int_Space     = (INT *)array_alloc(__FILE__, __LINE__, 1,
                                     (7*Proc_Info[0] + 1), sizeof(INT));
  INT *Int_Node_Num  = Int_Space     + 1;
  INT *Bor_Node_Num  = Int_Node_Num  +Proc_Info[0];
  INT *Ext_Node_Num  = Bor_Node_Num  +Proc_Info[0];
  INT *Int_Elem_Num  = Ext_Node_Num  +Proc_Info[0];
  INT *Bor_Elem_Num  = Int_Elem_Num  +Proc_Info[0];
  INT *Node_Comm_Num = Bor_Elem_Num  +Proc_Info[0];
  INT *Elem_Comm_Num = Node_Comm_Num +Proc_Info[0];

  /* Read the initial information contained in the load balance file */
  read_lb_init(lb_exoid, Int_Space, Int_Node_Num,
               Bor_Node_Num, Ext_Node_Num, Int_Elem_Num,
               Bor_Elem_Num, Node_Comm_Num, Elem_Comm_Num,
               Title);

  /* Allocate memory for the communication map arrays */
  globals.N_Comm_Map = (NODE_COMM_MAP<INT>**)malloc(Proc_Info[2]*sizeof(NODE_COMM_MAP<INT> *));
  globals.E_Comm_Map = (ELEM_COMM_MAP<INT>**)malloc(Proc_Info[2]*sizeof(ELEM_COMM_MAP<INT> *));
  if(!globals.N_Comm_Map || !globals.E_Comm_Map) {
    fprintf(stderr, "ERROR: Insufficient memory!\n");
    exit(1);
  }

  for (int iproc=0; iproc <Proc_Info[2]; iproc++) {

    /*
     * Error check:
     *	Currently a maximum of one nodal communication map and one
     * elemental communication map is supported.
     */
    if(globals.Num_N_Comm_Maps[iproc] > 1 || globals.Num_E_Comm_Maps[iproc] > 1) {
      fprintf(stderr, "%s: ERROR. Only 1 nodal and elemental comm map "
              "is supported\n", yo);
      exit(1);
    }
    else {

      /* Always allocate at least one and initialize the counts to 0 */
      globals.N_Comm_Map[iproc] = (NODE_COMM_MAP<INT> *)malloc(PEX_MAX(1,
								  globals.Num_N_Comm_Maps[iproc]) *
							  sizeof(NODE_COMM_MAP<INT>));
      if(globals.N_Comm_Map[iproc] == NULL && globals.Num_N_Comm_Maps[iproc] > 0) {
        fprintf(stderr,
                "%s: ERROR. Insufficient memory for nodal comm. map!\n",
                yo);
        exit(1);
      }

      for(size_t ijump=0; ijump < PEX_MAX(1, globals.Num_N_Comm_Maps[iproc]); ijump++)
        ((globals.N_Comm_Map[iproc])+ijump)->node_cnt = 0;

      globals.E_Comm_Map[iproc] = (ELEM_COMM_MAP<INT> *)malloc(PEX_MAX(1,
								  globals.Num_E_Comm_Maps[iproc]) *
							  sizeof(ELEM_COMM_MAP<INT>));
      if(globals.E_Comm_Map[iproc] == NULL && globals.Num_E_Comm_Maps[iproc] > 0) {
        fprintf(stderr,
                "%s: ERROR. Insufficient memory for elemental comm. map!\n",
                yo);
        exit(1);
      }

      for(size_t ijump=0; ijump < PEX_MAX(1, globals.Num_E_Comm_Maps[iproc]); ijump++)
        ((globals.E_Comm_Map[iproc])+ijump)->elem_cnt = 0;

    }

  } /* End "for (int iproc=0; iproc <Proc_Info[2]; iproc++)" */

  /* Set up each processor for the communication map parameters */
  read_cmap_params(lb_exoid, Node_Comm_Num, Elem_Comm_Num,
                   globals.Num_N_Comm_Maps, globals.Num_E_Comm_Maps,
                   globals.E_Comm_Map, globals.N_Comm_Map, &cmap_max_size,
                   &comm_vec);

  /* Allocate enough space to read the LB_data for one processor */
  INT *Integer_Vector = (INT *)array_alloc(__FILE__, __LINE__, 1,
                                      Int_Space[0] + cmap_max_size,
                                      sizeof (INT));

  /*
   * loop through the processors, one at a time, to read
   * their load balance information
   *
   * NOTE: From here on there are no provisions for multiple nodal
   *       or elemental communication maps.
   */

  size_t ijump = 0; /* keep track of where in comm_vec we are */
  for (int iproc = 0; iproc <Proc_Info[0]; iproc++) {

    /* Get the node map for processor "iproc" */
    if(ex_get_processor_node_maps(lb_exoid, &Integer_Vector[0],
				  &Integer_Vector[Int_Node_Num[iproc]],
				  &Integer_Vector[Int_Node_Num[iproc]+
						  Bor_Node_Num[iproc]],
				  iproc) < 0) {
      fprintf(stderr, "%s: ERROR, failed to get node map for Proc %d!\n",
	      yo, iproc);
      exit(1);
    }

    size_t vec_indx = Int_Node_Num[iproc] + Bor_Node_Num[iproc] +
      Ext_Node_Num[iproc];

    /* Get the element map for processor number "iproc" */
    if(ex_get_processor_elem_maps(lb_exoid, &Integer_Vector[vec_indx],
				  &Integer_Vector[vec_indx+Int_Elem_Num[iproc]],
				  iproc) < 0) {
      fprintf(stderr, "%s: ERROR, failed to get element map for Proc %d!\n",
	      yo, iproc);
      exit(1);
    }

    if(Node_Comm_Num[iproc] > 0) {
      vec_indx += Int_Elem_Num[iproc] + Bor_Elem_Num[iproc];

      if(ex_get_node_cmap(lb_exoid, comm_vec[ijump],
			  &Integer_Vector[vec_indx],
			  &Integer_Vector[vec_indx+comm_vec[ijump+1]],
			  iproc) < 0) {
	/*
	 * If there are disconnected mesh pieces, then it is
	 * possible that there is no comminication between the
	 * pieces and there will be no communication maps.  Normally
	 * this is a problem, so output a warning, but don't abort.
	 */
	fprintf(stderr,
		"%s: WARNING. Failed to get nodal comm map for Proc %d!\n",
		yo, iproc);
      }
    }

    if(Elem_Comm_Num[iproc] > 0) {
      vec_indx += 2*comm_vec[ijump+1];

      if(ex_get_elem_cmap(lb_exoid, comm_vec[ijump+2],
                          &Integer_Vector[vec_indx],
                          &Integer_Vector[vec_indx+comm_vec[ijump+3]],
                          &Integer_Vector[vec_indx+2*comm_vec[ijump+3]],
                          iproc) < 0) {
	fprintf(stderr,
		"%s: ERROR. Failed to get elemental comm map for Proc %d!\n",
		yo, iproc);
	exit(1);
      }
    }

    /*
     * Communicate load balance information to the correct processor
     *   - if iproc = Proc_Ids[*] then process the data instead.
     */
    assert(Proc_Ids[iproc] == iproc);
    process_lb_data (Integer_Vector, iproc);

    /*
     * now move ijump to the next communications map
     * make sure to check if there are any for this processor
     */
    if (Node_Comm_Num[iproc] > 0) ijump += 2;
    if (Elem_Comm_Num[iproc] > 0) ijump += 2;

  }

  /* Close the load balance file - we are finished with it */
  if(ex_close (lb_exoid) == -1) {
    fprintf (stderr, "%sERROR: Error in closing load balance file\n", yo);
    exit(1);
  }

  /************************* Cleanup and Printout Phase ***********************/

  /* Free temporary memory */
  safe_free((void **) &Integer_Vector);

  if(num_qa_rec > 0) {
    for(int i = 0; i < length_qa; i++) safe_free((void **) &(qa_record_ptr[i]));
    safe_free((void **) &qa_record_ptr);
  }

  if(num_inf_rec > 0) {
    for(int i = 0; i < num_inf_rec; i++) safe_free((void **) &(inf_record_ptr[i]));
    safe_free((void **) &inf_record_ptr);
  }

  safe_free((void **) &Int_Space);

  safe_free((void **) &comm_vec);

  /*========================================================================*/

  if(Debug_Flag)
    printf ("\nFinished distributing load balance info\n");

  /* Output Detailed timing information for the progam */
  /*
   * Print out a Large table of Load Balance Information if the debug_flag
   * setting is large enough
   */

  if(Debug_Flag >= 7) {
    printf ("\n\n");
    print_line ("=", 79);
    for (int iproc=0; iproc <Proc_Info[2]; iproc++) {
      printf("\n\t***For Processor %d***\n", Proc_Ids[iproc]);
      printf("\tInternal nodes owned by the current processor\n\t");
      for(INT i = 0; i < globals.Num_Internal_Nodes[iproc]; i++)
        printf(" %lu", (size_t)globals.GNodes[iproc][i]);

      printf("\n");

      printf("\tBorder nodes owned by the current processor\n\t");
      for(INT i = 0; i < globals.Num_Border_Nodes[iproc]; i++)
        printf(" %lu", (size_t)globals.GNodes[iproc][i + globals.Num_Internal_Nodes[iproc]]);

      printf("\n");

      if(globals.Num_External_Nodes[iproc] > 0) {
        printf("\tExternal nodes needed by the current processor\n\t");
        for(INT i = 0; i < globals.Num_External_Nodes[iproc]; i++)
          printf(" %lu", (size_t)globals.GNodes[iproc][i + globals.Num_Internal_Nodes[iproc] +
						       globals.Num_Border_Nodes[iproc]]);

        printf("\n");
      }

      printf("\tInternal elements owned by the current processor\n\t");
      for(INT i = 0; i < globals.Num_Internal_Elems[iproc]; i++)
        printf(" %lu", (size_t)globals.GElems[iproc][i]);

      printf("\n");

      if(globals.Num_Border_Elems[iproc] > 0) {
        printf("\tBorder elements owned by the current processor\n\t");
        for(INT i=0; i < globals.Num_Border_Elems[iproc]; i++)
          printf(" %lu", (size_t)globals.GElems[iproc][i+globals.Num_Internal_Elems[iproc]]);

        printf("\n");
      }

      if(globals.Num_N_Comm_Maps[iproc] > 0) {
        printf("\tNodal Comm Map for the current processor\n");
        printf("\t\tnode IDs:");
        for(size_t i=0; i < globals.N_Comm_Map[iproc]->node_cnt; i++)
          printf(" %lu", (size_t)globals.N_Comm_Map[iproc]->node_ids[i]);
        printf("\n\t\tproc IDs:");
        for(size_t i=0; i < globals.N_Comm_Map[iproc]->node_cnt; i++)
          printf(" %lu", (size_t)globals.N_Comm_Map[iproc]->proc_ids[i]);

        printf("\n");
      }

      if(globals.Num_E_Comm_Maps[iproc] > 0) {
        printf("\tElemental Comm Map for the current processor\n");
        printf("\t\telement IDs:");
        for(size_t i=0; i < globals.E_Comm_Map[iproc]->elem_cnt; i++)
          printf(" %lu", (size_t)globals.E_Comm_Map[iproc]->elem_ids[i]);
        printf("\n\t\tside IDs:");
        for(size_t i=0; i < globals.E_Comm_Map[iproc]->elem_cnt; i++)
          printf(" %lu", (size_t)globals.E_Comm_Map[iproc]->side_ids[i]);
        printf("\n\t\tproc IDs:");
        for(size_t i=0; i < globals.E_Comm_Map[iproc]->elem_cnt; i++)
          printf(" %lu", (size_t)globals.E_Comm_Map[iproc]->proc_ids[i]);

        printf("\n");
      }
    }
    printf("\n");
    print_line ("=", 79);
  }

} /* END of routine load_lb_info () ******************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T,INT>::process_lb_data (INT *Integer_Vector, int indx)

/*
*       Function which reads the load balance information from the processor's
*       communication buffer.
*
*       Author:          Scott Hutchinson (1421)
*       Date:            11 March 1993
*       Revised:         11 March 1993
*
* int *Integer_Vector;   Data structure for sending int data to processors
*               (ptr to vector of length, length_int)
*
*/
{
  /* Local variables */
  int icount = 0, itotal_nodes, itotal_elems, ig_count = 0;

  /* Function declarations */
  extern void sort_int  ( int n, int ra[] );

/***************************** execution begins ******************************/

  /* Calculate the length of the globals.GNodes and globals.GElems global variable */
  itotal_nodes = globals.Num_Internal_Nodes[indx] + globals.Num_Border_Nodes[indx] +
                 globals.Num_External_Nodes[indx];
  itotal_elems = globals.Num_Internal_Elems[indx] + globals.Num_Border_Elems[indx];

  /* Allocate Permament Arrays on the current processor */
  globals.GNodes[indx] = (INT *)
    array_alloc(__FILE__, __LINE__, 1,
               itotal_nodes + 2*itotal_elems + 2*(globals.N_Comm_Map[indx]->node_cnt) +
               3*(globals.E_Comm_Map[indx]->elem_cnt), sizeof(INT));

  globals.GElems[indx]   = globals.GNodes[indx] + itotal_nodes;
  globals.Elem_Map[indx] = globals.GElems[indx] + itotal_elems;
  globals.N_Comm_Map[indx]->node_ids = globals.Elem_Map[indx] + itotal_elems;
  globals.N_Comm_Map[indx]->proc_ids = globals.N_Comm_Map[indx]->node_ids +
    globals.N_Comm_Map[indx]->node_cnt;
  globals.E_Comm_Map[indx]->elem_ids = globals.N_Comm_Map[indx]->proc_ids +
    globals.N_Comm_Map[indx]->node_cnt;
  globals.E_Comm_Map[indx]->side_ids = globals.E_Comm_Map[indx]->elem_ids +
    globals.E_Comm_Map[indx]->elem_cnt;
  globals.E_Comm_Map[indx]->proc_ids = globals.E_Comm_Map[indx]->side_ids +
    globals.E_Comm_Map[indx]->elem_cnt;

  /*
   *            Extract the load balance information, and store it in
   *            permanent vectors:
   *              -   globals.GNodes[0]          - Internal_Nodes
   *                  globals.GNodes[+]          - Border_Nodes
   *                  globals.GNodes[+]          - External_Nodes
   *                  globals.GElems[0]          - Internal_Elems
   *                  globals.GElems[+]          - Border_Elems
   */

  for(INT i = 0; i < globals.Num_Internal_Nodes[indx]; i++)
    globals.GNodes[indx][ig_count++] = Integer_Vector[icount++];
  for(INT i = 0; i < globals.Num_Border_Nodes[indx]; i++)
    globals.GNodes[indx][ig_count++] = Integer_Vector[icount++];
  for(INT i = 0; i < globals.Num_External_Nodes[indx]; i++)
    globals.GNodes[indx][ig_count++] = Integer_Vector[icount++];

  ig_count = 0;
  for(INT i = 0; i < globals.Num_Internal_Elems[indx]; i++) {
    globals.GElems[indx][ig_count] = Integer_Vector[icount++];
    globals.Elem_Map[indx][ig_count] = globals.GElems[indx][ig_count];
    ig_count++;
  }
  for(INT i = 0; i < globals.Num_Border_Elems[indx]; i++) {
    globals.GElems[indx][ig_count] = Integer_Vector[icount++];
    globals.Elem_Map[indx][ig_count] = globals.GElems[indx][ig_count];
    ig_count++;
  }

  for(size_t i = 0; i < globals.N_Comm_Map[indx]->node_cnt; i++)
    (globals.N_Comm_Map[indx]->node_ids)[i] = Integer_Vector[icount++];
  for(size_t i = 0; i < globals.N_Comm_Map[indx]->node_cnt; i++)
    (globals.N_Comm_Map[indx]->proc_ids)[i] = Integer_Vector[icount++];
  for(size_t i = 0; i < globals.E_Comm_Map[indx]->elem_cnt; i++)
    (globals.E_Comm_Map[indx]->elem_ids)[i] = Integer_Vector[icount++];
  for(size_t i = 0; i < globals.E_Comm_Map[indx]->elem_cnt; i++)
    (globals.E_Comm_Map[indx]->side_ids)[i] = Integer_Vector[icount++];
  for(size_t i = 0; i < globals.E_Comm_Map[indx]->elem_cnt; i++)
    (globals.E_Comm_Map[indx]->proc_ids)[i] = Integer_Vector[icount++];

  /*
   * Sort the local element numbers in ascending global element numbers.
   * This means that globals.GElems will be monotonic.
   */
  gds_qsort(globals.GElems[indx],   globals.Num_Internal_Elems[indx]);
  gds_qsort(globals.Elem_Map[indx], globals.Num_Internal_Elems[indx]);

  /* Check that globals.GNodes is monotonic, from i = 0 to Num_Internal_Nodes */
#ifdef DEBUG
  
  assert(check_monot (globals.GNodes[indx], globals.Num_Internal_Nodes[indx]));

  /*
   * Check that globals.GNodes is monotonic, from i = Num_Internal_Nodes to
   *    (Num_Internal_Nodes + Num_Border_Nodes)
   */
  assert(check_monot(&(globals.GNodes[indx][globals.Num_Internal_Nodes[indx]]),
		     globals.Num_Border_Nodes[indx]));
#endif

} /* END of process_lb_data () ***********************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T,INT>::read_proc_init(int lb_exoid, int proc_info[], int **proc_ids_ptr)

/*----------------------------------------------------------------------------
 *  read_proc_init:
 *
 *	This function reads information about the processor configuration
 * which the load balance was generated for and makes assignments for each
 * processor.
 *----------------------------------------------------------------------------
 * Variable Glossary (after return):
 *
 *		proc_info[0] = # procs, from load balance file
 *		proc_info[1] = # procs for, from load balance file
 *		proc_info[2] = # procs this processor is responsible for
 *		proc_info[3] = # of extra procs
 *
 */
{
  const char *yo="read_proc_init";

  char  ftype[2];
  if(ex_get_init_info(lb_exoid, &proc_info[0], &proc_info[1], ftype) < 0) {
    fprintf(stderr, "%s: ERROR, could not get init info!\n",
	    yo);
    exit(1);
  }

  /* Calculate which processor is responsible for what */
  proc_info[2] = proc_info[0];
  int *proc_ids = (int *)array_alloc(__FILE__, __LINE__, 1, proc_info[2],
				     sizeof(int));

  for(int i1=0; i1 < proc_info[2]; i1++)
    proc_ids[i1] = i1;

  *proc_ids_ptr = proc_ids;

  return;

} /* End of read_proc_init() *************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T,INT>::read_lb_init(int lb_exoid,
				    INT *Int_Space,
				    INT *Int_Node_Num,
				    INT *Bor_Node_Num,
				    INT *Ext_Node_Num,
				    INT *Int_Elem_Num,
				    INT *Bor_Elem_Num,
				    INT *Node_Comm_Num,
				    INT *Elem_Comm_Num,
				    char *Title
				    )

/*
 *    read_lb_init:
 *
 *         This function reads the initial information contained in the
 * load balance file
 *
 *
 */
{
  const char *yo = "read_lb_init";

  /*********************BEGIN EXECUTABLE STATEMENTS****************************/

  /* Read Set-up information from the load-balance file on Proc 0 */

  /*
   * If debugging is not on go ahead and report errors from init. This
   * will show version mismatch information by default.
   */
  if(Debug_Flag == 0)
    ex_opts(EX_VERBOSE);

  /* Read the title of the LB File and about the size of the mesh */
  INT   num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;
  int error = ex_get_init_global(lb_exoid, &num_nodes, &num_elem,
				 &num_elem_blk, &num_node_sets, &num_side_sets);

  check_exodus_error (error, "ex_get_init");

  if(Debug_Flag == 0)
    ex_opts(!EX_VERBOSE);

#ifdef DEBUG
  if(Debug_Flag >= 2) {
    printf("---------------------------------------------------------\n");
    printf("\t\tLoad balance file global information\n");
    printf("---------------------------------------------------------\n");
    printf("\tNumber of nodes: %d\n", num_nodes);
    printf("\tNumber of elements: %d\n", num_elem);
    printf("\tNumber of element blocks: %d\n", num_elem_blk);
    printf("---------------------------------------------------------\n");
  }
#endif

  /* Cross-check the load balance file info against the mesh info */
  if(((size_t)num_nodes != globals.Num_Node) || ((size_t)num_elem != globals.Num_Elem) ||
     (num_elem_blk != globals.Num_Elem_Blk)) {
    fprintf(stderr,
	    "%s: ERROR: Problem dimensions in the LB File don't match with those \
in mesh file",yo);
    exit(1);
  }

  /* Read the QA Records */
  error = ex_inquire(lb_exoid, EX_INQ_QA, &num_qa_rec, NULL, NULL);
  check_exodus_error(error, "ex_inquire");
  if(num_qa_rec > 0) {
    length_qa = 4 * num_qa_rec;
    qa_record_ptr = (char **)array_alloc(__FILE__, __LINE__, 1, length_qa,
					 sizeof(char *));
    for(int i = 0; i < length_qa; i++) {
      qa_record_ptr[i] = (char *)
	array_alloc (__FILE__, __LINE__, 1, (MAX_STR_LENGTH+1),
		     sizeof(char));
    }
    error = ex_get_qa(lb_exoid, (char *(*)[4]) &qa_record_ptr[0]);
    check_exodus_error(error, "ex_get_qa");
  }

  /* Read the Info Records */
  error = ex_inquire(lb_exoid, EX_INQ_INFO, &num_inf_rec, NULL,  NULL);
  check_exodus_error(error, "ex_inquire");
  if(num_inf_rec > 0) {
    inf_record_ptr = (char **)array_alloc(__FILE__, __LINE__, 1, num_inf_rec,
					  sizeof(char *));
    for(int i = 0; i < num_inf_rec; i++) {
      inf_record_ptr[i] = (char *)
	array_alloc(__FILE__, __LINE__, 1, (MAX_LINE_LENGTH+2),
		    sizeof(char));
    }
    error = ex_get_info(lb_exoid, inf_record_ptr);
    check_exodus_error(error, "ex_get_info");
  }

  Int_Space[0] = 0;

  for(int i=0; i <Proc_Info[0]; i++) {

    if(ex_get_loadbal_param(lb_exoid,
			    &Int_Node_Num[i], &Bor_Node_Num[i],
			    &Ext_Node_Num[i], &Int_Elem_Num[i],
			    &Bor_Elem_Num[i], &Node_Comm_Num[i],
			    &Elem_Comm_Num[i], i) < 0) {
      fprintf(stderr, "%s: ERROR, could not get load balance params!\n",
	      yo);
      exit(1);
    }

#ifdef DEBUG
    if(Debug_Flag >= 5) {
      if(i == 0) {
	printf("--------------------------------------------------------\n");
	printf("\t\tLoad balance parameters as read by Processor 0\n");
	printf("--------------------------------------------------------\n");
      }
      printf("Read on processor 0 for processor %d:\n", i);
      printf("\tNumber internal nodes: %d\n", Int_Node_Num[i]);
      printf("\tNumber border nodes: %d\n", Bor_Node_Num[i]);
      printf("\tNumber external nodes: %d\n", Ext_Node_Num[i]);
      printf("\tNumber internal elements: %d\n", Int_Elem_Num[i]);
      printf("\tNumber border elements: %d\n", Bor_Elem_Num[i]);
      printf("\tNumber of nodal comm maps: %d\n", Node_Comm_Num[i]);
      printf("\tNumber of elemental comm maps: %d\n", Elem_Comm_Num[i]);
      printf("--------------------------------------------------------\n");
    }
#endif /* DEBUG */

    Int_Space[0] = PEX_MAX(Int_Space[0], Int_Node_Num[i]+Bor_Node_Num[i]+
			   Ext_Node_Num[i]+Int_Elem_Num[i]+
			   Bor_Elem_Num[i]);

  }

  /* Each processor extracts the information that it needs */
  for(int i=0; i <Proc_Info[2]; i++) {
    globals.Num_Internal_Nodes[i]  = Int_Space[1 + Proc_Ids[i]];
    globals.Num_Border_Nodes[i]    = Int_Space[1 + Proc_Ids[i] +  Proc_Info[0]];
    globals.Num_External_Nodes[i]  = Int_Space[1 + Proc_Ids[i] + 2*Proc_Info[0]];
    globals.Num_Internal_Elems[i]  = Int_Space[1 + Proc_Ids[i] + 3*Proc_Info[0]];
    globals.Num_Border_Elems[i]    = Int_Space[1 + Proc_Ids[i] + 4*Proc_Info[0]];
    globals.Num_N_Comm_Maps[i]     = Int_Space[1 + Proc_Ids[i] + 5*Proc_Info[0]];
    globals.Num_E_Comm_Maps[i]     = Int_Space[1 + Proc_Ids[i] + 6*Proc_Info[0]];
  }

  /*
   * Print Out a Summary of the Load Balance Information, as distributed
   * across the processors, for the appropriate value of Debug_Flag
   */
  if(Debug_Flag >= 3) {
    print_line("=", 79);
    printf("\n\t\tTABLE OF LOAD BALANCE STATISTICS\n\n");
    printf("%s%s\n", "globals. Int_Nodes Bor_Nodes Ext_Nodes",
	   " Int_Elems Bor_Elems N_Comm_Maps E_Comm_Maps");
    print_line("-", 79); printf ("\n");
    for(int i=0; i <Proc_Info[2]; i++) {
      printf("%6lu  %6lu  %6lu   %6lu    %6lu    %6lu     %6lu     %6lu\n",
             (size_t)Proc_Ids[i], (size_t)globals.Num_Internal_Nodes[i], (size_t)globals.Num_Border_Nodes[i],
             (size_t)globals.Num_External_Nodes[i], (size_t)globals.Num_Internal_Elems[i],
             (size_t)globals.Num_Border_Elems[i], (size_t)globals.Num_N_Comm_Maps[i], (size_t)globals.Num_E_Comm_Maps[i]);
    }
    print_line("=", 79);
    printf ("\n\n");
  }


} /* END of read_lb_init () **************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T,INT>::read_cmap_params(int lb_exoid,
					INT *Node_Comm_Num,
					INT *Elem_Comm_Num,
					INT *Num_N_Comm_Maps,
					INT *Num_E_Comm_Maps,
					ELEM_COMM_MAP<INT> **E_Comm_Map,
					NODE_COMM_MAP<INT> **N_Comm_Map,
					INT *cmap_max_size,
					INT **comm_vec)
/*
 * read_cmap_params:
 *
 *	This function reads the parameters for the communication maps.
 * Processor 0 reads each processors paramters and then a broadcast of
 * this information is performed and each processor extracts it's
 * information.
 */
{
  const char *yo="read_cmap_params";

  /*
   * Calculate the length of the longest vector needed. There is a
   * factor of 2 here, one for the counts, and one for the IDs.
   */
  size_t read_len = 0;
  for (int iproc=0; iproc <Proc_Info[0]; iproc++)
    read_len += 2*Node_Comm_Num[iproc] + 2*Elem_Comm_Num[iproc];

  /* Allocate a buffer */
  if(read_len != 0)
    *comm_vec = (INT *)array_alloc(__FILE__, __LINE__, 1, read_len,
				   sizeof(INT));
  else
    *comm_vec = NULL;

  size_t vec_start = 0;
  for (int iproc=0; iproc <Proc_Info[0]; iproc++) {

    /* Reset pointer for this processors information */
    INT* node_cm_ids  = *comm_vec    + vec_start;
    INT* node_cm_cnts = node_cm_ids  + Node_Comm_Num[iproc];
    INT* elem_cm_ids  = node_cm_cnts + Node_Comm_Num[iproc];
    INT*elem_cm_cnts = elem_cm_ids  + Elem_Comm_Num[iproc];

    /* Read the communication map IDs for processor "iproc" */
    if(ex_get_cmap_params(lb_exoid, node_cm_ids, node_cm_cnts,
			  elem_cm_ids, elem_cm_cnts, iproc) < 0) {
      fprintf(stderr,
	      "%s: ERROR, unable to read communication map params\n",
	      yo);
      exit(1);
    }

    /* Increment starting pointer */
    vec_start += 2*Node_Comm_Num[iproc] + 2*Elem_Comm_Num[iproc];

    /* Calculate the length of the broadcast information */
    read_len = 2*Node_Comm_Num[iproc] + 2*Elem_Comm_Num[iproc];

    /* For processor and node IDs */
    INT psum = 0;
    for(int i1=0; i1 < Node_Comm_Num[iproc]; i1++)
      psum += 2*node_cm_cnts[i1];

    /* For processor, element and side IDs */
    for(int i1=0; i1 < Elem_Comm_Num[iproc]; i1++)
      psum += 3*elem_cm_cnts[i1];

    *cmap_max_size = PEX_MAX(*cmap_max_size, psum);

    /*
     * NOTE: This code all assumes that the information for each
     * processor is read sequentially by processor 0.
     */
    int i1;
    for(i1=0; i1 <Proc_Info[2]; i1++)
      if(Proc_Ids[i1] == iproc)break;

    /*
     * This would need to be changed when multiple node/element
     * maps per processor are supported.
     */
    if(Node_Comm_Num[iproc] > 0) {
      N_Comm_Map[i1]->map_id   = node_cm_ids[0];
      N_Comm_Map[i1]->node_cnt = node_cm_cnts[0];
    }
    if(Elem_Comm_Num[iproc] > 0) {
      E_Comm_Map[i1]->map_id   = elem_cm_ids[0];
      E_Comm_Map[i1]->elem_cnt = elem_cm_cnts[0];
    }

  }

  if(Debug_Flag >= 4) {
    print_line("=", 79);
    printf("\t\tCOMMUNICATION MAP INFORMATION\n");
    printf("\t\t   largest cmap = %lu integers\n", (size_t)*cmap_max_size);
    print_line("=", 79);

    int bprnt_n=0;
    int bprnt_e=0;

    for(int i1=0; i1 <Proc_Info[2]; i1++) {
      if(globals.Num_N_Comm_Maps[i1] > 0)bprnt_n = 1;
      if(globals.Num_E_Comm_Maps[i1] > 0)bprnt_e = 1;
    }

    if(bprnt_n > 0) {
      printf("\tFor Proc\tNode Map ID\tNode Count\n");
      printf("\t------------------------------------------------\n");
      for (int iproc=0; iproc <Proc_Info[2]; iproc++) {
        if(globals.Num_N_Comm_Maps[iproc] > 0) {
          for(int i1=0; i1 < globals.Num_N_Comm_Maps[iproc]; i1++) {
            printf("\t     %lu\t\t    %lu\t\t    %lu\n",
		   (size_t)Proc_Ids[iproc],
                   (size_t)(N_Comm_Map[iproc]+i1)->map_id,
                   (size_t)(N_Comm_Map[iproc]+i1)->node_cnt);
          }
        }
      }
    }

    if(bprnt_e > 0) {
      printf("\tFor Proc\tElem Map ID\tElem Count\n");
      printf("\t------------------------------------------------\n");
      for (int iproc=0; iproc <Proc_Info[2]; iproc++) {
        if(globals.Num_E_Comm_Maps[iproc] > 0) {
          for(int i1=0; i1 < globals.Num_E_Comm_Maps[iproc]; i1++) {
            printf("\t     %lu\t\t    %lu\t\t    %lu\n",
                   (size_t)Proc_Ids[iproc],
                   (size_t)(E_Comm_Map[iproc]+i1)->map_id,
                   (size_t)(E_Comm_Map[iproc]+i1)->elem_cnt);
          }
        }
      }
    }

    print_line("=", 79);

  }

} /* End of read_cmap_params() ***********************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
