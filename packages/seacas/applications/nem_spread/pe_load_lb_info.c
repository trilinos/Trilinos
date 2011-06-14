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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "netcdf.h"
#include "exodusII.h"

#include "sort_utils.h"
#include "rf_salsa.h"
#include "pe_common.h"
#include "el_geom_const.h"
#include "rf_message.h"
#include "rf_util.h"

#include "rf_allo.h"
#include "rf_io_const.h"
#include "rf_mp_const.h"

#include "ps_pario_const.h"

#include "ne_nemesisI.h"

/************ R O U T I N E S   I N   T H I S   F I L E ***********************
*
*  Name_of_Routine		type		     Called by
*  ---------------------------------------------------------------
*
*  load_lb_info ()			        main:rf_salsa.c
*     set_lb_info_for_1Proc ()		    	load_lb_info
*     comm_lb_data ()				load_lb_info
*     recv_lb_data ()				load_lb_info
*     process_lb_data ()			load_lb_info
*     read_lb_init ()			        load_lb_info
*     read_cmap_params()			load_lb_info
*
******************************************************************************/

static void comm_lb_data (int iproc, int proc_for,
                          int *Integer_Vector_Ptr, int read_length);
static void recv_lb_data (int *Integer_Vector, int indx);
static void process_lb_data (int *Integer_Vector, int index);
static void read_proc_init(int lb_exoid, int proc_info[], int **proc_ids_ptr);
static void read_lb_init (int exoid, int *Int_Space, int Proc_Info[],
                          int *Int_Node_Num, int *Bor_Node_Num,
                          int *Ext_Node_Num, int *Int_Elem_Num,
                          int *Bor_Elem_Num, int *Node_Comm_Num,
                          int *Elem_Comm_Num, char *Title);
static void read_cmap_params(int exoid, int *Node_Comm_Num,
                             int *Elem_Comm_Num, int *Num_N_Comm_Maps,
                             int *Num_E_Comm_Maps, ELEM_COMM_MAP **E_Comm_Map,
                             NODE_COMM_MAP **N_Comm_Map, int *cmap_max_size,
                             int **comm_vec);

char **qa_record_ptr, **inf_record_ptr;
int num_inf_rec=0, num_qa_rec=0, length_qa=0;
char LB_File_Version[MAX_STR_LENGTH+1];

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
void load_lb_info(void)

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
 *        Num_Internal_Elems
 *        Num_Border_Elems
 *
 *    and the following vectors, also specific to each processor:
 *
 *        GNodes [Num_Internal_Nodes+Num_Border_Nodes+Num_External_Nodes]
 *        GElems [Num_Internal_Elems+Num_Border_Elems]
 *
 *    For the 1 processor case, this routine fills in appropriate values
 *    for these quantities.
 */
{

  int    iproc, lb_exoid, i, vec_indx, ijump, do_proc, send_to;
  int   *Integer_Vector = NULL;   /* Data structure for sending int data to  *
			           * processors                              */
  int   *Int_Space = NULL, *Int_Node_Num, *Bor_Node_Num, *Ext_Node_Num;
  int   *Node_Comm_Num, *Elem_Comm_Num;
  int   *Int_Elem_Num, *Bor_Elem_Num;
  int    read_length, cmap_max_size=0, *comm_vec, psum;
  char  *yo = "load_lb_info";
  char   Title[MAX_LINE_LENGTH+1];
  float  version;

  int    cpu_ws=0, io_ws=0;

#ifdef DEBUG_TIME
  double        start_time, Time_Read_lb_info = 0.0, Time_Comm_lb_info = 0.0,
    Time_File_Ops = 0.0, Time_Proc_lb_info = 0.0,
    Time_Load_lb_init = 0.0;
#endif

/******************************** START EXECUTION ****************************/

  /* The rest of the procedure handles multiple processor jobs */
  if(Proc == 0) {

    if(Debug_Flag)
      printf ("\nStart to read in and distribute the load balance info\n");

    /* Open the Load Balance exoII file for reading */
#ifdef DEBUG_TIME
    start_time = second ();
#endif

    printf ("EXODUS II load-balance file: %s\n", Exo_LB_File);
    if((lb_exoid = ex_open(Exo_LB_File, EX_READ, &cpu_ws,
			   &io_ws, &version)) == -1) {
      fprintf(stderr, "%sERROR: Couldn\'t open lb file, %s\n", yo,
              Exo_LB_File);
      exit(1);
    }

#ifdef DEBUG_TIME
    Time_File_Ops += second () - start_time;
#endif

  }

  /* Read information about the processor configuration */
  read_proc_init(lb_exoid, Proc_Info, &Proc_Ids);

  if(Debug_Flag >= 2) {
    if(Proc == 0) {
      print_line("=", 79);
      printf("\t\tTABLE OF PROCESSOR IDS\n\n");
      printf("Responsible Processor ID\t| Processors being handled\n");
      print_line("-", 79);
    }
    print_sync_start(Proc, Num_Proc, FALSE);
    printf("\t%d\t\t\t|", Proc);

    for(iproc=0; iproc < Proc_Info[2]; iproc++)
      printf(" %d", Proc_Ids[iproc]);

    printf("\n");

    print_sync_end(Proc, Num_Proc, FALSE);

    if(Proc == 0) {
      print_line("=", 79);
      printf("\n");
    }
  }

  /* Allocate space for the counts */
  Num_Internal_Nodes = (int *)array_alloc(__FILE__, __LINE__, 1,
                                          7*Proc_Info[2], sizeof(int));
  Num_Border_Nodes   = Num_Internal_Nodes + Proc_Info[2];
  Num_External_Nodes = Num_Border_Nodes   + Proc_Info[2];
  Num_Internal_Elems = Num_External_Nodes + Proc_Info[2];
  Num_Border_Elems   = Num_Internal_Elems + Proc_Info[2];
  Num_N_Comm_Maps    = Num_Border_Elems   + Proc_Info[2];
  Num_E_Comm_Maps    = Num_N_Comm_Maps    + Proc_Info[2];

  /* Allocate space for each processor entity */
  GNodes   = (int **)array_alloc(__FILE__, __LINE__, 1, 3*Proc_Info[2],
                                 sizeof(int *));
  GElems   = GNodes + Proc_Info[2];
  Elem_Map = GElems + Proc_Info[2];

  /* Allocate contiguous space for the pointer vectors on all processors */
  Int_Space     = (int *)array_alloc(__FILE__, __LINE__, 1,
                                     (7*Proc_Info[0] + 1), sizeof(int));
  Int_Node_Num  = Int_Space     + 1;
  Bor_Node_Num  = Int_Node_Num  + Proc_Info[0];
  Ext_Node_Num  = Bor_Node_Num  + Proc_Info[0];
  Int_Elem_Num  = Ext_Node_Num  + Proc_Info[0];
  Bor_Elem_Num  = Int_Elem_Num  + Proc_Info[0];
  Node_Comm_Num = Bor_Elem_Num  + Proc_Info[0];
  Elem_Comm_Num = Node_Comm_Num + Proc_Info[0];

#ifdef DEBUG_TIME
  start_time = second ();
#endif

  /* Read the initial information contained in the load balance file */
  read_lb_init(lb_exoid, Int_Space, Proc_Info, Int_Node_Num,
               Bor_Node_Num, Ext_Node_Num, Int_Elem_Num,
               Bor_Elem_Num, Node_Comm_Num, Elem_Comm_Num,
               Title);

  /* Allocate memory for the communication map arrays */
  N_Comm_Map = malloc(Proc_Info[2]*sizeof(NODE_COMM_MAP *));
  E_Comm_Map = malloc(Proc_Info[2]*sizeof(ELEM_COMM_MAP *));
  if(!N_Comm_Map || !E_Comm_Map) {
    fprintf(stderr, "ERROR: Insufficient memory!\n");
    exit(1);
  }

  for(iproc=0; iproc < Proc_Info[2]; iproc++) {

    /*
     * Error check:
     *	Currently a maximum of one nodal communication map and one
     * elemental communication map is supported.
     */
    if(Num_N_Comm_Maps[iproc] > 1 || Num_E_Comm_Maps[iproc] > 1) {
      fprintf(stderr, "%s: ERROR. Only 1 nodal and elemental comm map "
              "is supported\n", yo);
      exit(1);
    }
    else {

      /* Always allocate at least one and initialize the counts to 0 */
      N_Comm_Map[iproc] = (NODE_COMM_MAP *)malloc(PEX_MAX(1,
                                                  Num_N_Comm_Maps[iproc]) *
                                                  sizeof(NODE_COMM_MAP));
      if(N_Comm_Map[iproc] == NULL && Num_N_Comm_Maps[iproc] > 0) {
        fprintf(stderr,
                "%s: ERROR. Insufficient memory for nodal comm. map!\n",
                yo);
        exit(1);
      }

      for(ijump=0; ijump < PEX_MAX(1, Num_N_Comm_Maps[iproc]); ijump++)
        ((N_Comm_Map[iproc])+ijump)->node_cnt = 0;

      E_Comm_Map[iproc] = (ELEM_COMM_MAP *)malloc(PEX_MAX(1,
                                                      Num_E_Comm_Maps[iproc]) *
                                                  sizeof(ELEM_COMM_MAP));
      if(E_Comm_Map[iproc] == NULL && Num_E_Comm_Maps[iproc] > 0) {
        fprintf(stderr,
                "%s: ERROR. Insufficient memory for elemental comm. map!\n",
                yo);
        exit(1);
      }

      for(ijump=0; ijump < PEX_MAX(1, Num_E_Comm_Maps[iproc]); ijump++)
        ((E_Comm_Map[iproc])+ijump)->elem_cnt = 0;

    }

  } /* End "for(iproc=0; iproc < Proc_Info[2]; iproc++)" */

#ifdef DEBUG_TIME
  Time_Load_lb_init += second () - start_time;
#endif

  /* Set up each processor for the communication map parameters */
  read_cmap_params(lb_exoid, Node_Comm_Num, Elem_Comm_Num,
                   Num_N_Comm_Maps, Num_E_Comm_Maps,
                   E_Comm_Map, N_Comm_Map, &cmap_max_size,
                   &comm_vec);

  /* Allocate enough space to read the LB_data for one processor */
  Integer_Vector = (int *)array_alloc(__FILE__, __LINE__, 1,
                                      Int_Space[0] + cmap_max_size,
                                      sizeof (int));

  /*
   * On Proc == 0, loop through the processors, one at a time, to read
   * their load balance information
   *
   * NOTE: From here on there are no provisions for multiple nodal
   *       or elemental communication maps.
   */

  if(Proc == 0) {
    ijump = 0; /* keep track of where in comm_vec we are */
    for(iproc = 0; iproc < Proc_Info[0]; iproc++) {

      psum = 0;

      /*
       * Calculate the size of the vector needed to store the communication
       * map information. 2*number of nodes in nodal map (1 for node IDs and
       * 1 for processor IDs associated with each node) + 3*number of
       * elements in the elemental map (1 for each element ID, 1 for each
       * side ID associated with each element and 1 for each processor ID
       * associated with each element).
       */
      if(Node_Comm_Num[iproc] > 0)
        psum = 2*comm_vec[ijump+1] ;
      if(Elem_Comm_Num[iproc] > 0)
        psum += 3*comm_vec[ijump+3];

      /*
       * Calculate the length of the read for the current processor
       *   - cross-check pointer vectors
       */
#ifdef DEBUG_TIME
      start_time = second ();
#endif

      read_length = Int_Node_Num[iproc] + Bor_Node_Num[iproc] +
                    Ext_Node_Num[iproc] + Int_Elem_Num[iproc] +
                    Bor_Elem_Num[iproc] + psum;

      /* Get the node map for processor "iproc" */
      if(ne_get_node_map(lb_exoid, &Integer_Vector[0],
                         &Integer_Vector[Int_Node_Num[iproc]],
                         &Integer_Vector[Int_Node_Num[iproc]+
                                        Bor_Node_Num[iproc]],
                         iproc) < 0) {
        fprintf(stderr, "%s: ERROR, failed to get node map for Proc %d!\n",
                yo, iproc);
        exit(1);
      }

      vec_indx = Int_Node_Num[iproc] + Bor_Node_Num[iproc] +
                 Ext_Node_Num[iproc];

      /* Get the element map for processor number "iproc" */
      if(ne_get_elem_map(lb_exoid, &Integer_Vector[vec_indx],
                         &Integer_Vector[vec_indx+Int_Elem_Num[iproc]],
                         iproc) < 0) {
        fprintf(stderr, "%s: ERROR, failed to get element map for Proc %d!\n",
                yo, iproc);
        exit(1);
      }

      if(Node_Comm_Num[iproc] > 0) {
        vec_indx += Int_Elem_Num[iproc] + Bor_Elem_Num[iproc];

        if(ne_get_node_cmap(lb_exoid, comm_vec[ijump],
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
                  yo, Proc);
        }
      }

      if(Elem_Comm_Num[iproc] > 0) {
        vec_indx += 2*comm_vec[ijump+1];

        if(ne_get_elem_cmap(lb_exoid, comm_vec[ijump+2],
                          &Integer_Vector[vec_indx],
                          &Integer_Vector[vec_indx+comm_vec[ijump+3]],
                          &Integer_Vector[vec_indx+2*comm_vec[ijump+3]],
                          iproc) < 0) {
          fprintf(stderr,
                  "%s: ERROR. Failed to get elemental comm map for Proc %d!\n",
                  yo, Proc);
          exit(1);
        }
      }

#ifdef DEBUG_TIME
      Time_Read_lb_info += second () - start_time;
      start_time = second ();
#endif

      /*
       * Communicate load balance information to the correct processor
       *   - if iproc = Proc_Ids[*] then process the data instead.
       */
      do_proc = -1;
      for(i=0; i < Proc_Info[2]; i++) {
        if(iproc == Proc_Ids[i]) {
          do_proc = i;
        }
      }

      if(do_proc >= 0) {
        process_lb_data (Integer_Vector, do_proc);

#ifdef DEBUG_TIME
        Time_Proc_lb_info += second () - start_time;
#endif
      }
      else {

        send_to = iproc;
        while(send_to >= Num_Proc)
          send_to -= Num_Proc;

        comm_lb_data (send_to, iproc, Integer_Vector, read_length);

#ifdef DEBUG_TIME
        Time_Comm_lb_info += second () - start_time;
#endif
      }

      /*
       * now move ijump to the next communications map
       * make sure to check if there are any for this processor
       */
      if (Node_Comm_Num[iproc] > 0) ijump += 2;
      if (Elem_Comm_Num[iproc] > 0) ijump += 2;

    } /* End "for(iproc = 0; iproc < Num_Proc; iproc++)" */

#ifdef DEBUG_TIME
    start_time = second ();
#endif

    /* Close the load balance file - we are finished with it */
    if(ex_close (lb_exoid) == -1) {
      fprintf (stderr, "%sERROR: Error in closing load balance file\n", yo);
      exit(1);
    }

#ifdef DEBUG_TIME
    Time_File_Ops += second () - start_time;
#endif

  } /* END if(Proc == 0) */
  else {

    /*
     * Receive the load balance information, sent by Proc 0,
     * from the communication buffers
     */
    for(iproc=0; iproc < Proc_Info[2]; iproc++) {

      recv_lb_data(Integer_Vector, iproc);

      /* Process the load balance information */
      process_lb_data(Integer_Vector, iproc);
    }

  }

 /************************* Cleanup and Printout Phase ***********************/

  /* Free temporary memory */
  safe_free((void **) &Integer_Vector);

  if(num_qa_rec > 0) {
    for(i = 0; i < length_qa; i++) safe_free((void **) &(qa_record_ptr[i]));
    safe_free((void **) &qa_record_ptr);
  }

  if(num_inf_rec > 0) {
    for(i = 0; i < num_inf_rec; i++) safe_free((void **) &(inf_record_ptr[i]));
    safe_free((void **) &inf_record_ptr);
  }

  safe_free((void **) &Int_Space);

  if(Proc == 0)
    safe_free((void **) &comm_vec);

  /*========================================================================*/

  if(Debug_Flag && Proc == 0)
    printf ("\nFinished distributing load balance info\n");

  /* Output Detailed timing information for the progam */
#ifdef DEBUG_TIME
  if(Debug_Flag && Proc == 0) {
    printf ("\n\nload_lb_info timing information:\n");
    printf ("\tTime to Process Initial load balance info       = %f seconds\n",
	    Time_Load_lb_init);
    printf ("\tTime to Read load balance data                  = %f seconds\n",
	    Time_Read_lb_info);
    printf ("\tTime to distribute load balance data            = %f seconds\n",
	    Time_Comm_lb_info);
    printf ("\tTime to Process load balance data on proc 0     = %f seconds\n",
	    Time_Proc_lb_info);
    printf ("\tTime to do file operations (open and close)     = %f seconds\n",
	    Time_File_Ops);
    printf
      ("\t---------------------------------------------------------------\n");
    printf ("\t TOTAL TIME ACCOUNTED FOR                       = %f seconds\n",
	    (Time_Load_lb_init + Time_Read_lb_info + Time_Comm_lb_info +
             Time_Proc_lb_info + Time_File_Ops)  );
    printf ("\n");
  }
#endif

  /*
   * Print out a Large table of Load Balance Information if the debug_flag
   * setting is large enough
   */

  if(Debug_Flag >= 7) {
    print_sync_start (Proc, Num_Proc, TRUE);
    printf ("\n\n");
    print_line ("=", 79);
    printf("Processor %d:\n", Proc);

    for(iproc=0; iproc < Proc_Info[2]; iproc++) {
      printf("\n\t***For Processor %d***\n", Proc_Ids[iproc]);
      printf("\tInternal nodes owned by the current processor\n\t");
      for(i = 0; i < Num_Internal_Nodes[iproc]; i++)
        printf(" %d", GNodes[iproc][i]);

      printf("\n");

      printf("\tBorder nodes owned by the current processor\n\t");
      for(i = 0; i < Num_Border_Nodes[iproc]; i++)
        printf(" %d", GNodes[iproc][i + Num_Internal_Nodes[iproc]]);

      printf("\n");

      if(Num_External_Nodes[iproc] > 0) {
        printf("\tExternal nodes needed by the current processor\n\t");
        for(i = 0; i < Num_External_Nodes[iproc]; i++)
          printf(" %d", GNodes[iproc][i + Num_Internal_Nodes[iproc] +
                                     Num_Border_Nodes[iproc]]);

        printf("\n");
      }

      printf("\tInternal elements owned by the current processor\n\t");
      for(i = 0; i < Num_Internal_Elems[iproc]; i++)
        printf(" %d", GElems[iproc][i]);

      printf("\n");

      if(Num_Border_Elems[iproc] > 0) {
        printf("\tBorder elements owned by the current processor\n\t");
        for(i=0; i < Num_Border_Elems[iproc]; i++)
          printf(" %d", GElems[iproc][i+Num_Internal_Elems[iproc]]);

        printf("\n");
      }

      if(Num_N_Comm_Maps[iproc] > 0) {
        printf("\tNodal Comm Map for the current processor\n");
        printf("\t\tnode IDs:");
        for(i=0; i < N_Comm_Map[iproc]->node_cnt; i++)
          printf(" %d", N_Comm_Map[iproc]->node_ids[i]);
        printf("\n\t\tproc IDs:");
        for(i=0; i < N_Comm_Map[iproc]->node_cnt; i++)
          printf(" %d", N_Comm_Map[iproc]->proc_ids[i]);

        printf("\n");
      }

      if(Num_E_Comm_Maps[iproc] > 0) {
        printf("\tElemental Comm Map for the current processor\n");
        printf("\t\telement IDs:");
        for(i=0; i < E_Comm_Map[iproc]->elem_cnt; i++)
          printf(" %d", E_Comm_Map[iproc]->elem_ids[i]);
        printf("\n\t\tside IDs:");
        for(i=0; i < E_Comm_Map[iproc]->elem_cnt; i++)
          printf(" %d", E_Comm_Map[iproc]->side_ids[i]);
        printf("\n\t\tproc IDs:");
        for(i=0; i < E_Comm_Map[iproc]->elem_cnt; i++)
          printf(" %d", E_Comm_Map[iproc]->proc_ids[i]);

        printf("\n");
      }
    }
    printf("\n");
    print_line ("=", 79);
    print_sync_end (Proc, Num_Proc, TRUE);
  }

} /* END of routine load_lb_info () ******************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void
comm_lb_data (int iproc, int proc_for, int *Integer_Vector, int read_length)

/*
*       Function which communicates the load balance information to a processor
*       (iproc).  This function is run on processor 0, only.
*
*/

{
  int byte_length_int, st;

  /* Send the vector */
  if(iproc != 0) {
    byte_length_int = read_length*sizeof(int);
    if(byte_length_int > 0) {
      if(nwrite_big ((char *)Integer_Vector, byte_length_int, iproc,
                     EXOII_MESH_INT+proc_for, &st) != 0 ) {
        fprintf (stderr, "comm_lb_data: ERROR on node %d\n", Proc);
        fprintf (stderr, "              nwrite_big failed, message type %d\n",
                 EXOII_MESH_INT);
        exit(1);
      }
    }
  }

} /* END of comm_lb_data () **************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void recv_lb_data (int *Integer_Vector, int indx)

/*
*  Function which receives the load_balance data from Proc 0
*
*/
{
  int iproc, itype, read_length, byte_length_int, st;

  read_length = Num_Internal_Nodes[indx] + Num_Border_Nodes[indx] +
                Num_External_Nodes[indx] + Num_Internal_Elems[indx] +
                Num_Border_Elems[indx]   +
                2*N_Comm_Map[indx]->node_cnt + 3*E_Comm_Map[indx]->elem_cnt;

  byte_length_int = read_length * sizeof (int);
  iproc = 0;
  itype = EXOII_MESH_INT+Proc_Ids[indx];
  if(byte_length_int > 0) {
    if(nread_big ((char *)Integer_Vector, byte_length_int, &iproc, &itype, &st)
       != byte_length_int ) {
      fprintf (stderr, "recv_lb_data: ERROR on node %d\n", Proc);
      fprintf (stderr, "              nread failed, message type %d/n", itype);
      exit(1);
    }
  }

} /* END of routine recv_lb_data () ******************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void process_lb_data (int *Integer_Vector, int indx)

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
  int i, icount = 0, itotal_nodes, itotal_elems, ig_count = 0;

  /* External Variables */
  extern int  Proc;

  /* Function declarations */
  extern void sort_int  ( int n, int ra[] );
  extern int  check_monot (int vector[], int length);

/***************************** execution begins ******************************/

  /* Calculate the length of the GNodes and GElems global variable */
  itotal_nodes = Num_Internal_Nodes[indx] + Num_Border_Nodes[indx] +
                 Num_External_Nodes[indx];
  itotal_elems = Num_Internal_Elems[indx] + Num_Border_Elems[indx];

  /* Allocate Permament Arrays on the current processor */
  GNodes[indx] = (int *)
    array_alloc(__FILE__, __LINE__, 1,
               itotal_nodes + 2*itotal_elems + 2*(N_Comm_Map[indx]->node_cnt) +
               3*(E_Comm_Map[indx]->elem_cnt), sizeof(int));

  GElems[indx]   = GNodes[indx] + itotal_nodes;
  Elem_Map[indx] = GElems[indx] + itotal_elems;
  N_Comm_Map[indx]->node_ids = Elem_Map[indx] + itotal_elems;
  N_Comm_Map[indx]->proc_ids = N_Comm_Map[indx]->node_ids +
                               N_Comm_Map[indx]->node_cnt;
  E_Comm_Map[indx]->elem_ids = N_Comm_Map[indx]->proc_ids +
                               N_Comm_Map[indx]->node_cnt;
  E_Comm_Map[indx]->side_ids = E_Comm_Map[indx]->elem_ids +
                               E_Comm_Map[indx]->elem_cnt;
  E_Comm_Map[indx]->proc_ids = E_Comm_Map[indx]->side_ids +
                               E_Comm_Map[indx]->elem_cnt;

  /*
   *            Extract the load balance information, and store it in
   *            permanent vectors:
   *              -   GNodes[0]          - Internal_Nodes
   *                  GNodes[+]          - Border_Nodes
   *                  GNodes[+]          - External_Nodes
   *                  GElems[0]          - Internal_Elems
   *                  GElems[+]          - Border_Elems
   */

  for(i = 0; i < Num_Internal_Nodes[indx]; i++)
    GNodes[indx][ig_count++] = Integer_Vector[icount++];
  for(i = 0; i < Num_Border_Nodes[indx]; i++)
    GNodes[indx][ig_count++] = Integer_Vector[icount++];
  for(i = 0; i < Num_External_Nodes[indx]; i++)
    GNodes[indx][ig_count++] = Integer_Vector[icount++];

  ig_count = 0;
  for(i = 0; i < Num_Internal_Elems[indx]; i++) {
    GElems[indx][ig_count] = Integer_Vector[icount++];
    Elem_Map[indx][ig_count] = GElems[indx][ig_count];
    ig_count++;
  }
  for(i = 0; i < Num_Border_Elems[indx]; i++) {
    GElems[indx][ig_count] = Integer_Vector[icount++];
    Elem_Map[indx][ig_count] = GElems[indx][ig_count];
    ig_count++;
  }

  for(i = 0; i < N_Comm_Map[indx]->node_cnt; i++)
    (N_Comm_Map[indx]->node_ids)[i] = Integer_Vector[icount++];
  for(i = 0; i < N_Comm_Map[indx]->node_cnt; i++)
    (N_Comm_Map[indx]->proc_ids)[i] = Integer_Vector[icount++];
  for(i = 0; i < E_Comm_Map[indx]->elem_cnt; i++)
    (E_Comm_Map[indx]->elem_ids)[i] = Integer_Vector[icount++];
  for(i = 0; i < E_Comm_Map[indx]->elem_cnt; i++)
    (E_Comm_Map[indx]->side_ids)[i] = Integer_Vector[icount++];
  for(i = 0; i < E_Comm_Map[indx]->elem_cnt; i++)
    (E_Comm_Map[indx]->proc_ids)[i] = Integer_Vector[icount++];

  /*
   * Sort the local element numbers in ascending global element numbers.
   * This means that GElems will be monotonic.
   */
  gds_qsort(GElems[indx],   Num_Internal_Elems[indx]);
  gds_qsort(Elem_Map[indx], Num_Internal_Elems[indx]);

  /* Check that GNodes is monotonic, from i = 0 to Num_Internal_Nodes */
#ifdef DEBUG
  
  assert(check_monot (GNodes[indx], Num_Internal_Nodes[indx]));

  /*
   * Check that GNodes is monotonic, from i = Num_Internal_Nodes to
   *    (Num_Internal_Nodes + Num_Border_Nodes)
   */
  assert(check_monot(&(GNodes[indx][Num_Internal_Nodes[indx]]),
		     Num_Border_Nodes[indx]));
#endif

} /* END of process_lb_data () ***********************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_proc_init(int lb_exoid, int proc_info[], int **proc_ids_ptr)

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
  char *yo="read_proc_init";

  int  *proc_ids, i1;
  char  ftype[2];

  /* Get information about the size of the run */

  if(Proc == 0) {

    /*
     * If debugging is not on go ahead and report errors from init. This
     * will show version mismatch information by default.
     */
    if(Debug_Flag == 0)
      ex_opts(EX_VERBOSE);

    if(ne_get_init_info(lb_exoid, &proc_info[0], &proc_info[1], ftype) < 0) {
      fprintf(stderr, "%s: ERROR, could not get init info!\n",
              yo);
      exit(1);
    }

    if(Debug_Flag == 0)
      ex_opts(!EX_VERBOSE);

    if(Num_Proc > proc_info[0]) {
      fprintf(stderr, "ERROR: Must run on a number of processors <= %d\n",
              proc_info[0]);
      exit(1);
    }

    /*
     * If this is a single processor generating a file for another
     * processor then assign the global Num_Proc_For.
     */
    if(Proc_For >= 0)
      Num_Proc_For = proc_info[0];
  }

  /* Broadcast the information about the processor configuration */
  brdcst(Proc, Num_Proc, (char *) proc_info, 2*sizeof(int), 0);

  /* Calculate which processor is responsible for what */
  proc_info[2] = proc_info[0] / Num_Proc;
  proc_info[3]  = proc_info[0] % Num_Proc;
  if(Proc <= (proc_info[3]-1))
    proc_info[2]++;

  proc_ids = (int *)array_alloc(__FILE__, __LINE__, 1, proc_info[2],
                                sizeof(int));

  for(i1=0; i1 < proc_info[2]; i1++)
    proc_ids[i1] = Proc + i1*Num_Proc;

  *proc_ids_ptr = proc_ids;

  return;

} /* End of read_proc_init() *************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_lb_init(int lb_exoid,
                         int *Int_Space,
                         int  Proc_Info[],
                         int *Int_Node_Num,
                         int *Bor_Node_Num,
                         int *Ext_Node_Num,
                         int *Int_Elem_Num,
                         int *Bor_Elem_Num,
                         int *Node_Comm_Num,
                         int *Elem_Comm_Num,
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
  int   i, error;
  int   num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;
  int   int_brd_len;
  char *yo = "read_lb_init";

 /*********************BEGIN EXECUTABLE STATEMENTS****************************/

  /* Read Set-up information from the load-balance file on Proc 0 */
  if(Proc == 0) {

    /*
     * If debugging is not on go ahead and report errors from init. This
     * will show version mismatch information by default.
     */
    if(Debug_Flag == 0)
      ex_opts(EX_VERBOSE);

    /* Read the title of the LB File and about the size of the mesh */
    error = ne_get_init_global(lb_exoid, &num_nodes, &num_elem,
                               &num_elem_blk, &num_node_sets, &num_side_sets);

    check_exodus_error (error, "ne_get_init");

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
    if((num_nodes != Num_Node) || (num_elem != Num_Elem) ||
        (num_elem_blk != Num_Elem_Blk)) {
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
      for(i = 0; i < length_qa; i++) {
        qa_record_ptr[i] = (char *)
          array_alloc (__FILE__, __LINE__, 1, (MAX_STR_LENGTH+1),
                       sizeof(char));
      }
      error = ex_get_qa(lb_exoid, (char *(*)[]) &qa_record_ptr[0]);
      check_exodus_error(error, "ex_get_qa");
    }

    /* Read the Info Records */
    error = ex_inquire(lb_exoid, EX_INQ_INFO, &num_inf_rec, NULL,  NULL);
    check_exodus_error(error, "ex_inquire");
    if(num_inf_rec > 0) {
      inf_record_ptr = (char **)array_alloc(__FILE__, __LINE__, 1, num_inf_rec,
                                            sizeof(char *));
      for(i = 0; i < num_inf_rec; i++) {
        inf_record_ptr[i] = (char *)
          array_alloc(__FILE__, __LINE__, 1, (MAX_LINE_LENGTH+2),
                      sizeof(char));
      }
      error = ex_get_info(lb_exoid, inf_record_ptr);
      check_exodus_error(error, "ex_get_info");
    }

    Int_Space[0] = 0;

    for(i=0; i < Proc_Info[0]; i++) {

      if(ne_get_loadbal_param(lb_exoid,
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

  } /* END if(Proc == 0) */

  /* Broadcast all of the set-up information to all of the processors */
  int_brd_len  = 7*Proc_Info[0] + 1;

  brdcst_maxlen(Proc, Num_Proc, (char *) Int_Space,
                (sizeof(int)*int_brd_len), 0);

  /* Each processor extracts the information that it needs */
  for(i=0; i < Proc_Info[2]; i++) {
    Num_Internal_Nodes[i]  = Int_Space[1 + Proc_Ids[i]];
    Num_Border_Nodes[i]    = Int_Space[1 + Proc_Ids[i] +   Proc_Info[0]];
    Num_External_Nodes[i]  = Int_Space[1 + Proc_Ids[i] + 2*Proc_Info[0]];
    Num_Internal_Elems[i]  = Int_Space[1 + Proc_Ids[i] + 3*Proc_Info[0]];
    Num_Border_Elems[i]    = Int_Space[1 + Proc_Ids[i] + 4*Proc_Info[0]];
    Num_N_Comm_Maps[i]     = Int_Space[1 + Proc_Ids[i] + 5*Proc_Info[0]];
    Num_E_Comm_Maps[i]     = Int_Space[1 + Proc_Ids[i] + 6*Proc_Info[0]];
  }

  /*
   * Print Out a Summary of the Load Balance Information, as distributed
   * across the processors, for the appropriate value of Debug_Flag
   */
  if(Debug_Flag >= 3) {
    print_sync_start(Proc, Num_Proc, FALSE);
    if(Proc == 0) {
      print_line("=", 79);
      printf("\n\t\tTABLE OF LOAD BALANCE STATISTICS\n\n");
      printf("%s%s\n", "Proc_# Int_Nodes Bor_Nodes Ext_Nodes",
             " Int_Elems Bor_Elems N_Comm_Maps E_Comm_Maps");
      print_line("-", 79); printf ("\n");
    }
    printf("Processor %d storing for:\n", Proc);
    for(i=0; i < Proc_Info[2]; i++) {
      printf("%6d  %6d  %6d   %6d    %6d    %6d     %6d     %6d\n",
             Proc_Ids[i], Num_Internal_Nodes[i], Num_Border_Nodes[i],
             Num_External_Nodes[i], Num_Internal_Elems[i],
             Num_Border_Elems[i], Num_N_Comm_Maps[i], Num_E_Comm_Maps[i]);
    }

    print_sync_end(Proc, Num_Proc, FALSE);
    if(Proc == 0) {
      print_line("=", 79);
      printf ("\n\n");
    }

  }


} /* END of read_lb_init () **************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void read_cmap_params(int lb_exoid, int *Node_Comm_Num,
                             int *Elem_Comm_Num, int *Num_N_Comm_Maps,
                             int *Num_E_Comm_Maps, ELEM_COMM_MAP **E_Comm_Map,
                             NODE_COMM_MAP **N_Comm_Map, int *cmap_max_size,
                             int **comm_vec)
/*
 * read_cmap_params:
 *
 *	This function reads the parameters for the communication maps.
 * Processor 0 reads each processors paramters and then a broadcast of
 * this information is performed and each processor extracts it's
 * information.
 */
{
  char *yo="read_cmap_params";

  int   bprnt_n=0, bprnt_e=0;

  int   i1, iproc, read_len, itype, st, psum, vec_start, send_to;
  int   *node_cm_ids, *node_cm_cnts, *elem_cm_ids, *elem_cm_cnts;

/******************************** START EXECUTION ****************************/

  if(Proc == 0) {

    /*
     * Calculate the length of the longest vector needed. There is a
     * factor of 2 here, one for the counts, and one for the IDs.
     */
    read_len = 0;
    for(iproc=0; iproc < Proc_Info[0]; iproc++)
      read_len += 2*Node_Comm_Num[iproc] + 2*Elem_Comm_Num[iproc];

    /* Allocate a buffer */
    if(read_len != 0)
      *comm_vec = (int *)array_alloc(__FILE__, __LINE__, 1, read_len,
                                     sizeof(int));
    else
      *comm_vec = NULL;

    vec_start = 0;
    for(iproc=0; iproc < Proc_Info[0]; iproc++) {

      /* Reset pointer for this processors information */
      node_cm_ids  = *comm_vec    + vec_start;
      node_cm_cnts = node_cm_ids  + Node_Comm_Num[iproc];
      elem_cm_ids  = node_cm_cnts + Node_Comm_Num[iproc];
      elem_cm_cnts = elem_cm_ids  + Elem_Comm_Num[iproc];

      /* Read the communication map IDs for processor "iproc" */
      if(ne_get_cmap_params(lb_exoid, node_cm_ids, node_cm_cnts,
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
      psum = 0;
      for(i1=0; i1 < Node_Comm_Num[iproc]; i1++)
        psum += 2*node_cm_cnts[i1];

      /* For processor, element and side IDs */
      for(i1=0; i1 < Elem_Comm_Num[iproc]; i1++)
        psum += 3*elem_cm_cnts[i1];

      *cmap_max_size = PEX_MAX(*cmap_max_size, psum);

      /* Calculate the processor ID to send to */
      send_to = iproc;
      while(send_to >= Num_Proc)
        send_to -= Num_Proc;

      /*
       * NOTE: This code all assumes that the information for each
       * processor is read sequentially by processor 0.
       */
      if(send_to == Proc) {

        for(i1=0; i1 < Proc_Info[2]; i1++)
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
      else {
        if(read_len > 0) {
          if(nwrite_big((char *)node_cm_ids, read_len * sizeof(int), send_to,
                        EXOII_MESH_INT+iproc, &st) != 0) {
            fprintf(stderr, "%s: ERROR, nwrite_big failed, message type %d\n",
                    yo, EXOII_MESH_INT+10);
            exit(1);
          }
        }
      }

    } /* End "for(iproc=0; iproc < Num_Proc; iproc++)" */

  }
  else {

    *cmap_max_size = 0;
    for(iproc=0; iproc < Proc_Info[2]; iproc++) {
      psum = 0;

      if((Num_N_Comm_Maps[iproc]+Num_E_Comm_Maps[iproc]) > 0) {
        node_cm_ids  = (int *)array_alloc(__FILE__, __LINE__, 1,
                                          2*Num_N_Comm_Maps[iproc] +
                                          2*Num_E_Comm_Maps[iproc],
                                          sizeof(int));
        node_cm_cnts = node_cm_ids  + Num_N_Comm_Maps[iproc];
        elem_cm_ids  = node_cm_cnts + Num_N_Comm_Maps[iproc];
        elem_cm_cnts = elem_cm_ids  + Num_E_Comm_Maps[iproc];
      }
      else {
        node_cm_ids  = NULL;
        node_cm_cnts = NULL;
        elem_cm_ids  = NULL;
        elem_cm_cnts = NULL;
      }

      /* Receive the comm map parameters */
      read_len = 2*Num_N_Comm_Maps[iproc] + 2*Num_E_Comm_Maps[iproc];
      i1 = 0;
      itype = EXOII_MESH_INT+Proc_Ids[iproc];
      if(read_len > 0) {
        if(nread_big((char *)node_cm_ids, read_len * sizeof(int), &i1,
                     &itype, &st) != read_len*sizeof(int)) {
          fprintf(stderr,
                  "%s: ERROR, nread_big failed on node %d, message type %d\n",
                  yo, Proc, itype);
          exit(1);
        }
      }
      if(Num_N_Comm_Maps[iproc] > 0) {
        N_Comm_Map[iproc]->map_id   = node_cm_ids[0];
        N_Comm_Map[iproc]->node_cnt = node_cm_cnts[0];
        psum += 2*node_cm_cnts[0];
        *cmap_max_size = PEX_MAX(psum, *cmap_max_size);
      }
      if(Num_E_Comm_Maps[iproc] > 0) {
        E_Comm_Map[iproc]->map_id   = elem_cm_ids[0];
        E_Comm_Map[iproc]->elem_cnt = elem_cm_cnts[0];
        psum += 3*elem_cm_cnts[0];
        *cmap_max_size = PEX_MAX(psum, *cmap_max_size);
      }
      safe_free((void **) &node_cm_ids);

    } /* End "for(iproc=0; iproc < Proc_Info[2]; iproc++)" */

  }

  if(Debug_Flag >= 4) {
    print_sync_start(Proc, Num_Proc, FALSE);
    if(Proc == 0) {
      print_line("=", 79);
      printf("\t\tCOMMUNICATION MAP INFORMATION\n");
      printf("\t\t   largest cmap = %d integers\n", *cmap_max_size);
      print_line("=", 79);
    }

    for(i1=0; i1 < Proc_Info[2]; i1++) {
      if(Num_N_Comm_Maps[i1] > 0)bprnt_n = 1;
      if(Num_E_Comm_Maps[i1] > 0)bprnt_e = 1;
    }

    if(bprnt_n > 0 || bprnt_e > 0)
      printf("Processor %d:\n", Proc);

    if(bprnt_n > 0) {
      printf("\tFor Proc\tNode Map ID\tNode Count\n");
      printf("\t------------------------------------------------\n");
      for(iproc=0; iproc < Proc_Info[2]; iproc++) {
        if(Num_N_Comm_Maps[iproc] > 0) {
          for(i1=0; i1 < Num_N_Comm_Maps[iproc]; i1++) {
            printf("\t     %d\t\t    %d\t\t    %d\n",
                 Proc_Ids[iproc],
                   (N_Comm_Map[iproc]+i1)->map_id,
                   (N_Comm_Map[iproc]+i1)->node_cnt);
          }
        }
      }
    }

    if(bprnt_e > 0) {
      printf("\tFor Proc\tElem Map ID\tElem Count\n");
      printf("\t------------------------------------------------\n");
      for(iproc=0; iproc < Proc_Info[2]; iproc++) {
        if(Num_E_Comm_Maps[iproc] > 0) {
          for(i1=0; i1 < Num_E_Comm_Maps[iproc]; i1++) {
            printf("\t     %d\t\t    %d\t\t    %d\n",
                   Proc_Ids[iproc],
                   (E_Comm_Map[iproc]+i1)->map_id,
                   (E_Comm_Map[iproc]+i1)->elem_cnt);
          }
        }
      }
    }

    print_sync_end(Proc, Num_Proc, FALSE);

    if(Proc == 0)
      print_line("=", 79);

  }

} /* End of read_cmap_params() ***********************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
