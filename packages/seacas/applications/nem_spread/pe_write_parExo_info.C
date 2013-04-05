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
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for fprintf, printf, NULL, etc
#include <stdlib.h>                     // for exit, free, malloc
#include <string.h>                     // for strcpy, strlen, memset, etc
#include <sys/select.h>                 // for time_t
#include <time.h>                       // for asctime, localtime, time
#include <vector>                       // for vector
#include "exodusII.h"                   // for ex_close, etc
#include "nem_spread.h"                 // for NemSpread, second, etc
#include "pe_common.h"                  // for PEX_MAX
#include "ps_pario_const.h"             // for PIO_Time_Array
#include "rf_allo.h"                    // for safe_free, array_alloc
#include "rf_io_const.h"                // for Debug_Flag
#include "sort_utils.h"                 // for gds_iqsort

template <typename INT> struct ELEM_COMM_MAP;
template <typename INT> struct NODE_COMM_MAP;

#define TOPTR(x) (x.empty() ? NULL : &x[0])

namespace {
  template <typename INT>
  void reverse_map(INT *global, int p01, size_t gsize,
		   INT *glmap, INT *index, size_t msize,
		   INT *mapout);
}
/*
 * need this variable for the 0 processor to hold on to the correct
 * Exodus II database title
 */
extern char GeomTitle[];

/****************************************************************************/
/* This function writes parallel specific mesh information out to the       */
/* parallel disk(s) for each processor in the salsa run.                    */
/*                                                                          */
/* Author(s): Gary L. Hennigan (1421)                                       */
/*--------------------------------------------------------------------------*/
/* Revision History                                                         */
/*                                                                          */
/*   Gary Hennigan:                                                         */
/*      11 November 1993 - Date of first working version on the nCube.      */
/****************************************************************************/

template void NemSpread<float,int>::write_parExo_data(int mesh_exoid, int max_name_length,
							int iproc,
							int *Num_Nodes_In_NS,
							int *Num_Elems_In_SS,
							int *Num_Elems_In_EB);

template void NemSpread<double,int>::write_parExo_data(int mesh_exoid, int max_name_length,
							int iproc,
							int *Num_Nodes_In_NS,
							int *Num_Elems_In_SS,
							int *Num_Elems_In_EB);

template void NemSpread<double,int64_t>::write_parExo_data(int mesh_exoid, int max_name_length,
							   int iproc,
							   int64_t *Num_Nodes_In_NS,
							   int64_t *Num_Elems_In_SS,
							   int64_t *Num_Elems_In_EB);

template void NemSpread<float,int64_t>::write_parExo_data(int mesh_exoid, int max_name_length,
							  int iproc,
							  int64_t *Num_Nodes_In_NS,
							  int64_t *Num_Elems_In_SS,
							  int64_t *Num_Elems_In_EB);

template <typename T, typename INT>
void NemSpread<T,INT>::write_parExo_data(int mesh_exoid, int max_name_length,
					 int iproc,
					 INT *Num_Nodes_In_NS,
					 INT *Num_Elems_In_SS,
					 INT *Num_Elems_In_EB)
{
  static   char yo[] = "write_parExo_data";

  /* Performance metrics. */
  unsigned long    bytes_out=0;
  double           total_out_time=0.0, tt1;

  int error;
  /****************************BEGIN EXECUTION*********************************/

  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*                    PARALLEL EXODUSII SECTION                            */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  int proc_for     = Proc_Ids[iproc];
  int num_proc_for =Proc_Info[0];

  size_t itotal_nodes = globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc]+
    globals.Num_External_Nodes[iproc];
  size_t itotal_elems = globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc];

  bytes_out += 5*sizeof(INT);
  tt1 = second();

  if(Debug_Flag >= 4) {
    printf("Putting init global info in file id: %d\n", mesh_exoid);
    printf("\tNumber Global Nodes: %lu\n", (size_t)globals.Num_Node);
    printf("\tNumber Global Elements: %lu\n", (size_t)globals.Num_Elem);
    printf("\tNumber Global Element Blocks: %lu\n", (size_t)globals.Num_Elem_Blk);
    printf("\tNumber Global Node Sets: %lu\n", (size_t)globals.Num_Node_Set);
    printf("\tNumber Global Side Sets: %lu\n", (size_t)globals.Num_Side_Set);
  }

  if(ex_put_init_global(mesh_exoid, globals.Num_Node, globals.Num_Elem,
                        globals.Num_Elem_Blk, globals.Num_Node_Set,
                        globals.Num_Side_Set) < 0) {
    fprintf(stderr, "[%s]: ERROR, Unable to put global initial "
            "information in parallel mesh file!\n", yo);
    exit(1);
  }

  PIO_Time_Array[0]  = (second() - tt1);
  total_out_time    += PIO_Time_Array[0];

  /*
   * Find the number of distinct nodal communication maps that will
   * be needed. This assumes that there is only a single input
   * nodal communication map that stores information for all
   * processors.
   */
  int ncomm_cnt = 0;
  std::vector<NODE_COMM_MAP<INT> > n_comm_map;
  std::vector<INT> n_comm_ids;
  std::vector<INT> n_comm_ncnts;
  
  if(globals.Num_N_Comm_Maps[iproc] > 0 && globals.N_Comm_Map[iproc]->node_cnt != 0 ) {
    ncomm_cnt = 1;
    ex_entity_id itemp = globals.N_Comm_Map[iproc]->proc_ids[0];

    /* First find the count */
    for(size_t i1=1; i1 < globals.N_Comm_Map[iproc]->node_cnt; i1++) {
      if(globals.N_Comm_Map[iproc]->proc_ids[i1] != itemp) {
        ncomm_cnt++;
        itemp = globals.N_Comm_Map[iproc]->proc_ids[i1];
      }
    }

    /* Allocate memory for the nodal communication maps */
    n_comm_map.resize(ncomm_cnt);
    n_comm_ids.resize(ncomm_cnt);
    n_comm_ncnts.resize(ncomm_cnt);

    /* Find the size of each map */
    ncomm_cnt = 0;
    itemp = globals.N_Comm_Map[iproc]->proc_ids[0];
    n_comm_map[0].node_cnt = 1;
    for(size_t i1=1; i1 < globals.N_Comm_Map[iproc]->node_cnt; i1++) {
      if(globals.N_Comm_Map[iproc]->proc_ids[i1] != itemp) {
        itemp = globals.N_Comm_Map[iproc]->proc_ids[i1];
        ncomm_cnt++;
        n_comm_map[ncomm_cnt].node_cnt = 1;
      }
      else
        (n_comm_map[ncomm_cnt].node_cnt)++;
    }

    ncomm_cnt++;
    size_t temp = 0;

    /* Allocate memory for the maps */
    for(int i1=0; i1 < ncomm_cnt; i1++) {
      n_comm_map[i1].proc_ids = (INT*)malloc(2*(n_comm_map[i1].node_cnt)*
					     sizeof(INT));
      if(!(n_comm_map[i1].proc_ids)) {
        fprintf(stderr, "[%s]: ERROR, insufficient memory!\n",
                yo);
        exit(1);
      }

      n_comm_map[i1].node_ids = n_comm_map[i1].proc_ids +
	n_comm_map[i1].node_cnt;

      for(size_t i2=0; i2 < n_comm_map[i1].node_cnt; i2++) {
        n_comm_map[i1].proc_ids[i2] = globals.N_Comm_Map[iproc]->proc_ids[temp];
        n_comm_map[i1].node_ids[i2] = globals.N_Comm_Map[iproc]->node_ids[temp++];
      }

      n_comm_ncnts[i1] = n_comm_map[i1].node_cnt;
      n_comm_ids[i1]   = n_comm_map[i1].proc_ids[0];
    }
  }
  else {
    ncomm_cnt = 0;
  }

  /*
   * Find the number of distinct elemental communication maps that will
   * be needed. This assumes that there is only a single input
   * elemental communication map that stores information for all
   * processors.
   */
  int ecomm_cnt = 0;
  std::vector<INT> e_comm_ids;
  std::vector<INT> e_comm_ecnts;
  std::vector<ELEM_COMM_MAP<INT> > e_comm_map;

  if(globals.Num_E_Comm_Maps[iproc] > 0) {
    ecomm_cnt = 1;
    int itemp = globals.E_Comm_Map[iproc]->proc_ids[0];

    /* First find the count */
    for(size_t i1=1; i1 < globals.E_Comm_Map[iproc]->elem_cnt; i1++) {
      if(globals.E_Comm_Map[iproc]->proc_ids[i1] != itemp) {
        ecomm_cnt++;
        itemp = globals.E_Comm_Map[iproc]->proc_ids[i1];
      }
    }

    /* Allocate memory for the elemental communication maps */
    e_comm_map.resize(ecomm_cnt);
    e_comm_ids.resize(ecomm_cnt);
    e_comm_ecnts.resize(ecomm_cnt);

    /* Find the size of each map */
    ecomm_cnt = 0;
    itemp = globals.E_Comm_Map[iproc]->proc_ids[0];
    e_comm_map[0].elem_cnt = 1;
    for(size_t i1=1; i1 < globals.E_Comm_Map[iproc]->elem_cnt; i1++) {
      if(globals.E_Comm_Map[iproc]->proc_ids[i1] != itemp) {
        itemp = globals.E_Comm_Map[iproc]->proc_ids[i1];
        ecomm_cnt++;
        e_comm_map[ecomm_cnt].elem_cnt = 1;
      }
      else
        (e_comm_map[ecomm_cnt].elem_cnt)++;
    }

    ecomm_cnt++;
    itemp = 0;

    /* Allocate memory for the maps */
    for(int i1=0; i1 < ecomm_cnt; i1++) {
      e_comm_map[i1].proc_ids = (INT*)malloc(3*(e_comm_map[i1].elem_cnt)*
					     sizeof(INT));
      if(!(e_comm_map[i1].proc_ids)) {
        fprintf(stderr, "[%s]: ERROR, insufficient memory!\n",
                yo);
        exit(1);
      }

      e_comm_map[i1].elem_ids = e_comm_map[i1].proc_ids + e_comm_map[i1].elem_cnt;
      e_comm_map[i1].side_ids = e_comm_map[i1].elem_ids + e_comm_map[i1].elem_cnt;

      for(size_t i2=0; i2 < e_comm_map[i1].elem_cnt; i2++) {
        e_comm_map[i1].proc_ids[i2] = globals.E_Comm_Map[iproc]->proc_ids[itemp];
        e_comm_map[i1].elem_ids[i2] = globals.E_Comm_Map[iproc]->elem_ids[itemp];
        e_comm_map[i1].side_ids[i2] = globals.E_Comm_Map[iproc]->side_ids[itemp++];
      }

      e_comm_ecnts[i1] = e_comm_map[i1].elem_cnt;
      e_comm_ids[i1]   = e_comm_map[i1].proc_ids[0];
    }
  }
  else {
    ecomm_cnt = 0;
  }

  /* Output load balance information */
  bytes_out += 9*sizeof(INT) + 2*sizeof(char);
  tt1 = second();

  if(Debug_Flag >= 4) {
    printf("Putting init Nemesis info in file id: %d\n", mesh_exoid);
    printf("\tNumber of Proccesor for: %d\n", num_proc_for);
  }

  if(ex_put_init_info(mesh_exoid, num_proc_for, 1, (char*)"p") < 0) {
    fprintf(stderr, "[%s]: ERROR, unable to output init info!\n",
            yo);
    exit(1);
  }

  if(Debug_Flag >= 6) {
    printf("Putting init load balance info in file id: %d\n",
           mesh_exoid);
    printf("\tNumber Internal Nodes: %lu\n", (size_t)globals.Num_Internal_Nodes[iproc]);
    printf("\tNumber Border Nodes: %lu\n", (size_t)globals.Num_Border_Nodes[iproc]);
    printf("\tNumber External Nodes: %lu\n", (size_t)globals.Num_External_Nodes[iproc]);
    printf("\tNumber Internal Elements: %lu\n", (size_t)globals.Num_Internal_Elems[iproc]);
    printf("\tNumber Border Elements: %lu\n", (size_t)globals.Num_Border_Elems[iproc]);
    printf("\tNumber Nodal Cmaps: %lu\n", (size_t)ncomm_cnt);
    printf("\tNumber Elemental Cmaps: %lu\n", (size_t)ecomm_cnt);
    printf("\tProccesor For: %d\n", proc_for);
  }

  if(ex_put_loadbal_param(mesh_exoid, globals.Num_Internal_Nodes[iproc],
                          globals.Num_Border_Nodes[iproc], globals.Num_External_Nodes[iproc],
                          globals.Num_Internal_Elems[iproc], globals.Num_Border_Elems[iproc],
                          ncomm_cnt, ecomm_cnt,
                          proc_for) < 0) {
    fprintf(stderr, "[%s]: ERROR, unable to output load balance info\n",
            yo);
    ex_close(mesh_exoid);
    exit(1);
  }

  PIO_Time_Array[1]  = (second() - tt1);
  total_out_time    += PIO_Time_Array[1];

  /* Output the communication map */
  bytes_out += 4*sizeof(INT);
  tt1 = second();

  if(ex_put_cmap_params(mesh_exoid,
                        TOPTR(n_comm_ids),
                        TOPTR(n_comm_ncnts),
			TOPTR(e_comm_ids),
			TOPTR(e_comm_ecnts),
                        proc_for) < 0) {
    fprintf(stderr, "[%s]: ERROR, unable to output comm map params!\n",
            yo);
    ex_close(mesh_exoid);
    exit(1);
  }

  PIO_Time_Array[2]  = (second() - tt1);
  total_out_time    += PIO_Time_Array[2];

  /*
   * The Nemesis node maps are lists of internal, border and external
   * FEM node numbers. These are output as local node numbers.
   */
  INT *nem_node_mapi = (INT *)array_alloc(__FILE__, __LINE__, 1,
					  itotal_nodes + itotal_elems,
					  sizeof(INT));
  INT *nem_node_mapb = nem_node_mapi + globals.Num_Internal_Nodes[iproc];
  INT *nem_node_mape = nem_node_mapb + globals.Num_Border_Nodes[iproc];

  for(size_t i1=0; i1 < itotal_nodes; i1++)
    nem_node_mapi[i1] = i1+1;

  /* Convert Elem_Map to local element numbering */
  reverse_map(globals.Elem_Map[iproc], 0, itotal_elems,
	      &globals.GElems[iproc][0], (INT*)NULL, itotal_elems,
	      &globals.Elem_Map[iproc][0]);
  /* Convert element IDs in the comm map to local numbering */
  for(int i0=0; i0 < ecomm_cnt; i0++) {
    reverse_map(e_comm_map[i0].elem_ids, 0, e_comm_map[i0].elem_cnt,
		&globals.GElems[iproc][0], (INT*)NULL, itotal_elems,
		e_comm_map[i0].elem_ids);
  }

  /* Sort the globals.GNodes array using the index array 'loc_index' */
  INT *loc_index = (INT *) array_alloc(__FILE__, __LINE__, 1, itotal_nodes,
				       sizeof(INT));

  /* Initialize index array */
  for (size_t i2=0; i2 < itotal_nodes; i2++)
    loc_index[i2] = i2;
      
  /*
   * Sort the globals.GNodes[iproc] array via the index array
   * 'loc_index' 
   */
  gds_iqsort(globals.GNodes[iproc], loc_index, itotal_nodes);
  
  /* Convert nodal IDs in the comm map to local numbering */
  for(int i0=0; i0 < ncomm_cnt; i0++) {
    reverse_map(n_comm_map[i0].node_ids, 0, n_comm_map[i0].node_cnt,
		&globals.GNodes[iproc][0], loc_index, itotal_nodes,
		n_comm_map[i0].node_ids);
  }

  PIO_Time_Array[3] = 0.0;

  if(globals.Num_N_Comm_Maps[iproc] > 0 && globals.N_Comm_Map[iproc]->node_cnt != 0 ) {

    bytes_out += 2*globals.Num_External_Nodes[iproc]*sizeof(INT);
    tt1 = second();

    for(int i1=0; i1 < ncomm_cnt; i1++) {
      if(ex_put_node_cmap(mesh_exoid, n_comm_ids[i1], n_comm_map[i1].node_ids,
                          n_comm_map[i1].proc_ids, proc_for) < 0) {
        fprintf(stderr, "[%s]: ERROR, unable to output nodal comm map!\n",
                yo);
        ex_close(mesh_exoid);
        exit(1);
      }

      free(n_comm_map[i1].proc_ids);
    }

    PIO_Time_Array[3] += (second() - tt1);

  }
  if(globals.Num_E_Comm_Maps[iproc] > 0) {
    bytes_out += 3*(globals.E_Comm_Map[iproc]->elem_cnt)*sizeof(INT);
    tt1 = second();

    for(int i1=0; i1 < ecomm_cnt; i1++) {
      if(ex_put_elem_cmap(mesh_exoid, e_comm_ids[i1], e_comm_map[i1].elem_ids,
			  e_comm_map[i1].side_ids, e_comm_map[i1].proc_ids,
			  proc_for) < 0) {
        fprintf(stderr,
                "[%s]: ERROR, unable to output elemental comm map!\n",
                yo);
        ex_close(mesh_exoid);
        exit(1);
      }

      free(e_comm_map[i1].proc_ids);
    }
    PIO_Time_Array[3] += (second() - tt1);
  }

  total_out_time    += PIO_Time_Array[3];

  /* Output the global node set parameters */
  if(globals.Num_Node_Set > 0) {

    std::vector<INT> glob_ns_df_cnts(globals.Num_Node_Set);

    bytes_out += 3*globals.Num_Node_Set*sizeof(INT);
    tt1 = second();

    if(ex_put_ns_param_global(mesh_exoid, Node_Set_Ids,
                              Num_Nodes_In_NS, TOPTR(glob_ns_df_cnts)) < 0) {
      fprintf(stderr,
              "[%s]: ERROR, unable to output global node-set params\n",
              yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    PIO_Time_Array[4]  = (second() - tt1);
    total_out_time    += PIO_Time_Array[4];
  }

  /* Output the global side set parameters */
  if(globals.Num_Side_Set > 0) {

    std::vector<INT> glob_ss_df_cnts(globals.Num_Side_Set);

    bytes_out += 3*globals.Num_Side_Set*sizeof(INT);
    tt1 = second();

    if(ex_put_ss_param_global(mesh_exoid, Side_Set_Ids,
                              Num_Elems_In_SS, TOPTR(glob_ss_df_cnts)) < 0) {
      fprintf(stderr,
              "[%s]: ERROR, unable to output global side-set params\n",
              yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    PIO_Time_Array[5]  = (second() - tt1);
    total_out_time    += PIO_Time_Array[5];
  }

  bytes_out += globals.Num_Elem_Blk*sizeof(INT);
  tt1 = second();

  if(ex_put_eb_info_global(mesh_exoid, Elem_Blk_Ids, Num_Elems_In_EB) < 0) {
    fprintf(stderr,
            "[%s]: ERROR, unable to output global elem blk IDs\n",
            yo);
    ex_close(mesh_exoid);
    exit(1);
  }

  PIO_Time_Array[6]  = (second() - tt1);
  total_out_time    += PIO_Time_Array[6];
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  if (iproc == 0) {
    /* Generate a QA record for the utility */
    time_t  time_val = time(NULL);
    char *ct_ptr   = asctime(localtime(&time_val));
    char tm_date[30];
    strcpy(tm_date, ct_ptr);

    /* Break string with null characters */
    tm_date[3]  = '\0';
    tm_date[7]  = '\0';
    tm_date[10] = '\0';
    tm_date[19] = '\0';

    char qa_date[MAX_STR_LENGTH+1];
    char qa_time[MAX_STR_LENGTH+1];
    char qa_name[MAX_STR_LENGTH+1];
    char qa_vers[MAX_STR_LENGTH+1];
    sprintf(qa_date, "%s %s %s", &tm_date[8], &tm_date[4], &tm_date[20]);
    sprintf(qa_time, "%s", &tm_date[11]);
    strcpy(qa_name, UTIL_NAME);
    strcpy(qa_vers, VER_STR);

    if(qa_date[strlen(qa_date)-1] == '\n')
      qa_date[strlen(qa_date)-1] = '\0';

    if (globals.Num_QA_Recs > 0) {
      strcpy(globals.QA_Record[4*(globals.Num_QA_Recs-1)],     qa_name);
      strcpy(globals.QA_Record[(4*(globals.Num_QA_Recs-1))+1], qa_vers);
      strcpy(globals.QA_Record[(4*(globals.Num_QA_Recs-1))+2], qa_date);
      strcpy(globals.QA_Record[(4*(globals.Num_QA_Recs-1))+3], qa_time);

      /* Output QA records to screen */
      if(Debug_Flag >= 4) {
	printf("Number of QA records: %d\n", globals.Num_QA_Recs);
	if(Debug_Flag >= 6) {
	  printf("QA Records:\n");
	  for(int i1=0; i1 < 4*(globals.Num_QA_Recs); i1++) {
	    printf("\t%s\n", globals.QA_Record[i1]);
	  }
	}
      }
    }

    /* Output the QA and Info records */
    for(int i1=0; i1 < 4*globals.Num_QA_Recs; i1++)
      bytes_out += (MAX_STR_LENGTH+MAX_LINE_LENGTH)*sizeof(char);

    tt1 = second();

    if(ex_put_qa(mesh_exoid, globals.Num_QA_Recs,
		 (char *(*)[4]) &globals.QA_Record[0]) < 0) {
      fprintf(stderr, "[%s]: ERROR Could not put QA records\n",
	      yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    if(globals.Num_Info_Recs > 0) {
      if(Debug_Flag >= 4) {
	printf("Number of info records: %d\n", globals.Num_Info_Recs);
      }

      if(ex_put_info(mesh_exoid, globals.Num_Info_Recs, globals.Info_Record) < 0) {
	fprintf(stderr, "[%s]: ERROR Could not put Info records\n",
		yo);
	ex_close(mesh_exoid);
	exit(1);
      }
    }
  }
  PIO_Time_Array[8]  = (second() - tt1);
  total_out_time    += PIO_Time_Array[8];

  /* Output the coordinate frame information (if any). This puts the
     file in/out of define mode, so should be early in the write stage
  */
  if(Debug_Flag >= 4) {
    printf("Number of Coordinate Frames: %d\n", globals.Num_Coordinate_Frames);
  }

  if (globals.Num_Coordinate_Frames > 0) {
    T *Coordinate_Frame_Coordinates = NULL;
    Coordinate_Frame_Coordinates = globals.Coordinate_Frame_Coordinates;
    if (ex_put_coordinate_frames(mesh_exoid, globals.Num_Coordinate_Frames,
				 globals.Coordinate_Frame_Ids,
				 Coordinate_Frame_Coordinates,
				 globals.Coordinate_Frame_Tags) < 0) {
      fprintf(stderr,
	      "[%s]: ERROR, Unable to put coordinate frame data in parallel mesh file\n",
	      yo);
      ex_close(mesh_exoid);
      exit(1);
    }
  }

  char    cTitle[MAX_LINE_LENGTH+1];

  if (proc_for == 0)
    strcpy(cTitle, GeomTitle);
  else
    sprintf(cTitle, "Parallel Mesh File for Processor %d", proc_for);

  /* Output the initial information to the parallel Exodus file(s) */
  bytes_out += strlen(cTitle)*sizeof(char) + 6*sizeof(INT);
  tt1 = second();

  if(Debug_Flag >= 4) {
    printf("Putting init info in file id: %d\n", mesh_exoid);
    printf("\tTitle: %s\n", cTitle);
    printf("\tNumber Dimensions: %d\n", globals.Num_Dim);
    printf("\tNumber Nodes: %lu\n", itotal_nodes);
    printf("\tNumber Elements: %lu\n",
           (size_t)globals.Num_Internal_Elems[iproc]+(size_t)globals.Num_Border_Elems[iproc]);
    printf("\tNumber Element Blocks: %lu\n", (size_t)globals.Num_Elem_Blk);
    printf("\tNumber Node Sets: %lu\n", (size_t)globals.Num_Node_Set);
    printf("\tNumber Side Sets: %lu\n", (size_t)globals.Num_Side_Set);
  }

  if(ex_put_init(mesh_exoid, cTitle, globals.Num_Dim, itotal_nodes,
                 globals.Num_Internal_Elems[iproc]+globals.Num_Border_Elems[iproc],
                 globals.Num_Elem_Blk,
                 globals.Num_Node_Set,
                 globals.Num_Side_Set) < 0) {
    fprintf(stderr,
	    "[%s]: ERROR, Unable to put initial info in parallel mesh file\n",
            yo);
    ex_close(mesh_exoid);
    exit(1);
  }

  /* Output names... */
  if (globals.Num_Elem_Blk > 0) {
    error = ex_put_names(mesh_exoid, EX_ELEM_BLOCK, &Elem_Blk_Names[0]);
    check_exodus_error(error, "ex_put_names");
  }
  if (globals.Num_Node_Set > 0) {
    error = ex_put_names(mesh_exoid, EX_NODE_SET, &Node_Set_Names[0]);
    check_exodus_error(error, "ex_put_names");
  }
  if (globals.Num_Side_Set > 0) {
    error = ex_put_names(mesh_exoid, EX_SIDE_SET, &Side_Set_Names[0]);
    check_exodus_error(error, "ex_put_names");
  }
  
  PIO_Time_Array[7]  = (second() - tt1);
  total_out_time    += PIO_Time_Array[7];

  /* Assign the coordinates to the coord_vector */
  T *x_coord=NULL;
  T *y_coord=NULL;
  T *z_coord=NULL;
  switch(globals.Num_Dim) {
  case 3:
    z_coord = globals.Coor[iproc][2];
    /* FALLTHROUGH */
  case 2:
    y_coord = globals.Coor[iproc][1];
    /* FALLTHROUGH */
  case 1:
    x_coord = globals.Coor[iproc][0];

    break;
  }

  /* Output the coordinates to the parallel Exodus file */
  bytes_out += globals.Num_Dim * itotal_nodes * io_ws;
  tt1 = second();

  if(Debug_Flag >= 4) {
    printf("Putting coordinate info in file id: %d\n",
           mesh_exoid);
  }
  if(ex_put_coord(mesh_exoid, x_coord, y_coord, z_coord) < 0) {
    fprintf(stderr,
            "[%s]: ERROR, could not write out nodal coordinates\n",
            yo);
    ex_close(mesh_exoid);
    exit(1);
  }

  PIO_Time_Array[9]   = (second() - tt1);
  total_out_time     += PIO_Time_Array[9];

  bytes_out += itotal_nodes*sizeof(INT);
  tt1 = second();

  if(ex_put_processor_node_maps(mesh_exoid, nem_node_mapi, nem_node_mapb,
				nem_node_mape, proc_for) < 0) {
    fprintf(stderr,
            "[%s]: ERROR, could not write Nemesis nodal number map!\n",
            yo);
    check_exodus_error(ex_close(mesh_exoid), "ex_close");
    ex_close(mesh_exoid);
    exit(1);
  }

  /* If non-NULL, output the global node id map which preserves
     the global node ids in the original mesh */
  if (globals.Proc_Global_Node_Id_Map[iproc] != NULL) {
    bytes_out += itotal_nodes * sizeof(INT);
    if (ex_put_map_param(mesh_exoid, 1, 0) < 0) {
      fprintf(stderr, "[%s]: ERROR, unable to define global node map parameters!\n",
	      yo);
      ex_close(mesh_exoid);
      exit(1);
    }
    
    if (ex_put_num_map(mesh_exoid, EX_NODE_MAP, 1, globals.Proc_Global_Node_Id_Map[iproc]) < 0) {
      fprintf(stderr, "[%s]: ERROR, unable to output global node id map!\n",
	      yo);
      ex_close(mesh_exoid);
      exit(1);
    }
    if (ex_put_name(mesh_exoid, EX_NODE_MAP, 1, "original_global_id_map") < 0) {
      fprintf(stderr, "[%s]: ERROR, unable to define global node map name!\n",
	      yo);
      ex_close(mesh_exoid);
      exit(1);
    }
  }
    
  PIO_Time_Array[11]  = (second() - tt1);
  total_out_time     += PIO_Time_Array[11];

  safe_free((void **) &nem_node_mapi);

  PIO_Time_Array[13] = 0.0;
  PIO_Time_Array[14] = 0.0;
  PIO_Time_Array[15] = 0.0;

  /*
   * Generate a list of the global element blocks, some of which may be
   * NULL on a given processor.
   */
  {
    /*
     * Start a local block so we can define some variables locally
     * here instead of a few thousand lines up from here...
     */
    std::vector<INT> EB_Ids(globals.Num_Elem_Blk);
    std::vector<INT> EB_Cnts(globals.Num_Elem_Blk);
    std::vector<INT> EB_NperE(globals.Num_Elem_Blk);
    std::vector<INT> EB_Nattr(globals.Num_Elem_Blk);

    char **EB_Types = (char**)array_alloc(__FILE__, __LINE__, 2, globals.Num_Elem_Blk,
					  MAX_STR_LENGTH+1, sizeof(char));
    if (!EB_Types) {
      fprintf(stderr, "%s: fatal: insufficient memory\n", yo);
      exit(1);
    }
    
    for (int i1=0; i1 < globals.Num_Elem_Blk; i1++) {
      memset(EB_Types[i1], '\0', (MAX_STR_LENGTH+1));
    }

    int cnt = 0;
    for(int i1=0; i1 < globals.Num_Elem_Blk; i1++) {
      bool ifound = false;
      for(int i2=0; i2 < globals.Proc_Num_Elem_Blk[iproc]; i2++) {
	if(globals.Proc_Elem_Blk_Ids[iproc][i2] == Elem_Blk_Ids[i1]) {
	  ifound = true;
	  break;
	}
      }
      
      /*
       * If it's not found then this element block is null on the current
       * processor.
       */
      if(!ifound) {
	int iblk = globals.Proc_Num_Elem_Blk[iproc]+cnt;
	globals.Proc_Elem_Blk_Ids[iproc][iblk] =  Elem_Blk_Ids[i1];
	globals.Proc_Num_Elem_In_Blk[iproc][iblk] = 0;
	globals.Proc_Nodes_Per_Elem[iproc][iblk]  = 0;
	globals.Proc_Num_Attr[iproc][iblk]        = 0;
	cnt++;
      }
    }

    /* Output the elemental block(s). */
    for(int i1=0; i1 < globals.Num_Elem_Blk; i1++) {

      ex_entity_id iglobal_blk = Elem_Blk_Ids[i1];

      /* Find the local element block index */
      int ilocal;
      for(ilocal=0; ilocal < globals.Num_Elem_Blk; ilocal++) {
	if(globals.Proc_Elem_Blk_Ids[iproc][ilocal] == iglobal_blk)
	  break;
      }

      /* Error check */
      if(ilocal >= globals.Num_Elem_Blk) {
	fprintf(stderr, "[%s]: Error finding local element block ID\n",
		yo);
	exit(1);
      }
    
      /* If it's a non-null block output all information */
      EB_Ids[i1]   = iglobal_blk;
      if (ilocal < globals.Proc_Num_Elem_Blk[iproc]) {

	/* Generate the ExodusII element name */
	strncpy(EB_Types[i1], Elem_Blk_Types[i1], MAX_STR_LENGTH);
	EB_Types[i1][MAX_STR_LENGTH] = '\0';

	EB_Cnts[i1]  = globals.Proc_Num_Elem_In_Blk[iproc][ilocal];
	EB_NperE[i1] = globals.Proc_Nodes_Per_Elem[iproc][ilocal];
	EB_Nattr[i1] = globals.Proc_Num_Attr[iproc][ilocal];
      }
    }
    if(Debug_Flag >= 4) {
      printf("Putting concat_elem_block info in file id: %d\n",
	     mesh_exoid);
    }
    error = ex_put_concat_elem_block(mesh_exoid,
				     &EB_Ids[0], &EB_Types[0],
				     &EB_Cnts[0], &EB_NperE[0],
				     &EB_Nattr[0], 1);
    check_exodus_error(error, "ex_put_concat_elem_block");

    safe_free((void **) &EB_Types);

    /* Output attribute names for each element block */
    for(int i1=0; i1 < globals.Num_Elem_Blk; i1++) {

      ex_entity_id iglobal_blk = Elem_Blk_Ids[i1];

      /* Find the local element block index */
      int ilocal;
      for (ilocal=0; ilocal < globals.Num_Elem_Blk; ilocal++) {
	if(globals.Proc_Elem_Blk_Ids[iproc][ilocal] == iglobal_blk)
	  break;
      }

      /* If it's a non-null block output attribute name information */
      if (ilocal < globals.Proc_Num_Elem_Blk[iproc]) {
	if (globals.Proc_Num_Attr[iproc][ilocal] > 0) {
	  if (ex_put_attr_names(mesh_exoid, EX_ELEM_BLOCK, Elem_Blk_Ids[i1],
				Elem_Blk_Attr_Names[i1]) < 0) {
	    fprintf(stderr,
		    "[%s]: ERROR, could not write Exodus attribute names!\n",
		    yo);
	    ex_close(mesh_exoid);
	    exit(1);
	  }
	}
      }
    }

  
    /* Reset globals.GNodes to start at 1 instead of 0 */
    for(size_t i1=0; i1 < itotal_nodes; (globals.GNodes[iproc][i1++])++);

    /* Output the Exodus node number map */
    bytes_out += itotal_nodes*sizeof(INT);
    tt1 = second();

    if(Debug_Flag >= 4) {
      printf("Putting node_num_map in file id: %d\n",
	     mesh_exoid);
    }
    if(ex_put_id_map(mesh_exoid, EX_NODE_MAP, globals.GNodes[iproc]) < 0) {
      fprintf(stderr,
	      "[%s]: ERROR, could not write Exodus node number map!\n",
	      yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    PIO_Time_Array[10]  = (second() - tt1);
    total_out_time     += PIO_Time_Array[10];

    /*
     * Allocate memory for the elemental map. Currently this map is assigned
     * as a linear array since it is not really used.
     */
    INT *iElem_Map = (INT *)array_alloc(__FILE__, __LINE__, 1,
					globals.Num_Internal_Elems[iproc] +
					globals.Num_Border_Elems[iproc],
					sizeof(INT));
    for(INT i1=0; i1 < globals.Num_Internal_Elems[iproc]+globals.Num_Border_Elems[iproc]; i1++)
      iElem_Map[i1] = globals.GElems[iproc][i1] + 1;

    bytes_out += 2 * globals.Num_Internal_Elems[iproc] * globals.Num_Border_Elems[iproc] *
      sizeof(INT);
    tt1 = second();

    if(Debug_Flag >= 4) {
      printf("Putting elem_num_map info in file id: %d\n",
	     mesh_exoid);
    }
    if(ex_put_id_map(mesh_exoid, EX_ELEM_MAP, iElem_Map) < 0) {
      fprintf(stderr, "[%s]: ERROR, unable to output element map\n",
	      yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    /* If non-NULL, output the global element id map which preserves
       the global element ids in the original mesh */
    if (globals.Proc_Global_Elem_Id_Map[iproc] != NULL) {
      bytes_out += globals.Num_Internal_Elems[iproc] * globals.Num_Border_Elems[iproc] * sizeof(INT);
      if (ex_put_map_param(mesh_exoid, 0, 1) < 0) {
	fprintf(stderr, "[%s]: ERROR, unable to define global map parameters!\n",
		yo);
	ex_close(mesh_exoid);
	exit(1);
      }
      
      if (ex_put_num_map(mesh_exoid, EX_ELEM_MAP, 1, globals.Proc_Global_Elem_Id_Map[iproc]) < 0) {
	fprintf(stderr, "[%s]: ERROR, unable to output global id map!\n",
		yo);
	ex_close(mesh_exoid);
	exit(1);
      }
      if (ex_put_name(mesh_exoid, EX_ELEM_MAP, 1, "original_global_id_map") < 0) {
	fprintf(stderr, "[%s]: ERROR, unable to define global map name!\n",
		yo);
	ex_close(mesh_exoid);
	exit(1);
      }

    }
      
    /* Also output the Nemesis element map */
    if(ex_put_processor_elem_maps(mesh_exoid, globals.Elem_Map[iproc],
				  globals.Elem_Map[iproc] + globals.Num_Internal_Elems[iproc],
				  proc_for) < 0) {
      fprintf(stderr, "[%s]: ERROR, unable to output nemesis element map!\n",
	      yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    PIO_Time_Array[12]  = (second() - tt1);
    total_out_time     += PIO_Time_Array[12];
    safe_free((void **) &iElem_Map);

    for(int i1=0; i1 < globals.Num_Elem_Blk; i1++) {

      ex_entity_id iglobal_blk = Elem_Blk_Ids[i1];

      /* Find the local element block index */
      int ilocal;
      for(ilocal=0; ilocal < globals.Num_Elem_Blk; ilocal++) {
	if(globals.Proc_Elem_Blk_Ids[iproc][ilocal] == iglobal_blk)
	  break;
      }

      /* Error check */
      if(ilocal >= globals.Num_Elem_Blk) {
	fprintf(stderr, "[%s]: Error finding local element block ID\n",
		yo);
	exit(1);
      }

      /* If it's a non-null block output all information */
      if(ilocal < globals.Proc_Num_Elem_Blk[iproc]) {
	/* Find the first index into the connectivity for this block */
	size_t iIndex0 = 0;
	for(int i2=0; i2 < ilocal; i2++) {
	  iIndex0 += globals.Proc_Num_Elem_In_Blk[iproc][i2] *
	    globals.Proc_Nodes_Per_Elem[iproc][i2];
	}

	PIO_Time_Array[13] += (second() - tt1);

	size_t tmp_cnt = globals.Proc_Num_Elem_In_Blk[iproc][ilocal] *
	  globals.Proc_Nodes_Per_Elem[iproc][ilocal];
      
	/* Generate the connectivity array for local node numbering */
	INT *proc_local_conn = (INT *) array_alloc(__FILE__, __LINE__, 1, tmp_cnt,
						   sizeof(INT));

	reverse_map(&globals.Proc_Elem_Connect[iproc][iIndex0], 1, tmp_cnt,
		    &globals.GNodes[iproc][0], loc_index, itotal_nodes,
		    proc_local_conn);
      
	bytes_out += globals.Proc_Nodes_Per_Elem[iproc][ilocal] *
	  globals.Proc_Num_Elem_In_Blk[iproc][ilocal]*sizeof(INT);
	tt1 = second();

	if(Debug_Flag >= 4) {
	  printf("Putting element_connectivity info in file id: %d\n",
		 mesh_exoid);
	}
	if(ex_put_conn(mesh_exoid, EX_ELEM_BLOCK, globals.Proc_Elem_Blk_Ids[iproc][ilocal],
		       proc_local_conn, NULL, NULL) < 0) {
	  fprintf(stderr, "[%s]: ERROR, unable to output connectivity\n",
		  yo);
	  ex_close(mesh_exoid);
	  exit(1);
	}

	PIO_Time_Array[14] += (second() - tt1);

	/* Free up the locally numbered connectivity */
	safe_free((void **) &proc_local_conn);

	if(globals.Proc_Num_Attr[iproc][ilocal] > 0) {

	  /* Find the first index into the attribute list for this block */
	  size_t iIndex1 = 0;
	  for(int i2=0; i2 < ilocal; i2++) {
	    iIndex1 += globals.Proc_Num_Attr[iproc][i2] *
	      globals.Proc_Num_Elem_In_Blk[iproc][i2];
	  }

	  bytes_out += globals.Proc_Num_Elem_In_Blk[iproc][ilocal] * io_ws;
	  tt1 = second();

	  T *ptr = &(globals.Proc_Elem_Attr[iproc][iIndex1]);

	  if(ex_put_attr(mesh_exoid, EX_ELEM_BLOCK, globals.Proc_Elem_Blk_Ids[iproc][ilocal],
			 ptr) < 0) {
	    fprintf(stderr,
		    "[%s]: ERROR, unable to output element attributes\n",
		    yo);
	    exit(1);
	  }

	  PIO_Time_Array[15] += (second() - tt1);
	  total_out_time     += PIO_Time_Array[15];
	}

      } /* End "if(ilocal < globals.Num_Elem_Blk[iproc])" */

    } /* End "for(i1=0; i1 < globals.Num_Elem_Block; i1++)" */
  }
  total_out_time += (PIO_Time_Array[13] + PIO_Time_Array[14] +
                     PIO_Time_Array[15]);

  /*
   * Write out the node-set information. Note that the value of the
   * node-set distribution factor is not currently used so only a
   * dummy set is output for this value.
   */
  INT iMaxLen = 0;
  for(int i1=0; i1 < globals.Proc_Num_Node_Sets[iproc]; i1++)
    iMaxLen = PEX_MAX(globals.Proc_NS_Count[iproc][i1], iMaxLen);

  /* Renumber Node set node lists to use local node numbers */
  INT *proc_local_ns = NULL;
  if(globals.Proc_Num_Node_Sets[iproc] > 0) {
    proc_local_ns = (INT *)array_alloc(__FILE__, __LINE__, 1,
                                       globals.Proc_NS_List_Length[iproc],
                                       sizeof(INT));

    reverse_map(&globals.Proc_NS_List[iproc][0], 1,
		globals.Proc_NS_List_Length[iproc],
		&globals.GNodes[iproc][0], loc_index, itotal_nodes,
		proc_local_ns);
  }

  safe_free((void **) &loc_index);

  PIO_Time_Array[16] = 0.0;
  PIO_Time_Array[17] = 0.0;

  /* Fill in the information for the NULL node sets */
  size_t cnt = 0;
  for(int i1=0; i1 < globals.Num_Node_Set; i1++) {
    bool ifound = false;
    for(int i2=0; i2 < globals.Proc_Num_Node_Sets[iproc]; i2++) {
      if(globals.Proc_NS_Ids[iproc][i2] == Node_Set_Ids[i1]) {
        ifound = true;
        break;
      }
    }

    if (!ifound) {
      globals.Proc_NS_Ids[iproc][globals.Proc_Num_Node_Sets[iproc]+cnt] = Node_Set_Ids[i1];
      globals.Proc_NS_Count[iproc][globals.Proc_Num_Node_Sets[iproc]+cnt] = 0;
      globals.Proc_NS_DF_Count[iproc][globals.Proc_Num_Node_Sets[iproc]+cnt] = 0;
      cnt++;
    }
  }

  tt1 = second();
  if (globals.Num_Node_Set > 0) {
    size_t dcount = 0;
    size_t ncount = 0;
    for (int i1=0; i1 < globals.Num_Node_Set; i1++) {
      dcount += globals.Proc_NS_DF_Count[iproc][i1];
      ncount += globals.Proc_NS_Count[iproc][i1];
    }

    std::vector<INT> conc_ids(globals.Num_Node_Set);
    std::vector<INT> conc_nodes(globals.Num_Node_Set);
    std::vector<INT> conc_df(globals.Num_Node_Set);
    std::vector<INT> conc_nind(globals.Num_Node_Set);
    std::vector<INT> conc_dind(globals.Num_Node_Set);
    std::vector<INT> conc_nlist(ncount);
    std::vector<T>   conc_sdf(dcount);
    
    ncount = 0;
    dcount = 0;
    for(int i1=0; i1 < globals.Num_Node_Set; i1++) {

      /* Find the local ID */
      int i2 = 0;
      for(i2 = 0; i2 < globals.Num_Node_Set; i2++) {
	if(globals.Proc_NS_Ids[iproc][i2] == Node_Set_Ids[i1])
	  break;
      }

      conc_ids[i1] = globals.Proc_NS_Ids[iproc][i2];
      conc_nodes[i1] = globals.Proc_NS_Count[iproc][i2];
      conc_df[i1] = globals.Proc_NS_DF_Count[iproc][i2];
      
      conc_nind[i1] = ncount;
      for (int i3=0; i3 < globals.Proc_NS_Count[iproc][i2]; i3++) {
	conc_nlist[ncount++] = proc_local_ns[globals.Proc_NS_Pointers[iproc][i2]+i3];
      }

      conc_dind[i1] = dcount;
      for (INT i3=0; i3 < globals.Proc_NS_DF_Count[iproc][i2]; i3++) {
	conc_sdf[dcount++] = globals.Proc_NS_Dist_Fact[iproc][globals.Proc_NS_Pointers[iproc][i2]+i3];
      }
    }

    ex_put_concat_node_sets(mesh_exoid, TOPTR(conc_ids), TOPTR(conc_nodes), TOPTR(conc_df),
			    TOPTR(conc_nind), TOPTR(conc_dind), TOPTR(conc_nlist), TOPTR(conc_sdf));
  }
  total_out_time += second() - tt1;

  /* Free local number array */
  if(globals.Proc_Num_Node_Sets[iproc] > 0)
    safe_free((void **) &proc_local_ns);

  /* Renumber element SS to use local element numbers */
  INT *proc_local_ss = NULL;
  if(globals.Proc_Num_Side_Sets[iproc] > 0) {
    proc_local_ss = (INT *)array_alloc(__FILE__, __LINE__, 1,
                                       globals.Proc_SS_Elem_List_Length[iproc],
                                       sizeof(INT));
    reverse_map(&globals.Proc_SS_Elem_List[iproc][0], 0,
		globals.Proc_SS_Elem_List_Length[iproc],
		&globals.GElems[iproc][0], (INT*)NULL,
		globals.Num_Internal_Elems[iproc]+globals.Num_Border_Elems[iproc],
		proc_local_ss);
  }

  PIO_Time_Array[18] = 0.0;
  PIO_Time_Array[19] = 0.0;

  /* Set up the null side sets */
  cnt = 0;
  for(int i1=0; i1 < globals.Num_Side_Set; i1++) {
    bool ifound = false;
    for(int i2=0; i2 < globals.Proc_Num_Side_Sets[iproc]; i2++) {
      if(globals.Proc_SS_Ids[iproc][i2] == Side_Set_Ids[i1]) {
        ifound = true;
        break;
      }
    }

    if (!ifound) {
      globals.Proc_SS_Ids[iproc][globals.Proc_Num_Side_Sets[iproc]+cnt] = Side_Set_Ids[i1];
      globals.Proc_SS_Elem_Count[iproc][globals.Proc_Num_Side_Sets[iproc]+cnt] = 0;
      globals.Proc_SS_DF_Count[iproc][globals.Proc_Num_Side_Sets[iproc]+cnt]   = 0;
      cnt++;
    }
  }

  /* Output concatenated sidesets.  For each processor need:
   * side_set ids          (size globals.Num_Side_Set)
   * num_side_per_set      (size globals.Num_Side_Set)
   * num_dist_per_set      (size globals.Num_Side_Set)
   * side_sets_elem_index  (size globals.Num_Side_Set)
   * side_sets_dist_index  (size globals.Num_Side_Set)
   * side_sets_elem_list
   * side_sets_side_list
   * side_sets_dist_fact
   */
  
  tt1 = second();
  if (globals.Num_Side_Set > 0) {
    int df_count = 0;
    int el_count = 0;
    for (int i1=0; i1 < globals.Num_Side_Set; i1++) {
      df_count += globals.Proc_SS_DF_Count[iproc][i1];
      el_count += globals.Proc_SS_Elem_Count[iproc][i1];
    }
    
    std::vector<INT> conc_ids(globals.Num_Side_Set);
    std::vector<INT> conc_sides(globals.Num_Side_Set);
    std::vector<INT> conc_dist(globals.Num_Side_Set);
    std::vector<INT> conc_eind(globals.Num_Side_Set);
    std::vector<INT> conc_dind(globals.Num_Side_Set);
    std::vector<INT> conc_elist(el_count);
    std::vector<INT> conc_slist(el_count);
    std::vector<T>   conc_sdflist(df_count);
    
    /* Fill in the arrays ... */
    df_count = 0;
    el_count = 0;
    for(int i1=0; i1 < globals.Num_Side_Set; i1++) {
      
      /* Find the local ID of this side set */
      int i2 = 0;
      for(i2=0; i2 < globals.Num_Side_Set; i2++) {
	if(globals.Proc_SS_Ids[iproc][i2] == Side_Set_Ids[i1])
	  break;
      }
      
      conc_ids[i1]   = globals.Proc_SS_Ids[iproc][i2];
      conc_sides[i1] = globals.Proc_SS_Elem_Count[iproc][i2];
      conc_dist[i1]  = globals.Proc_SS_DF_Count[iproc][i2];
      
      conc_eind[i1]  = el_count;
      for (INT i3 = 0; i3 < globals.Proc_SS_Elem_Count[iproc][i2]; i3++) {
	conc_elist[el_count] = proc_local_ss[globals.Proc_SS_Elem_Pointers[iproc][i2]+i3];
	conc_slist[el_count] = globals.Proc_SS_Side_List[iproc][globals.Proc_SS_Elem_Pointers[iproc][i2]+i3];
	el_count++;
      }

      conc_dind[i1] = df_count;
      for (INT i3 = 0; i3 < globals.Proc_SS_DF_Count[iproc][i2]; i3++) {
	conc_sdflist[df_count++] =
	  globals.Proc_SS_Dist_Fact[iproc][globals.Proc_SS_DF_Pointers[iproc][i2]+i3];
      }	  
    }
    ex_put_concat_side_sets(mesh_exoid, TOPTR(conc_ids), TOPTR(conc_sides), TOPTR(conc_dist),
			    TOPTR(conc_eind), TOPTR(conc_dind), TOPTR(conc_elist), TOPTR(conc_slist),
			    TOPTR(conc_sdflist));
  }
  PIO_Time_Array[19] += (second() - tt1);
  total_out_time += (PIO_Time_Array[18] + PIO_Time_Array[19]);

  /* Free unneeded memory */
  if(globals.Proc_Num_Side_Sets[iproc] > 0)
    safe_free((void **) &proc_local_ss);

  /*
   * Write out the name of the coordinate axes to the parallel ExodusII
   * files.
   */
  bytes_out += globals.Num_Dim * 8 * sizeof(char);
  tt1 = second();
  if(ex_put_coord_names(mesh_exoid, Coord_Name) < 0) {
    fprintf(stderr, "[%s]: ERROR, could not output coordinate names\n",
            yo);
    ex_close(mesh_exoid);
    exit(1);
  }
  PIO_Time_Array[20]  = (second() - tt1);
  total_out_time     += PIO_Time_Array[20];

  if (Restart_Info.Flag > 0) {

    tt1 = second();
    bytes_out += write_var_param(mesh_exoid, max_name_length,
				 Restart_Info.NVar_Glob, Restart_Info.GV_Name,
				 Restart_Info.NVar_Node, Restart_Info.NV_Name, 
				 Restart_Info.NVar_Elem, Restart_Info.EV_Name,  TOPTR(Restart_Info.GElem_TT),
				 Restart_Info.NVar_Nset, Restart_Info.NSV_Name, TOPTR(Restart_Info.GNset_TT),
				 Restart_Info.NVar_Sset, Restart_Info.SSV_Name, TOPTR(Restart_Info.GSset_TT));

    PIO_Time_Array[21] = (second() - tt1);
    total_out_time    += PIO_Time_Array[21];
  }

  /* Calculate the overall performance */
  tt1 = bytes_out;

  if( total_out_time == 0 )
    tt1 = 0;
  else
    tt1 = tt1/total_out_time;

  tt1 /= 1024.0;

  PIO_Time_Array[25] = tt1;

  return;

} /* END write_parExo_data() */

template <typename T, typename INT>
int NemSpread<T,INT>::write_var_param(int mesh_exoid, int max_name_length,
			       int num_glob, char **gv_names,
			       int num_node, char **nv_names, 
			       int num_elem, char **ev_names, int *local_ebtt,
			       int num_nset, char **ns_names, int *local_nstt,
			       int num_sset, char **ss_names, int *local_sstt)
{
  size_t bytes_out=0;
  int error;

  bytes_out += (5 + globals.Num_Elem_Blk * num_elem +
		    globals.Num_Side_Set * num_sset +
		    globals.Num_Node_Set * num_nset) * sizeof(INT);
  error = ex_put_all_var_param(mesh_exoid, num_glob, num_node,
			       num_elem, local_ebtt,
			       num_nset, local_nstt,
			       num_sset, local_sstt); 
  check_exodus_error(error, "ex_put_all_var_param");

  if (gv_names != NULL) {
    bytes_out += Restart_Info.NVar_Glob * max_name_length;
    error = ex_put_variable_names(mesh_exoid, EX_GLOBAL, num_glob, gv_names);
    check_exodus_error(error, "ex_put_var_names");
  }
  if (nv_names != NULL) {
    bytes_out += num_node * max_name_length;
    error = ex_put_variable_names(mesh_exoid, EX_NODAL, num_node, nv_names);
    check_exodus_error(error, "ex_put_var_names");
  }
  if (ev_names != NULL) {
    bytes_out += Restart_Info.NVar_Elem * max_name_length;
    error = ex_put_variable_names(mesh_exoid, EX_ELEM_BLOCK, num_elem, ev_names);
    check_exodus_error(error, "ex_put_var_names");
  }
  if (ns_names != NULL) {
    bytes_out += Restart_Info.NVar_Nset * max_name_length;
    error = ex_put_variable_names(mesh_exoid, EX_NODE_SET, num_nset, ns_names);
    check_exodus_error(error, "ex_put_var_names");
  }
  if (ss_names != NULL) {
    bytes_out += Restart_Info.NVar_Sset * max_name_length;
    error = ex_put_variable_names(mesh_exoid, EX_SIDE_SET, num_sset, ss_names);
    check_exodus_error(error, "ex_put_var_names");
  }
  return (bytes_out);
}

template void NemSpread<double,int>::write_var_timestep(int exoid, int proc, int time_step, 
							int *eb_ids_global, int *ss_ids_global, int *ns_ids_global);
template void NemSpread<float,int>::write_var_timestep(int exoid, int proc, int time_step, 
						       int *eb_ids_global, int *ss_ids_global, int *ns_ids_global);
template void NemSpread<double,int64_t>::write_var_timestep(int exoid, int proc, int time_step, 
							    int64_t *eb_ids_global, int64_t *ss_ids_global, int64_t *ns_ids_global);
template void NemSpread<float,int64_t>::write_var_timestep(int exoid, int proc, int time_step, 
							   int64_t *eb_ids_global, int64_t *ss_ids_global, int64_t *ns_ids_global);
  

template <typename T, typename INT>
void NemSpread<T,INT>::write_var_timestep(int exoid, int proc, int time_step, 
					  INT *eb_ids_global, INT *ss_ids_global, INT *ns_ids_global)
{
  int error;

  /* output the time */
  {
    T *var_ptr = (T *) &(Restart_Info.Time);
    error = ex_put_time(exoid, time_step, var_ptr);
    check_exodus_error(error, "ex_put_time");
  }

  /* start by outputting the global variables */
  if (Restart_Info.NVar_Glob > 0) {

    T *var_ptr = &Restart_Info.Glob_Vals[0];

    error = ex_put_glob_vars(exoid, time_step, Restart_Info.NVar_Glob,
                             var_ptr);

    check_exodus_error(error, "ex_put_glob_vars");
  }

  if (Restart_Info.NVar_Node > 0) {
    size_t num_nodes = globals.Num_Internal_Nodes[proc] + globals.Num_Border_Nodes[proc] +
      globals.Num_External_Nodes[proc];

    for (int var_num=0; var_num < Restart_Info.NVar_Node; var_num++) {

      size_t var_offset = var_num * num_nodes;


      T *var_ptr = &(Restart_Info.Node_Vals[proc][var_offset]);

      error = ex_put_nodal_var(exoid, time_step, (var_num+1), num_nodes,
                               var_ptr);

      check_exodus_error(error, "ex_put_nodal_var");
    }
  }

  if (Restart_Info.NVar_Elem > 0) {

    size_t num_elem = globals.Num_Internal_Elems[proc] + globals.Num_Border_Elems[proc];

    for (int var_num=0; var_num < Restart_Info.NVar_Elem; var_num++) {
      int eb_num_g = 0;

      size_t var_offset = var_num * num_elem;
      T *var_ptr = &(Restart_Info.Elem_Vals[proc][var_offset]);

      for (int eb_num=0; eb_num < globals.Proc_Num_Elem_Blk[proc]; eb_num++) {

        /* now I have to find the appropriate entry in the truth table */
	/* can always assume this eb num is greater than the last one */
	for (int cnt1=eb_num_g; cnt1<globals.Num_Elem_Blk; cnt1++) {
	  if (globals.Proc_Elem_Blk_Ids[proc][eb_num] == eb_ids_global[cnt1]) {
	    eb_num_g = cnt1;
	    break;
	  }
	}

        if (Restart_Info.GElem_TT[eb_num_g*Restart_Info.NVar_Elem+var_num]) {
	  
          error = ex_put_var(exoid, time_step, EX_ELEM_BLOCK, (var_num+1),
                                  globals.Proc_Elem_Blk_Ids[proc][eb_num],
                                  globals.Proc_Num_Elem_In_Blk[proc][eb_num], var_ptr);
	  
          check_exodus_error(error, "ex_put_elem_var");
	}
	/* and now move the variable pointer */

	/* Note that the offsetting here must match the 'var_offset'
	 * treatment in ps_restart.c function read_elem_vars.
	 * Currently, the offset is applied even if the variable does
	 * not exist on a particular block.
	 */
	var_ptr += globals.Proc_Num_Elem_In_Blk[proc][eb_num];
      }
    }
  }

  if (Restart_Info.NVar_Sset > 0) {
    int ss_num_g = 0;
    size_t num_elem = globals.Proc_SS_Elem_List_Length[proc];
    for (int var_num=0; var_num < Restart_Info.NVar_Sset; var_num++) {

      size_t var_offset = var_num * num_elem;
      T *var_ptr = &(Restart_Info.Sset_Vals[proc][var_offset]);

      for (int ss_num=0; ss_num < globals.Proc_Num_Side_Sets[proc]; ss_num++) {

        /* now I have to find the appropriate entry in the truth table */
	for (int cnt1=0; cnt1<globals.Num_Side_Set; cnt1++) {
	  if (globals.Proc_SS_Ids[proc][ss_num] == ss_ids_global[cnt1]) {
	    ss_num_g = cnt1;
	    break;
	  }
	}
	assert(globals.Proc_SS_Ids[proc][ss_num] == ss_ids_global[ss_num_g]);

        if (Restart_Info.GSset_TT[ss_num_g*Restart_Info.NVar_Sset+var_num]) {
	  
          error = ex_put_var(exoid, time_step, EX_SIDE_SET, (var_num+1),
                                  globals.Proc_SS_Ids[proc][ss_num],
                                  globals.Proc_SS_Elem_Count[proc][ss_num], var_ptr);
	  
          check_exodus_error(error, "ex_put_sset_var");
	}
	/* and now move the variable pointer */

	/* Note that the offsetting here must match the 'var_offset'
	 * treatment in ps_restart.c function read_elem_vars.
	 * Currently, the offset is applied even if the variable does
	 * not exist on a particular block.
	 */
	var_ptr += globals.Proc_SS_Elem_Count[proc][ss_num];
      }
    }
  }

  if (Restart_Info.NVar_Nset > 0) {
    int ns_num_g = 0;
    size_t num_elem = globals.Proc_NS_List_Length[proc];
    for (int var_num=0; var_num < Restart_Info.NVar_Nset; var_num++) {

      size_t var_offset = var_num * num_elem;
      T *var_ptr = &(Restart_Info.Nset_Vals[proc][var_offset]);

      for (int ns_num=0; ns_num < globals.Proc_Num_Node_Sets[proc]; ns_num++) {

        /* now I have to find the appropriate entry in the truth table */
	for (int cnt1=0; cnt1<globals.Num_Node_Set; cnt1++) {
	  if (globals.Proc_NS_Ids[proc][ns_num] == ns_ids_global[cnt1]) {
	    ns_num_g = cnt1;
	    break;
	  }
	}
	assert(globals.Proc_NS_Ids[proc][ns_num] == ns_ids_global[ns_num_g]);

        if (Restart_Info.GNset_TT[ns_num_g*Restart_Info.NVar_Nset+var_num]) {
	  
          error = ex_put_var(exoid, time_step, EX_NODE_SET, (var_num+1),
                                  globals.Proc_NS_Ids[proc][ns_num],
                                  globals.Proc_NS_Count[proc][ns_num], var_ptr);
	  
          check_exodus_error(error, "ex_put_nset_var");
	}
	/* and now move the variable pointer */

	/* Note that the offsetting here must match the 'var_offset'
	 * treatment in ps_restart.c function read_elem_vars.
	 * Currently, the offset is applied even if the variable does
	 * not exist on a particular block.
	 */
	var_ptr += globals.Proc_NS_Count[proc][ns_num];
      }
    }
  }
}

namespace {
template <typename INT>
void reverse_map(INT *global, int p01, size_t gsize,
		 INT *glmap, INT *index, size_t msize,
		 INT *mapout)
{
  /*
   * The 'global' array is an array of node or element numbers
   * in the global id space.  It needs to be converted to local
   * numbers via the 'glmap' array.  The glmap array is sorted
   * by the 'index' array.  The map from global to local is
   * glmap[local_id] = global_id
   *
   * The 'p01' is either 0 or 1 and is an offset to the values
   * in 'global'
   */
  
  /*
   * The algorithm used is to sort the 'global' array via the
   * 'tmp_index' array (global[tmp_index[0..gsize]] is sorted)
   * Then, progress through the 'global' array in sorted order
   * and find the location in 'glmap'.  Note that since both are
   * sorted, it should be easy to progress sequentially through
   * both arrays.
   */
  
  INT *tmp_index = (INT *) array_alloc(__FILE__, __LINE__, 1, gsize,
				       sizeof(INT));

  /* Initialize index array */
  for (size_t i2=0; i2 < gsize; i2++)
    tmp_index[i2] = i2;
      
  /* Sort the 'global' array via the index array 'tmp_index' */
  gds_iqsort(global, tmp_index, gsize);

  size_t i3 = 0;
  if (index != NULL) {
    for(size_t i2 = 0; i2 < gsize; i2++) {
      INT gval = global[tmp_index[i2]] + p01;
      
      while (glmap[index[i3]] < gval)
	i3++;

      assert(glmap[index[i3]] == gval);

      mapout[tmp_index[i2]] = index[i3] + 1;
    }
  } else {
    for(size_t i2 = 0; i2 < gsize; i2++) {
      INT gval = global[tmp_index[i2]] + p01;

      while (glmap[i3] < gval)
	i3++;

      assert(glmap[i3] == gval);

      mapout[tmp_index[i2]] = i3 + 1;
    }
  }
  safe_free((void **) &tmp_index);
}
}
