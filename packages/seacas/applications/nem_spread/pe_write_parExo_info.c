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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "netcdf.h"
#include "exodusII.h"
#include "ne_nemesisI.h"

#include "sort_utils.h"
#include "rf_salsa.h"

#include "pe_common.h"

#include "el_geom_const.h"
#include "el_elm.h"
#include "rf_message.h"

#include "rf_allo.h"
#include "rf_io_const.h"
#include "rf_mp_const.h"

#include "ps_pario_const.h"

static int write_var_param(int mesh_exoid, int max_name_length, 
			   int num_glob, char **gv_names,
			   int num_node, char **nv_names, 
                           int num_elem, char **ev_names, int *local_ebtt,
			   int num_nset, char **ns_names, int *local_nstt,
			   int num_sset, char **ss_names, int *local_sstt);

static void reverse_map(int *global, int p01, int gsize,
			int *glmap, int *index, int msize,
			int *mapout);
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

void write_parExo_data(int mesh_exoid, int max_name_length,
                       int iproc,
                       int io_ws,
                       int num_qa_records,
                       char **qa_record,
                       int num_info_records,
                       char **info_record,
                       char **Elem_Blk_Types,
                       int *Node_Set_Ids,
                       int *Side_Set_Ids,
                       int *Elem_Blk_Ids,
                       int *Num_Nodes_In_NS,
                       int *Num_Elems_In_SS,
                       int *Num_Elems_In_EB,
                       ELEM_COMM_MAP *E_Comm_Map,
                       NODE_COMM_MAP *N_Comm_Map,
		       char **Node_Set_Names,
		       char **Side_Set_Names,
		       char **Elem_Blk_Names,
		       char ***Elem_Blk_Attr_Names)
{
  static   char yo[] = "write_parExo_data";

  /* For QA record */
  time_t  time_val;
  char   *ct_ptr, tm_date[30];
  char    qa_date[MAX_STR_LENGTH+1], qa_time[MAX_STR_LENGTH+1], qa_name[MAX_STR_LENGTH+1];
  char    qa_vers[MAX_STR_LENGTH+1];

  extern int Num_Side_Set, Num_Node_Set, Num_Elem_Blk;

  int     itotal_nodes, itotal_elems, *iElem_Map, iglobal_blk, iMaxLen;
  int     ilocal, proc_for;
  int     i0, i1, i2, iIndex[2]={0,0};
  int     num_proc_for;
  int    *nem_node_mapi=NULL, *nem_node_mapb=NULL, *nem_node_mape=NULL;
  int    *nem_elem_mapi=NULL, *nem_elem_mapb=NULL;

  void   *x_coord=NULL, *y_coord=NULL, *z_coord=NULL;

  void   *ptr;

  char    cTitle[MAX_LINE_LENGTH+1];

  int    *proc_local_conn, *glob_ns_df_cnts, *glob_ss_df_cnts;
  int    *loc_index;
  int    *proc_local_ns, *proc_local_ss, itemp;
  int    *local_tt, *local_nstt, *local_sstt;

  int     ifound, cnt, tmp_cnt;

  /* Local communication maps */
  int     ncomm_cnt;
  int    *n_comm_ids, *n_comm_ncnts;

  int     ecomm_cnt;
  int    *e_comm_ids, *e_comm_ecnts;

  NODE_COMM_MAP *n_comm_map;
  ELEM_COMM_MAP *e_comm_map;

  /* Performance metrics. */
  unsigned long    bytes_out=0;
  double           total_out_time=0.0, tt1;

  int error;
  /****************************BEGIN EXECUTION*********************************/

  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*                    PARALLEL EXODUSII SECTION                            */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


  if (io_ws < sizeof(float)) io_ws = sizeof(float);

  proc_for     = Proc_Ids[iproc];
  num_proc_for = Proc_Info[0];

  itotal_nodes = Num_Internal_Nodes[iproc] + Num_Border_Nodes[iproc]+
    Num_External_Nodes[iproc];
  itotal_elems = Num_Internal_Elems[iproc] + Num_Border_Elems[iproc];

  bytes_out += 5*sizeof(int);
  tt1 = second();

  if(Debug_Flag >= 4) {
    printf("[%d]: Putting init global info in file id: %d\n", Proc,
           mesh_exoid);
    printf("\tNumber Global Nodes: %d\n", Num_Node);
    printf("\tNumber Global Elements: %d\n", Num_Elem);
    printf("\tNumber Global Element Blocks: %d\n", Num_Elem_Blk);
    printf("\tNumber Global Node Sets: %d\n", Num_Node_Set);
    printf("\tNumber Global Side Sets: %d\n", Num_Side_Set);
  }

  if(ne_put_init_global(mesh_exoid, Num_Node, Num_Elem,
                        Num_Elem_Blk, Num_Node_Set,
                        Num_Side_Set) < 0) {
    fprintf(stderr, "[%d, %s]: ERROR, Unable to put global initial "
            "information in parallel mesh file!\n", Proc, yo);
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
  if(Num_N_Comm_Maps[iproc] > 0 && N_Comm_Map->node_cnt != 0 ) {
    ncomm_cnt = 1;
    itemp = N_Comm_Map->proc_ids[0];

    /* First find the count */
    for(i1=1; i1 < N_Comm_Map->node_cnt; i1++) {
      if(N_Comm_Map->proc_ids[i1] != itemp) {
        ncomm_cnt++;
        itemp = N_Comm_Map->proc_ids[i1];
      }
    }

    /* Allocate memory for the nodal communication maps */
    n_comm_map = malloc(ncomm_cnt*sizeof(NODE_COMM_MAP));
    if(!n_comm_map) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n",
              Proc, yo);
      exit(1);
    }

    n_comm_ids = malloc(2*ncomm_cnt*sizeof(int));
    if(!n_comm_ids) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n",
              Proc, yo);
      exit(1);
    }
    n_comm_ncnts = n_comm_ids + ncomm_cnt;

    /* Find the size of each map */
    ncomm_cnt = 0;
    itemp = N_Comm_Map->proc_ids[0];
    n_comm_map[0].node_cnt = 1;
    for(i1=1; i1 < N_Comm_Map->node_cnt; i1++) {
      if(N_Comm_Map->proc_ids[i1] != itemp) {
        itemp = N_Comm_Map->proc_ids[i1];
        ncomm_cnt++;
        n_comm_map[ncomm_cnt].node_cnt = 1;
      }
      else
        (n_comm_map[ncomm_cnt].node_cnt)++;
    }

    ncomm_cnt++;
    itemp = 0;

    /* Allocate memory for the maps */
    for(i1=0; i1 < ncomm_cnt; i1++) {
      n_comm_map[i1].proc_ids = malloc(2*(n_comm_map[i1].node_cnt)*
                                       sizeof(int));
      if(!(n_comm_map[i1].proc_ids)) {
        fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n",
                Proc, yo);
        exit(1);
      }

      n_comm_map[i1].node_ids = n_comm_map[i1].proc_ids +
	n_comm_map[i1].node_cnt;

      for(i2=0; i2 < n_comm_map[i1].node_cnt; i2++) {
        n_comm_map[i1].proc_ids[i2] = N_Comm_Map->proc_ids[itemp];
        n_comm_map[i1].node_ids[i2] = N_Comm_Map->node_ids[itemp++];
      }

      n_comm_ncnts[i1] = n_comm_map[i1].node_cnt;
      n_comm_ids[i1]   = n_comm_map[i1].proc_ids[0];
    }
  }
  else {
    ncomm_cnt = 0;
    n_comm_map = NULL;
    n_comm_ids = NULL;
    n_comm_ncnts = NULL;
  }

  /*
   * Find the number of distinct elemental communication maps that will
   * be needed. This assumes that there is only a single input
   * elemental communication map that stores information for all
   * processors.
   */
  if(Num_E_Comm_Maps[iproc] > 0) {
    ecomm_cnt = 1;
    itemp = E_Comm_Map->proc_ids[0];

    /* First find the count */
    for(i1=1; i1 < E_Comm_Map->elem_cnt; i1++) {
      if(E_Comm_Map->proc_ids[i1] != itemp) {
        ecomm_cnt++;
        itemp = E_Comm_Map->proc_ids[i1];
      }
    }

    /* Allocate memory for the elemental communication maps */
    e_comm_map = malloc(ecomm_cnt*sizeof(ELEM_COMM_MAP));
    if(!e_comm_map) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n",
              Proc, yo);
      exit(1);
    }

    e_comm_ids = malloc(2*ecomm_cnt*sizeof(int));
    if(!e_comm_ids) {
      fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n",
              Proc, yo);
      exit(1);
    }
    e_comm_ecnts = e_comm_ids + ecomm_cnt;

    /* Find the size of each map */
    ecomm_cnt = 0;
    itemp = E_Comm_Map->proc_ids[0];
    e_comm_map[0].elem_cnt = 1;
    for(i1=1; i1 < E_Comm_Map->elem_cnt; i1++) {
      if(E_Comm_Map->proc_ids[i1] != itemp) {
        itemp = E_Comm_Map->proc_ids[i1];
        ecomm_cnt++;
        e_comm_map[ecomm_cnt].elem_cnt = 1;
      }
      else
        (e_comm_map[ecomm_cnt].elem_cnt)++;
    }

    ecomm_cnt++;
    itemp = 0;

    /* Allocate memory for the maps */
    for(i1=0; i1 < ecomm_cnt; i1++) {
      e_comm_map[i1].proc_ids = malloc(3*(e_comm_map[i1].elem_cnt)*
                                       sizeof(int));
      if(!(e_comm_map[i1].proc_ids)) {
        fprintf(stderr, "[%d, %s]: ERROR, insufficient memory!\n",
                Proc, yo);
        exit(1);
      }

      e_comm_map[i1].elem_ids = e_comm_map[i1].proc_ids +
	e_comm_map[i1].elem_cnt;
      e_comm_map[i1].side_ids = e_comm_map[i1].elem_ids +
	e_comm_map[i1].elem_cnt;

      for(i2=0; i2 < e_comm_map[i1].elem_cnt; i2++) {
        e_comm_map[i1].proc_ids[i2] = E_Comm_Map->proc_ids[itemp];
        e_comm_map[i1].elem_ids[i2] = E_Comm_Map->elem_ids[itemp];
        e_comm_map[i1].side_ids[i2] = E_Comm_Map->side_ids[itemp++];
      }

      e_comm_ecnts[i1] = e_comm_map[i1].elem_cnt;
      e_comm_ids[i1]   = e_comm_map[i1].proc_ids[0];
    }
  }
  else {
    ecomm_cnt = 0;
    e_comm_map = NULL;
    e_comm_ids = NULL;
    e_comm_ecnts = NULL;
  }

  /* Output load balance information */
  bytes_out += 9*sizeof(int) + 2*sizeof(char);
  tt1 = second();

  if(Debug_Flag >= 4) {
    printf("[%d]: Putting init Nemesis info in file id: %d\n", Proc,
           mesh_exoid);
    printf("\tNumber of Proccesor for: %d\n", num_proc_for);
  }

  if(ne_put_init_info(mesh_exoid, num_proc_for, 1, "p") < 0) {
    fprintf(stderr, "[%d, %s]: ERROR, unable to output init info!\n",
            Proc, yo);
    exit(1);
  }

  if(Debug_Flag >= 6) {
    printf("[%d]: Putting init load balance info in file id: %d\n", Proc,
           mesh_exoid);
    printf("\tNumber Internal Nodes: %d\n", Num_Internal_Nodes[iproc]);
    printf("\tNumber Border Nodes: %d\n", Num_Border_Nodes[iproc]);
    printf("\tNumber External Nodes: %d\n", Num_External_Nodes[iproc]);
    printf("\tNumber Internal Elements: %d\n", Num_Internal_Elems[iproc]);
    printf("\tNumber Border Elements: %d\n", Num_Border_Elems[iproc]);
    printf("\tNumber Nodal Cmaps: %d\n", ncomm_cnt);
    printf("\tNumber Elemental Cmaps: %d\n", ecomm_cnt);
    printf("\tProccesor For: %d\n", proc_for);
  }

  if(ne_put_loadbal_param(mesh_exoid, Num_Internal_Nodes[iproc],
                          Num_Border_Nodes[iproc], Num_External_Nodes[iproc],
                          Num_Internal_Elems[iproc], Num_Border_Elems[iproc],
                          ncomm_cnt, ecomm_cnt,
                          proc_for) < 0) {
    fprintf(stderr, "[%d, %s]: ERROR, unable to output load balance info\n",
            Proc, yo);
    ex_close(mesh_exoid);
    exit(1);
  }

  PIO_Time_Array[1]  = (second() - tt1);
  total_out_time    += PIO_Time_Array[1];

  /* Output the communication map */
  bytes_out += 4*sizeof(int);
  tt1 = second();

  if(ne_put_cmap_params(mesh_exoid,
                        n_comm_ids,
                        n_comm_ncnts,
                        e_comm_ids,
                        e_comm_ecnts,
                        proc_for) < 0) {
    fprintf(stderr, "[%d, %s]: ERROR, unable to output comm map params!\n",
            Proc, yo);
    ex_close(mesh_exoid);
    exit(1);
  }

  PIO_Time_Array[2]  = (second() - tt1);
  total_out_time    += PIO_Time_Array[2];

  /*
   * The Nemesis node maps are lists of internal, border and external
   * FEM node numbers. These are output as local node numbers.
   */
  nem_node_mapi = (int *)array_alloc(__FILE__, __LINE__, 1,
                                     itotal_nodes + itotal_elems,
                                     sizeof(int));
  nem_node_mapb = nem_node_mapi + Num_Internal_Nodes[iproc];
  nem_node_mape = nem_node_mapb + Num_Border_Nodes[iproc];
  nem_elem_mapi = nem_node_mape + Num_External_Nodes[iproc];
  nem_elem_mapb = nem_elem_mapi + Num_Border_Elems[iproc];

  for(i1=0; i1 < itotal_nodes; i1++)
    nem_node_mapi[i1] = i1+1;

  /* Convert Elem_Map to local element numbering */
  reverse_map(Elem_Map[iproc], 0, itotal_elems,
	      &GElems[iproc][0], NULL, itotal_elems,
	      &Elem_Map[iproc][0]);
  /* Convert element IDs in the comm map to local numbering */
  for(i0=0; i0 < ecomm_cnt; i0++) {
    reverse_map(e_comm_map[i0].elem_ids, 0, e_comm_map[i0].elem_cnt,
		&GElems[iproc][0], NULL, itotal_elems,
		e_comm_map[i0].elem_ids);
  }

  /* Sort the GNodes array using the index array 'loc_index' */
  loc_index = (int *) array_alloc(__FILE__, __LINE__, 1, itotal_nodes,
				  sizeof(int));

  /* Initialize index array */
  for (i2=0; i2 < itotal_nodes; i2++)
    loc_index[i2] = i2;
      
  /*
   * Sort the GNodes[iproc] array via the index array
   * 'loc_index' 
   */
  gds_iqsort(GNodes[iproc], loc_index, itotal_nodes);
  
  /* Convert nodal IDs in the comm map to local numbering */
  for(i0=0; i0 < ncomm_cnt; i0++) {
    reverse_map(n_comm_map[i0].node_ids, 0, n_comm_map[i0].node_cnt,
		&GNodes[iproc][0], loc_index, itotal_nodes,
		n_comm_map[i0].node_ids);
  }

  PIO_Time_Array[3] = 0.0;

  if(Num_N_Comm_Maps[iproc] > 0 && N_Comm_Map->node_cnt != 0 ) {

    bytes_out += 2*Num_External_Nodes[iproc]*sizeof(int);
    tt1 = second();

    for(i1=0; i1 < ncomm_cnt; i1++) {
      if(ne_put_node_cmap(mesh_exoid, n_comm_ids[i1], n_comm_map[i1].node_ids,
                          n_comm_map[i1].proc_ids, proc_for) < 0) {
        fprintf(stderr, "[%d, %s]: ERROR, unable to output nodal comm map!\n",
                Proc, yo);
        ex_close(mesh_exoid);
        exit(1);
      }

      free(n_comm_map[i1].proc_ids);
    }

    free(n_comm_ids);
    free(n_comm_map);
    PIO_Time_Array[3] += (second() - tt1);

  }
  if(Num_E_Comm_Maps[iproc] > 0) {
    bytes_out += 3*(E_Comm_Map->elem_cnt)*sizeof(int);
    tt1 = second();

    for(i1=0; i1 < ecomm_cnt; i1++) {
      if(ne_put_elem_cmap(mesh_exoid, e_comm_ids[i1], e_comm_map[i1].elem_ids,
			  e_comm_map[i1].side_ids, e_comm_map[i1].proc_ids,
			  proc_for) < 0) {
        fprintf(stderr,
                "[%d, %s]: ERROR, unable to output elemental comm map!\n",
                Proc, yo);
        ex_close(mesh_exoid);
        exit(1);
      }

      free(e_comm_map[i1].proc_ids);
    }

    free(e_comm_ids);
    free(e_comm_map);

    PIO_Time_Array[3] += (second() - tt1);
  }

  total_out_time    += PIO_Time_Array[3];

  /* Output the global node set parameters */
  if(Num_Node_Set > 0) {

    glob_ns_df_cnts = (int *)array_alloc(__FILE__, __LINE__, 1, Num_Node_Set,
                                         sizeof(int));
    for(i1=0; i1 < Num_Node_Set; i1++)
      glob_ns_df_cnts[i1] = 0;

    bytes_out += 3*Num_Node_Set*sizeof(int);
    tt1 = second();

    if(ne_put_ns_param_global(mesh_exoid, Node_Set_Ids,
                              Num_Nodes_In_NS, glob_ns_df_cnts) < 0) {
      fprintf(stderr,
              "[%d, %s]: ERROR, unable to output global node-set params\n",
              Proc, yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    PIO_Time_Array[4]  = (second() - tt1);
    total_out_time    += PIO_Time_Array[4];

    safe_free((void **) &glob_ns_df_cnts);
  }

  /* Output the global side set parameters */
  if(Num_Side_Set > 0) {

    glob_ss_df_cnts = (int *)array_alloc(__FILE__, __LINE__, 1, Num_Side_Set,
                                         sizeof(int));
    for(i1=0; i1 < Num_Side_Set; i1++)
      glob_ss_df_cnts[i1] = 0;

    bytes_out += 3*Num_Side_Set*sizeof(int);
    tt1 = second();

    if(ne_put_ss_param_global(mesh_exoid, Side_Set_Ids,
                              Num_Elems_In_SS, glob_ss_df_cnts) < 0) {
      fprintf(stderr,
              "[%d, %s]: ERROR, unable to output global side-set params\n",
              Proc, yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    PIO_Time_Array[5]  = (second() - tt1);
    total_out_time    += PIO_Time_Array[5];

    safe_free((void **) &glob_ss_df_cnts);
  }

  bytes_out += Num_Elem_Blk*sizeof(int);
  tt1 = second();

  if(ne_put_eb_info_global(mesh_exoid, Elem_Blk_Ids, Num_Elems_In_EB) < 0) {
    fprintf(stderr,
            "[%d, %s]: ERROR, unable to output global elem blk IDs\n",
            Proc, yo);
    ex_close(mesh_exoid);
    exit(1);
  }

  PIO_Time_Array[6]  = (second() - tt1);
  total_out_time    += PIO_Time_Array[6];
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  if (iproc == 0) {
    /* Generate a QA record for the utility */
    time_val = time(NULL);
    ct_ptr   = asctime(localtime(&time_val));
    strcpy(tm_date, ct_ptr);

    /* Break string with null characters */
    tm_date[3]  = '\0';
    tm_date[7]  = '\0';
    tm_date[10] = '\0';
    tm_date[19] = '\0';

    sprintf(qa_date, "%s %s %s", &tm_date[8], &tm_date[4], &tm_date[20]);
    sprintf(qa_time, "%s", &tm_date[11]);
    strcpy(qa_name, UTIL_NAME);
    strcpy(qa_vers, VER_STR);

    if(qa_date[strlen(qa_date)-1] == '\n')
      qa_date[strlen(qa_date)-1] = '\0';

    strcpy(qa_record[4*(num_qa_records-1)],     qa_name);
    strcpy(qa_record[(4*(num_qa_records-1))+1], qa_vers);
    strcpy(qa_record[(4*(num_qa_records-1))+2], qa_date);
    strcpy(qa_record[(4*(num_qa_records-1))+3], qa_time);

    /* Output QA records to screen */
    if(Debug_Flag >= 4) {
      printf("[%d]: Number of QA records: %d\n", Proc, num_qa_records);
      if(Debug_Flag >= 6) {
	printf("QA Records:\n");
	for(i1=0; i1 < 4*(num_qa_records); i1++) {
	  printf("\t%s\n", qa_record[i1]);
	}
      }
    }

    /* Output the QA and Info records */
    for(i1=0; i1 < 4*num_qa_records; i1++)
      bytes_out += (MAX_STR_LENGTH+MAX_LINE_LENGTH)*sizeof(char);

    tt1 = second();

    if(ex_put_qa(mesh_exoid, num_qa_records,
		 (char *(*)[]) &qa_record[0]) < 0) {
      fprintf(stderr, "[%d, %s]: ERROR Could not put QA records\n",
	      Proc, yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    if(num_info_records > 0) {
      if(Debug_Flag >= 4) {
	printf("[%d]: Number of info records: %d\n", Proc, num_info_records);
      }

      if(ex_put_info(mesh_exoid, num_info_records, info_record) < 0) {
	fprintf(stderr, "[%d, %s]: ERROR Could not put Info records\n",
		Proc, yo);
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
    printf("[%d]: Number of Coordinate Frames: %d\n", Proc, Num_Coordinate_Frames);
  }

  if (Num_Coordinate_Frames > 0) {
    void *Coordinate_Frame_Coordinates = NULL;
    if (io_ws <= sizeof(float)) {
      Coordinate_Frame_Coordinates = Coordinate_Frame_Coordinates_sp;
    } else {
      Coordinate_Frame_Coordinates = Coordinate_Frame_Coordinates_dp;
    }
    if (ex_put_coordinate_frames(mesh_exoid, Num_Coordinate_Frames,
				 Coordinate_Frame_Ids,
				 Coordinate_Frame_Coordinates,
				 Coordinate_Frame_Tags) < 0) {
      fprintf(stderr,
	      "[%d, %s]: ERROR, Unable to put coordinate frame data in parallel mesh file\n",
	      Proc, yo);
      ex_close(mesh_exoid);
      exit(1);
    }
  }

  if (proc_for == 0)
    strcpy(cTitle, GeomTitle);
  else
    sprintf(cTitle, "Parallel Mesh File for Processor %d", proc_for);

  /* Output the initial information to the parallel Exodus file(s) */
  bytes_out += strlen(cTitle)*sizeof(char) + 6*sizeof(int);
  tt1 = second();

  if(Debug_Flag >= 4) {
    printf("[%d]: Putting init info in file id: %d\n", Proc, mesh_exoid);
    printf("\tTitle: %s\n", cTitle);
    printf("\tNumber Dimensions: %d\n", Num_Dim);
    printf("\tNumber Nodes: %d\n", itotal_nodes);
    printf("\tNumber Elements: %d\n",
           Num_Internal_Elems[iproc]+Num_Border_Elems[iproc]);
    printf("\tNumber Element Blocks: %d\n", Num_Elem_Blk);
    printf("\tNumber Node Sets: %d\n", Num_Node_Set);
    printf("\tNumber Side Sets: %d\n", Num_Side_Set);
  }

  if(ex_put_init(mesh_exoid, cTitle, Num_Dim, itotal_nodes,
                 Num_Internal_Elems[iproc]+Num_Border_Elems[iproc],
                 Num_Elem_Blk,
                 Num_Node_Set,
                 Num_Side_Set) < 0) {
    fprintf(stderr,
	    "[%d, %s]: ERROR, Unable to put initial info in parallel mesh file\n",
            Proc, yo);
    ex_close(mesh_exoid);
    exit(1);
  }

  /* Output names... */
  if (Num_Elem_Blk > 0) {
    error = ex_put_names(mesh_exoid, EX_ELEM_BLOCK, &Elem_Blk_Names[0]);
    check_exodus_error(error, "ex_put_names");
  }
  if (Num_Node_Set > 0) {
    error = ex_put_names(mesh_exoid, EX_NODE_SET, &Node_Set_Names[0]);
    check_exodus_error(error, "ex_put_names");
  }
  if (Num_Side_Set > 0) {
    error = ex_put_names(mesh_exoid, EX_SIDE_SET, &Side_Set_Names[0]);
    check_exodus_error(error, "ex_put_names");
  }
  

  PIO_Time_Array[7]  = (second() - tt1);
  total_out_time    += PIO_Time_Array[7];

  /* Assign the coordinates to the coord_vector */
  if (io_ws == sizeof(float)) {
    switch(Num_Dim) {
    case 3:
      z_coord = (void *) Coor_sp[iproc][2];
      /* FALLTHROUGH */
    case 2:
      y_coord = (void *) Coor_sp[iproc][1];
      /* FALLTHROUGH */
    case 1:
      x_coord = (void *) Coor_sp[iproc][0];

      break;
    }
  }
  else {
    switch(Num_Dim) {
    case 3:
      z_coord = (void *) Coor_dp[iproc][2];
      /* FALLTHROUGH */
    case 2:
      y_coord = (void *) Coor_dp[iproc][1];
      /* FALLTHROUGH */
    case 1:
      x_coord = (void *) Coor_dp[iproc][0];

      break;
    }
  }

  /* Output the coordinates to the parallel Exodus file */
  bytes_out += Num_Dim * itotal_nodes * io_ws;
  tt1 = second();

  if(Debug_Flag >= 4) {
    printf("[%d]: Putting coordinate info in file id: %d\n", Proc,
           mesh_exoid);
  }
  if(ex_put_coord(mesh_exoid, x_coord, y_coord, z_coord) < 0) {
    fprintf(stderr,
            "[%d, %s]: ERROR, could not write out nodal coordinates\n",
            Proc, yo);
    ex_close(mesh_exoid);
    exit(1);
  }

  PIO_Time_Array[9]   = (second() - tt1);
  total_out_time     += PIO_Time_Array[9];

  bytes_out += itotal_nodes*sizeof(int);
  tt1 = second();

  if(ne_put_node_map(mesh_exoid, nem_node_mapi, nem_node_mapb,
		     nem_node_mape, proc_for) < 0) {
    fprintf(stderr,
            "[%d, %s]: ERROR, could not write Nemesis nodal number map!\n",
            Proc, yo);
    check_exodus_error(ex_close(mesh_exoid), "ex_close");
    ex_close(mesh_exoid);
    exit(1);
  }

  /* If non-NULL, output the global node id map which preserves
     the global node ids in the original mesh */
  if (Proc_Global_Node_Id_Map[iproc] != NULL) {
    bytes_out += itotal_nodes * sizeof(int);
    if (ex_put_map_param(mesh_exoid, 1, 0) < 0) {
      fprintf(stderr, "[%d, %s]: ERROR, unable to define global node map parameters!\n",
	      Proc, yo);
      ex_close(mesh_exoid);
      exit(1);
    }
    
    if (ex_put_num_map(mesh_exoid, EX_NODE_MAP, 1, Proc_Global_Node_Id_Map[iproc]) < 0) {
      fprintf(stderr, "[%d, %s]: ERROR, unable to output global node id map!\n",
	      Proc, yo);
      ex_close(mesh_exoid);
      exit(1);
    }
    if (ex_put_name(mesh_exoid, EX_NODE_MAP, 1, "original_global_id_map") < 0) {
      fprintf(stderr, "[%d, %s]: ERROR, unable to define global node map name!\n",
	      Proc, yo);
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
    int *EB_Ids;
    int *EB_Cnts;
    int *EB_NperE;
    int *EB_Nattr;
    char **EB_Types;
    EB_Ids = (int*)array_alloc(__FILE__, __LINE__, 1, (4*Num_Elem_Blk),
			       sizeof(int));
    if (!EB_Ids) {
      fprintf(stderr, "%s: fatal: insufficient memory\n", yo);
      exit(1);
    }

    EB_Cnts  = EB_Ids + 1 * Num_Elem_Blk;
    EB_NperE = EB_Ids + 2 * Num_Elem_Blk;
    EB_Nattr = EB_Ids + 3 * Num_Elem_Blk;

    EB_Types = (char**)array_alloc(__FILE__, __LINE__, 2, Num_Elem_Blk,
				   MAX_STR_LENGTH+1, sizeof(char));
    if (!EB_Types) {
      fprintf(stderr, "%s: fatal: insufficient memory\n", yo);
      exit(1);
    }
    
    /* Initialize some of the arrays */
    for (i1=0; i1 < Num_Elem_Blk; i1++) {
      EB_Ids[i1]   = 0;
      EB_Cnts[i1]  = 0;
      EB_NperE[i1] = 0;
      EB_Nattr[i1] = 0;
      memset(EB_Types[i1], '\0', (MAX_STR_LENGTH+1));
    }

    cnt = 0;
    for(i1=0; i1 < Num_Elem_Blk; i1++) {
      ifound = FALSE;
      for(i2=0; i2 < Proc_Num_Elem_Blk[iproc]; i2++) {
	if(Proc_Elem_Blk_Ids[iproc][i2] == Elem_Blk_Ids[i1]) {
	  ifound = TRUE;
	  break;
	}
      }
      
      /*
       * If it's not found then this element block is null on the current
       * processor.
       */
      if(ifound == FALSE) {
	int iblk = Proc_Num_Elem_Blk[iproc]+cnt;
	Proc_Elem_Blk_Ids[iproc][iblk] =  Elem_Blk_Ids[i1];
	Proc_Num_Elem_In_Blk[iproc][iblk] = 0;
	Proc_Nodes_Per_Elem[iproc][iblk]  = 0;
	Proc_Num_Attr[iproc][iblk]        = 0;
	cnt++;
	
	/*
	 * if there are elemental variables in the restart data,
	 * take care of the truth table for this processor here
	 */
	if (Restart_Info.Flag > 0)
	  for (i2=0; i2<Restart_Info.NVar_Elem; i2++)
	    Restart_Info.Elem_TT[iproc][i1*Restart_Info.NVar_Elem+i2] = 0;
      }
    }

    /* Output the elemental block(s). */
    for(i1=0; i1 < Num_Elem_Blk; i1++) {

      iglobal_blk = Elem_Blk_Ids[i1];

      /* Find the local element block index */
      for(ilocal=0; ilocal < Num_Elem_Blk; ilocal++)
	{
	  if(Proc_Elem_Blk_Ids[iproc][ilocal] == iglobal_blk)
	    break;
	}

      /* Error check */
      if(ilocal >= Num_Elem_Blk) {
	fprintf(stderr, "[%d, %s]: Error finding local element block ID\n",
		Proc, yo);
	exit(1);
      }
    
      /* If it's a non-null block output all information */
      EB_Ids[i1]   = iglobal_blk;
      if (ilocal < Proc_Num_Elem_Blk[iproc]) {

	/* Generate the ExodusII element name */
	strncpy(EB_Types[i1], Elem_Blk_Types[i1], MAX_STR_LENGTH);
	EB_Types[i1][MAX_STR_LENGTH] = '\0';

	EB_Cnts[i1]  = Proc_Num_Elem_In_Blk[iproc][ilocal];
	EB_NperE[i1] = Proc_Nodes_Per_Elem[iproc][ilocal];
	EB_Nattr[i1] = Proc_Num_Attr[iproc][ilocal];
      }
    }
    if(Debug_Flag >= 4) {
      printf("[%d]: Putting concat_elem_block info in file id: %d\n", Proc,
	     mesh_exoid);
    }
    error = ex_put_concat_elem_block(mesh_exoid,
				     &EB_Ids[0], &EB_Types[0],
				     &EB_Cnts[0], &EB_NperE[0],
				     &EB_Nattr[0], 1);
    check_exodus_error(error, "ex_put_concat_elem_block");

    safe_free((void **) &EB_Ids);
    safe_free((void **) &EB_Types);

    /* Output attribute names for each element block */
    for(i1=0; i1 < Num_Elem_Blk; i1++) {

      iglobal_blk = Elem_Blk_Ids[i1];

      /* Find the local element block index */
      for(ilocal=0; ilocal < Num_Elem_Blk; ilocal++)
	{
	  if(Proc_Elem_Blk_Ids[iproc][ilocal] == iglobal_blk)
	    break;
	}

      /* If it's a non-null block output attribute name information */
      if (ilocal < Proc_Num_Elem_Blk[iproc]) {

	if (Proc_Num_Attr[iproc][ilocal] > 0) {
	  if (ex_put_attr_names(mesh_exoid, EX_ELEM_BLOCK, Elem_Blk_Ids[i1],
				Elem_Blk_Attr_Names[i1]) < 0) {
	    fprintf(stderr,
		    "[%d, %s]: ERROR, could not write Exodus attribute names!\n",
		    Proc, yo);
	    ex_close(mesh_exoid);
	    exit(1);
	  }
	}
      }
    }

  
    /* Reset GNodes to start at 1 instead of 0 */
    for(i1=0; i1 < itotal_nodes; (GNodes[iproc][i1++])++);

    /* Output the Exodus node number map */
    bytes_out += itotal_nodes*sizeof(int);
    tt1 = second();

    if(Debug_Flag >= 4) {
      printf("[%d]: Putting node_num_map in file id: %d\n", Proc,
	     mesh_exoid);
    }
    if(ex_put_node_num_map(mesh_exoid, GNodes[iproc]) < 0) {
      fprintf(stderr,
	      "[%d, %s]: ERROR, could not write Exodus node number map!\n",
	      Proc, yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    PIO_Time_Array[10]  = (second() - tt1);
    total_out_time     += PIO_Time_Array[10];

    /*
     * Allocate memory for the elemental map. Currently this map is assigned
     * as a linear array since it is not really used.
     */
    iElem_Map = (int *)array_alloc(__FILE__, __LINE__, 1,
				   Num_Internal_Elems[iproc] +
				   Num_Border_Elems[iproc],
				   sizeof(int));
    for(i1=0; i1 < Num_Internal_Elems[iproc]+Num_Border_Elems[iproc]; i1++)
      iElem_Map[i1] = GElems[iproc][i1] + 1;

    bytes_out += 2 * Num_Internal_Elems[iproc] * Num_Border_Elems[iproc] *
      sizeof(int);
    tt1 = second();

    if(Debug_Flag >= 4) {
      printf("[%d]: Putting elem_num_map info in file id: %d\n", Proc,
	     mesh_exoid);
    }
    if(ex_put_elem_num_map(mesh_exoid, iElem_Map) < 0) {
      fprintf(stderr, "[%d, %s]: ERROR, unable to output element map\n",
	      Proc, yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    /* If non-NULL, output the global element id map which preserves
       the global element ids in the original mesh */
    if (Proc_Global_Elem_Id_Map[iproc] != NULL) {
      bytes_out += Num_Internal_Elems[iproc] * Num_Border_Elems[iproc] * sizeof(int);
      if (ex_put_map_param(mesh_exoid, 0, 1) < 0) {
	fprintf(stderr, "[%d, %s]: ERROR, unable to define global map parameters!\n",
		Proc, yo);
	ex_close(mesh_exoid);
	exit(1);
      }
      
      if (ex_put_num_map(mesh_exoid, EX_ELEM_MAP, 1, Proc_Global_Elem_Id_Map[iproc]) < 0) {
	fprintf(stderr, "[%d, %s]: ERROR, unable to output global id map!\n",
		Proc, yo);
	ex_close(mesh_exoid);
	exit(1);
      }
      if (ex_put_name(mesh_exoid, EX_ELEM_MAP, 1, "original_global_id_map") < 0) {
	fprintf(stderr, "[%d, %s]: ERROR, unable to define global map name!\n",
		Proc, yo);
	ex_close(mesh_exoid);
	exit(1);
      }

    }
      
    /* Also output the Nemesis element map */
    if(ne_put_elem_map(mesh_exoid, Elem_Map[iproc],
		       Elem_Map[iproc] + Num_Internal_Elems[iproc],
		       proc_for) < 0) {
      fprintf(stderr, "[%d, %s]: ERROR, unable to output nemesis element map!\n",
	      Proc, yo);
      ex_close(mesh_exoid);
      exit(1);
    }

    PIO_Time_Array[12]  = (second() - tt1);
    total_out_time     += PIO_Time_Array[12];
    safe_free((void **) &iElem_Map);

    for(i1=0; i1 < Num_Elem_Blk; i1++) {

      iglobal_blk = Elem_Blk_Ids[i1];

      /* Find the local element block index */
      for(ilocal=0; ilocal < Num_Elem_Blk; ilocal++)
	{
	  if(Proc_Elem_Blk_Ids[iproc][ilocal] == iglobal_blk)
	    break;
	}

      /* Error check */
      if(ilocal >= Num_Elem_Blk)
	{
	  fprintf(stderr, "[%d, %s]: Error finding local element block ID\n",
		  Proc, yo);
	  exit(1);
	}

      /* If it's a non-null block output all information */
      if(ilocal < Proc_Num_Elem_Blk[iproc])
	{
	  /* Find the first index into the connectivity for this block */
	  iIndex[0] = 0;
	  for(i2=0; i2 < ilocal; i2++) {
	    iIndex[0] += Proc_Num_Elem_In_Blk[iproc][i2] *
	      Proc_Nodes_Per_Elem[iproc][i2];
	  }

	  PIO_Time_Array[13] += (second() - tt1);

	  tmp_cnt = Proc_Num_Elem_In_Blk[iproc][ilocal] *
	    Proc_Nodes_Per_Elem[iproc][ilocal];
      
	  /* Generate the connectivity array for local node numbering */
	  proc_local_conn = (int *) array_alloc(__FILE__, __LINE__, 1, tmp_cnt,
						sizeof(int));

	  reverse_map(&Proc_Elem_Connect[iproc][iIndex[0]], 1, tmp_cnt,
		      &GNodes[iproc][0], loc_index, itotal_nodes,
		      proc_local_conn);
      
	  bytes_out += Proc_Nodes_Per_Elem[iproc][ilocal] *
	    Proc_Num_Elem_In_Blk[iproc][ilocal]*sizeof(int);
	  tt1 = second();

	  if(Debug_Flag >= 4) {
	    printf("[%d]: Putting element_connectivity info in file id: %d\n", Proc,
		   mesh_exoid);
	  }
	  if(ex_put_elem_conn(mesh_exoid, Proc_Elem_Blk_Ids[iproc][ilocal],
			      proc_local_conn) < 0) {
	    fprintf(stderr, "[%d, %s]: ERROR, unable to output connectivity\n",
		    Proc, yo);
	    ex_close(mesh_exoid);
	    exit(1);
	  }

	  PIO_Time_Array[14] += (second() - tt1);

	  /* Free up the locally numbered connectivity */
	  safe_free((void **) &proc_local_conn);

	  if(Proc_Num_Attr[iproc][ilocal] > 0) {

	    /* Find the first index into the attribute list for this block */
	    iIndex[1] = 0;
	    for(i2=0; i2 < ilocal; i2++) {
	      iIndex[1] += Proc_Num_Attr[iproc][i2] *
		Proc_Num_Elem_In_Blk[iproc][i2];
	    }

	    bytes_out += Proc_Num_Elem_In_Blk[iproc][ilocal] * io_ws;
	    tt1 = second();

	    if (io_ws == sizeof(float))
	      ptr = (void *) &(Proc_Elem_Attr_sp[iproc][iIndex[1]]);
	    else
	      ptr = (void *) &(Proc_Elem_Attr_dp[iproc][iIndex[1]]);

	    if(ex_put_elem_attr(mesh_exoid, Proc_Elem_Blk_Ids[iproc][ilocal],
				ptr) < 0) {
	      fprintf(stderr,
		      "[%d, %s]: ERROR, unable to output element attributes\n",
		      Proc, yo);
	      exit(1);
	    }

	    PIO_Time_Array[15] += (second() - tt1);
	    total_out_time     += PIO_Time_Array[15];
	  }

	} /* End "if(ilocal < Proc_Num_Elem_Blk[iproc])" */

    } /* End "for(i1=0; i1 < Num_Elem_Block; i1++)" */
  }
  total_out_time += (PIO_Time_Array[13] + PIO_Time_Array[14] +
                     PIO_Time_Array[15]);

  /*
   * Write out the node-set information. Note that the value of the
   * node-set distribution factor is not currently used so only a
   * dummy set is output for this value.
   */
  iMaxLen = 0;
  for(i1=0; i1 < Proc_Num_Node_Sets[iproc]; i1++)
    iMaxLen = PEX_MAX(Proc_NS_Count[iproc][i1], iMaxLen);

  /* Renumber Node set node lists to use local node numbers */
  if(Proc_Num_Node_Sets[iproc] > 0) {
    proc_local_ns = (int *)array_alloc(__FILE__, __LINE__, 1,
                                       Proc_NS_List_Length[iproc],
                                       sizeof(int));

    reverse_map(&Proc_NS_List[iproc][0], 1,
		Proc_NS_List_Length[iproc],
		&GNodes[iproc][0], loc_index, itotal_nodes,
		proc_local_ns);
  }

  safe_free((void **) &loc_index);

  PIO_Time_Array[16] = 0.0;
  PIO_Time_Array[17] = 0.0;

  /* Fill in the information for the NULL node sets */
  cnt = 0;
  for(i1=0; i1 < Num_Node_Set; i1++) {
    ifound = FALSE;
    for(i2=0; i2 < Proc_Num_Node_Sets[iproc]; i2++) {
      if(Proc_NS_Ids[iproc][i2] == Node_Set_Ids[i1]) {
        ifound = TRUE;
        break;
      }
    }

    if(ifound == FALSE) {
      Proc_NS_Ids[iproc][Proc_Num_Node_Sets[iproc]+cnt] = Node_Set_Ids[i1];
      Proc_NS_Count[iproc][Proc_Num_Node_Sets[iproc]+cnt] = 0;
      Proc_NS_DF_Count[iproc][Proc_Num_Node_Sets[iproc]+cnt] = 0;
      cnt++;
    }
  }

  tt1 = second();
  if (Num_Node_Set > 0) {
    int i3;
    int *conc_ids;
    int *conc_nodes;
    int *conc_df;
    int *conc_nind;
    int *conc_dind;
    int *conc_nlist;
    float  *conc_sdf;
    double *conc_ddf;

    int dcount;
    int ncount;

    /* Count nodes and distribution factors */
    dcount = 0;
    ncount = 0;
    for (i1=0; i1 < Num_Node_Set; i1++) {
      dcount += Proc_NS_DF_Count[iproc][i1];
      ncount += Proc_NS_Count[iproc][i1];
    }

    conc_ids = malloc(Num_Node_Set * sizeof(int));
    conc_nodes = malloc(Num_Node_Set * sizeof(int));
    conc_df = malloc(Num_Node_Set * sizeof(int));
    conc_nind = malloc(Num_Node_Set * sizeof(int));
    conc_dind = malloc(Num_Node_Set * sizeof(int));
    conc_nlist = malloc(ncount * sizeof(int));
    if (io_ws == sizeof(float))
      conc_sdf = malloc(dcount * io_ws);
    else
      conc_ddf = malloc(dcount * io_ws);
    
    ncount = 0;
    dcount = 0;
    for(i1=0; i1 < Num_Node_Set; i1++) {

      /* Find the local ID */
      for(i2=0; i2 < Num_Node_Set; i2++) {
	if(Proc_NS_Ids[iproc][i2] == Node_Set_Ids[i1])
	  break;
      }

      conc_ids[i1] = Proc_NS_Ids[iproc][i2];
      conc_nodes[i1] = Proc_NS_Count[iproc][i2];
      conc_df[i1] = Proc_NS_DF_Count[iproc][i2];
      
      conc_nind[i1] = ncount;
      for (i3=0; i3 < Proc_NS_Count[iproc][i2]; i3++) {
	conc_nlist[ncount++] = proc_local_ns[Proc_NS_Pointers[iproc][i2]+i3];
      }

      conc_dind[i1] = dcount;
      if (io_ws == sizeof(float)) {
	for (i3=0; i3 < Proc_NS_DF_Count[iproc][i2]; i3++) {
	  conc_sdf[dcount++] = Proc_NS_Dist_Fact_sp[iproc][Proc_NS_Pointers[iproc][i2]+i3];
	}
      } else {
	for (i3=0; i3 < Proc_NS_DF_Count[iproc][i2]; i3++) {
	  conc_ddf[dcount++] = Proc_NS_Dist_Fact_dp[iproc][Proc_NS_Pointers[iproc][i2]+i3];
	}
      }
    }

    if (io_ws == sizeof(float)) {
      ex_put_concat_node_sets(mesh_exoid, conc_ids, conc_nodes, conc_df, conc_nind,
			      conc_dind, conc_nlist, conc_sdf);
    } else {
      ex_put_concat_node_sets(mesh_exoid, conc_ids, conc_nodes, conc_df, conc_nind,
			      conc_dind, conc_nlist, conc_ddf);
    }

    if (conc_ids   != NULL) free(conc_ids);
    if (conc_nodes != NULL) free(conc_nodes);
    if (conc_df    != NULL) free(conc_df);
    if (conc_nind  != NULL) free(conc_nind);
    if (conc_dind  != NULL) free(conc_dind);
    if (conc_nlist != NULL) free(conc_nlist);
    if (io_ws == sizeof(float)) {
      if (conc_sdf != NULL) free(conc_sdf);
    } else {
      if (conc_ddf != NULL) free(conc_ddf);
    }
  }
  total_out_time += second() - tt1;

  /* Free local number array */
  if(Proc_Num_Node_Sets[iproc] > 0)
    safe_free((void **) &proc_local_ns);

  /* Renumber element SS to use local element numbers */
  if(Proc_Num_Side_Sets[iproc] > 0) {
    proc_local_ss = (int *)array_alloc(__FILE__, __LINE__, 1,
                                       Proc_SS_Elem_List_Length[iproc],
                                       sizeof(int));
    reverse_map(&Proc_SS_Elem_List[iproc][0], 0,
		Proc_SS_Elem_List_Length[iproc],
		&GElems[iproc][0], NULL,
		Num_Internal_Elems[iproc]+Num_Border_Elems[iproc],
		proc_local_ss);
  }

  PIO_Time_Array[18] = 0.0;
  PIO_Time_Array[19] = 0.0;

  /* Set up the null side sets */
  cnt = 0;
  for(i1=0; i1 < Num_Side_Set; i1++) {
    ifound = FALSE;
    for(i2=0; i2 < Proc_Num_Side_Sets[iproc]; i2++) {
      if(Proc_SS_Ids[iproc][i2] == Side_Set_Ids[i1]) {
        ifound = TRUE;
        break;
      }
    }

    if(ifound == FALSE) {
      Proc_SS_Ids[iproc][Proc_Num_Side_Sets[iproc]+cnt] = Side_Set_Ids[i1];
      Proc_SS_Elem_Count[iproc][Proc_Num_Side_Sets[iproc]+cnt] = 0;
      Proc_SS_DF_Count[iproc][Proc_Num_Side_Sets[iproc]+cnt]   = 0;
      cnt++;
    }
  }

  /* Output concatenated sidesets.  For each processor need:
   * side_set ids          (size Num_Side_Set)
   * num_side_per_set      (size Num_Side_Set)
   * num_dist_per_set      (size Num_Side_Set)
   * side_sets_elem_index  (size Num_Side_Set)
   * side_sets_dist_index  (size Num_Side_Set)
   * side_sets_elem_list
   * side_sets_side_list
   * side_sets_dist_fact
   */
  
  tt1 = second();
  if (Num_Side_Set > 0) {
    int i3;
    int df_count;
    int el_count;
    int *conc_ids;
    int *conc_sides;
    int *conc_dist;
    int *conc_eind;
    int *conc_dind;
    int *conc_elist;
    int *conc_slist;
    float  *conc_sdflist;
    double *conc_ddflist;
    
    df_count = 0;
    el_count = 0;
    for (i1=0; i1 < Num_Side_Set; i1++) {
      df_count += Proc_SS_DF_Count[iproc][i1];
      el_count += Proc_SS_Elem_Count[iproc][i1];
    }
    
    conc_ids   = malloc(Num_Side_Set * sizeof(int));
    conc_sides = malloc(Num_Side_Set * sizeof(int));
    conc_dist  = malloc(Num_Side_Set * sizeof(int));
    conc_eind  = malloc(Num_Side_Set * sizeof(int));
    conc_dind  = malloc(Num_Side_Set * sizeof(int));

    conc_elist = malloc(el_count * sizeof(int));
    conc_slist = malloc(el_count * sizeof(int));
    if (io_ws == sizeof(float)) 
      conc_sdflist = malloc(df_count * io_ws);
    else
      conc_ddflist = malloc(df_count * io_ws);
	
    /* Fill in the arrays ... */
    df_count = 0;
    el_count = 0;
    for(i1=0; i1 < Num_Side_Set; i1++) {
      
      /* Find the local ID of this side set */
      for(i2=0; i2 < Num_Side_Set; i2++) {
	if(Proc_SS_Ids[iproc][i2] == Side_Set_Ids[i1])
	  break;
      }
      
      conc_ids[i1]   = Proc_SS_Ids[iproc][i2];
      conc_sides[i1] = Proc_SS_Elem_Count[iproc][i2];
      conc_dist[i1]  = Proc_SS_DF_Count[iproc][i2];
      
      conc_eind[i1]  = el_count;
      for (i3 = 0; i3 < Proc_SS_Elem_Count[iproc][i2]; i3++) {
	conc_elist[el_count] = proc_local_ss[Proc_SS_Elem_Pointers[iproc][i2]+i3];
	conc_slist[el_count] = Proc_SS_Side_List[iproc][Proc_SS_Elem_Pointers[iproc][i2]+i3];
	el_count++;
      }

      conc_dind[i1] = df_count;
      if (io_ws == sizeof(float)) {
	for (i3 = 0; i3 < Proc_SS_DF_Count[iproc][i2]; i3++) {
	  conc_sdflist[df_count++] =
	    Proc_SS_Dist_Fact_sp[iproc][Proc_SS_DF_Pointers[iproc][i2]+i3];
	}	  
      } else {
	for (i3 = 0; i3 < Proc_SS_DF_Count[iproc][i2]; i3++) {
	  conc_ddflist[df_count++] =
	    Proc_SS_Dist_Fact_dp[iproc][Proc_SS_DF_Pointers[iproc][i2]+i3];
	}	  
      }
    }
    if (io_ws == sizeof(float)) 
      ex_put_concat_side_sets(mesh_exoid, conc_ids, conc_sides, conc_dist, conc_eind,
			      conc_dind, conc_elist, conc_slist, conc_sdflist);
    else
      ex_put_concat_side_sets(mesh_exoid, conc_ids, conc_sides, conc_dist, conc_eind,
			      conc_dind, conc_elist, conc_slist, conc_ddflist);
    
    if (conc_ids   != NULL) free(conc_ids);
    if (conc_sides != NULL) free(conc_sides);
    if (conc_dist  != NULL) free(conc_dist);
    if (conc_eind  != NULL) free(conc_eind);
    if (conc_dind  != NULL) free(conc_dind);
    if (conc_elist != NULL) free(conc_elist);
    if (conc_slist != NULL) free(conc_slist);
    if (io_ws == sizeof(float)) {
      if (conc_sdflist != NULL) free(conc_sdflist);
    } else {
      if (conc_ddflist != NULL) free(conc_ddflist);
    }
  }
  PIO_Time_Array[19] += (second() - tt1);
  total_out_time += (PIO_Time_Array[18] + PIO_Time_Array[19]);

  /* Free unneeded memory */
  if(Proc_Num_Side_Sets[iproc] > 0)
    safe_free((void **) &proc_local_ss);

  /*
   * Write out the name of the coordinate axes to the parallel ExodusII
   * files.
   */
  bytes_out += Num_Dim * 8 * sizeof(char);
  tt1 = second();
  if(ex_put_coord_names(mesh_exoid, Coord_Name) < 0) {
    fprintf(stderr, "[%d, %s]: ERROR, could not output coordinate names\n",
            Proc, yo);
    ex_close(mesh_exoid);
    exit(1);
  }
  PIO_Time_Array[20]  = (second() - tt1);
  total_out_time     += PIO_Time_Array[20];

  if (Restart_Info.Flag > 0) {

    tt1 = second();

    if (Restart_Info.NVar_Elem > 0)
      local_tt = Restart_Info.Elem_TT[iproc];
    else
      local_tt = NULL;

    if (Restart_Info.NVar_Nset > 0)
      local_nstt = Restart_Info.Nset_TT[iproc];
    else
      local_nstt = NULL;

    if (Restart_Info.NVar_Sset > 0)
      local_sstt = Restart_Info.Sset_TT[iproc];
    else
      local_sstt = NULL;

    bytes_out += write_var_param(mesh_exoid, max_name_length,
				 Restart_Info.NVar_Glob, Restart_Info.GV_Name,
				 Restart_Info.NVar_Node, Restart_Info.NV_Name, 
				 Restart_Info.NVar_Elem, Restart_Info.EV_Name,  local_tt,
				 Restart_Info.NVar_Nset, Restart_Info.NSV_Name, local_nstt,
				 Restart_Info.NVar_Sset, Restart_Info.SSV_Name, local_sstt);

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

static int write_var_param(int mesh_exoid, int max_name_length,
			   int num_glob, char **gv_names,
			   int num_node, char **nv_names, 
                           int num_elem, char **ev_names, int *local_ebtt,
			   int num_nset, char **ns_names, int *local_nstt,
			   int num_sset, char **ss_names, int *local_sstt)
{
  int     bytes_out=0, error;

  bytes_out += (5 + Num_Elem_Blk * num_elem +
		    Num_Side_Set * num_sset +
		    Num_Node_Set * num_nset) * sizeof(int);
  error = ex_put_all_var_param(mesh_exoid, num_glob, num_node,
			       num_elem, local_ebtt,
			       num_nset, local_nstt,
			       num_sset, local_sstt); 
  check_exodus_error(error, "ex_put_all_var_param");

  if (gv_names != NULL) {
    bytes_out += Restart_Info.NVar_Glob * max_name_length;
    error = ex_put_var_names(mesh_exoid, "g", num_glob, gv_names);
    check_exodus_error(error, "ex_put_var_names");
  }
  if (nv_names != NULL) {
    bytes_out += num_node * max_name_length;
    error = ex_put_var_names(mesh_exoid, "n", num_node, nv_names);
    check_exodus_error(error, "ex_put_var_names");
  }
  if (ev_names != NULL) {
    bytes_out += Restart_Info.NVar_Elem * max_name_length;
    error = ex_put_var_names(mesh_exoid, "e", num_elem, ev_names);
    check_exodus_error(error, "ex_put_var_names");
  }
  if (ns_names != NULL) {
    bytes_out += Restart_Info.NVar_Nset * max_name_length;
    error = ex_put_var_names(mesh_exoid, "m", num_nset, ns_names);
    check_exodus_error(error, "ex_put_var_names");
  }
  if (ss_names != NULL) {
    bytes_out += Restart_Info.NVar_Sset * max_name_length;
    error = ex_put_var_names(mesh_exoid, "s", num_sset, ss_names);
    check_exodus_error(error, "ex_put_var_names");
  }
  return (bytes_out);
}

void write_var_timestep(int exoid, int proc, int time_step, int blk_cnt,
                        int *eb_ids_global, int *ss_ids_global, int *ns_ids_global,
			int io_ws)
{
  int     cnt1, var_num, error, bytes_out=0;
  int     num_nodes, num_elem, eb_num, eb_num_g, var_offset;
  float  *sp_ptr;
  double *dp_ptr;
  void   *var_ptr;

  /* check to see if the io_ws is smaller than the machine precision */
  if (io_ws < sizeof(float)) io_ws = sizeof(float);

  /* output the time */
  if (io_ws == sizeof(float))
    var_ptr = (void *) &(Restart_Info.Time_sp[blk_cnt]);
  else
    var_ptr = (void *) &(Restart_Info.Time_dp[blk_cnt]);

  error = ex_put_time(exoid, time_step, var_ptr);

  check_exodus_error(error, "ex_put_time");

  /* start by outputting the global variables */
  if (Restart_Info.NVar_Glob > 0) {
    bytes_out += Restart_Info.NVar_Glob * io_ws;

    if (io_ws == sizeof(float))
      var_ptr = (void *) Restart_Info.Glob_Vals_sp[blk_cnt];
    else
      var_ptr = (void *) Restart_Info.Glob_Vals_dp[blk_cnt];

    error = ex_put_glob_vars(exoid, time_step, Restart_Info.NVar_Glob,
                             var_ptr);

    check_exodus_error(error, "ex_put_glob_vars");
  }

  if (Restart_Info.NVar_Node > 0) {
    num_nodes = Num_Internal_Nodes[proc] + Num_Border_Nodes[proc] +
                Num_External_Nodes[proc];

    for (var_num=0; var_num < Restart_Info.NVar_Node; var_num++) {

      var_offset = var_num * num_nodes;

      bytes_out += num_nodes * io_ws;

      if (io_ws == sizeof(float))
        var_ptr =
          (void *) &(Restart_Info.Node_Vals_sp[proc][blk_cnt][var_offset]);
      else
        var_ptr =
          (void *) &(Restart_Info.Node_Vals_dp[proc][blk_cnt][var_offset]);

      error = ex_put_nodal_var(exoid, time_step, (var_num+1), num_nodes,
                               var_ptr);

      check_exodus_error(error, "ex_put_nodal_var");
    }
  }

  if (Restart_Info.NVar_Elem > 0) {

    num_elem = Num_Internal_Elems[proc] + Num_Border_Elems[proc];

    for (var_num=0; var_num < Restart_Info.NVar_Elem; var_num++) {
      eb_num_g = 0;

      var_offset = var_num * num_elem;
      if (io_ws == sizeof(float))
        sp_ptr = &(Restart_Info.Elem_Vals_sp[proc][blk_cnt][var_offset]);
      else
        dp_ptr = &(Restart_Info.Elem_Vals_dp[proc][blk_cnt][var_offset]);

      for (eb_num=0; eb_num < Proc_Num_Elem_Blk[proc]; eb_num++) {

        if (io_ws == sizeof(float)) var_ptr = (void *) sp_ptr;
        else                        var_ptr = (void *) dp_ptr;

        /* now I have to find the appropriate entry in the truth table */
	/* can always assume this eb num is greater than the last one */
	for (cnt1=eb_num_g; cnt1<Num_Elem_Blk; cnt1++) {
	  if (Proc_Elem_Blk_Ids[proc][eb_num] == eb_ids_global[cnt1]) {
	    eb_num_g = cnt1;
	    break;
	  }
	}

        if (Restart_Info.Elem_TT[proc]
	    [eb_num_g*Restart_Info.NVar_Elem+var_num]) {
	  
          bytes_out += Proc_Num_Elem_In_Blk[proc][eb_num] * io_ws;
	  
          error = ex_put_elem_var(exoid, time_step, (var_num+1),
                                  Proc_Elem_Blk_Ids[proc][eb_num],
                                  Proc_Num_Elem_In_Blk[proc][eb_num], var_ptr);
	  
          check_exodus_error(error, "ex_put_elem_var");
	}
	/* and now move the variable pointer */

	/* Note that the offsetting here must match the 'var_offset'
	 * treatment in ps_restart.c function read_elem_vars.
	 * Currently, the offset is applied even if the variable does
	 * not exist on a particular block.
	 */
	if (io_ws == sizeof(float))
	  sp_ptr += Proc_Num_Elem_In_Blk[proc][eb_num];
	else
	  dp_ptr += Proc_Num_Elem_In_Blk[proc][eb_num];
      }
    }
  }

  if (Restart_Info.NVar_Sset > 0) {
    int ss_num;
    int ss_num_g = 0;
    num_elem = Proc_SS_Elem_List_Length[proc];
    for (var_num=0; var_num < Restart_Info.NVar_Sset; var_num++) {

      var_offset = var_num * num_elem;
      if (io_ws == sizeof(float))
        sp_ptr = &(Restart_Info.Sset_Vals_sp[proc][blk_cnt][var_offset]);
      else
        dp_ptr = &(Restart_Info.Sset_Vals_dp[proc][blk_cnt][var_offset]);

      for (ss_num=0; ss_num < Proc_Num_Side_Sets[proc]; ss_num++) {

        if (io_ws == sizeof(float)) var_ptr = (void *) sp_ptr;
        else                        var_ptr = (void *) dp_ptr;

        /* now I have to find the appropriate entry in the truth table */
	for (cnt1=0; cnt1<Num_Side_Set; cnt1++) {
	  if (Proc_SS_Ids[proc][ss_num] == ss_ids_global[cnt1]) {
	    ss_num_g = cnt1;
	    break;
	  }
	}
	assert(Proc_SS_Ids[proc][ss_num] == ss_ids_global[ss_num_g]);

        if (Restart_Info.Sset_TT[proc][ss_num_g*Restart_Info.NVar_Sset+var_num]) {
          bytes_out += Proc_SS_Elem_Count[proc][ss_num] * io_ws;
	  
          error = ex_put_sset_var(exoid, time_step, (var_num+1),
                                  Proc_SS_Ids[proc][ss_num],
                                  Proc_SS_Elem_Count[proc][ss_num], var_ptr);
	  
          check_exodus_error(error, "ex_put_sset_var");
	}
	/* and now move the variable pointer */

	/* Note that the offsetting here must match the 'var_offset'
	 * treatment in ps_restart.c function read_elem_vars.
	 * Currently, the offset is applied even if the variable does
	 * not exist on a particular block.
	 */
	if (io_ws == sizeof(float))
	  sp_ptr += Proc_SS_Elem_Count[proc][ss_num];
	else
	  dp_ptr += Proc_SS_Elem_Count[proc][ss_num];
      }
    }
  }

  if (Restart_Info.NVar_Nset > 0) {
    int ns_num;
    int ns_num_g = 0;
    num_elem = Proc_NS_List_Length[proc];
    for (var_num=0; var_num < Restart_Info.NVar_Nset; var_num++) {

      var_offset = var_num * num_elem;
      if (io_ws == sizeof(float))
        sp_ptr = &(Restart_Info.Nset_Vals_sp[proc][blk_cnt][var_offset]);
      else
        dp_ptr = &(Restart_Info.Nset_Vals_dp[proc][blk_cnt][var_offset]);

      for (ns_num=0; ns_num < Proc_Num_Node_Sets[proc]; ns_num++) {

        if (io_ws == sizeof(float)) var_ptr = (void *) sp_ptr;
        else                        var_ptr = (void *) dp_ptr;

        /* now I have to find the appropriate entry in the truth table */
	for (cnt1=0; cnt1<Num_Node_Set; cnt1++) {
	  if (Proc_NS_Ids[proc][ns_num] == ns_ids_global[cnt1]) {
	    ns_num_g = cnt1;
	    break;
	  }
	}
	assert(Proc_NS_Ids[proc][ns_num] == ns_ids_global[ns_num_g]);

        if (Restart_Info.Nset_TT[proc][ns_num_g*Restart_Info.NVar_Nset+var_num]) {
          bytes_out += Proc_NS_Count[proc][ns_num] * io_ws;
	  
          error = ex_put_nset_var(exoid, time_step, (var_num+1),
                                  Proc_NS_Ids[proc][ns_num],
                                  Proc_NS_Count[proc][ns_num], var_ptr);
	  
          check_exodus_error(error, "ex_put_nset_var");
	}
	/* and now move the variable pointer */

	/* Note that the offsetting here must match the 'var_offset'
	 * treatment in ps_restart.c function read_elem_vars.
	 * Currently, the offset is applied even if the variable does
	 * not exist on a particular block.
	 */
	if (io_ws == sizeof(float))
	  sp_ptr += Proc_NS_Count[proc][ns_num];
	else
	  dp_ptr += Proc_NS_Count[proc][ns_num];
      }
    }
  }
}

void reverse_map(int *global, int p01, int gsize,
		 int *glmap, int *index, int msize,
		 int *mapout)
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
  

  int i2, i3;
  int *tmp_index = (int *) array_alloc(__FILE__, __LINE__, 1, gsize,
				       sizeof(int));

  /* Initialize index array */
  for (i2=0; i2 < gsize; i2++)
    tmp_index[i2] = i2;
      
  /*
   * Sort the 'global' array via the index array 'tmp_index' 
   */
  gds_iqsort(global, tmp_index, gsize);
  

  i3 = 0;
  if (index != NULL) {
    for(i2 = 0; i2 < gsize; i2++) {
      int gval = global[tmp_index[i2]] + p01;

      while (glmap[index[i3]] < gval)
	i3++;

      assert(glmap[index[i3]] == gval);

      mapout[tmp_index[i2]] = index[i3] + 1;
    }
  } else {
    for(i2 = 0; i2 < gsize; i2++) {
      int gval = global[tmp_index[i2]] + p01;

      while (glmap[i3] < gval)
	i3++;

      assert(glmap[i3] == gval);

      mapout[tmp_index[i2]] = i3 + 1;
    }
  }
  safe_free((void **) &tmp_index);
}
