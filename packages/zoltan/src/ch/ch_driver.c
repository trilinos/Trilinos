/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_ch_driver_id = "$Id$";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "lb_const.h"
#include "lbi_const.h"
#include "all_allo_const.h"
#include "ch_input_const.h"
#include "ch_input.h"
#include "ch_driver_const.h"

/*
 * Simple driver that reads a Chaco file and calls Zoltan,
 * Currently ParMetis is used.
 *
 */

main(int argc, char **argv)
{
  LB  *lb;
  int i, nprocs, myproc, nvtxs, imax, changes;
  int *start, *adjncy, *vwgts, *vtxdist;
  int *imp_procs, *exp_procs, num_imp, num_exp;
  float version, *ewgts;
  double params[5];
  char s[80];
  LB_GID *imp_gids, *exp_gids;
  LB_LID *imp_lids, *exp_lids;
  Graph graph;

  /* Initialize MPI */
  MPI_Init( &argc, &argv );
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs );
  MPI_Comm_rank (MPI_COMM_WORLD, &myproc );


  /* Read Chaco file */
  if (myproc==0) {
    if (DEBUG_TRACE) 
      printf("Debug: Reading Chaco file...\n");
	    chaco_input_graph(stdin, "STDIN", &start, &adjncy, &nvtxs,
                              &vwgts, &ewgts);
  }

  if (DEBUG_TRACE) 
    printf("[%1d] Debug: distributing graph...\n", myproc);

  /* Distribute graph */
  vtxdist = NULL;
  chaco_dist_graph(MPI_COMM_WORLD, 0, &nvtxs, &vtxdist, &start, &adjncy, 
                   &vwgts, &ewgts);

  if (DEBUG_INPUT){
    printf("[%1d] Debug: nvtxs=%d, nedges=%d\n", myproc, nvtxs, start[nvtxs]);
    imax = (nvtxs>10) ? 10 : nvtxs;
    for (i=0; i<imax; i++) sprintf(&s[3*i], "%3d", start[i]);
    printf("[%1d] Debug: start = %s\n", myproc, s);
    imax = (start[nvtxs]>10) ? 10 : start[nvtxs];
    for (i=0; i<imax; i++) sprintf(&s[3*i], "%3d", adjncy[i]);
    printf("[%1d] Debug: adjncy = %s\n", myproc, s);
  }

  graph.lnvtxs = nvtxs;
  graph.vtxdist = vtxdist;
  graph.start = start;
  graph.adjncy = adjncy;
  graph.vwgts = vwgts;
  graph.ewgts = ewgts;
  graph.nprocs = nprocs;
  graph.myproc = myproc;

  /* Make sure all processors are sync'ed before we call Zoltan */
  MPI_Barrier(MPI_COMM_WORLD);

  if (DEBUG_TRACE) 
    printf("[%1d] Debug: Initializing Zoltan...\n", myproc);

  /* Initialize  Zoltan */
  LB_Initialize(argc, argv, &version);
/* BAH: need to put in new parameter initialization. */
  /*LB_Initialize_Params_Array(params);*/

  if (DEBUG_TRACE) 
    printf("[%1d] Debug: Creating Zoltan LB object...\n",myproc);
  lb = LB_Create_Object(MPI_COMM_WORLD);

  /* Register Zoltan functions */
  if (DEBUG_TRACE) 
    printf("[%1d] Debug: Registering LB functions...\n", myproc);
  LB_Set_Method(lb, "ParMetis_Part", params);
  LB_Set_Fn(lb, LB_NUM_OBJ_FN_TYPE, (void *)num_obj_fn, &graph);
  LB_Set_Fn(lb, LB_OBJ_LIST_FN_TYPE, (void *)list_obj_fn, &graph);
  LB_Set_Fn(lb, LB_OBJ_WEIGHT_FN_TYPE, (void *)obj_weight_fn, &graph);
  LB_Set_Fn(lb, LB_NUM_EDGES_FN_TYPE, (void *)num_edges_fn, &graph);
  LB_Set_Fn(lb, LB_EDGE_LIST_FN_TYPE, (void *)edge_list_fn, &graph);

  if (DEBUG_TRACE) 
    printf("[%1d] Debug: Calling Zoltan LB_Balance ...\n", myproc);

  /* Call Zoltan */
  LB_Balance(lb, &changes, &num_imp, &imp_gids, &imp_lids, &imp_procs,
    &num_exp, &exp_gids, &exp_lids, &exp_procs);

  if (DEBUG_TRACE) 
    printf("[%1d] Debug: Finished Zoltan LB_Balance! \n", myproc);

  /* Check results */
  if (myproc==0){
    printf("Proc %1d: #import = %d, #export = %d\n", myproc, 
      num_imp, num_exp);
  }

  /* End MPI */
  MPI_Finalize();

  if (DEBUG_TRACE) 
    printf("[%1d] Debug: Driver program finished! \n", myproc);

  return LB_OK;
}

