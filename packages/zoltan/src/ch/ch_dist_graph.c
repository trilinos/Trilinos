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
static char *cvs_chaco_dist_graph_id = "$Id$";
#endif

#include <stdio.h>
#include <mpi.h>
#include "lbi_const.h"
#include "all_allo_const.h"
#include "ch_input_const.h"

/*
 * Distribute a graph from one processor to all processors.
 * The ParMetis format is used for the distributed graph.
 * The memory for the graph on the host node is freed
 * and fresh memory is allocated for the distr. graph.
 *
 * Geometry information is not supported in this version.
 */

int chaco_dist_graph(
  MPI_Comm Comm,		/* MPI Communicator */
  int     host_proc,		/* processor where all the data is initially */
  int     *nvtxs,		/* number of vertices in graph */
  int     **vtxdist, 		/* vertex distribution data */
  int     **xadj,		/* start of edge list for each vertex */
  int     **adjncy,		/* edge list data */
  int     **vwgts,		/* vertex weight list data */
  float   **ewgts 		/* edge weight list data */
)
{
  extern int CHECK_INPUT;	/* print warnings or not? */
  extern int DEBUG_INPUT;	/* echo that input file read successful? */
  extern int DEBUG_TRACE;	/* trace main execution path */

  int nprocs, myproc, i, n, p, nedges, nsend, rest, flag;
  int offset, use_vwgts, use_ewgts;
  int *old_xadj, *old_adjncy, *old_vwgts, *size;
  float *old_ewgts;
  MPI_Status status;

  /* Determine number of processors and my rank. */
  MPI_Comm_size (Comm, &nprocs );
  MPI_Comm_rank (Comm, &myproc );

  if (DEBUG_TRACE > 0) {
    if (myproc==0) printf("<Entering ch_dist_graph>\n");
  }

  /* Initialize */
  use_ewgts = !(*ewgts == NULL);
  use_vwgts = !(*vwgts == NULL);

  /* Set up vtxdist data */
  flag = 0;
  if (*vtxdist == NULL){
    if (myproc==host_proc) flag=1;
    *vtxdist = (int *)LB_SMALLOC((nprocs+1) * sizeof(int));
    if (*vtxdist == NULL) return 1;
  }
  if (myproc==host_proc){
    if (flag){
      /* Calculate uniform vertex distribution */
      (*vtxdist)[0] = 0;
      rest = *nvtxs;
      for (i=0; i<nprocs; i++){
        n = rest/(nprocs-i);
        (*vtxdist)[i+1] = (*vtxdist)[i] + n;
        rest -= n;
      }
    }
  }
  /* Broadcast vtxdist to all procs */
  MPI_Bcast( *vtxdist, nprocs+1, MPI_INT, host_proc, Comm);
  MPI_Bcast( &use_vwgts, 1, MPI_INT, host_proc, Comm);
  MPI_Bcast( &use_ewgts, 1, MPI_INT, host_proc, Comm);
  
  /* Store pointers to original data */
  if (myproc == host_proc){
    old_xadj   = *xadj;
    old_adjncy = *adjncy;
  }

  /* Allocate space for new distributed graph data */
  n = (*vtxdist)[myproc+1]- (*vtxdist)[myproc]; /* local # of nodes */
  *nvtxs = n;
  *xadj = (int *)LB_SMALLOC((n+1)*sizeof(int));
  if (*xadj == NULL) return 1;
  if (use_vwgts){
    old_vwgts = *vwgts;
    *vwgts = (int *)LB_SMALLOC((n+1)*sizeof(int));
    if (*vwgts == NULL) return 1;
  }

  /* Distribute graph data to all procs */

  /* First send xadj data */
  if (myproc== host_proc){
    size = (int *)LB_SMALLOC(nprocs*sizeof(int));
    if (size == NULL) return 1;
    for (p=0; p<nprocs; p++){
      size[p] = old_xadj[(*vtxdist)[myproc+1]]-old_xadj[(*vtxdist)[myproc]];
      offset = (*vtxdist)[myproc];
      nsend = (*vtxdist)[myproc+1] - offset;
      if (p == myproc){
        /* do a local copy */
        for (i=0; i<=nsend; i++)
          (*xadj)[i] = old_xadj[offset+i];
        if (use_vwgts)
          for (i=0; i<nsend; i++)
            (*vwgts)[i] = old_vwgts[offset+i];
      }
      else {
        MPI_Send( &old_xadj[offset], nsend+1, MPI_INT, p, 1, Comm);
        if (use_vwgts)
          MPI_Send( &old_vwgts[offset], nsend, MPI_INT, p, 2, Comm);
      }
    }
  }
  else {
    MPI_Recv (*xadj, (*nvtxs)+1, MPI_INT, host_proc, 1, Comm, &status);
    if (use_vwgts)
      MPI_Recv (*vwgts, *nvtxs, MPI_INT, host_proc, 2, Comm, &status);
  }

  /* Adjust xadj values */
  offset = (*xadj)[0];
  for (i=0; i<= *nvtxs; i++)
    (*xadj)[i] -= offset;
  nedges = (*xadj)[ *nvtxs];
  *adjncy = (int *)LB_SMALLOC(nedges * sizeof (int));
  if (*adjncy == NULL) return 1;
  if (use_ewgts){
    old_ewgts = *ewgts;
    *ewgts = (float *)LB_SMALLOC(nedges * sizeof (float));
    if (*ewgts == NULL) return 1;
  }

  /* Next send adjncy data */
  if (myproc== host_proc){
    for (p=0; p<nprocs; p++){
      offset = old_xadj[(*vtxdist)[p]];
      if (p == myproc){
        /* do a local copy */
        for (i=0; i<size[p]; i++)
          (*adjncy)[i] = old_adjncy[offset+i];
        if (use_ewgts)
          for (i=0; i<size[p]; i++)
            (*ewgts)[i] = old_ewgts[offset+i];
      }
      else {
        MPI_Send( &old_adjncy[offset], size[p], MPI_INT, p, 3, Comm);
        if (use_ewgts)
          MPI_Send( &old_ewgts[offset], size[p], MPI_FLOAT, p, 4, Comm);
      }
    }
  }
  else {
    MPI_Recv (*adjncy, nedges, MPI_INT, host_proc, 3, Comm, &status);
    if (use_ewgts)
      MPI_Recv (*ewgts, nedges, MPI_FLOAT, host_proc, 4, Comm, &status);
  }

  /* Free space on host proc */
  if (myproc == host_proc){
    LB_safe_free((void **) &old_xadj);
    LB_safe_free((void **) &old_adjncy);
    if (use_vwgts)
      LB_safe_free((void **) &old_vwgts);
    if (use_ewgts)
      LB_safe_free((void **) &old_ewgts);
  }
   
  if (DEBUG_INPUT > 0) {
    if (myproc==0) printf("Done distributing graph \n");
  }
  return 0;
}
