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
#include <stdlib.h>
#include <mpi.h>
#include "lbi_const.h"
#include "all_allo_const.h"
#include "ch_input_const.h"
#include "dr_err_const.h"

/*
 * Distribute a graph from one processor to all processors.
 * The ParMetis format is used for the distributed graph.
 * The memory for the graph on the host node is freed
 * and fresh memory is allocated for the distr. graph.
 */

int chaco_dist_graph(
  MPI_Comm comm,		/* MPI Communicator */
  int     host_proc,		/* processor where all the data is initially */
  int     *nvtxs,		/* number of vertices in graph */
  int     **vtxdist, 		/* vertex distribution data */
  int     **xadj,		/* start of edge list for each vertex */
  int     **adjncy,		/* edge list data */
  int     **vwgts,		/* vertex weight list data */
  float   **ewgts,		/* edge weight list data */
  int     *ndim,                /* dimension of the geometry */
  float   **x,                  /* x-coordinates of the vertices */
  float   **y,                  /* y-coordinates of the vertices */
  float   **z                   /* z-coordinates of the vertices */
)
{
  int nprocs, myproc, i, n, p, nedges, nsend, rest;
  int offset, use_vwgts, use_ewgts, use_graph;
  int *old_xadj = NULL, *old_adjncy = NULL, *old_vwgts = NULL, *size = NULL;
  float *old_x = NULL, *old_y = NULL, *old_z = NULL;
  float *old_ewgts = NULL;
  MPI_Status status;

  /* Determine number of processors and my rank. */
  MPI_Comm_size (comm, &nprocs );
  MPI_Comm_rank (comm, &myproc );

  if (DEBUG_TRACE > 0) {
    if (myproc==0) printf("<Entering ch_dist_graph>\n");
  }

  /* Initialize */
  use_ewgts = (*ewgts != NULL);
  use_vwgts = (*vwgts != NULL);
  use_graph = (*xadj  != NULL);
 
  /* Broadcast to all procs */
  MPI_Bcast( &use_vwgts, 1, MPI_INT, host_proc, comm);
  MPI_Bcast( &use_ewgts, 1, MPI_INT, host_proc, comm);
  MPI_Bcast( &use_graph, 1, MPI_INT, host_proc, comm);
  MPI_Bcast( ndim, 1, MPI_INT, host_proc, comm);
  MPI_Bcast( nvtxs, 1, MPI_INT, host_proc, comm);
  
  /* Set up vtxdist data */
  if (*vtxdist == NULL){
    *vtxdist = (int *) malloc((nprocs+1) * sizeof(int));
    if (*vtxdist == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
  }
  /* Calculate uniform vertex distribution */
  (*vtxdist)[0] = 0;
  rest = *nvtxs;
  for (i=0; i<nprocs; i++){
    n = rest/(nprocs-i);
    (*vtxdist)[i+1] = (*vtxdist)[i] + n;
    rest -= n;
  }

  /* Store pointers to original data */
  if (myproc == host_proc) {
    old_xadj   = *xadj;
    old_adjncy = *adjncy;
    old_x      = *x;
    old_y      = *y;
    old_z      = *z;
  }

  /* Allocate space for new distributed graph data */
  n = (*vtxdist)[myproc+1]- (*vtxdist)[myproc]; /* local # of nodes */
  *nvtxs = n;
  if (use_graph) {
    *xadj = (int *) malloc((n+1)*sizeof(int));
    if (*xadj == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
  }
  if (use_vwgts){
    old_vwgts = *vwgts;
    *vwgts = (int *) malloc(n*sizeof(int));
    if (*vwgts == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
  }
  if (*ndim > 0) {
    *x = (float *)  malloc(n*sizeof(float));
    if (*ndim > 1) {
      *y = (float *)  malloc(n*sizeof(float));
      if (*ndim > 2) {
        *z = (float *)  malloc(n*sizeof(float));
      }
    }
  }


  /* Distribute graph data to all procs */
  /* Send xadj and coordinates, if appropriate */

  if (myproc == host_proc){
    if (use_graph) {
      size = (int *) malloc(nprocs*sizeof(int));
      if (size == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
    }

    for (p = 0; p < nprocs; p++){
      if (use_graph)
        size[p] = old_xadj[(*vtxdist)[p+1]]-old_xadj[(*vtxdist)[p]];
      offset = (*vtxdist)[p];
      nsend = (*vtxdist)[p+1] - offset;
      if (p == myproc){
        /* do a local copy */
        if (use_graph)
          for (i=0; i<=nsend; i++)
            (*xadj)[i] = old_xadj[offset+i];
        for (i=0; i<nsend; i++) {
          if (use_vwgts)
            (*vwgts)[i] = old_vwgts[offset+i];
          if (*ndim > 0) {
            (*x)[i] = old_x[offset+i];
            if (*ndim > 1) {
              (*y)[i] = old_y[offset+i];
              if (*ndim > 2) {
                (*z)[i] = old_z[offset+i];
              }
            }
          }
        }
      }
      else {
        if (use_graph)
          MPI_Send( &old_xadj[offset], nsend+1, MPI_INT, p, 1, comm);
        if (use_vwgts)
          MPI_Send( &old_vwgts[offset], nsend, MPI_INT, p, 2, comm);
        if (*ndim > 0) {
          MPI_Send(&old_x[offset], nsend, MPI_FLOAT, p, 3, comm);
          if (*ndim > 1) {
            MPI_Send(&old_y[offset], nsend, MPI_FLOAT, p, 4, comm);
            if (*ndim > 2) {
              MPI_Send(&old_z[offset], nsend, MPI_FLOAT, p, 5, comm);
            }
          }
        }
      }
    }
  }
  else {
    if (use_graph)
      MPI_Recv (*xadj, (*nvtxs)+1, MPI_INT, host_proc, 1, comm, &status);
    if (use_vwgts)
      MPI_Recv (*vwgts, *nvtxs, MPI_INT, host_proc, 2, comm, &status);
    if (*ndim > 0) {
      MPI_Recv(*x, *nvtxs,  MPI_FLOAT, host_proc, 3, comm, &status);
      if (*ndim > 1) {
        MPI_Recv(*y, *nvtxs, MPI_FLOAT, host_proc, 4, comm, &status);
        if (*ndim > 2) {
          MPI_Recv(*z, *nvtxs, MPI_FLOAT, host_proc, 5, comm, &status);
        }
      }
    }
  }

  if (use_graph) {
    /* Adjust xadj values */
    offset = (*xadj)[0];
    for (i=0; i<= *nvtxs; i++)
      (*xadj)[i] -= offset;

    /* Distribute adjacency info */
    nedges = (*xadj)[ *nvtxs];
    if (nedges > 0) {
      *adjncy = (int *) malloc(nedges * sizeof (int));
      if (*adjncy == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
      if (use_ewgts){
        old_ewgts = *ewgts;
        *ewgts = (float *) malloc(nedges * sizeof (float));
        if (*ewgts == NULL) {
          Gen_Error(0, "fatal: insufficient memory");
          return 0;
        }
      }
    }

    /* Next send adjncy data */
    if (myproc== host_proc){
      for (p=0; p<nprocs; p++){
        if (size[p] == 0) continue;
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
          MPI_Send(&old_adjncy[offset], size[p], MPI_INT, p, 6, comm);
          if (use_ewgts)
            MPI_Send(&old_ewgts[offset], size[p], MPI_FLOAT, p, 7, comm);
        }
      }
    }
    else {
      if (nedges > 0) {
        MPI_Recv (*adjncy, nedges, MPI_INT, host_proc, 6, comm, &status);
        if (use_ewgts)
          MPI_Recv (*ewgts, nedges, MPI_FLOAT, host_proc, 7, comm, &status);
      }
    }
  }

  /* Free space on host proc */
  if (myproc == host_proc){
    if (old_xadj != NULL)
      free(old_xadj);
    if (old_adjncy != NULL)
      free(old_adjncy);
    if (use_vwgts)
      free(old_vwgts);
    if (use_ewgts)
      free(old_ewgts);

    if (*ndim > 0) {
      free(old_x);
      if (*ndim > 1) {
        free(old_y);
        if (*ndim > 2) {
          free(old_z);
        }
      }
    }
  }
   
  if (DEBUG_INPUT > 0) {
    if (myproc==0) printf("Done distributing graph \n");
  }
  return 1;
}
