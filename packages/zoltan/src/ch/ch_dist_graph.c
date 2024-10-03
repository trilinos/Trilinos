// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include "ch_input_const.h"
#include "dr_err_const.h"
#include "dr_const.h"
#include "dr_externs.h"
#include "dr_util_const.h"
#include "ch_init_dist_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/*
 * Distribute a graph from one processor to all processors.
 * The ParMetis format is used for the distributed graph.
 * The memory for the graph on the host node is freed
 * and fresh memory is allocated for the distributed graph.
 */

int chaco_dist_graph(
  MPI_Comm comm,		/* MPI Communicator */
  PARIO_INFO_PTR pio_info,      /* Parallel IO info */
  int     host_proc,		/* processor where all the data is initially */
  int     *gnvtxs,		/* number of vertices in global graph */
  int     *nvtxs,		/* number of vertices in local graph */
  int     **xadj,		/* start of edge list for each vertex */
  int     **adjncy,		/* edge list data */
  int     *vwgt_dim,            /* number of weights per vertex */
  float   **vwgts,		/* vertex weight list data */
  int     *ewgt_dim,            /* number of weights per edge */
  float   **ewgts,		/* edge weight list data */
  int     *ndim,                /* dimension of the geometry */
  float   **x,                  /* x-coordinates of the vertices */
  float   **y,                  /* y-coordinates of the vertices */
  float   **z,                  /* z-coordinates of the vertices */
  short   **assignments         /* assignments from Chaco file; may be NULL */
)
{
  const char *yo = "chaco_dist_graph";
  int nprocs, myproc, i, j, k, n, p, nedges, nsend, max_nvtxs, v, adj_cnt;
  int offset, use_graph, nvtx_edges;
  int *old_xadj = NULL, *old_adjncy = NULL, *size = NULL;
  int *send_xadj = NULL, *send_adjncy = NULL;
  int *vtx_list = NULL;
  float *old_x = NULL, *old_y = NULL, *old_z = NULL;
  float *send_x = NULL, *send_y = NULL, *send_z = NULL;
  float *old_vwgts = NULL, *old_ewgts = NULL;
  float *send_vwgts = NULL, *send_ewgts = NULL;
  MPI_Status status;

  /* Determine number of processors and my rank. */
  MPI_Comm_size (comm, &nprocs );
  MPI_Comm_rank (comm, &myproc );

  DEBUG_TRACE_START(myproc, yo);

  /* Initialize */
  use_graph = (*xadj  != NULL);
 
  /* Handle serial case and return. */
  if (nprocs == 1) {
    /* Set values expected to be returned by this function. */
    /* All array pointers are unchanged.                    */
    *gnvtxs = *nvtxs;
    /* Initialize distribution, so other routines using it work. */
    ch_dist_init(nprocs, *gnvtxs, pio_info, assignments, host_proc, comm);

    return 1;
  }

  /* Broadcast to all procs */
  MPI_Bcast( vwgt_dim, 1, MPI_INT, host_proc, comm);
  MPI_Bcast( ewgt_dim, 1, MPI_INT, host_proc, comm);
  MPI_Bcast( &use_graph, 1, MPI_INT, host_proc, comm);
  MPI_Bcast( ndim, 1, MPI_INT, host_proc, comm);
  MPI_Bcast( nvtxs, 1, MPI_INT, host_proc, comm);
  *gnvtxs = *nvtxs;

  /* Initialize the chaco distribution on all processors */
  ch_dist_init(nprocs, *gnvtxs, pio_info, assignments, host_proc, comm);
  
  /* Store pointers to original data */
  if (myproc == host_proc) {
    old_xadj   = *xadj;
    *xadj      = NULL;
    old_adjncy = *adjncy;
    *adjncy    = NULL;
    old_x      = *x;
    old_y      = *y;
    old_z      = *z;
  }

  /* Allocate space for new distributed graph data */
  n = *nvtxs = ch_dist_num_vtx(myproc, *assignments);

  if (use_graph) {
    *xadj = (int *) malloc((n+1)*sizeof(int));
    if (*xadj == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
  }
  if (*vwgt_dim){
    old_vwgts = *vwgts;
    *vwgts = NULL;
    if (n > 0) {
      *vwgts = (float *) malloc(n*(*vwgt_dim)*sizeof(float));
      if (*vwgts == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
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


  /* 
   * Distribute vertex data (xadj, coordinates, etc ) to all procs 
   */

  if (myproc == host_proc){

    /* Allocate space for send buffers  (size = max num vtx per proc ) */
    max_nvtxs = ch_dist_max_num_vtx(*assignments);
    if (use_graph) {
      send_xadj = (int *) malloc((max_nvtxs+1)*sizeof(int));
      if (send_xadj == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
    }
    if (*vwgt_dim) {
      send_vwgts = (float *) malloc(max_nvtxs*(*vwgt_dim)*sizeof(float));
      if (send_vwgts == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
    }
    if (*ndim > 0) {
      send_x = (float *) malloc(max_nvtxs*sizeof(float));
      if (send_x == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
      if (*ndim > 1) {
        send_y = (float *) malloc(max_nvtxs*sizeof(float));
        if (send_y == NULL) {
          Gen_Error(0, "fatal: insufficient memory");
          return 0;
        }
        if (*ndim > 2) {
          send_z = (float *) malloc(max_nvtxs*sizeof(float));
          if (send_z == NULL) {
            Gen_Error(0, "fatal: insufficient memory");
            return 0;
          }
        }
      }
    }

    /* Allocate space for list of vertices on a given processor */
    vtx_list = (int *) malloc(max_nvtxs*sizeof(int));
    if (vtx_list == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }

    /* Allocate array to accumulate number of edges to be sent to each proc. */
    if (use_graph) {
      size = (int *) malloc(nprocs*sizeof(int));
      if (size == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
    }

    /* For each processor, gather its vertex information and send it. */
    for (p = 0; p < nprocs; p++){
      if (use_graph) size[p] = 0;

      /* Get list of vertices to be assigned to processor p */
      ch_dist_vtx_list(vtx_list, &nsend, p, *assignments);

      if (p == myproc){

        /* Loop over vertices assigned to myproc; copy the vertex */
        /* data into local arrays.                                */

        if (use_graph) (*xadj)[0] = 0;
        for (i = 0; i < nsend; i++) {
          v = vtx_list[i];
          if (use_graph) {
            size[p] += old_xadj[v+1]-old_xadj[v];
            (*xadj)[i+1] = (*xadj)[i] + old_xadj[v+1] - old_xadj[v];
          }

          if (*vwgt_dim){
            for (j=0; j<*vwgt_dim; j++)
              (*vwgts)[i*(*vwgt_dim)+j] = old_vwgts[v*(*vwgt_dim)+j];
          }

          if (*ndim > 0) {
            (*x)[i] = old_x[v];
            if (*ndim > 1) {
              (*y)[i] = old_y[v];
              if (*ndim > 2) {
                (*z)[i] = old_z[v];
              }
            }
          }
        }
      }
      else {

        /* Loop over vertices assigned to proc p to gather */
        /* vertex data into send buffers                   */

        if (use_graph) send_xadj[0] = 0;
        for (i = 0; i < nsend; i++) {
          v = vtx_list[i];
          if (use_graph) {
            size[p] += old_xadj[v+1]-old_xadj[v];
            send_xadj[i+1] = send_xadj[i] + old_xadj[v+1] - old_xadj[v];
          }
          if (*vwgt_dim){
            for (j=0; j<*vwgt_dim; j++)
              send_vwgts[i*(*vwgt_dim)+j] = old_vwgts[v*(*vwgt_dim)+j];
          }
          if (*ndim > 0) {
            send_x[i] = old_x[v];
            if (*ndim > 1) {
              send_y[i] = old_y[v];
              if (*ndim > 2) {
                send_z[i] = old_z[v];
              }
            }
          }
        }

        /* Send vertex data to proc p. */
        if (use_graph)
          MPI_Send(send_xadj, nsend+1, MPI_INT, p, 1, comm);
        if (*vwgt_dim)
          MPI_Send(send_vwgts, nsend*(*vwgt_dim), MPI_FLOAT, p, 2, comm);
        if (*ndim > 0) {
          MPI_Send(send_x, nsend, MPI_FLOAT, p, 3, comm);
          if (*ndim > 1) {
            MPI_Send(send_y, nsend, MPI_FLOAT, p, 4, comm);
            if (*ndim > 2) {
              MPI_Send(send_z, nsend, MPI_FLOAT, p, 5, comm);
            }
          }
        }
      }
    }
    safe_free((void **)(void *) &send_xadj);
    safe_free((void **)(void *) &send_vwgts);
    safe_free((void **)(void *) &send_x);
    safe_free((void **)(void *) &send_y);
    safe_free((void **)(void *) &send_z);
  }
  else {
    /* host_proc != myproc; receive vertex data from host_proc */
    if (use_graph)
      MPI_Recv (*xadj, (*nvtxs)+1, MPI_INT, host_proc, 1, comm, &status);
    if (*vwgt_dim)
      MPI_Recv (*vwgts, (*nvtxs)*(*vwgt_dim), MPI_FLOAT, host_proc, 2, comm, &status);
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


  /* 
   * Distribute edge data to all procs 
   */

  if (use_graph) {

    if (*ewgt_dim) {
      old_ewgts = *ewgts;
      *ewgts = NULL;
    }

    /* Allocate space for edge data */
    nedges = (*xadj)[ *nvtxs];
    if (nedges > 0) {
      *adjncy = (int *) malloc(nedges * sizeof (int));
      if (*adjncy == NULL) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
      if (*ewgt_dim){
        *ewgts = (float *) malloc(nedges*(*ewgt_dim) * sizeof (float));
        if (*ewgts == NULL) {
          Gen_Error(0, "fatal: insufficient memory");
          return 0;
        }
      }
    }

    /* Gather and send/receive edge data */

    if (myproc == host_proc){

      /* For each processor, gather its edge data and send it. */
      for (p = 0; p < nprocs; p++){
        if (size[p] == 0) continue;

        /* Get list of vertices to be assigned to processor p */
        ch_dist_vtx_list(vtx_list, &nsend, p, *assignments);

        adj_cnt = 0;
        if (p == myproc) {

          /* Loop over vertices assigned to myproc copy the edge */
          /* data into local arrays.                             */

          for (i = 0; i < nsend; i++) {
            v = vtx_list[i];
            offset = old_xadj[v];
            nvtx_edges = old_xadj[v+1] - old_xadj[v];
            for (j = 0; j < nvtx_edges; j++) {
              (*adjncy)[adj_cnt] = old_adjncy[offset+j];
              if (*ewgt_dim){ 
                for (k=0; k<*ewgt_dim; k++)
                  (*ewgts)[adj_cnt*(*ewgt_dim)+k] = old_ewgts[(offset+j)*(*ewgt_dim)+k];
              }
              adj_cnt++;
            }
          }
        }
        else { /* p != myproc */

          /* allocate send buffers; size = num edges to send to proc p */
          nvtx_edges = 0;
          for (i = 0; i < nsend; i++) {
            v = vtx_list[i];
            nvtx_edges += old_xadj[v+1] - old_xadj[v];
          }
          send_adjncy = (int *) malloc(nvtx_edges * sizeof(int));
          if (send_adjncy == NULL) {
            Gen_Error(0, "fatal: insufficient memory");
            return 0;
          }
          if (*ewgt_dim) {
            send_ewgts = (float *) malloc(nvtx_edges*(*ewgt_dim) * sizeof(float));
            if (send_ewgts == NULL) {
              Gen_Error(0, "fatal: insufficient memory");
              return 0;
            }
          }

          /* Loop over vertices assigned to proc p to gather */
          /* edge data into send buffers                     */

          for (i = 0; i < nsend; i++) {
            v = vtx_list[i];
            offset = old_xadj[v];
            nvtx_edges = old_xadj[v+1] - old_xadj[v];
            for (j = 0; j < nvtx_edges; j++) {
              send_adjncy[adj_cnt] = old_adjncy[offset+j];
              if (*ewgt_dim){
                for (k=0; k<*ewgt_dim; k++)
                  send_ewgts[adj_cnt*(*ewgt_dim)+k] = old_ewgts[(offset+j)*(*ewgt_dim)+k];
              }
              adj_cnt++;
            }
          }
          /* Send edge data to proc p. */
          MPI_Send(send_adjncy, size[p], MPI_INT, p, 6, comm);
          if (*ewgt_dim)
            MPI_Send(send_ewgts, size[p]*(*ewgt_dim), MPI_FLOAT, p, 7, comm);
          safe_free((void **)(void *) &send_adjncy);
          safe_free((void **)(void *) &send_ewgts);
        }
      }
    }
    else {
      /* host_proc != myproc; receive edge data from host_proc */
      if (nedges > 0) {
        MPI_Recv (*adjncy, nedges, MPI_INT, host_proc, 6, comm, &status);
        if (*ewgt_dim)
          MPI_Recv (*ewgts, nedges*(*ewgt_dim), MPI_FLOAT, host_proc, 7, comm, &status);
      }
    }
  }

  /* Free space on host proc */
  if (myproc == host_proc){
    safe_free((void **)(void *) &old_xadj);
    safe_free((void **)(void *) &old_adjncy);
    safe_free((void **)(void *) &old_vwgts);
    safe_free((void **)(void *) &old_ewgts);

    safe_free((void **)(void *) &old_x);
    safe_free((void **)(void *) &old_y);
    safe_free((void **)(void *) &old_z);
    safe_free((void **)(void *) &vtx_list);
  }
  if (size != NULL) safe_free((void **)(void *) &size);
   
  DEBUG_TRACE_END(myproc, yo);
  return 1;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
