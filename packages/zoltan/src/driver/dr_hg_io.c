/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_util_const.h"
#include "dr_par_util_const.h"
#include "dr_err_const.h"
#include "dr_output_const.h"
#include "dr_elem_util_const.h"
#include "dr_maps_const.h"
#include "ch_input_const.h"
#include "ch_init_dist_const.h"
#include "dr_hg_readfile.h"

#ifndef MAX_STR_LENGTH
#define MAX_STR_LENGTH 80
#endif

static int dist_hyperedges(MPI_Comm comm, PARIO_INFO_PTR, int, int, int *,
                           int *, int **, int **, int *, float **, short **);

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int read_hypergraph_file(
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh
)
{
  /* Local declarations. */
  char  *yo = "read_hypergraph_file";
  char   cmesg[256];

  int    i, gnvtxs; 
  int    nvtxs = 0, gnhedges = 0, nhedges = 0, npins = 0;
  int    vwgt_dim=0, hewgt_dim=0;
  int   *hindex = NULL, *hvertex = NULL;
  float *hewgts = NULL, *vwgts = NULL;
  FILE  *fp;
  int base = 0;   /* Smallest vertex number; usually zero or one. */
  char filename[256];

  /* Variables that allow graph-based functions to be reused. */
  /* If no chaco.graph or chaco.coords files exist, values are NULL or 0, 
   * since graph is not being built. If chaco.graph and/or chaco.coords
   * exist, these arrays are filled and values stored in mesh. 
   * Including these files allows for comparison of HG methods with other
   * methods, along with visualization of results and comparison of 
   * LB_Eval results.
   */
  int    ch_nvtxs = 0;        /* Temporary values for chaco_read_graph.   */
  int    ch_vwgt_dim = 0;     /* Their values are ignored, as vertex      */
  float *ch_vwgts = NULL;     /* info is provided by hypergraph file.     */
  int   *ch_start = NULL, *ch_adj = NULL, ch_ewgt_dim = 0;
  short *ch_assignments = NULL;
  float *ch_ewgts = NULL;
  int    ch_ndim = 0;
  float *ch_x = NULL, *ch_y = NULL, *ch_z = NULL;
  int    ch_no_geom = TRUE;   /* Assume no geometry info is given; reset if
                                 it is provided. */

/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  if (Proc == 0) {

    /* Open and read the hypergraph file. */
    sprintf(filename, "%s.hg", pio_info->pexo_fname);
    fp = fopen(filename, "r");
    if (fp == NULL) {
      sprintf(cmesg, "fatal:  Could not open hypergraph file %s", filename);
      Gen_Error(0, cmesg);
      return 0;
    }

    /* read the array in on processor 0 */
    if (Zoltan_HG_Readfile(Proc, fp, &nvtxs, &nhedges, &npins,
                           &hindex, &hvertex,
                           &vwgt_dim, &vwgts, &hewgt_dim, &hewgts, &base) != 0){
      Gen_Error(0, "fatal: Error returned from Zoltan_HG_Readfile");
      return 0;
    }

    fclose(fp);

    /* If CHACO graph file is available, read it. */
    sprintf(filename, "%s.graph", pio_info->pexo_fname);
    fp = fopen(filename, "r");
    if (fp != NULL) {
      /* CHACO graph file is available. */
      /* Assuming hypergraph vertices are same as chaco vertices. */
      /* Chaco vertices and their weights are ignored in rest of function. */
      if (chaco_input_graph(fp, filename, &ch_start, &ch_adj, &ch_nvtxs,
                      &ch_vwgt_dim, &ch_vwgts, &ch_ewgt_dim, &ch_ewgts) != 0) {
        Gen_Error(0, "fatal: Error returned from chaco_input_graph");
        return 0;
      }
      
      fclose(fp);
    }


    /* If coordinate file is available, read it. */
    sprintf(filename, "%s.coords", pio_info->pexo_fname);
    fp = fopen(filename, "r");
    if (fp != NULL) {
      /* CHACO coordinates file is available. */
      ch_no_geom = FALSE;
      if (chaco_input_geom(fp, filename, ch_nvtxs, &ch_ndim, 
                           &ch_x, &ch_y, &ch_z) != 0) {
        Gen_Error(0, "fatal: Error returned from chaco_input_geom");
        return 0;
      }

      fclose(fp);
    }

  }
  
  MPI_Bcast(&base, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Distribute hypergraph graph */
  /* Use hypergraph vertex information and chaco edge information. */
  if (!chaco_dist_graph(MPI_COMM_WORLD, pio_info, 0, &gnvtxs, &nvtxs, 
             &ch_start, &ch_adj, &vwgt_dim, &vwgts, &ch_ewgt_dim, &ch_ewgts,
             &ch_ndim, &ch_x, &ch_y, &ch_z, &ch_assignments) != 0) {
    Gen_Error(0, "fatal: Error returned from chaco_dist_graph");
    return 0;
  }

  if (!dist_hyperedges(MPI_COMM_WORLD, pio_info, 0, base, &gnhedges,
                       &nhedges, &hindex, &hvertex,
                       &hewgt_dim, &hewgts, &ch_assignments)) {
    Gen_Error(0, "fatal: Error returned from dist_hyperedges");
    return 0;
  }
                       

  /* Initialize mesh structure for Hypergraph. */
  mesh->data_type = HYPERGRAPH;
  mesh->num_elems = nvtxs;
  mesh->vwgt_dim = vwgt_dim;
  mesh->ewgt_dim = ch_ewgt_dim;
  mesh->elem_array_len = mesh->num_elems + 5;
  mesh->num_dims = ch_ndim;
  mesh->num_el_blks = 1;

  mesh->gnhedges = gnhedges;
  mesh->nhedges = nhedges;
  mesh->hewgt_dim = hewgt_dim;

  mesh->hindex = hindex;
  mesh->hvertex = hvertex;
  mesh->hewgts = hewgts;

  mesh->eb_etypes = (int *) malloc (5 * mesh->num_el_blks * sizeof(int));
  if (!mesh->eb_etypes) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  mesh->eb_ids = mesh->eb_etypes + mesh->num_el_blks;
  mesh->eb_cnts = mesh->eb_ids + mesh->num_el_blks;
  mesh->eb_nnodes = mesh->eb_cnts + mesh->num_el_blks;
  mesh->eb_nattrs = mesh->eb_nnodes + mesh->num_el_blks;

  mesh->eb_names = (char **) malloc (mesh->num_el_blks * sizeof(char *));
  if (!mesh->eb_names) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  mesh->eb_etypes[0] = -1;
  mesh->eb_ids[0] = 1;
  mesh->eb_cnts[0] = nvtxs;
  mesh->eb_nattrs[0] = 0;
  /*
   * Each element has one set of coordinates (i.e., node) if a coords file
   * was provided; zero otherwise.
   */
  MPI_Bcast( &ch_no_geom, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ch_no_geom)
    mesh->eb_nnodes[0] = 0;
  else
    mesh->eb_nnodes[0] = 1;

  /* allocate space for name */
  mesh->eb_names[0] = (char *) malloc((MAX_STR_LENGTH+1) * sizeof(char));
  if (!mesh->eb_names[0]) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  strcpy(mesh->eb_names[0], "hypergraph");

  /* allocate the element structure array */
  mesh->elements = (ELEM_INFO_PTR) malloc (mesh->elem_array_len 
                                         * sizeof(ELEM_INFO));
  if (!(mesh->elements)) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  /*
   * intialize all of the element structs as unused by
   * setting the globalID to -1
   */
  for (i = 0; i < mesh->elem_array_len; i++) 
    initialize_element(&(mesh->elements[i]));

  /*
   * now fill the element structure array with the
   * information from the Chaco file
   * Use hypergraph vertex information and chaco edge information.
   */
  if (!chaco_fill_elements(Proc, Num_Proc, prob, mesh, gnvtxs, nvtxs,
                     ch_start, ch_adj, vwgt_dim, vwgts, ch_ewgt_dim, ch_ewgts, 
                     ch_ndim, ch_x, ch_y, ch_z, ch_assignments, base)) {
    Gen_Error(0, "fatal: Error returned from chaco_fill_elements");
    return 0;
  }

  ZOLTAN_FREE((void **) &vwgts);
  safe_free((void **) &ch_ewgts);
  safe_free((void **) &ch_vwgts);
  safe_free((void **) &ch_x);
  safe_free((void **) &ch_y);
  safe_free((void **) &ch_z);
  safe_free((void **) &ch_start);
  safe_free((void **) &ch_adj);
  safe_free((void **) &ch_assignments);

#ifdef KDDKDD_DO_NOT_COMMIT
 if (Debug_Driver > 3)
#endif /* KDDKDD_DO_NOT_COMMIT */
    print_distributed_mesh(Proc, Num_Proc, mesh);

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

/*****************************************************************************/

static int dist_hyperedges(
  MPI_Comm comm,		/* MPI Communicator */
  PARIO_INFO_PTR pio_info,      /* Parallel IO info */
  int     host_proc,		/* processor where all the data is initially */
  int     base,                 /* indicates whether input is 0-based or
                                   1-based (i.e., is lowest vertex number 
                                   0 or 1?). */
  int     *gnhedges,		/* global number of hyperedges */
  int     *nhedges,		/* local number of hyperedges */
  int     **hindex,		/* Starting hvertex entry for hyperedges */
  int     **hvertex,		/* Array of vertices in hyperedges */
  int     *hewgt_dim,           /* number of weights per hyperedge */
  float   **hewgts,		/* hyperedge weight list data */
  short   **assignments         /* assignments from Chaco file; may be NULL */
)
{
/*
 * Distribute hyperedges from one processor to all processors.
 * Vertex distribution is assumed already done through chaco_dist_graph.
 * The memory for the hyperedges on the host node is freed
 * and fresh memory is allocated for the distributed hyperedges.
 */

char *yo = "dist_hypergraph";
int nprocs, myproc, i, h, p;
int *old_hindex = NULL, *old_hvertex = NULL; 
int *size = NULL, *num_send = NULL;
int **send = NULL;
int *send_hindex = NULL, *send_hvertex = NULL;
int max_size, max_num_send;
int hcnt[2];
int hecnt, hvcnt;
float *old_hewgts = NULL;
float *send_hewgts = NULL;
MPI_Status status;

  /* Determine number of processors and my rank. */
  MPI_Comm_size (comm, &nprocs );
  MPI_Comm_rank (comm, &myproc );

  DEBUG_TRACE_START(myproc, yo);

  if (nprocs == 1) {
    /* Set values expected to be returned by this function. */
    /* All array pointers are unchanged.                    */
    *gnhedges = *nhedges;
    return 1;
  }

  /* Broadcast to all procs */
  MPI_Bcast( hewgt_dim, 1, MPI_INT, host_proc, comm);
  MPI_Bcast( nhedges, 1, MPI_INT, host_proc, comm);
  *gnhedges = *nhedges;

  /* Initialize */
  if (*gnhedges == 0) {
    *hindex = NULL;
    *hvertex = NULL;
    return 1;
  }
 
  /* Store pointers to original data */
  if (myproc == host_proc) {
    old_hindex  = *hindex;
    old_hvertex = *hvertex;

    /* Allocate space for size and send flags */
    size = (int *) calloc(2 * nprocs, sizeof(int));
    num_send = size + nprocs;
    send = (int **) malloc(nprocs * sizeof(int *));
    for (i = 0; i < nprocs; i++)
      send[i] = (int *) calloc(*gnhedges, sizeof(int));

    /* Determine to which processors hyperedges should be sent */
    for (h = 0; h < *gnhedges; h++) {
      for (i = old_hindex[h]; i < old_hindex[h+1]; i++) {
        p = ch_dist_proc(old_hvertex[i], *assignments, base);
        if (send[p][h] == 0) {
          size[p] += (old_hindex[h+1] - old_hindex[h]);
          num_send[p]++;
          send[p][h] = 1;
        }
      }
    }

    /* Determine size of send buffers and allocate them. */
    max_size = 0;
    max_num_send = 0;
    for (p = 0; p < nprocs; p++) {
      if (size[p] > max_size) max_size = size[p];
      if (num_send[p] > max_num_send) max_num_send = num_send[p];
    }

    send_hindex = (int *) malloc((max_num_send+1) * sizeof(int));
    send_hvertex = (int *) malloc(max_size * sizeof(int));
    send_hewgts = (float *) malloc(max_num_send * (*hewgt_dim) * sizeof(float));

    /* Load and send data */
    for (p = 0; p < nprocs; p++) {

      if (p == myproc) continue;

      hecnt = 0;
      hvcnt = 0;
      send_hindex[0] = 0;
      for (h = 0; h < *gnhedges; h++) {
        if (send[p][h]) {
          send_hindex[hecnt+1] = send_hindex[hecnt] 
                               + (old_hindex[h+1] - old_hindex[h]);
          for (i = 0; i < *hewgt_dim; i++)
            send_hewgts[hecnt*(*hewgt_dim)+i] = old_hewgts[h*(*hewgt_dim)+i];
          hecnt++;
          for (i = old_hindex[h]; i < old_hindex[h+1]; i++) {
            send_hvertex[hvcnt++] = old_hvertex[i];
          }
        }
      }

      hcnt[0] = hecnt;
      hcnt[1] = hvcnt;

      MPI_Send(hcnt, 2, MPI_INT, p, 1, comm);
      MPI_Send(send_hindex, hecnt+1, MPI_INT, p, 2, comm);
      MPI_Send(send_hvertex, hvcnt, MPI_INT, p, 3, comm);
      MPI_Send(send_hewgts, hecnt*(*hewgt_dim), MPI_FLOAT, p, 4, comm);
    }

    safe_free((void **) &send_hindex);
    safe_free((void **) &send_hvertex);
    safe_free((void **) &send_hewgts);

    /* Copy data owned by myproc into new local storage */
    *nhedges = num_send[myproc];
    *hindex = (int *) malloc((*nhedges+1) * sizeof(int));
    *hvertex = (int *) malloc(size[myproc] * sizeof(int));
    *hewgts = (float *) malloc(*nhedges * *hewgt_dim * sizeof(float));

    hecnt = 0;
    hvcnt = 0;
    (*hindex)[0] = 0;

    for (h = 0; h < *gnhedges; h++) {
      if (send[myproc][h]) {
        (*hindex)[hecnt+1] = (*hindex)[hecnt]
                           + (old_hindex[h+1] - old_hindex[h]);
        for (i = 0; i < *hewgt_dim; i++)
          (*hewgts)[hecnt*(*hewgt_dim)+i] = old_hewgts[h*(*hewgt_dim)+i];
        hecnt++;
        for (i = old_hindex[h]; i < old_hindex[h+1]; i++) {
          (*hvertex)[hvcnt++] = old_hvertex[i];
        }
      }
    }

    for (p = 0; p < nprocs; p++) safe_free((void **) &(send[p]));
    safe_free((void **) &send);
    safe_free((void **) &size);
    safe_free((void **) &old_hindex);
    safe_free((void **) &old_hvertex);
    safe_free((void **) &old_hewgts);
  }
  else {  
    /* host_proc != myproc; receive hedge data from host_proc */
    MPI_Recv(hcnt, 2, MPI_INT, host_proc, 1, comm, &status);
    *nhedges = hcnt[0];
    *hindex = (int *) malloc((hcnt[0]+1) * sizeof(int));
    *hvertex = (int *) malloc(hcnt[1] * sizeof(int));
    *hewgts = (float *) malloc(hcnt[0] * *hewgt_dim * sizeof(float));
    MPI_Recv(*hindex, hcnt[0]+1, MPI_INT, host_proc, 2, comm, &status);
    MPI_Recv(*hvertex, hcnt[1], MPI_INT, host_proc, 3, comm, &status);
    MPI_Recv(*hewgts, hcnt[0]* *hewgt_dim, MPI_FLOAT, host_proc, 4, comm, 
             &status);
  }

  DEBUG_TRACE_END(myproc, yo);
  return 1;
}

/*****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
