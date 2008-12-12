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

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#ifndef MAX_STR_LENGTH
#define MAX_STR_LENGTH 80
#endif

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/* Don't generate coords and graph file if nvtxs > this number */
#define OUTPUT_FILES_MAX_NVTXS 1000

/* Use "size" as number of triangles - generate three points for each,
 * this gives us some adjacencies in the .graph file - also can
 * visualize this with vtk_view, which can't be used to visualize
 * points alone.
 */

int create_random_triangles(
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh)
{
  /* Local declarations. */
  const char  *yo = "create_random_input";

  int    i, j, w, nvtxs, gnvtxs, ntri;
  int    vwgt_dim=0, ewgt_dim=0;
  int    ndim = 0;
  int   *start = NULL, *adj = NULL;
  int    no_geom = FALSE;

  float *ewgts = NULL, *vwgts = NULL;
  float *x = NULL, *y = NULL, *z = NULL;
  float wgt, diff;

  short *assignments = NULL;

  char filename[256];
  FILE *fpg = NULL, *fpc = NULL;  /* Files to echo random input */

/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  if (Proc == 0) {

    /* Adapted from read_chaco_graph; build same values as read_chaco_graph
     * and then let random input be distributed as a Chaco graph would be. */

    /* read the array in on processor 0 */
    ntri = pio_info->init_size;
    nvtxs = ntri * 3;
    ndim = pio_info->init_dim;
    vwgt_dim = pio_info->init_vwgt_dim;
    if (vwgt_dim<1) vwgt_dim=1; /* For now, insist on 1 or more weights. */

    if (ntri <= OUTPUT_FILES_MAX_NVTXS) {
      sprintf(filename, "%s.graph", pio_info->pexo_fname);
      fpg = fopen(filename, "w");
      sprintf(filename, "%s.coords", pio_info->pexo_fname);
      fpc = fopen(filename, "w");
    }

    if (nvtxs > 0) {
      /* Allocate space for vertex weights and coordinates. */
      x = (float *) malloc(nvtxs * sizeof(double));
      if (ndim > 1) y = (float *) malloc(nvtxs * sizeof(double));
      if (ndim > 2) z = (float *) malloc(nvtxs * sizeof(double));
      vwgts = (float *) malloc(vwgt_dim*nvtxs * sizeof(float));
      if (!x || (ndim > 1 && !y) || (ndim > 2 && !z) || !vwgts) {
        Gen_Error(0, "fatal: Error allocating memory.");
        return 0;
      }

      /* Generate random triangles and vertex weights. */
      srand(0);
      diff = 1.0/50.0; 

      switch (ndim) {
      case 1:
        for (i = 0; i < nvtxs; i+=3)  {
          x[i] = ((float) rand())/RAND_MAX;
          x[i+1] = x[i] - diff;
          x[i+2] = x[i] + diff;
          if (fpc != NULL) fprintf(fpc, "%e\n%e\n%e\n", x[i],x[i+1],x[i+2]);
        }
        break;
      case 2:
        for (i = 0; i < nvtxs; i+=3)  {
          x[i] = ((float) rand())/RAND_MAX;
          y[i] = ((float) rand())/RAND_MAX;
          x[i+1] = x[i] - diff;
          y[i+1] = y[i];
          x[i+2] = x[i];
          y[i+2] = y[i] + diff;
          if (fpc != NULL) fprintf(fpc, "%e %e\n%e %e\n%e %e\n", 
                           x[i], y[i],x[i+1], y[i+1],x[i+2], y[i+2]);
        }
        break;
      case 3:
        for (i = 0; i < nvtxs; i+=3)  {
          x[i] = ((float) rand())/RAND_MAX;
          y[i] = ((float) rand())/RAND_MAX;
          z[i] = ((float) rand())/RAND_MAX;
          x[i+1] = x[i] - diff;
          y[i+1] = y[i];
          z[i+1] = z[i];
          x[i+2] = x[i];
          y[i+2] = y[i] + diff;
          z[i+2] = z[i];
          if (fpc != NULL) fprintf(fpc, "%e %e %e\n%e %e %e\n%e %e %e\n", 
                   x[i], y[i], z[i],x[i+1], y[i+1], z[i+1],x[i+2], y[i+2], z[i+2]);
        }
        break;
      }

      if (pio_info->init_vwgt_dim == 0) {
        for (i = 0; i < nvtxs; i++)  {
          /* Unit weights if no weights were requested. */
          vwgts[i] = 1.0;
        }
      }
      else{
        memset(vwgts, 0, nvtxs * sizeof(float));
        for (i=0,w=0; i < ntri; i++)  {
          for (j = 0; j < vwgt_dim; j++)  {
            /* Each point of triangle gets the same weight. */
            /* Only assign one of the weight dimensions a weight>0. */
            /* Modify to get more complicated test cases. */
            if (j == i%vwgt_dim){
              wgt = ((float) rand())/RAND_MAX;
              vwgts[w+j] = wgt;
              w += vwgt_dim;
              vwgts[w+j] = wgt;
              w += vwgt_dim;
              vwgts[w+j] = wgt;
              w += vwgt_dim;
            }
          }
        }
      }

      /* Each vertex has two neighbors */

      ewgt_dim = 0;
      start = (int *) malloc((nvtxs+1) * sizeof(int));
      adj = (int *)malloc(nvtxs*2*sizeof(int));
      ewgts = NULL;
      for (i = 0; i <= nvtxs; i++) {
        start[i] = i*2;
      }

      for (i=0, w=0, j=1; i < ntri; i++,j+=3){
        /* j is first vertex (1-based) in triangle */
        adj[w++] = j+1;   /* adjacencies of vertex j */
        adj[w++] = j+2;
        adj[w++] = j+2;   /* adjacencies of vertex j+1 */
        adj[w++] = j;
        adj[w++] = j;     /* adjacencies of vertex j+2 */
        adj[w++] = j+1;
      }
      if (fpg != NULL) {
        if (vwgt_dim==1)
          fprintf(fpg, "%d %d 010\n", nvtxs, nvtxs);
        else
          fprintf(fpg, "%d %d 010 %d\n", nvtxs, nvtxs, vwgt_dim);
      
        for (i = 0, w=0; i < nvtxs; i++, w += 2) {
          for (j = 0; j < vwgt_dim; j++) 
            fprintf(fpg, "%e ", vwgts[i*vwgt_dim+j]);
          fprintf(fpg, "%d %d",adj[w],adj[w+1]);
          fprintf(fpg, "\n");
        }
      }
    }
    if (fpg != NULL) fclose(fpg);
    if (fpc != NULL) fclose(fpc);
  }

  /* Distribute graph */
  if (!chaco_dist_graph(MPI_COMM_WORLD, pio_info, 0, &gnvtxs, &nvtxs, 
             &start, &adj, &vwgt_dim, &vwgts, &ewgt_dim, &ewgts, 
             &ndim, &x, &y, &z, &assignments)) {
    Gen_Error(0, "fatal: Error returned from chaco_dist_graph");
    return 0;
  }

  if (!chaco_setup_mesh_struct(Proc, Num_Proc, prob, mesh, gnvtxs, nvtxs,
                     start, adj, vwgt_dim, vwgts, ewgt_dim, ewgts,
                     ndim, x, y, z, assignments, 1, no_geom)) {
    Gen_Error(0, "fatal: Error returned from chaco_setup_mesh_struct");
    return 0;
  }

  safe_free((void **)(void *) &adj);
  safe_free((void **)(void *) &vwgts);
  safe_free((void **)(void *) &ewgts);
  safe_free((void **)(void *) &start);
  safe_free((void **)(void *) &x);
  safe_free((void **)(void *) &y);
  safe_free((void **)(void *) &z);
  safe_free((void **)(void *) &assignments);

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int create_random_input(
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh)
{
  /* Local declarations. */
  const char  *yo = "create_random_input";

  int    i, j, nvtxs, gnvtxs;
  int    vwgt_dim=0, ewgt_dim=0;
  int    ndim = 0;
  int   *start = NULL, *adj = NULL;
  int    no_geom = FALSE;

  float *ewgts = NULL, *vwgts = NULL;
  float *x = NULL, *y = NULL, *z = NULL;

  short *assignments = NULL;

  char filename[256];
  FILE *fpg = NULL, *fpc = NULL;  /* Files to echo random input */

/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  if (Proc == 0) {

    /* Adapted from read_chaco_graph; build same values as read_chaco_graph
     * and then let random input be distributed as a Chaco graph would be. */

    /* read the array in on processor 0 */
    nvtxs = pio_info->init_size;
    ndim = pio_info->init_dim;
    vwgt_dim = pio_info->init_vwgt_dim;
    if (vwgt_dim<1) vwgt_dim=1; /* For now, insist on 1 or more weights. */

    if (nvtxs <= OUTPUT_FILES_MAX_NVTXS) {
      sprintf(filename, "%s.graph", pio_info->pexo_fname);
      fpg = fopen(filename, "w");
      sprintf(filename, "%s.coords", pio_info->pexo_fname);
      fpc = fopen(filename, "w");
    }

    if (nvtxs > 0) {
      /* Allocate space for vertex weights and coordinates. */
      x = (float *) malloc(nvtxs * sizeof(double));
      if (ndim > 1) y = (float *) malloc(nvtxs * sizeof(double));
      if (ndim > 2) z = (float *) malloc(nvtxs * sizeof(double));
      vwgts = (float *) malloc(vwgt_dim*nvtxs * sizeof(float));
      if (!x || (ndim > 1 && !y) || (ndim > 2 && !z) || !vwgts) {
        Gen_Error(0, "fatal: Error allocating memory.");
        return 0;
      }

      /* Generate random coordinates and weights. */
      srand(0);
      switch (ndim) {
      case 1:
        for (i = 0; i < nvtxs; i++)  {
          x[i] = ((float) rand())/RAND_MAX;
          if (fpc != NULL) fprintf(fpc, "%e\n", x[i]);
        }
        break;
      case 2:
        for (i = 0; i < nvtxs; i++)  {
          x[i] = ((float) rand())/RAND_MAX;
          y[i] = ((float) rand())/RAND_MAX;
          if (fpc != NULL) fprintf(fpc, "%e %e\n", x[i], y[i]);
        }
        break;
      case 3:
        for (i = 0; i < nvtxs; i++)  {
          x[i] = ((float) rand())/RAND_MAX;
          y[i] = ((float) rand())/RAND_MAX;
          z[i] = ((float) rand())/RAND_MAX;
          if (fpc != NULL) fprintf(fpc, "%e %e %e\n", x[i], y[i], z[i]);
        }
        break;
      }

      for (i = 0; i < nvtxs; i++)  {
        if (pio_info->init_vwgt_dim == 0) 
          /* Unit weights if no weights were requested. */
          vwgts[i] = 1.0;
        else
          for (j = 0; j < vwgt_dim; j++)  {
            /* Only assign one of the weight dimensions a weight>0. */
            /* Modify to get more complicated test cases. */
            if (j == i%vwgt_dim)
              vwgts[i*vwgt_dim+j] = ((float) rand())/RAND_MAX;
            else
              vwgts[i*vwgt_dim+j] = 0.0;
          }
      }

      /* KDDSJP:  Need to add edge info here.  Later.  For now, set no edges */
      ewgt_dim = 0;
      start = (int *) malloc((nvtxs+1) * sizeof(int));
      adj = NULL;
      ewgts = NULL;
      for (i = 0; i < nvtxs; i++) {
        start[i] = 0;
      }
      start[nvtxs] = 0;
      if (fpg != NULL) {
        if (vwgt_dim==1)
          fprintf(fpg, "%d %d 010\n", nvtxs, start[nvtxs]);
        else
          fprintf(fpg, "%d %d 010 %d\n", nvtxs, start[nvtxs], vwgt_dim);
      }
      for (i = 0; i < nvtxs; i++) {
        if (fpg != NULL) 
          for (j = 0; j < vwgt_dim; j++) 
            fprintf(fpg, "%e ", vwgts[i*vwgt_dim+j]);
        /* KDDSJP Print edges here */
        if (fpg != NULL) fprintf(fpg, "\n");
      }
    }
    if (fpg != NULL) fclose(fpg);
    if (fpc != NULL) fclose(fpc);
  }

  /* Distribute graph */
  if (!chaco_dist_graph(MPI_COMM_WORLD, pio_info, 0, &gnvtxs, &nvtxs, 
             &start, &adj, &vwgt_dim, &vwgts, &ewgt_dim, &ewgts, 
             &ndim, &x, &y, &z, &assignments)) {
    Gen_Error(0, "fatal: Error returned from chaco_dist_graph");
    return 0;
  }

  if (!chaco_setup_mesh_struct(Proc, Num_Proc, prob, mesh, gnvtxs, nvtxs,
                     start, adj, vwgt_dim, vwgts, ewgt_dim, ewgts,
                     ndim, x, y, z, assignments, 1, no_geom)) {
    Gen_Error(0, "fatal: Error returned from chaco_setup_mesh_struct");
    return 0;
  }

  safe_free((void **)(void *) &adj);
  safe_free((void **)(void *) &vwgts);
  safe_free((void **)(void *) &ewgts);
  safe_free((void **)(void *) &start);
  safe_free((void **)(void *) &x);
  safe_free((void **)(void *) &y);
  safe_free((void **)(void *) &z);
  safe_free((void **)(void *) &assignments);

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
