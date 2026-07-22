// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif
#ifndef M_PI_2
# define M_PI_2         1.57079632679489661923  /* pi/2 */
#endif

#include "dr_const.h"
#include "dr_externs.h"
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

#ifdef M_PIl
#error "HAVE M_PI1"
#endif

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/* Don't generate coords and graph file if nvtxs > this number */
#define OUTPUT_FILES_MAX_NVTXS 1

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
  const char  *yo = "create_random_triangles";

  int    i, j, w, nvtxs, gnvtxs, ntri;
  int    vwgt_dim=0, ewgt_dim=0;
  int    ndim = 0;
  int   *start = NULL, *adj = NULL;
  int    no_geom = FALSE;

  float *ewgts = NULL, *vwgts = NULL;
  float *x = NULL, *y = NULL, *z = NULL;
  float wgt, diff;

  short *assignments = NULL;

  char filename[FILENAME_MAX+9];
  FILE *fpg = NULL, *fpc = NULL;  /* Files to echo random input */

/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  if (Proc == 0) {

    /* Adapted from read_chaco_graph; build same values as read_chaco_graph
     * and then let random input be distributed as a Chaco graph would be. */

    /* read the array in on processor 0 */
    ntri = (int)pio_info->init_size;
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
      x = (float *) malloc(nvtxs * sizeof(float));
      if (ndim > 1) y = (float *) malloc(nvtxs * sizeof(float));
      if (ndim > 2) z = (float *) malloc(nvtxs * sizeof(float));
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
          x[i] = ((float) rand())/(float)RAND_MAX;
          x[i+1] = x[i] - diff;
          x[i+2] = x[i] + diff;
          if (fpc != NULL) fprintf(fpc, "%e\n%e\n%e\n", x[i],x[i+1],x[i+2]);
        }
        break;
      case 2:
        for (i = 0; i < nvtxs; i+=3)  {
          x[i] = ((float) rand())/(float)RAND_MAX;
          y[i] = ((float) rand())/(float)RAND_MAX;
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
          x[i] = ((float) rand())/(float)RAND_MAX;
          y[i] = ((float) rand())/(float)RAND_MAX;
          z[i] = ((float) rand())/(float)RAND_MAX;
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
              wgt = ((float) rand())/(float)RAND_MAX;
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

  if (!chaco_dist_graph(zoltan_get_global_comm(), pio_info, 0, &gnvtxs, &nvtxs, 
             &start, &adj, &vwgt_dim, &vwgts, &ewgt_dim, &ewgts, 
             &ndim, &x, &y, &z, &assignments)) {
    Gen_Error(0, "fatal: Error returned from chaco_dist_graph");
    return 0;
  }

  if (!chaco_setup_mesh_struct(Proc, Num_Proc, prob, mesh, pio_info, gnvtxs, nvtxs,
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

  char filename[FILENAME_MAX+9];
  FILE *fpg = NULL, *fpc = NULL;  /* Files to echo random input */

/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  if (Proc == 0) {

    /* Adapted from read_chaco_graph; build same values as read_chaco_graph
     * and then let random input be distributed as a Chaco graph would be. */

    /* read the array in on processor 0 */
    nvtxs = pio_info->init_size*Num_Proc;
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
      x = (float *) malloc(nvtxs * sizeof(float));
      if (ndim > 1) y = (float *) malloc(nvtxs * sizeof(float));
      if (ndim > 2) z = (float *) malloc(nvtxs * sizeof(float));
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
          x[i] = ((float) rand())/(float)RAND_MAX;
          if (fpc != NULL) fprintf(fpc, "%e\n", x[i]);
        }
        break;
      case 2:
        for (i = 0; i < nvtxs; i++)  {
          x[i] = ((float) rand())/(float)RAND_MAX;
          y[i] = ((float) rand())/(float)RAND_MAX;
          if (fpc != NULL) fprintf(fpc, "%e %e\n", x[i], y[i]);
        }
        break;
      case 3:
        for (i = 0; i < nvtxs; i++)  {
          x[i] = ((float) rand())/(float)RAND_MAX;
          y[i] = ((float) rand())/(float)RAND_MAX;
          z[i] = ((float) rand())/(float)RAND_MAX;
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
              vwgts[i*vwgt_dim+j] = ((float) rand())/(float)RAND_MAX;
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

  if (!chaco_dist_graph(zoltan_get_global_comm(), pio_info, 0, &gnvtxs, &nvtxs, 
             &start, &adj, &vwgt_dim, &vwgts, &ewgt_dim, &ewgts, 
             &ndim, &x, &y, &z, &assignments)) {
    Gen_Error(0, "fatal: Error returned from chaco_dist_graph");
    return 0;
  }

  if (!chaco_setup_mesh_struct(Proc, Num_Proc, prob, mesh, pio_info, gnvtxs, nvtxs,
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
static int setup_mesh_struct(int, int,
  PROB_INFO_PTR, MESH_INFO_PTR, PARIO_INFO_PTR,
  ZOLTAN_ID_TYPE, int, int*, ZOLTAN_ID_TYPE *, int, float *, int,
  float *, int, float *, float *, float *);

static ZOLTAN_ID_TYPE *vertex_global_ids=NULL;

static ZOLTAN_ID_TYPE local_to_global_id_map(int local_id, int rank)
{
  return vertex_global_ids[rank] + local_id;
}

static int global_to_proc_owner_map(ZOLTAN_ID_TYPE global_id, int numProc, int rank)
{
  int i;
  for (i=numProc-1; i >= 0; i--){
    if (global_id >= vertex_global_ids[i]){
      return i;
    }
  }
  return -1;
}

static void initialize_vertex_global_id_info(ZOLTAN_ID_TYPE numMyGIDs, int numProc)
{
  int i;
  vertex_global_ids = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * (numProc + 1));
  vertex_global_ids[0] = 0;

  for (i=1; i <= numProc; i++){
    vertex_global_ids[i] = vertex_global_ids[i-1] + numMyGIDs;
  }
}

int create_a_graph(
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh)
{
  const char  *yo = "create_a_graph";

  /* The graph (and geometry) is created in parallel by each process, as opposed to being 
   * created by process 0 and then dealt out to the other processes.  This allows us to
   * create graphs where the number of vertices is larger than a number which would
   * fit in the memory of one process.
   *
   * Geometrically the graph is a cylinder extending in the z-direction.
   *
   * Each process creates points along a circle in an x-y plane, and knows which process has the 
   * plane adjacent to it and what global ID has been assigned to the neighbors with the same 
   * x and y value as its points.  So adjacency information is easily created.
   */

  ZOLTAN_ID_TYPE    i, nvtxs, gnvtxs, num4;
  ZOLTAN_ID_TYPE    gid;
  long left=0, right=0;
  int    vwgt_dim=0, ewgt_dim=0;
  int    ndim = 0, next;

  int *start;
  ZOLTAN_ID_TYPE *adj = NULL;

  float *ewgts = NULL, *vwgts = NULL;
  float *x = NULL, *y = NULL, *z = NULL;

  double theta, delta, radius, m, length, step; 

/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  gnvtxs = pio_info->init_size;
  ndim = pio_info->init_dim;
  vwgt_dim = pio_info->init_vwgt_dim;
  if (vwgt_dim<1) vwgt_dim=1; /* For now, insist on 1 or more weights. */

  /* for simplicity coerce number of vertices on a process to a multiple of 4 */

  nvtxs = gnvtxs / Num_Proc;

  if (nvtxs > 4){ 
    num4 = nvtxs / 4;
    nvtxs = num4 * 4;
  }
  else{
    num4 = 1;
    nvtxs = 4;
  }

  gnvtxs = (ZOLTAN_ID_TYPE)nvtxs * Num_Proc;

  if (Proc == 0){
    printf("create_a_graph: Graph will have " ZOLTAN_ID_SPEC " vertices, " ZOLTAN_ID_SPEC " on each process\n",
               gnvtxs, nvtxs);
  }

  /* Each process has the same number of vertices.  Let's determine their global IDs */

  initialize_vertex_global_id_info(nvtxs, Num_Proc);

  /* Calculate vertex coordinates and adjacencies */

  x = (float *) malloc(nvtxs * sizeof(float));
  if (ndim > 1) y = (float *) malloc(nvtxs * sizeof(float));
  if (ndim > 2) z = (float *) malloc(nvtxs * sizeof(float));
  vwgts = (float *) malloc(vwgt_dim*nvtxs * sizeof(float));
  if (!x || (ndim > 1 && !y) || (ndim > 2 && !z) || !vwgts) {
    Gen_Error(0, "fatal: Error allocating memory.");
    return 0;
  }

  if (ndim == 1){
    /* a line */

    step = 1.0 / 500.0;
    length = (double)nvtxs * step;
    x[0] = length * (float)Proc;
    
    for (i=1; i < nvtxs; i++){
      x[i] = x[i+1] + step;
    }
  }
  else if (ndim == 2){
    /* a circle */
    radius = (double)nvtxs/500.0;
    theta = (2 * M_PI ) / (double)Num_Proc;
    delta = theta / (double)nvtxs;
    m = (theta * Proc);
    
    for (i=0; i < nvtxs; i++, m += delta){
      x[i] = radius * cos(m);
      y[i] = radius * sin(m);
    } 
  }
  else if (ndim == 3){
    /* a cylinder */

    radius = (double)nvtxs/500.0;
    delta = M_PI_2 / (double)(num4 + 1);
    theta = delta;
    i = 0;

    while (theta < M_PI_2){
      /* points along first quadrant of a circle in the plane z=Proc */
      x[i] = radius * cos(theta);
      y[i] = radius * sin(theta);
      z[i] = (float)Proc;
      theta += delta;
      i++;
    }

    for (i=0; i < num4; i++){
      /* second quadrant */
      x[num4+i] = -x[num4 - i - 1];
      y[num4+i] = y[num4 - i - 1];
      z[num4+i] = (float)Proc;

      /* third quadrant */
      x[2*num4+i] = -x[i];
      y[2*num4+i] = -y[i];
      z[2*num4+i] = (float)Proc;

      /* third quadrant */
      x[3*num4+i] = x[num4 - i - 1];
      y[3*num4+i] = -y[num4 - i - 1];
      z[3*num4+i] = (float)Proc;
    }
  }

  for (i = 0; i < nvtxs; i++)  {
    if (pio_info->init_vwgt_dim == 0) 
      /* Unit weights if no weights were requested. */
      vwgts[i] = 1.0;
    else {
      int jj;
      srand(0);
      for (jj = 0; jj < vwgt_dim; jj++)  {
        /* Only assign one of the weight dimensions a weight>0. */
        /* Modify to get more complicated test cases. */
        if (jj == (int)(i%vwgt_dim))
          vwgts[i*vwgt_dim+jj] = ((float) rand())/(float)RAND_MAX;
        else
          vwgts[i*vwgt_dim+jj] = 0.0;
      }
    }
  }

  start = (int *)malloc(sizeof(int) * (nvtxs + 1));

  if (ndim < 3){

    /* each vertex has one or two neighbors */

    adj = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * 2 * nvtxs);

    start[0] = 0;
    next = 0;

    for (i=0; i < nvtxs; i++){
      gid = local_to_global_id_map(i, Proc);

      if (ndim == 1){
        left = gid - 1;
        right = gid + 1;
      }
      else if (ndim == 2){
        left = (gid == 0) ? (gnvtxs - 1) : (gid - 1);
        right =  (gid == gnvtxs - 1) ? 0 : (gid + 1);
      }

      start[i+1] = start[i];

      if (left >= 0){
        adj[next++] = left;
        start[i+1]++;
      }

      if (right < gnvtxs){
        adj[next++] = left;
        start[i+1]++;
      }
    }
  }
  else{  

    /* each vertex has 2 neighbors on this process, and one or two off process */

    adj = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * 4 * nvtxs);

    start[0] = 0;
    next = 0;

    for (i=0; i < nvtxs; i++){
      gid = local_to_global_id_map(i, Proc);

      left = (i == 0) ? local_to_global_id_map(nvtxs-1, Proc) : local_to_global_id_map(i-1, Proc);
      right = (i == nvtxs-1) ? local_to_global_id_map(0, Proc) : local_to_global_id_map(i+1, Proc);

      start[i+1] = start[i];

      adj[next++] = left;
      start[i+1]++;

      adj[next++] = right;
      start[i+1]++;

      left = gid - nvtxs;
      right = gid + nvtxs;

      if (left >= 0){
        adj[next++] = left;
        start[i+1]++;
      }
      if (right < gnvtxs){
        adj[next++] = right;
        start[i+1]++;
      }
    }
  }

  /* TODO - edge weights? */

  if (!setup_mesh_struct(Proc, Num_Proc, prob, mesh, pio_info, gnvtxs, nvtxs,
                     start, adj, vwgt_dim, vwgts, ewgt_dim, ewgts,
                     ndim, x, y, z)){
    Gen_Error(0, "fatal: Error returned from chaco_setup_mesh_struct");
    return 0;
  }

  safe_free((void **) &adj);
  safe_free((void **) &vwgts);
  safe_free((void **) &ewgts);
  safe_free((void **) &start);
  safe_free((void **) &x);
  safe_free((void **) &y);
  safe_free((void **) &z);

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

/* Building zdrive's mesh (taken mostly from chaco_setup_mesh_struct) */

static int setup_mesh_struct(
  int        Proc,
  int        Num_Proc,
  PROB_INFO_PTR prob,            /* problem description */
  MESH_INFO_PTR mesh,            /* mesh information for the problem */
  PARIO_INFO_PTR pio_info,       /* element distribution info*/
  ZOLTAN_ID_TYPE   gnvtxs,             /* global number of vertices across all procs*/
  int        nvtxs,              /* number of vertices in local graph */
  int       *start,              /* start of edge list for each vertex */
  ZOLTAN_ID_TYPE  *adj,                /* edge list data */
  int        vwgt_dim,           /* # of weights per vertex */
  float     *vwgts,              /* vertex weight list data */
  int        ewgt_dim,           /* # of weights per edge */
  float     *ewgts,              /* edge weight list data */
  int        ndim,               /* dimension of the geometry */
  float     *x,                  /* x-coordinates of the vertices */
  float     *y,                  /* y-coordinates of the vertices */
  float     *z                   /* z-coordinates of the vertices */
)
{
const char *yo = "setup_mesh_struct";
int i, j, k;
ZOLTAN_ID_TYPE elem_id;
ZOLTAN_ID_TYPE min_vtx;

  DEBUG_TRACE_START(Proc, yo);

  /* Initialize mesh structure for Chaco mesh. */
  mesh->data_type = ZOLTAN_GRAPH;
  mesh->vwgt_dim = vwgt_dim;
  mesh->ewgt_dim = ewgt_dim;
  mesh->num_elems = nvtxs;
  mesh->elem_array_len = mesh->num_elems + 5;
  mesh->num_dims = ndim;
  mesh->num_el_blks = 1;

  mesh->eb_etypes = (int *) malloc (4 * mesh->num_el_blks * sizeof(int));
  if (!mesh->eb_etypes) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  mesh->eb_ids = mesh->eb_etypes + mesh->num_el_blks;
  mesh->eb_nnodes = mesh->eb_ids + mesh->num_el_blks;
  mesh->eb_nattrs = mesh->eb_nnodes + mesh->num_el_blks;

  mesh->eb_cnts = (ZOLTAN_ID_TYPE *) malloc (mesh->num_el_blks * sizeof(ZOLTAN_ID_TYPE));
  if (!mesh->eb_cnts) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  mesh->eb_names = (char **) malloc (mesh->num_el_blks * sizeof(char *));
  if (!mesh->eb_names) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  mesh->eb_etypes[0] = -1;
  mesh->eb_ids[0] = 1;
  mesh->eb_cnts[0] = (ZOLTAN_ID_TYPE)nvtxs;
  mesh->eb_nattrs[0] = 0;

  mesh->hindex = (int *) malloc(sizeof(int));
  mesh->hindex[0] = 0;
  mesh->eb_nnodes[0] = 1;

  /* allocate space for name */
  mesh->eb_names[0] = (char *) malloc(16* sizeof(char));
  if (!mesh->eb_names[0]) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  strcpy(mesh->eb_names[0], "random-graph");

  /* allocate the element structure array */
  mesh->elements = (ELEM_INFO_PTR) malloc (mesh->elem_array_len * sizeof(ELEM_INFO));
  if (!(mesh->elements)) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  /* write element data */

  for (i = 0; i < mesh->elem_array_len; i++) 
    initialize_element(&(mesh->elements[i]));

  min_vtx = local_to_global_id_map(0, Proc);

  for (i = 0; i < nvtxs; i++) {
    mesh->elements[i].globalID = local_to_global_id_map(i, Proc);
                                                       
    if (vwgts != NULL){
      for (j=0; j<vwgt_dim; j++) {
        mesh->elements[i].cpu_wgt[j] = vwgts[i*vwgt_dim+j];
      }
    }
    else
      mesh->elements[i].cpu_wgt[0] = 1.0;

    mesh->elements[i].elem_blk = 0;
    mesh->elements[i].my_part = Proc;

    if (mesh->num_dims > 0) {
      /* One set of coords per element. */
      mesh->elements[i].connect = (ZOLTAN_ID_TYPE *) malloc(sizeof(ZOLTAN_ID_TYPE));
      mesh->elements[i].connect[0] = mesh->elements[i].globalID;
      mesh->elements[i].coord = (float **) malloc(sizeof(float *));
      mesh->elements[i].coord[0] = (float *) calloc(mesh->num_dims, sizeof(float));  
      mesh->elements[i].coord[0][0] = x[i];
      mesh->elements[i].avg_coord[0] = x[i];
      if (mesh->num_dims > 1) {
        mesh->elements[i].coord[0][1] = y[i];
        mesh->elements[i].avg_coord[1] = y[i];
        if (mesh->num_dims > 2) {
          mesh->elements[i].coord[0][2] = z[i];
          mesh->elements[i].avg_coord[2] = z[i];
        }
      }
    }
  }

  for (i = 0; i < nvtxs; i++) {
    /* now start with the adjacencies */
    if (start != NULL)
      mesh->elements[i].nadj = start[i+1] - start[i];
    else
      mesh->elements[i].nadj = 0;
    if (mesh->elements[i].nadj > 0) {
      mesh->elements[i].adj_len = mesh->elements[i].nadj;
      mesh->elements[i].adj = (ZOLTAN_ID_TYPE *) malloc (mesh->elements[i].nadj * sizeof(ZOLTAN_ID_TYPE));
      mesh->elements[i].adj_proc = (int *) malloc (mesh->elements[i].nadj * sizeof(int));
      if (!(mesh->elements[i].adj) || !(mesh->elements[i].adj_proc)) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
      if (ewgts != NULL) {
        mesh->elements[i].edge_wgt = (float *) malloc (mesh->elements[i].nadj * sizeof(float));
        if (!(mesh->elements[i].edge_wgt)) {
          Gen_Error(0, "fatal: insufficient memory");
          return 0;
        }
      }
      else
        mesh->elements[i].edge_wgt = NULL;

      for (j = 0; j < mesh->elements[i].nadj; j++) {
        elem_id = adj[start[i] + j];

        k = global_to_proc_owner_map(elem_id, Num_Proc, Proc);

        /*
         * if the adjacent element is on this processor
         * then find the local id for that element
         */
        if (k == Proc) 
          mesh->elements[i].adj[j] = elem_id-min_vtx;
        else /* use the global id */
          mesh->elements[i].adj[j] = elem_id;

        mesh->elements[i].adj_proc[j] = k;

        if (ewgts != NULL)
          mesh->elements[i].edge_wgt[j] = ewgts[start[i] + j];
      }
    } /* End: "if (mesh->elements[i].nadj > 0)" */
  } /* End: "for (i = 0; i < mesh->num_elems; i++)" */

  if (!build_elem_comm_maps(Proc, mesh)) {
    Gen_Error(0, "Fatal: error building initial elem comm maps");
    return 0;
  }

  if (Debug_Driver > 3)
    print_distributed_mesh(Proc, Num_Proc, mesh);

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
