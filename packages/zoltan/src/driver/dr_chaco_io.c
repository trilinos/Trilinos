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
#include <limits.h>

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
#include "dr_compress_const.h"

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

int read_chaco_file(int Proc,
                    int Num_Proc,
                    PROB_INFO_PTR prob,
                    PARIO_INFO_PTR pio_info,
                    MESH_INFO_PTR mesh)
{
  /* Local declarations. */
  const char  *yo = "read_chaco_mesh";
  char   cmesg[256];
  char   chaco_fname[FILENAME_MAX + 8];

  int    nvtxs,base;
  int    vwgt_dim=0, ewgt_dim=0;
  int    ndim = 0;
  int   *start = NULL, *adj = NULL;
  int    no_geom = FALSE;
  int    gnvtxs;

  float *ewgts = NULL, *vwgts = NULL;
  float *x = NULL, *y = NULL, *z = NULL;

  short *assignments = NULL;

  ZOLTAN_FILE *fp = NULL;
  int file_error;
/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  /* Set Debug_Chaco_Input */
  if (Debug_Driver > 2) Debug_Chaco_Input = 1;

  if (Proc == 0) {

    /* Open and read the Chaco graph file. */
    sprintf(chaco_fname, "%s.graph", pio_info->pexo_fname);

    fp = ZOLTAN_FILE_open(chaco_fname, "r", pio_info->file_comp);
    file_error = (fp == NULL);
  }

  MPI_Bcast(&file_error, 1, MPI_INT, 0, zoltan_get_global_comm());

  if (file_error) {
    sprintf(cmesg, "fatal:  Could not open Chaco graph file %s",
            chaco_fname);
    Gen_Error(0, cmesg);
    return 0;
  }

  if (Proc == 0) {
    /* read the array in on processor 0 */
    if (chaco_input_graph(fp, chaco_fname, &start, &adj, &nvtxs,
                           &vwgt_dim, &vwgts, &ewgt_dim, &ewgts) != 0) {
      Gen_Error(0, "fatal: Error returned from chaco_input_graph");
      return 0;
    }

    base = 1;  // chaco files are one-based

    /* Read Chaco geometry file, if provided. */
    sprintf(chaco_fname, "%s.coords", pio_info->pexo_fname);

    fp = ZOLTAN_FILE_open(chaco_fname, "r", pio_info->file_comp);
    if (fp == NULL) {
      no_geom = TRUE;
      sprintf(cmesg, "warning:  Could not open Chaco geometry file %s; "
              "no geometry data will be read",
              chaco_fname);
      Gen_Error(0, cmesg);
    }
    else {
      /* read the coordinates in on processor 0 */
      if (chaco_input_geom(fp, chaco_fname, nvtxs, &ndim, &x, &y, &z) != 0) {
        Gen_Error(0, "fatal: Error returned from chaco_input_geom");
        return 0;
      }
    }

#ifdef MESS_UP_POINTS
/* Try to create a geometry that isn't nicely symmetric */
{
int i, j;
double min[3], max[3], a[3], b[3];
min[0] = max[0] = x[0];
if (ndim > 1){
  min[1] = max[1] = y[0];
  if (ndim > 2)
    min[2] = max[2] = z[0];
}
for (i=1; i<nvtxs; i++) {
  if (x[i] < min[0]) min[0] = x[i];
  else if (x[i] > max[0]) max[0] = x[i];
  if (ndim > 1){
    if (y[i] < min[1]) min[1] = y[i];
    else if (y[i] > max[1]) max[1] = y[i];
    if (ndim > 2){
      if (z[i] < min[2]) min[2] = z[i];
      else if (z[i] > max[2]) max[2] = z[i];
    }
  }
}
for (i=0; i<ndim; i++)  /* point inside but near edge of geometry */
  a[i] = ((max[i] - min[i]) * .1) + min[i];

for (i=0; i<nvtxs; i++) { /* move 2/3 of points much closer to "a" */
  if (i%3 == 0) continue;
  b[0] = x[i];
  if (ndim > 1){
    b[1] = y[i];
    if (ndim > 2){
      b[2] = z[i];
    }
  }
  for (j=0; j<ndim; j++){
    b[j] = a[j] + ((b[j] - a[j])*.1);
  }
  x[i] = b[0];
  if (ndim > 1){
    y[i] = b[1];
    if (ndim > 2){
      z[i] = b[2];
    }
  }
}
}
#endif

    /* Read Chaco assignment file, if requested */
    if (pio_info->init_dist_type == INITIAL_FILE) {
      sprintf(chaco_fname, "%s.assign", pio_info->pexo_fname);
      if (pio_info->file_comp == GZIP)
	sprintf(chaco_fname, "%s.gz", chaco_fname);

      fp = ZOLTAN_FILE_open(chaco_fname, "r", pio_info->file_comp);
      if (fp == NULL) {
        sprintf(cmesg, "Error:  Could not open Chaco assignment file %s; "
                "initial distribution cannot be read",
                chaco_fname);
        Gen_Error(0, cmesg);
        return 0;
      }
      else {
        /* read the coordinates in on processor 0 */
        assignments = (short *) malloc(nvtxs * sizeof(short));
        if (!assignments) {
          Gen_Error(0, "fatal: insufficient memory");
          return 0;
        }
        if (chaco_input_assign(fp, chaco_fname, nvtxs, assignments) != 0) {
          Gen_Error(0, "fatal: Error returned from chaco_input_assign");
          return 0;
        }
      }
    }
  }

  /* Distribute graph */

  if (!chaco_dist_graph(zoltan_get_global_comm(), pio_info, 0, &gnvtxs, &nvtxs, 
             &start, &adj, &vwgt_dim, &vwgts, &ewgt_dim, &ewgts, 
             &ndim, &x, &y, &z, &assignments)) {
    Gen_Error(0, "fatal: Error returned from chaco_dist_graph");
    return 0;
  }

  MPI_Bcast(&base, 1, MPI_INT, 0, zoltan_get_global_comm());

  if (!chaco_setup_mesh_struct(Proc, Num_Proc, prob, mesh, pio_info, gnvtxs, nvtxs,
                     start, adj, vwgt_dim, vwgts, ewgt_dim, ewgts, 
                     ndim, x, y, z, assignments, base, no_geom)) {
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

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int chaco_setup_mesh_struct(
  int        Proc,
  int        Num_Proc,
  PROB_INFO_PTR prob,            /* problem description */
  MESH_INFO_PTR mesh,            /* mesh information for the problem */
  PARIO_INFO_PTR pio_info,       /* element distribution info*/
  int        gnvtxs,             /* global number of vertices across all procs*/
  int        nvtxs,              /* number of vertices in local graph */
  int       *start,              /* start of edge list for each vertex */
  int       *adj,                /* edge list data */
  int        vwgt_dim,           /* # of weights per vertex */
  float     *vwgts,              /* vertex weight list data */
  int        ewgt_dim,           /* # of weights per edge */
  float     *ewgts,              /* edge weight list data */
  int        ndim,               /* dimension of the geometry */
  float     *x,                  /* x-coordinates of the vertices */
  float     *y,                  /* y-coordinates of the vertices */
  float     *z,                  /* z-coordinates of the vertices */
  short     *assignments,        /* assignments from Chaco file; may be NULL */
  int       base,                /* smallest vertex number to use; 
                                    base == 1 for Chaco; 
                                    may be 0 or 1 for HG files. */
  int       no_geom              /* flag indicating whether coords are avail. */
)
{
const char *yo = "chaco_setup_mesh_struct";
int i;

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

  /*
   * Each element has one set of coordinates (i.e., node) if a coords file
   * was provided; zero otherwise. 
   */
  MPI_Bcast( &no_geom, 1, MPI_INT, 0, zoltan_get_global_comm());
  if (no_geom)
    mesh->eb_nnodes[0] = 0;
  else
    mesh->eb_nnodes[0] = 1;

  /* allocate space for name */
  mesh->eb_names[0] = (char *) malloc((MAX_STR_LENGTH+1) * sizeof(char));
  if (!mesh->eb_names[0]) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  strcpy(mesh->eb_names[0], "chaco");

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
   */
  if (!chaco_fill_elements(Proc, Num_Proc, prob, mesh, pio_info, gnvtxs, nvtxs,
                     start, adj, vwgt_dim, vwgts, ewgt_dim, ewgts, 
                     ndim, x, y, z, assignments, base)) 
  {
    Gen_Error(0, "fatal: Error returned from chaco_fill_elements");
    return 0;
  }

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int chaco_fill_elements(
  int        Proc,
  int        Num_Proc,
  PROB_INFO_PTR prob,            /* problem description */
  MESH_INFO_PTR mesh,            /* mesh information for the problem */
  PARIO_INFO_PTR pio_info,       /* element distribution info*/
  int        gnvtxs,             /* global number of vertices across all procs*/
  int        nvtxs,              /* number of vertices in local graph */
  int       *start,              /* start of edge list for each vertex */
  int       *adj,                /* edge list data */
  int        vwgt_dim,           /* # of weights per vertex */
  float     *vwgts,              /* vertex weight list data */
  int        ewgt_dim,           /* # of weights per edge */
  float     *ewgts,              /* edge weight list data */
  int        ndim,               /* dimension of the geometry */
  float     *x,                  /* x-coordinates of the vertices */
  float     *y,                  /* y-coordinates of the vertices */
  float     *z,                  /* z-coordinates of the vertices */
  short     *assignments,        /* assignments from Chaco file; may be NULL */
  int       base                 /* smallest vertex number to use; */
)
{
  /* Local declarations. */
  int i, j, k, *local_ids = NULL;
  int num_vtx;
  int elem_id;
  int min_vtx, max_vtx; 
  int *vtx_list = NULL;
  const char *yo = "chaco_fill_elements";
/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  chaco_init_local_ids(&local_ids, &vtx_list, &min_vtx, &max_vtx, &num_vtx, 
                       assignments, base);

  for (i = 0; i < num_vtx; i++) {
    mesh->elements[i].globalID = vtx_list[i]+base;      /* GlobalIDs are 1-based
                                                       in Chaco; may be 0-based
                                                       or 1-based in HG files */
    if (vwgts != NULL){
      for (j=0; j<vwgt_dim; j++) {
        mesh->elements[i].cpu_wgt[j] = vwgts[i*vwgt_dim+j];
      }
    }
    else
      mesh->elements[i].cpu_wgt[0] = 1.0;
    mesh->elements[i].elem_blk = 0; /* only one elem block for all vertices */

    if (assignments)
      mesh->elements[i].my_part = (int)assignments[vtx_list[i]];
    else
      mesh->elements[i].my_part = Proc;  /* Init partition is starting proc.*/

    if (mesh->num_dims > 0) {
      /* One set of coords per element. */
      mesh->elements[i].connect = (ZOLTAN_ID_TYPE *) malloc(sizeof(ZOLTAN_ID_TYPE));
      mesh->elements[i].connect[0] = mesh->elements[i].globalID;
      mesh->elements[i].coord = (float **) malloc(sizeof(float *));
      mesh->elements[i].coord[0] = (float *) calloc(mesh->num_dims,
                                                    sizeof(float));  
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

  for (i = 0; i < num_vtx; i++) {
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
#if 0
        elem_id = adj[start[i] + j] - (1-base);  /* Chaco is 1-based;
                                                    HG may be 0 or 1 based. */
#else
        /* when called from hypergraph build, the line above gives the wrong
         * value for the adjacency.  (it is one too low). I'm guessing that
         * this is the fix.  Works for hypergraph & chaco files in testing.
         */
        elem_id = adj[start[i] + j];
#endif

        /* determine which processor the adjacent vertex is on */

        k = ch_dist_proc(elem_id, assignments, base);

        /*
         * if the adjacent element is on this processor
         * then find the local id for that element
         */
        if (k == Proc) 
          mesh->elements[i].adj[j] = (ZOLTAN_ID_TYPE)local_ids[elem_id-base-min_vtx];
        else /* use the global id */
          mesh->elements[i].adj[j] = (ZOLTAN_ID_TYPE)elem_id;

        mesh->elements[i].adj_proc[j] = k;

        if (ewgts != NULL)
          mesh->elements[i].edge_wgt[j] = ewgts[start[i] + j];
      }
    } /* End: "if (mesh->elements[i].nadj > 0)" */
  } /* End: "for (i = 0; i < mesh->num_elems; i++)" */

  safe_free((void **)(void *) &vtx_list);
  safe_free((void **)(void *) &local_ids);

  if (!build_elem_comm_maps(Proc, mesh)) {
    Gen_Error(0, "Fatal: error building initial elem comm maps");
    return 0;
  }

  if (Debug_Driver > 3)
    print_distributed_mesh(Proc, Num_Proc, mesh);

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

/****************************************************************************/
void chaco_init_local_ids(
  int  **local_ids, 
  int **vtx_list, 
  int *min_vtx, 
  int *max_vtx, 
  int *num_vtx, 
  short *assignments,
  int    base
)
{
/* Initialize an array of local IDs for vertices for quick search. */
int i;
int Proc;

  MPI_Comm_rank(zoltan_get_global_comm(), &Proc);

  *num_vtx = ch_dist_max_num_vtx(assignments);
  *vtx_list = (int *) malloc(((int)*num_vtx) * sizeof(int));
  ch_dist_vtx_list(*vtx_list, num_vtx, Proc, assignments);

  if (*num_vtx > 0) {
    *min_vtx = INT_MAX;
    *max_vtx = -1;
    for (i = 0; i < *num_vtx; i++) {
      if ((*vtx_list)[i] > *max_vtx) *max_vtx = (*vtx_list)[i];
      if ((*vtx_list)[i] < *min_vtx) *min_vtx = (*vtx_list)[i];
    }
    *local_ids = (int *) malloc((*max_vtx - *min_vtx + 1) * sizeof(int));
  }

  for (i = 0; i < *num_vtx; i++) 
    (*local_ids)[(*vtx_list)[i] - (*min_vtx)] = i;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
