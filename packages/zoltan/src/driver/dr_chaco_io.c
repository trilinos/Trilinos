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

#ifndef MAX_STR_LENGTH
#define MAX_STR_LENGTH 80
#endif

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int read_chaco_mesh(int Proc,
                    int Num_Proc,
                    PROB_INFO_PTR prob,
                    PARIO_INFO_PTR pio_info,
                    MESH_INFO_PTR mesh)
{
  /* Local declarations. */
  char  *yo = "read_chaco_mesh";
  char   cmesg[256];
  char   chaco_fname[FILENAME_MAX + 8];

  int    i, nvtxs, gnvtxs;
  int    vwgt_dim=0, ewgt_dim=0;
  int    ndim = 0;
  int   *start = NULL, *adj = NULL;
  int    no_geom = FALSE;

  float *ewgts = NULL, *vwgts = NULL;
  float *x = NULL, *y = NULL, *z = NULL;

  short *assignments = NULL;

  FILE  *fp;
/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  /* Set Debug_Chaco_Input */
  if (Debug_Driver > 2) Debug_Chaco_Input = 1;

  if (Proc == 0) {

    /* Open and read the Chaco graph file. */
    sprintf(chaco_fname, "%s.graph", pio_info->pexo_fname);   
    fp = fopen(chaco_fname, "r");
    if (fp == NULL) {
      sprintf(cmesg, "fatal:  Could not open Chaco graph file %s",
              chaco_fname);
      Gen_Error(0, cmesg);
      return 0;
    }

    /* read the array in on processor 0 */
    if (chaco_input_graph(fp, chaco_fname, &start, &adj, &nvtxs,
                           &vwgt_dim, &vwgts, &ewgt_dim, &ewgts) != 0) {
      Gen_Error(0, "fatal: Error returned from chaco_input_graph");
      return 0;
    }

    /* Read Chaco geometry file, if provided. */
    sprintf(chaco_fname, "%s.coords", pio_info->pexo_fname);
    fp = fopen(chaco_fname, "r");
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

    /* Read Chaco assignment file, if requested */
    if (pio_info->init_dist_type == INITIAL_FILE) {
      sprintf(chaco_fname, "%s.assign", pio_info->pexo_fname);
      fp = fopen(chaco_fname, "r");
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
  if (!chaco_dist_graph(MPI_COMM_WORLD, pio_info, 0, &gnvtxs, &nvtxs, 
             &start, &adj, &vwgt_dim, &vwgts, &ewgt_dim, &ewgts, 
             &ndim, &x, &y, &z, &assignments) != 0) {
      Gen_Error(0, "fatal: Error returned from chaco_dist_graph");
      return 0;
  }

  /* Initialize mesh structure for Chaco mesh. */
  mesh->data_type = GRAPH;
  mesh->vwgt_dim = vwgt_dim;
  mesh->ewgt_dim = ewgt_dim;
  mesh->num_elems = nvtxs;
  mesh->elem_array_len = mesh->num_elems + 5;
  mesh->num_dims = ndim;
  mesh->num_el_blks = 1;

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
  MPI_Bcast( &no_geom, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
  if (!chaco_fill_elements(Proc, Num_Proc, prob, mesh, gnvtxs, nvtxs,
                     start, adj, vwgt_dim, vwgts, ewgt_dim, ewgts, 
                     ndim, x, y, z, assignments, 1)) {
    Gen_Error(0, "fatal: Error returned from chaco_fill_elements");
    return 0;
  }

  if (adj != NULL) free(adj);
  if (vwgts != NULL) free(vwgts);
  if (ewgts != NULL) free(ewgts);
  if (start != NULL) free(start);
  if (x != NULL) free(x);
  if (y != NULL) free(y);
  if (z != NULL) free(z);
  if (assignments != NULL) free(assignments);

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
  int       base                 /* smallest vertex number to use; 
                                    base == 1 for Chaco; 
                                    may be 0 or 1 for HG files. */
)
{
  /* Local declarations. */
  int i, j, k, elem_id, local_id;
  int num_vtx; 
  int *vtx_list = NULL;
  char *yo = "chaco_fill_elements";
/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  num_vtx = ch_dist_max_num_vtx(assignments);
  vtx_list = (int *) malloc(num_vtx * sizeof(int));
  ch_dist_vtx_list(vtx_list, &num_vtx, Proc, assignments);

  for (i = 0; i < num_vtx; i++) {
    mesh->elements[i].globalID = vtx_list[i]+base;  /* GlobalIDs are 1-based
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
    mesh->elements[i].my_part = Proc;  /* Initial partition is starting proc.*/
    if (mesh->num_dims > 0) {
      /* One set of coords per element. */
      mesh->elements[i].connect = (int *) malloc(sizeof(int));
      mesh->elements[i].connect[0] = mesh->elements[i].globalID;
      mesh->elements[i].coord = (float **) malloc(sizeof(float *));
      mesh->elements[i].coord[0] = (float *) calloc(mesh->num_dims,
                                                    sizeof(float));  
      mesh->elements[i].coord[0][0] = x[i];
      if (mesh->num_dims > 1) {
        mesh->elements[i].coord[0][1] = y[i];
        if (mesh->num_dims > 2) {
          mesh->elements[i].coord[0][2] = z[i];
        }
      }
    }

    /* now start with the adjacencies */
    if (start != NULL)
      mesh->elements[i].nadj = start[i+1] - start[i];
    else
      mesh->elements[i].nadj = 0;
    if (mesh->elements[i].nadj > 0) {
      mesh->elements[i].adj_len = mesh->elements[i].nadj;
      mesh->elements[i].adj = (int *) malloc (mesh->elements[i].nadj * sizeof(int));
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

        /* determine which processor the adjacent vertex is on */
        k = ch_dist_proc(elem_id, assignments, base);

        /*
         * if the adjacent element is on this processor
         * then find the local id for that element
         */
        if (k == Proc) {
          local_id = in_list((elem_id-1), num_vtx, vtx_list);
          mesh->elements[i].adj[j] = local_id;
        }
        else /* use the global id */
          mesh->elements[i].adj[j] = elem_id;

        mesh->elements[i].adj_proc[j] = k;

        if (ewgts != NULL)
          mesh->elements[i].edge_wgt[j] = ewgts[start[i] + j];
      }
    } /* End: "if (mesh->elements[i].nadj > 0)" */
  } /* End: "for (i = 0; i < mesh->num_elems; i++)" */

  safe_free((void **) &vtx_list);

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
