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

#ifndef MAX_STR_LENGTH
#define MAX_STR_LENGTH 80
#endif

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
  int    nvtxs = 0, nhedges = 0, npins = 0;
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


    /* If CHACO coordinate file is available, read it. */
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

  /* Assuming serial implementation only. */
  /* KDDKDD Need to handle distribution of hyperedges here. */
  mesh->nhedges = nhedges;
  mesh->hewgt_dim = hewgt_dim;

  mesh->hindex = (int *) malloc((nhedges+1) * sizeof(int));
  memcpy(mesh->hindex, hindex, (nhedges+1) * sizeof(int));

  if (npins > 0) {
    mesh->hvertex = (int *) malloc(npins * sizeof(int));
    memcpy(mesh->hvertex, hvertex, npins * sizeof(int));
  }
  if (nhedges > 0 && hewgt_dim > 0) {
    mesh->hewgts = (float *) malloc(hewgt_dim * nhedges * sizeof(float));
    memcpy(mesh->hewgts, hewgts, hewgt_dim * nhedges * sizeof(float));
  }
  /* KDDKDD Done distribution of hyperedges for serial case. */

  /* Initialize mesh structure for Hypergraph. */
  mesh->data_type = HYPERGRAPH;
  mesh->num_elems = nvtxs;
  mesh->hewgt_dim = hewgt_dim;
  mesh->vwgt_dim = vwgt_dim;
  mesh->ewgt_dim = ch_ewgt_dim;
  mesh->elem_array_len = mesh->num_elems + 5;
  mesh->num_dims = ch_ndim;
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

  ZOLTAN_FREE((void **) &hindex);
  ZOLTAN_FREE((void **) &hvertex);
  ZOLTAN_FREE((void **) &vwgts);
  ZOLTAN_FREE((void **) &hewgts);
  safe_free((void **) &ch_ewgts);
  safe_free((void **) &ch_vwgts);
  safe_free((void **) &ch_x);
  safe_free((void **) &ch_y);
  safe_free((void **) &ch_z);
  safe_free((void **) &ch_start);
  safe_free((void **) &ch_adj);
  safe_free((void **) &ch_assignments);

 if (Debug_Driver > 3)
    print_distributed_mesh(Proc, Num_Proc, mesh);

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
