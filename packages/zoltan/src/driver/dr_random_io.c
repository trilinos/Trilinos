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

/* Don't generate coords and graph file if nvtxs > this number */
#define OUTPUT_FILES_MAX_NVTXS 1000

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
  char  *yo = "create_random_input";

  int    i, nvtxs, gnvtxs;
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
    vwgt_dim = 1;


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
      vwgts = (float *) malloc(nvtxs * sizeof(float));
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
        vwgts[i] = ((float) rand())/RAND_MAX;
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
      if (fpg != NULL) fprintf(fpg, "%d %d 010\n", nvtxs, start[nvtxs]);
      for (i = 0; i < nvtxs; i++) {
        if (fpg != NULL) fprintf(fpg, "%e ", vwgts[i]);
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

  safe_free((void **) &adj);
  safe_free((void **) &vwgts);
  safe_free((void **) &ewgts);
  safe_free((void **) &start);
  safe_free((void **) &x);
  safe_free((void **) &y);
  safe_free((void **) &z);
  safe_free((void **) &assignments);

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
