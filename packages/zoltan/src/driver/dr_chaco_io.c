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
static char *cvs_dr_chaco_io = "$Id$";
#endif

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*--------------------------------------------------------------------------*/
/* Author(s):  Matthew M. St.John (9226)                                    */
/*--------------------------------------------------------------------------*/
/* Revision History:                                                        */
/*                                                                          */
/*    24 May 1999:      Date of creation                                    */
/*--------------------------------------------------------------------------*/


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
#include "ch_input_const.h"
#include "ch_input.h"

static int fill_elements(int, int, PROB_INFO_PTR, ELEM_INFO *, int, int *, 
                         int *, int *, int *, float *, int, 
                         float *, float *, float *);

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int read_chaco_mesh(int Proc,
                    int Num_Proc,
                    PROB_INFO_PTR prob,
                    PARIO_INFO_PTR pio_info,
                    ELEM_INFO **elements)
{
  /* Local declarations. */
  char   cmesg[256];
  char   chaco_fname[FILENAME_MAX + 8];

  int    i, nvtxs;
  int    ndim = 0;
  int   *start = NULL, *adj = NULL, *vwgts = NULL, *vtxdist = NULL;

  float *ewgts = NULL;
  float *x = NULL, *y = NULL, *z = NULL;

  FILE  *fp;
/***************************** BEGIN EXECUTION ******************************/

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
                           &vwgts, &ewgts) != 0) {
      Gen_Error(0, "fatal: Error returned from chaco_input_graph");
      return 0;
    }

    fclose (fp);

    /* If method requires coordinate info, read Chaco geometry file. */
    if (prob->read_coord) {
      sprintf(chaco_fname, "%s.coords", pio_info->pexo_fname);
      fp = fopen(chaco_fname, "r");
      if (fp == NULL) {
        sprintf(cmesg, "fatal:  Could not open Chaco geometry file %s",
                chaco_fname);
        Gen_Error(0, cmesg);
        return 0;
      }

      /* read the coordinates in on processor 0 */
      if (chaco_input_geom(fp, chaco_fname, nvtxs, &ndim, &x, &y, &z) != 0) {
        Gen_Error(0, "fatal: Error returned from chaco_input_geom");
        return 0;
      }
    }
  }

  /* Distribute graph */
  if (!chaco_dist_graph(MPI_COMM_WORLD, 0, &nvtxs, &vtxdist, &start, &adj, 
                        &vwgts, &ewgts, &ndim, &x, &y, &z) != 0) {
      Gen_Error(0, "fatal: Error returned from chaco_dist_graph");
      return 0;
  }

  /* Initialize Mesh structure for Chaco mesh. */
  Mesh.num_elems = nvtxs;
  Mesh.elem_array_len = Mesh.num_elems + 5;
  Mesh.num_dims = ndim;
  Mesh.num_el_blks = 1;

  Mesh.eb_ids = (int *) malloc (4 * Mesh.num_el_blks * sizeof(int));
  if (!Mesh.eb_ids) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  Mesh.eb_cnts = Mesh.eb_ids + Mesh.num_el_blks;
  Mesh.eb_nnodes = Mesh.eb_cnts + Mesh.num_el_blks;
  Mesh.eb_nattrs = Mesh.eb_nnodes + Mesh.num_el_blks;

  Mesh.eb_names = (char **) malloc (Mesh.num_el_blks * sizeof(char *));
  if (!Mesh.eb_names) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  /* Each element has only one set of coordinates (i.e., node). */
  Mesh.eb_ids[0] = 1;
  Mesh.eb_cnts[0] = nvtxs;
  Mesh.eb_nnodes[0] = 1;
  Mesh.eb_nattrs[0] = 0;
  Mesh.eb_names[0] = "chaco";

  /* allocate the element structure array */
  *elements = (ELEM_INFO_PTR) malloc (Mesh.elem_array_len * sizeof(ELEM_INFO));
  if (!(*elements)) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  /*
   * intialize all of the element structs as unused by
   * setting the globalID to -1
   */
  for (i = 0; i < Mesh.elem_array_len; i++)
    (*elements)[i].globalID = -1;

  /*
   * now fill the element structure array with the
   * information from the Chaco file
   */
  if (!fill_elements(Proc, Num_Proc, prob, *elements, nvtxs, vtxdist, 
                     start, adj, vwgts, ewgts, ndim, x, y, z)) {
    Gen_Error(0, "fatal: Error returned from fill_elements");
    return 0;
  }

  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int fill_elements(
  int        Proc,
  int        Num_Proc,
  PROB_INFO_PTR prob,            /* problem description */
  ELEM_INFO  elem[],             /* array of element information */
  int        nvtxs,              /* number of vertices in graph */
  int       *vtxdist,            /* vertex distribution data */
  int       *start,              /* start of edge list for each vertex */
  int       *adj,                /* edge list data */
  int       *vwgts,              /* vertex weight list data */
  float     *ewgts,              /* edge weight list data */
  int        ndim,               /* dimension of the geometry */
  float     *x,                  /* x-coordinates of the vertices */
  float     *y,                  /* y-coordinates of the vertices */
  float     *z                   /* z-coordinates of the vertices */
)
{
  /* Local declarations. */
  int i, j, k, start_id, elem_id, local_id;
  char cmesg[256];
/***************************** BEGIN EXECUTION ******************************/

  start_id = vtxdist[Proc]+1;  /* global ids start at 1 */

  for (i = 0; i < Mesh.num_elems; i++) {
    elem[i].globalID = start_id + i;
    if (vwgts != NULL)
      elem[i].cpu_wgt = vwgts[i];
    else
      elem[i].cpu_wgt = 1.0;
    elem[i].elem_blk = 0;        /* only one element block for all vertices */
    if (Mesh.num_dims > 0) {
      /* One set of coords per element. */
      elem[i].connect = (int *) malloc(sizeof(int));
      elem[i].connect[0] = elem[i].globalID;
      elem[i].coord = (float **) malloc(sizeof(float *));
      elem[i].coord[0] = (float *) malloc(Mesh.num_dims * sizeof(float));  
      elem[i].coord[0][0] = x[i];
      if (Mesh.num_dims > 1) {
        elem[i].coord[0][1] = y[i];
        if (Mesh.num_dims > 2) {
          elem[i].coord[0][2] = z[i];
        }
      }
    }

    /* now start with the adjacencies */
    if (start != NULL)
      elem[i].nadj = start[i+1] - start[i];
    else
      elem[i].nadj = 0;
    if (elem[i].nadj > 0) {
      elem[i].adj_len = elem[i].nadj;
      elem[i].adj = (int *) malloc (elem[i].nadj * sizeof(int));
      elem[i].adj_proc = (int *) malloc (elem[i].nadj * sizeof(int));
      if (!(elem[i].adj) || !(elem[i].adj_proc)) {
        Gen_Error(0, "fatal: insufficient memory");
        return 0;
      }
      if (ewgts != NULL) {
        elem[i].edge_wgt = (float *) malloc (elem[i].nadj * sizeof(float));
        if (!(elem[i].edge_wgt)) {
          Gen_Error(0, "fatal: insufficient memory");
          return 0;
        }
      }
      else
        elem[i].edge_wgt = NULL;

      for (j = 0; j < elem[i].nadj; j++) {
        elem_id = adj[start[i] + j];

        /* determine which processor the adjacent vertex is on */
        for (k = 0; k < Num_Proc; k++)
          /* Compare with <= since elem_id is 1-based and vtxdist is 0-based. */
          if (elem_id <= vtxdist[k+1]) break;

        /* sanity check */
        if (k == Num_Proc) {
          sprintf(cmesg, 
                  "fatal:  adjacent element %d not in vtxdist array (%d,%d)",
                  elem_id, i, j);   
          Gen_Error(0, cmesg);
          return 0;
        } 

        /*
         * if the adjacent element is on this processor
         * then find the local id for that element
         */
        if (k == Proc) {
          local_id = elem_id - start_id;
          elem[i].adj[j] = local_id;
        }
        else /* use the global id */
          elem[i].adj[j] = elem_id;

        elem[i].adj_proc[j] = k;

        if (ewgts != NULL)
          elem[i].edge_wgt[j] = ewgts[start[i] + j];
      }
    } /* End: "if (elem[i].nadj > 0)" */
  } /* End: "for (i = 0; i < Mesh.num_elems; i++)" */

  print_distributed_mesh(Proc, Num_Proc, prob, elem);

  return 1;
}
