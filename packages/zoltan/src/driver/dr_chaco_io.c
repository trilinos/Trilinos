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
#include "ch_input_const.h"
#include "ch_input.h"

int fill_elements(int, int, ELEM_INFO *, int, int *, int *, int *, int *,
                  float *);


int read_chaco_mesh(int Proc,
                    int Num_Proc,
                    PROB_INFO_PTR prob,
                    PARIO_INFO_PTR pio_info,
                    ELEM_INFO **elements)
{
  /* Local declarations. */
  char   cmesg[256];

  int    i, start_id, nvtxs;
  int   *start, *adj, *vwgts, *vtxdist;

  float *ewgts;

  FILE  *fp;
/***************************** BEGIN EXECUTION ******************************/

  if (Proc == 0) {
    fp = fopen(pio_info->pexo_fname, "r");

    /* read the array in on processor 0 */
    if (chaco_input_graph(fp, pio_info->pexo_fname, &start, &adj, &nvtxs,
                           &vwgts, &ewgts) != 0) {
      Gen_Error(0, "fatal: Error returned from chaco_input_graph");
      return 0;
    }

    fclose (fp);
  }

  /* Distribute graph */
  vtxdist = NULL;
  if (chaco_dist_graph(MPI_COMM_WORLD, 0, &nvtxs, &vtxdist, &start, &adj, 
                   &vwgts, &ewgts) != 0) {
      Gen_Error(0, "fatal: Error returned from chaco_dist_graph");
      return 0;
  }

  Mesh.num_elems = nvtxs;
  Mesh.elem_array_len = Mesh.num_elems + 5;

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
  if (!fill_elements(Proc, Num_Proc, *elements, nvtxs, vtxdist, start, adj,
                     vwgts, ewgts)) {
    Gen_Error(0, "fatal: Error returned from fill_elements");
    return 0;
  }

  return 1;
}

int fill_elements(
  int        Proc,
  int        Num_Proc,
  ELEM_INFO  elem[],             /* array of element information */
  int        nvtxs,              /* number of vertices in graph */
  int       *vtxdist,            /* vertex distribution data */
  int       *start,              /* start of edge list for each vertex */
  int       *adj,                /* edge list data */
  int       *vwgts,              /* vertex weight list data */
  float     *ewgts               /* edge weight list data */
)
{
  /* Local declarations. */
  int i, j, k, start_id, elem_id, local_id;
/***************************** BEGIN EXECUTION ******************************/

  start_id = vtxdist[Proc];

  for (i = 0; i < Mesh.num_elems; i++) {
    elem[i].globalID = start_id + i + 1; /* global ids start at 1 */
    if (vwgts != NULL)
      elem[i].cpu_wgt = vwgts[i];
    elem[i].cpu_wgt = 0.0;
    elem[i].elem_blk = 1;        /* only one element block for all vertices */
    elem[i].coord = NULL;
    elem[i].connect = NULL;

    /* now start with the adjacencies */
    elem[i].nadj = start[i+1] - start[i];

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
          if (elem_id < vtxdist[k]) break;

        /* sanity check */
        if (k == Num_Proc) {
          Gen_Error(0, "fatal: adjacent element not in vtxdist array");
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

  return 1;
}
