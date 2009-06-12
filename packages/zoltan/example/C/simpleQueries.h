/***************************************************************
** Application defined query functions
***************************************************************/

#ifndef SIMPLEQUERIES_H
#define SIMPLEQUERIES_H

#include "zoltan.h"
#include "simpleGraph.h"

static int myRank;
static int numProcs;

/* 
 **************************************************************
 * Prototype: ZOLTAN_NUM_OBJ_FN
 * Return the number of objects I own.  
 **************************************************************
 * Zoltan partitions a collection of objects distributed
 * across processes. In this example objects are vertices, and 
 * they are dealt out like cards (cyclic) based on process rank.
 */
static int get_number_of_objects(void *data, int *ierr)
{
int i, numobj=0;

  for (i=0; i<simpleNumVertices; i++){
    if (i % numProcs == myRank) numobj++;
  }

  *ierr = ZOLTAN_OK;

  return numobj;
}
/* 
 **************************************************************
 * Prototype: ZOLTAN_OBJ_LIST_FN
 * Return list of my objects, with optional object weights.  
 **************************************************************
 * Zoltan requires that objects be identified with an application
 * global ID.  The ID for a single object can be an integer array 
 * of arbitrary length.  Each process can optionally provide a
 * local ID for each object.  In this case, Zoltan will provide
 * the local ID in addition to the application global ID when 
 * querying the process.  If the application indicated that it 
 * will be supplying object weights, the application provides
 * those as well with this query function.
 */
static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
int i, next;
int j, k;

  /* By setting parameters, I previously gave the Zoltan library 
   * the length of my global IDs, local IDs, and object weights.
   */

  if (sizeGID != 1){  /* My global IDs are 1 integer */
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (i=0, next=0, j=0; i<simpleNumVertices; i++){
    if (i % numProcs == myRank){
      globalID[next] = i+1;   /* application wide global ID */
      localID[next] = next;   /* process specific local ID  */
      if (wgt_dim>0)
      for (k=0; k < wgt_dim; k++){
        obj_wgts[j++] = (float)simpleNumEdges[i] * (k+1);  /* weight */
      }
      next++;
    }
  }

  *ierr = ZOLTAN_OK;

  return;
}
/* 
 **************************************************************
 * Prototype: ZOLTAN_NUM_GEOM_FN
 * Return the dimension of a vertex, for geometric methods
 **************************************************************
 */
static int get_num_geometry(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 2;
}
/* 
 **************************************************************
 * Prototype: ZOLTAN_GEOM_MULTI_FN
 * Return the coordinates of my objects (vertices), for
 * geometric methods.
 **************************************************************
 */
static void get_geometry_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr)
{
int i;
int row, col;

  if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 2)){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (i=0;  i < num_obj ; i++){
    row = (globalID[i] - 1) / 5;
    col = (globalID[i] - 1) % 5;

    geom_vec[2*i] = (double)col;
    geom_vec[2*i + 1] = (double)row;
  }

  *ierr = ZOLTAN_OK;
  return;
}
/* 
 **************************************************************
 * Prototype: ZOLTAN_NUM_EDGES_MULTI_FN
 * Return the number of edges for each vertex in the ID lists.
 * For graph methods.
 **************************************************************
 */
static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int *numEdges, int *ierr)
{
int i, idx;

  if ( (sizeGID != 1) || (sizeLID != 1)){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (i=0;  i < num_obj ; i++){
    idx = globalID[i] - 1;
    numEdges[i] = simpleNumEdges[idx];
  }

  *ierr = ZOLTAN_OK;
  return;
}
/* 
 **************************************************************
 * Prototype: ZOLTAN_EDGE_LIST_MULTI_FN
 * For every vertex in the ID list, return a list of all its
 * adjacent vertices, and the processes on which they reside.
 * Also include the edge weights if any.
 * For graph methods.  
 **************************************************************
 */
static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges, 
        ZOLTAN_ID_PTR nborGID, int *nborProc,
        int wgt_dim, float *ewgts, int *ierr)
{
int i, j, idx;
ZOLTAN_ID_PTR nextID; 
int *nextProc;
#ifdef USE_EDGE_WEIGHTS
float *nextWeight;
#endif

  if ( (sizeGID != 1) || (sizeLID != 1) ){
    *ierr = ZOLTAN_FATAL;
    return;
  }

#ifdef USE_EDGE_WEIGHTS
  if ( wgt_dim != 1){    
    *ierr = ZOLTAN_FATAL;
    return;
  }
#else
  if ( wgt_dim != 0){    
    *ierr = ZOLTAN_FATAL;
    return;
  }
#endif

  nextID = nborGID;
  nextProc = nborProc;
#ifdef USE_EDGE_WEIGHTS
  nextWeight = ewgts;
#endif

  for (i=0;  i < num_obj ; i++){
    idx = globalID[i] - 1;
    if (num_edges[i] != simpleNumEdges[idx]){
      *ierr = ZOLTAN_FATAL;
      return;
    }
    for (j=0; j<num_edges[i]; j++){
      *nextID++ = simpleEdges[idx][j];
      *nextProc++ = (simpleEdges[idx][j] - 1) % numProcs;
#ifdef USE_EDGE_WEIGHTS
      /* Test that Zoltan creates edge weights equal to 3 */
      if (globalID[i] < simpleEdges[idx][j])
        *nextWeight++ = 1.0;
      else
        *nextWeight++ = 2.0;
#endif
    }
  }

  *ierr = ZOLTAN_OK;

  return;
}
/* 
 **************************************************************
 * Prototype: ZOLTAN_HG_SIZE_CS_FN
 * This function tells the Zoltan library how many hyperedges
 * I will supply, and how many total pins (vertices) are in those
 * hyperedges, and my choice of compressed hyperedge format.
 *
 * For hypergraph methods.
 *
 * We will create a hypergraph out of the simple graph in the
 * following way.  For every vertex in the graph, we will 
 * create a hyperedge connecting that vertex and each of its
 * neighbors.
 *
 * We will supply this hypergraph to Zoltan in compressed
 * hyperedge format.  (Compressed vertices is other choice.  See
 * the Zoltan User's Guide section on this query function for
 * more information on compressed formats.)
 *
 * We will initially divide the hyperedges across the processes
 * in a round robin manner.  However Zoltan does not require 
 * that each process supply an entire hyperedge.  Hyperedges 
 * could be distributed across processes. 
 **************************************************************
 */

static struct _hg_data{
  int numEdges;
  int numPins;
}hg_data;

void get_hg_size(void *data, int *num_lists, int *num_pins, 
                 int *format, int *ierr)
{
int i;
struct _hg_data *hgd;

  hgd = (struct _hg_data *)data;

  hgd->numEdges = 0;
  hgd->numPins  = 0;

  for (i=0; i<simpleNumVertices; i++){
    if (i % numProcs == myRank){    
      hgd->numPins += simpleNumEdges[i] + 1;  /* my hyperedge */
      hgd->numEdges++;
    }
  }

  *num_lists = hgd->numEdges;
  *num_pins  = hgd->numPins;
  *format    = ZOLTAN_COMPRESSED_EDGE;
  *ierr      = ZOLTAN_OK;

  return;
}
/* 
 **************************************************************
 * Prototype: ZOLTAN_HG_CS_FN
 *
 * Return a compressed list describing the hypergraph sparse matrix.
 * We can think of the columns as vertices and the rows as
 * hyperedges.  For each row we have ones in the column representing
 * vertices in the hyperedge, and zeroes in the other entries.
 * 
 * For hypergraph methods.
 **************************************************************
 */
void get_hg(void *data, int sizeGID, int num_rows, int num_pins, 
            int format, ZOLTAN_ID_PTR edge_GID, int *edge_ptr, 
            ZOLTAN_ID_PTR pin_GID, int *ierr) 
{
int i, j, npins;
struct _hg_data *hgd;

  hgd    = (struct _hg_data *)data;

  if ((num_rows != hgd->numEdges) || (num_pins != hgd->numPins) ||
      (format != ZOLTAN_COMPRESSED_EDGE)){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  npins = 0;

  for (i=0; i<simpleNumVertices; i++){
    if (i % numProcs == myRank){    /* my hyperedge */
      *edge_ptr++ = npins;      /* index into start of pin list */
      *edge_GID++ = i+1;        /* hyperedge global ID */

      /* list global ID of each pin (vertex) in hyperedge */
      for (j=0; j<simpleNumEdges[i]; j++){
        *pin_GID++ = simpleEdges[i][j];
      }
      *pin_GID++ = i+1;
      npins += simpleNumEdges[i] + 1;
    }
  }

  *ierr = ZOLTAN_OK;

  return;
}
/* 
 **************************************************************
 * Prototype: ZOLTAN_HG_SIZE_EDGE_WTS_FN
 * This query function returns the number of hyperedges for which
 * the process will supply edge weights.  The edges for which
 * a process supplies edge weights need not be the same edges
 * for which it supplied pins.  Multiple processes can supply
 * weights for the same edge.  The multiple weights are resolved
 * according the the setting of the PGH_EDGE_WEIGHT_OPERATION 
 * parameter.  Edge weights are optional.
 *
 * For hypergraph methods.
 **************************************************************
 */
void get_hg_num_edge_weights(void *data, int *num_edges, int *ierr) 
{
struct _hg_data *hgd;

  hgd    = (struct _hg_data *)data;

  *num_edges = hgd->numEdges;
  *ierr = ZOLTAN_OK;

  return;
}
/* 
 **************************************************************
 * Prototype: ZOLTAN_HG_EDGE_WTS_FN
 * This query function supplies edge weights to the
 * Zoltan library.  The edge global ID corresponds to the
 * hyperedge global IDs supplied in the ZOLTAN_HG_CS_FN.
 * The edge local IDs are for consistency with the rest
 * of the Zoltan library, but at this point in time, there
 * is no interface that uses local IDs to refer to hyperedges.
 *
 * The edge weight dimension was supplied by the application 
 * as the value of the EDGE_WEIGHT_DIM parameter.  If the dimension 
 * is greater than one, list all weights for the first edge,
 * followed by all weights for the second edge, and so on.
 *
 * For hypergraph methods.
 **************************************************************
 */
void get_hyperedge_weights(void *data, int sizeGID, 
         int sizeLID, int num_edges, int edge_weight_dim, 
         ZOLTAN_ID_PTR edge_GID, ZOLTAN_ID_PTR edge_LID, 
         float  *edge_weight, int *ierr) 
{
int i;
struct _hg_data *hgd;

  hgd    = (struct _hg_data *)data;

  if ((sizeGID != 1) || (sizeLID != 0) || 
      (num_edges != hgd->numEdges) || (edge_weight_dim != 1)){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (i=0; i<simpleNumVertices; i++){
    if (i % numProcs == myRank){ 
      *edge_GID++ = i+1;        /* hyperedge global ID */
      *edge_weight++ = 1.0;  /* all hyperedges same weight */
    }
  }
  *ierr = ZOLTAN_OK;
  return;
}
/* 
 **************************************************************
 * Prototype: ZOLTAN_OBJ_LIST_FN
 * Return list of my objects.
 * This version assumes we are not supplying object local IDs or
 * object weights.  We'll use it for the hypergraph example.
 **************************************************************
 */
static void get_hg_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
int i, next;

  if (sizeGID != 1){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (i=0, next=0; i<simpleNumVertices; i++){
    if (i % numProcs == myRank){
      globalID[next] = i+1;   /* application wide global ID */
      next++;
    }
  }

  *ierr = ZOLTAN_OK;

  return;
}
#endif
