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
 * they are dealt out like cards based on process rank.
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

  /* By setting parameters, I previously gave the Zoltan library 
   * the length of my global IDs, local IDs, and object weights.
   */

  if ( (sizeGID != 1) ||  /* My global IDs are 1 integer */
       (sizeLID != 1) ||  /* My local IDs are 1 integer */
       (wgt_dim != 1)){   /* My vertex weights are 1 float */
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (i=0, next=0; i<simpleNumVertices; i++){
    if (i % numProcs == myRank){
      globalID[next] = i+1;   /* application wide global ID */
      localID[next] = next;   /* process specific local ID  */
      obj_wgts[next] = (float)simpleNumEdges[i];  /* weight */
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

  if ( (sizeGID != 1) || (sizeLID != 1) ||
       (wgt_dim != 0)){      /* we are not using edge weights */
    *ierr = ZOLTAN_FATAL;
    return;
  }

  nextID = nborGID;
  nextProc = nborProc;

  for (i=0;  i < num_obj ; i++){
    idx = globalID[i] - 1;
    if (num_edges[i] != simpleNumEdges[idx]){
      *ierr = ZOLTAN_FATAL;
      return;
    }
    for (j=0; j<num_edges[i]; j++){
      *nextID++ = simpleEdges[idx][j];
      *nextProc++ = (simpleEdges[idx][j] - 1) % numProcs;
    }
  }

  *ierr = ZOLTAN_OK;

  return;
}
#endif
