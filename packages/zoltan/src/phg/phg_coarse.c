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

#include "phypergraph.h"

/****************************************************************************/



/*Procedure to coarsen a hypergraph based on a matching. All vertices of one match
  are clustered to one vertex. Identical hyperedges are collapsed to a single
  hyperedge with combined weight. The array LevelMap is the mapping of the old
  vertices to the new vertices. It will be used to pass a partition of the coarse
  graph back to the original graph. Time O(|I|*log(|I|)), due to sorting. */
  


int Zoltan_PHG_Coarsening
( ZZ *zz,             /* the Zoltan data structure */
  PHGraph  *hg,
  int      *match,      /* Matching, Packing or Grouping array */
  PHGraph  *c_hg,       /* points to a working copy of hg structure */
  int      *LevelMap)
{
  int i, j, k, l, old, vertex, new_vertex, deleted_he, deleted_pins, *hsize=NULL;
  int *sum=NULL, *used_vertices=NULL, *sorted=NULL, *c_hindex=NULL; 
  int *c_hvertex=NULL;
  float *c_ewgt=NULL;
  char *yo = "Zoltan_PHG_Coarsening";

  ZOLTAN_TRACE_ENTER(zz, yo);

  Zoltan_PHG_HGraph_Init(c_hg);
  c_hg->info  = hg->info + 1;
  c_hg->ratio = hg->ratio;
  c_hg->redl  = hg->redl;

/* RTHRTH:  CHANGES FOR PARALLEL FOR THIS ROUTINE: */  
/* Assume that the entire column has the matching information at this point */
/* Process entire match array.  For each entry that is not local and I own, 
/* (bottom bit of sum of global matched vertices)c
/* create a message (requesting vertex information) to its processor.
/* Send all messages using unstructured communications
/* Process received messages:  create response message with GID, weight, edge list
/* Use hash table to associate external GIDs with the response information
/* Scan match array to create LevelMap assuming local matching, packing, grouping but
/* external vertices may only be matched.
/* Instead of processing edge list with vertex pins (serial code), process vertex list
/*  with hyperedge pins. Use same trick to insure each edge is added once.
/* If vertex is external, use hash lookup to get its info.  Add hyperedge pins and
/*  sum vertex weights.
/* Create mirror for opposite information.
  
  
  
  
  /* Calculate the number of coarse vertices. match[vertex] -> -match[vertex]-1 */
  c_hg->nVtx = 0;
  for (i = 0; i < hg->nVtx; i++)
    if (match[i] >= 0) {
      (c_hg->nVtx)++;
      vertex = i;
      while (match[vertex] >= 0) {
        old       =  vertex;
        vertex    =  match[old];
        match[old] = -match[old] - 1;
      }
  }
  if (!(c_hg->vwgt = (float*) ZOLTAN_CALLOC (c_hg->nVtx, sizeof(float)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT (zz, yo);
    return ZOLTAN_MEMERR;
  }

  /* Construct the LevelMap. match[vertex] is changed back to original value.
   * Also sum up the vertex weights. */
  new_vertex = 0;
  for (i = 0; i < hg->nVtx; i++)
    if (match[i] < 0) {
      vertex = i;
      while (match[vertex] < 0) {
        LevelMap[vertex] = new_vertex;
        c_hg->vwgt[new_vertex] += hg->vwgt ? hg->vwgt[vertex] : 1.0;
        match[vertex] = -match[vertex] - 1;
        vertex       =  match[vertex];
      }
      new_vertex++;
    }

  /* Coarsen the hyperedges and avoid edges with only one vertex */
  if (!(used_vertices = (int*)  ZOLTAN_CALLOC (c_hg->nVtx,    sizeof(int)))
   || !(c_ewgt        = (float*)ZOLTAN_MALLOC (hg->nEdge    * sizeof(float)))
   || !(c_hindex      = (int*)  ZOLTAN_MALLOC ((hg->nEdge+1)* sizeof(int)))
   || !(c_hvertex     = (int*)  ZOLTAN_MALLOC (hg->nInput   * sizeof(int))) ) {
      Zoltan_Multifree (__FILE__, __LINE__, 4, &used_vertices, &c_ewgt,
       &c_hindex, &c_hvertex);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }
  c_hindex[0] = c_hindex[1] = c_hg->nEdge = c_hg->nInput = 0;
  for (i = 0; i < hg->nEdge; i++) {
    for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++) {
      new_vertex = LevelMap[hg->hvertex[j]];
      if (used_vertices[new_vertex] <= i) {
        used_vertices[new_vertex] = i + 1;
        c_hvertex[(c_hg->nInput)++] = new_vertex;
      }
    }
    if (c_hg->nInput > c_hindex[c_hg->nEdge] + 1) {  /* a new hyperedge */
      c_ewgt[c_hg->nEdge++] = hg->ewgt ? hg->ewgt[i] : 1.0;
      c_hindex[c_hg->nEdge] = c_hg->nInput;
    }
    else /* no new hyperedge, because it only covers one vertex */
      c_hg->nInput = c_hindex[c_hg->nEdge];
  }
  ZOLTAN_FREE((void**) &used_vertices);

  /* Done if there are no remaining edges */
  if (c_hg->nEdge == 0) {
    c_hg->ewgt = NULL;
    if (!(c_hg->hindex = (int*) ZOLTAN_CALLOC (1, sizeof(int)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
    }
    c_hg->hvertex = NULL;
    Zoltan_Multifree (__FILE__, __LINE__, 3, &c_ewgt, &c_hindex, &c_hvertex);
    ZOLTAN_TRACE_EXIT (zz, yo);
    return Zoltan_PHG_Create_Mirror(zz, c_hg);
  }

  /* RTHRTH: NOTE removed code per Umit's speedup hack from serial version HERE*/  
  c_hg->ewgt    = c_ewgt;
  c_hg->hindex  = c_hindex;
  c_hg->hvertex = c_hvertex;

  ZOLTAN_TRACE_EXIT (zz, yo);
  return Zoltan_PHG_Create_Mirror(zz, c_hg);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
