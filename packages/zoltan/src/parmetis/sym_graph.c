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

#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "all_allo_const.h"
#include "third_library_const.h"
#include "third_library_tools.h"
#include "params_const.h"

#ifdef __GNUC__
#define inline static __inline__
#else
#ifndef __USE_ISOC99
#define inline static
#endif
#endif

typedef struct Zoltan_Arcs_ {
  indextype src;
  indextype tgt;
} Zoltan_Arcs;

typedef struct Zoltan_Weighted_Arcs_ {
  Zoltan_Arcs arc;
  float       wgt;
} Zoltan_Weighted_Arcs;


inline int
compar_unweighted(const Zoltan_Arcs* arc1, const Zoltan_Arcs* arc2)
{
  if (arc1->src == arc2->src)
    return (arc1->tgt - arc2->tgt);
  return (arc1->src - arc2->src);
}

inline int
compar_weighted(const Zoltan_Weighted_Arcs* arc1, const Zoltan_Weighted_Arcs* arc2)
{
  return compar_unweighted(&arc1->arc, &arc2->arc);
}


int
give_proc (indextype vertex, const indextype *vtxdist, int numProc, int *myproc)
{
  int currentproc;


  if ((((*myproc) >= 0) && (*myproc) < numProc) &&
    (vertex >= vtxdist[*myproc]) && (vertex < vtxdist[*myproc+1])) {
    return (*myproc);
  }

  currentproc = vertex / (vtxdist[1]-vtxdist[0]);  /* Assume that vertices are balanced */

  if (currentproc >= numProc)
    currentproc = numProc - 1;

  if ((vertex < vtxdist[0])||( vertex >= vtxdist[numProc])) {
    ZOLTAN_PRINT_WARN ((*myproc), "Symmetrize Graph problem (1)", "Unknown vertex");
    return (-1);
  }

  while (1) {
    if (vertex >= vtxdist[currentproc + 1]) {
      currentproc ++;
      continue;
    }
    if (vertex < vtxdist[currentproc]) {
      currentproc --;
      continue;
    }
    break;
  }

  *myproc =currentproc;
  return (currentproc);
}


  /*
    This function create a couple (src, tgt) in the local (rcv) array and a (tgt,src) in
    the local or global array
  */

inline void
add_edge(indextype src, indextype tgt, float wgt,
	 Zoltan_Arcs* sndarctab, Zoltan_Weighted_Arcs *sndarcwgttab, int* proctab,
	 Zoltan_Arcs* rcvarctab, Zoltan_Weighted_Arcs *rcvarcwgttab,
	 int *cursersnd, int* curserrcv, Zoltan_Arcs **hashtab, int nnz_loc,
	 const indextype *vtxdist, int procnum, int myproc)
{
    unsigned int hash;
    int couple[2];
    int i;

    couple[0] = src; couple [1] = tgt;

    for (i = 0; i < 2 ; ++i) {
      hash = Zoltan_Hash((ZOLTAN_ID_PTR)couple, 2, 2*nnz_loc);
      if (hash < 0) return;

      while((hashtab[hash] != NULL) && ((hashtab[hash]->src != couple[0])
	    || (hashtab[hash]->tgt != couple[1]))){
	hash++;
	if (hash >= 2*nnz_loc)
	  hash =0;
      }
      if (hashtab[hash] != NULL) /* If already in the hash table */
	return;

      if (rcvarcwgttab != NULL) {
	rcvarcwgttab[*curserrcv].wgt = wgt;
	hashtab[hash] = &(rcvarcwgttab[*curserrcv].arc);
      }
      else {
	hashtab[hash] = &(rcvarctab[*curserrcv]);
      }

      (*curserrcv)++;
      hashtab[hash]->src = couple[0];
      hashtab[hash]->tgt = couple[1];

      if ((tgt >= vtxdist[myproc]) && (tgt < vtxdist[myproc+1])) {
	couple[0] = tgt; couple [1] = src;
      }
      else { /* Non local vertex */
	if (sndarcwgttab != NULL) {
	  sndarcwgttab[*cursersnd].arc.src = tgt;
	  sndarcwgttab[*cursersnd].arc.tgt = src;
	  sndarcwgttab[*cursersnd].wgt = wgt;
	}
	else {
	  sndarctab[*cursersnd].src = tgt;
	  sndarctab[*cursersnd].tgt = src;
	}
	proctab[*cursersnd]=procnum;
	(*cursersnd)++;
	break;
      }
    }

}


/*
 * Symmetrize the 1D distributed graph
 *
 * We don't deal with pseudo-symmetric graph, ie symmetric for arcs but
 * with different edge weights.
 * Maybe we need to sum weigts to deal with such arcs ?
 */

int Zoltan_Symmetrize_Graph(
    const ZZ *zz, int graph_type, int check_graph, int num_obj,
    ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
    int obj_wgt_dim, int edge_wgt_dim,
    indextype ** vtxdist, indextype **xadj, indextype **adjncy,
    float **ewgts, int **adjproc)
{
  /* Local variables */
  indextype vertsrc;
  int edgenum;
  int nnz_loc;
  int numProc, myproc, currentproc;
  Zoltan_Arcs *sndarctab = NULL, *rcvarctab = NULL;
  Zoltan_Weighted_Arcs *sndarcwgttab = NULL, *rcvarcwgttab = NULL;
  int cursersnd, curserrcv;
  Zoltan_Arcs **hashtab;
  int *proctab;
  int rcvsize = 0;
  indextype prevsrc, prevtgt, currentvtx;
  int currentedge;
  int edge;


  MPI_Comm comm = zz->Communicator;
  ZOLTAN_COMM_OBJ *comm_plan;

  char *yo = "Zoltan_Symmetrize_Graph";
  int ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);

  numProc = zz->Num_Proc;
  myproc = zz->Proc;

  nnz_loc = (*xadj)[num_obj];

  hashtab = (Zoltan_Arcs**)ZOLTAN_MALLOC(2*nnz_loc*sizeof(Zoltan_Arcs*));
  if (edge_wgt_dim == 1){
    sndarcwgttab = (Zoltan_Weighted_Arcs*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(Zoltan_Weighted_Arcs));
    rcvarcwgttab = (Zoltan_Weighted_Arcs*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(Zoltan_Weighted_Arcs));
  }
  else {
    sndarctab = (Zoltan_Arcs*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(Zoltan_Arcs));
    rcvarctab = (Zoltan_Arcs*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(Zoltan_Arcs));
  }
  proctab = (int*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(int));

  if ((nnz_loc >0 ) && (hashtab == NULL ||
      (sndarcwgttab == NULL && sndarctab == NULL) || (rcvarcwgttab == NULL && rcvarctab == NULL)
      || proctab == NULL)) {
    ZOLTAN_FREE(&proctab);
    ZOLTAN_FREE(&rcvarctab);
    ZOLTAN_FREE(&rcvarcwgttab);
    ZOLTAN_FREE(&sndarctab);
    ZOLTAN_FREE(&sndarcwgttab);
    ZOLTAN_FREE(&hashtab);
    ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory");
    return (ZOLTAN_MEMERR);
  }

  for (edgenum = 0 ; edgenum < 2*nnz_loc ; edgenum++)
   hashtab[edgenum] = NULL;


  /* Construct a locally symmetrized graph */
  for (vertsrc = 0, curserrcv=0, cursersnd =0 ; vertsrc < num_obj ; vertsrc ++) { /* For all vertices */
    for (edgenum = (*xadj)[vertsrc] ; edgenum < (*xadj)[vertsrc+1] ; ++edgenum) { /* For all neighbors */
      float wgt = 0;

      if (edge_wgt_dim == 1) {
	wgt = *ewgts[edgenum];
      }
      add_edge(vertsrc + (*vtxdist)[myproc], (*adjncy)[edgenum], wgt,
	       sndarctab, sndarcwgttab, proctab, rcvarctab, rcvarcwgttab,
	       &cursersnd, &curserrcv, hashtab, nnz_loc, *vtxdist, (*adjproc)[edgenum], myproc);
    }
  }
  ZOLTAN_FREE(&hashtab);

  /* These arrays may be reconstructed */
  ZOLTAN_FREE(ewgts);
  ZOLTAN_FREE(adjncy);
  ZOLTAN_FREE(adjproc);

  /* Communicate the arcs */
  Zoltan_Comm_Create(&comm_plan, cursersnd, proctab, comm, 6241984, &rcvsize);
  rcvsize += curserrcv;
  if ((rcvsize > 0) && (rcvsize >= 2*nnz_loc)) { /* reception buffer is too small */
    if (edge_wgt_dim == 1) {
      rcvarcwgttab = (Zoltan_Weighted_Arcs*)
	ZOLTAN_REALLOC(rcvarcwgttab, rcvsize*sizeof(Zoltan_Weighted_Arcs));
    }
    else {
      rcvarctab = (Zoltan_Arcs*) ZOLTAN_REALLOC(rcvarctab, rcvsize*sizeof(Zoltan_Arcs));
    }
    if ((rcvarcwgttab == NULL && rcvarctab == NULL)) {
      Zoltan_Comm_Destroy (&comm_plan);
      ZOLTAN_FREE(&proctab);
      ZOLTAN_FREE(&rcvarctab);
      ZOLTAN_FREE(&rcvarcwgttab);
      ZOLTAN_FREE(&sndarctab);
      ZOLTAN_FREE(&sndarcwgttab);
      ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory (2)")
      return (ZOLTAN_MEMERR);
    }
  }

  if (edge_wgt_dim == 1) {
    Zoltan_Comm_Do (comm_plan, 8161982, (char*)sndarcwgttab, sizeof(Zoltan_Weighted_Arcs),
			 (char*)(rcvarcwgttab + curserrcv));
  }
  else {
    Zoltan_Comm_Do (comm_plan, 8161982, (char*)sndarctab, sizeof(Zoltan_Arcs),
			 (char*)(rcvarctab + curserrcv));
  }

  /* I use blocking communication plan because I can free the send buffer before doing allocations */
  ZOLTAN_FREE(&sndarctab);
  ZOLTAN_FREE(&sndarcwgttab);
  ZOLTAN_FREE(&proctab);
  Zoltan_Comm_Destroy (&comm_plan);

  if (rcvsize > 0) {
    *adjncy = (indextype*)ZOLTAN_MALLOC(rcvsize*sizeof(indextype));
    if ((*adjncy) == NULL) {
      ZOLTAN_FREE(&rcvarctab);
      ZOLTAN_FREE(&rcvarcwgttab);
      ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory (3)");
      return (ZOLTAN_MEMERR);
    }
    *adjproc = (indextype*)ZOLTAN_MALLOC(rcvsize*sizeof(indextype));
    if ((*adjproc) == NULL) {
      ZOLTAN_FREE(&rcvarctab);
      ZOLTAN_FREE(&rcvarcwgttab);
      ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory (4)");
      return (ZOLTAN_MEMERR);
    }

    if (edge_wgt_dim == 1) {
      *ewgts = (float*)ZOLTAN_MALLOC(rcvsize*sizeof(float));
      if ((*ewgts) == NULL) {
	ZOLTAN_FREE(&rcvarctab);
	ZOLTAN_FREE(&rcvarcwgttab);
	ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory (5)");
	return (ZOLTAN_MEMERR);
      }
    }
  }

  /* Sort Edges */
  if (edge_wgt_dim == 1) {
    qsort ((void*)rcvarcwgttab, rcvsize, sizeof(Zoltan_Weighted_Arcs),
	   (int (*)(const void*,const void*)) compar_weighted);
  }
  else {
    qsort ((void*)rcvarctab, rcvsize, sizeof(Zoltan_Arcs),
	   (int (*)(const void*,const void*))compar_unweighted);
  }

  /* Now I only need to copy these informations in the edge array and to recompute
     the indirection array */
  prevsrc = (*vtxdist)[myproc] ; prevtgt = -1;
  (*xadj)[0] = 0; currentedge = 0; currentvtx=1;
  currentproc = 0;
  for (edge = 0 ; edge < rcvsize ; ++ edge) {
    Zoltan_Arcs *arc;

    if (edge_wgt_dim == 1)
      arc = &(rcvarcwgttab[edge].arc);
    else
      arc = &(rcvarctab[edge]);

    if (arc->src != prevsrc) {
#ifdef CC_DEBUG
      if ((arc->src < (*vtxdist)[myproc])|| (arc->src >= (*vtxdist)[myproc+1]))
	*(int*)NULL = 0;
#endif /* CC_DEBUG */
      do {
	(*xadj)[currentvtx++] = currentedge;
      } while (currentvtx + (*vtxdist)[myproc] <= arc->src); /* currentvtx is +1 for ending */
      prevsrc = arc->src;
      currentproc = 0;                /* Reset procid for searching ownership */
    }
    else { /* If the same src vertex, check if it's the same edge */
      if (prevtgt == arc->tgt)
	continue;
    }

    prevtgt = arc->tgt;
    (*adjncy)[currentedge] = arc->tgt;
    (*adjproc)[currentedge] = give_proc(arc->tgt, *vtxdist, numProc, &currentproc);
    if (edge_wgt_dim == 1)
      (*ewgts)[currentedge] = rcvarcwgttab[edge].wgt;

    currentedge++;
  }
  (*xadj)[num_obj] = currentedge;


  /* TODO: Add Realloc to be memory more efficient */

  ZOLTAN_FREE(&rcvarctab);
  ZOLTAN_FREE(&rcvarcwgttab);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}

/**************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
