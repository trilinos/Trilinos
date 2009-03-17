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
    return (arc2->tgt - arc2->tgt);
  return (arc2->src - arc2->src);
}

inline int
compar_weighted(const Zoltan_Weighted_Arcs* arc1, const Zoltan_Weighted_Arcs* arc2)
{
  return compar_unweighted(&arc1->arc, &arc2->arc);
}


static int
give_proc (indextype vertex, const indextype *vtxdist, int  procnum, int myproc)
{
  int currentproc;

  if ((vertex >= vtxdist[myproc]) && (vertex < vtxdist[myproc+1])) {
    return (myproc);
  }

  currentproc = vertex / (vtxdist[1]-vtxdist[0]);  /* Assume that vertices are balanced */

  if (currentproc >= procnum)
    currentproc = procnum - 1;

  if ((vertex < vtxdist[0])||( vertex >= vtxdist[procnum])) {
    ZOLTAN_PRINT_WARN (myproc, "Symmetrize Graph problem (1)", "Unknown vertex");
    return (-1);
  }

  while ((vertex < vtxdist[currentproc]) || (vertex >= vtxdist[currentproc + 1])) {
    currentproc ++;
  }

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
	 int *cursersnd, int* curserrcv, int *hashtab, int nnz_loc,
	 const indextype *vtxdist, int procnum, int myproc)
{
  unsigned int hash;
  int couple[2];
  int i;

  couple[0] = src;
  couple[1] = tgt;

  for (i=0 ; i < 2 ; ++i) {
    Zoltan_Arcs* arctab;
    Zoltan_Weighted_Arcs* arcwgttab;
    int *current;

    hash = Zoltan_Hash((ZOLTAN_ID_PTR)couple, 2, 2*nnz_loc);
    /* Hash table avoid duplication for local edges */

    if ((hash < 0) || hashtab[hash] != ~0)
      return;

    if ((i == 0) ||                                  /* local arc */
	(couple[0] >= vtxdist[myproc] && (couple[0] < vtxdist[myproc+1]))) {
      arcwgttab = rcvarcwgttab;
      arctab = rcvarctab;
      current = curserrcv;
    }
    else {
      arcwgttab = sndarcwgttab;
      arctab = sndarctab;
      current = cursersnd;
      proctab[*current] = give_proc(src, vtxdist, procnum, myproc);;
    }

    if (arcwgttab != NULL) {
      arcwgttab[*current].arc.src = couple[0];
      arcwgttab[*current].arc.tgt = couple[1];
      arcwgttab[*current].wgt = wgt;
    }
    else {
      arctab[*current].src = couple[0];
      arctab[*current].tgt = couple[1];
    }
    hashtab[hash] = *current;
    *current++;

    couple[0] = tgt;
    couple[1] = src;
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
    const indextype * const * vtxdist, indextype **xadj, indextype **adjncy,
    float **ewgts, const int * const *adjproc)
{
  /* Local variables */
  indextype vertsrc;
  int edgenum;
  int nnz_loc;
  int procnum, myproc;
  Zoltan_Arcs *sndarctab = NULL, *rcvarctab = NULL;
  Zoltan_Weighted_Arcs *sndarcwgttab = NULL, *rcvarcwgttab = NULL;
  int cursersnd, curserrcv;
  int *hashtab;
  int *proctab;
  int rcvsize;
  indextype prevsrc, prevtgt, currentvtx;
  int currentedge;
  int edge;


  MPI_Comm comm = zz->Communicator;
  ZOLTAN_COMM_OBJ *comm_plan;

  char *yo = "Zoltan_Symmetrize_Graph";
  int ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);

  procnum = zz->Num_Proc;
  myproc = zz->Proc;

  nnz_loc = *xadj[num_obj];

  hashtab = (int*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(int));
  if (edge_wgt_dim == 1){
    sndarcwgttab = (Zoltan_Weighted_Arcs*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(Zoltan_Weighted_Arcs));
    rcvarcwgttab = (Zoltan_Weighted_Arcs*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(Zoltan_Weighted_Arcs));
  }
  else {
    sndarctab = (Zoltan_Arcs*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(Zoltan_Arcs));
    rcvarctab = (Zoltan_Arcs*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(Zoltan_Arcs));
  }
  proctab = (int*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(int));

  if (hashtab == NULL ||
      (sndarcwgttab == NULL && sndarctab == NULL) || (rcvarcwgttab == NULL && rcvarctab == NULL)
      || proctab == NULL) {
    ZOLTAN_FREE(&proctab);
    ZOLTAN_FREE(&rcvarctab);
    ZOLTAN_FREE(&rcvarcwgttab);
    ZOLTAN_FREE(&sndarctab);
    ZOLTAN_FREE(&sndarcwgttab);
    ZOLTAN_FREE(&hashtab);
    ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory");
    return (ZOLTAN_MEMERR);
  }

  memset (hashtab, ~0, 2*nnz_loc*sizeof(int));

  /* Construct a locally symmetrized graph */
  for (vertsrc = 0, curserrcv=0, cursersnd =0 ; vertsrc < num_obj ; vertsrc ++) { /* For all vertices */
    for (edgenum = *xadj[vertsrc] ; edgenum < *xadj[vertsrc+1] ; ++edgenum) { /* For all neighbors */
      float wgt = 0;

      if (edge_wgt_dim == 1) {
	wgt = *ewgts[edgenum];
      }
      add_edge(vertsrc + *vtxdist[myproc], *adjncy[edgenum], wgt,
	       sndarctab, sndarcwgttab, proctab, rcvarctab, rcvarcwgttab,
	       &cursersnd, &cursersnd, hashtab, nnz_loc, *vtxdist, procnum, myproc);
    }
  }
  ZOLTAN_FREE(&hashtab);

  /* These arrays may be reconstructed */
  ZOLTAN_FREE(ewgts);
  ZOLTAN_FREE(adjncy);

  /* Communicate the arcs */
  Zoltan_Comm_Create(&comm_plan, cursersnd, proctab, comm, 6241984, &rcvsize);
  rcvsize += curserrcv;
  if (rcvsize >= 2*nnz_loc) { /* reception buffer is too small */
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

  *adjncy = (indextype*)ZOLTAN_MALLOC(rcvsize*sizeof(indextype));
  if (*adjncy == NULL) {
    ZOLTAN_FREE(&rcvarctab);
    ZOLTAN_FREE(&rcvarcwgttab);
    ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory (3)");
    return (ZOLTAN_MEMERR);
  }
  if (edge_wgt_dim == 1) {
    *ewgts = (float*)ZOLTAN_MALLOC(rcvsize*sizeof(float));
    if (*ewgts == NULL) {
      ZOLTAN_FREE(&rcvarctab);
      ZOLTAN_FREE(&rcvarcwgttab);
      ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory (4)");
      return (ZOLTAN_MEMERR);
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
  prevsrc = *vtxdist[myproc] ; prevtgt = -1;
  (*xadj)[0] = 0; currentedge = 0; currentvtx=0;
  for (edge = 0 ; edge < rcvsize ; ++ edge) {
    Zoltan_Arcs arc;

    if (edge_wgt_dim == 1)
      arc = rcvarctab[edge];
    else
      arc = rcvarcwgttab[edge].arc;

    if (arc.src != prevsrc) {
      do {
	*xadj[currentvtx++] = currentedge;
      } while (currentvtx + *vtxdist[myproc] < arc.src);
      prevsrc = arc.src;
    }
    else { /* If the same src vertex, check if it's the same edge */
      if (prevtgt == arc.tgt)
	continue;
    }

    prevtgt = arc.tgt;
    *adjncy[currentedge] = arc.tgt;
    if (edge_wgt_dim == 1)
      *ewgts[currentedge] = rcvarcwgttab[edge].wgt;

    currentedge++;
  }
  (*xadj)[num_obj] = currentedge;

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}

/**************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
