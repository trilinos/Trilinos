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
  int wgtoff;             /* Used to have a reference in wgttab */
} Zoltan_Arcs;


  /* Definition of the merge wgt operator : arguments are accumulation wgt, partial wgt, wgt dim */
  /* Return !0 when an error is found */
typedef int(*WgtFctPtr)(float*, float*, int);


/* This function compare if the wgt are the same for both arcs*/
#ifdef NEEDED
static int
wgtFctCmp(float* current, float* new, int dim)
{
  int i; int diff;

  for (i = 0, diff=0 ; (i <dim)&&(!diff) ; ++i) {
    diff |= (new[i] != current[i]);
  }
  return (diff);
}
#endif

/* This function adds the weights */
static int
wgtFctAdd(float* current, float* new, int dim)
{
  int i;
  for (i = 0 ; i <dim ; ++i) {
    current[i] += new[i];
  }
  return (0);
}

#ifdef NEEDED
/* This function chooses the maximum weight */ 
static int
wgtFctMax(float* current, float* new, int dim)
{
  int i;
  for (i = 0 ; i <dim ; ++i) {
    current[i] = MAX(current[i],new[i]);
  }
  return (0);
}
#endif

inline int
compar_arc(const Zoltan_Arcs* arc1, const Zoltan_Arcs* arc2)
{
  if (arc1->src == arc2->src)
    return (arc1->tgt - arc2->tgt);
  return (arc1->src - arc2->src);
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
add_edge(indextype src, indextype tgt, float *wgt, int wgt_dim, WgtFctPtr wgtfct,
	 Zoltan_Arcs* sndarctab, float *sndarcwgttab,
	 int* proctab, Zoltan_Arcs* rcvarctab, float *rcvarcwgttab,
	 int *cursersnd, int* curserrcv, int *hashtab, int nnz_loc,
	 const indextype *vtxdist, int procnum, int myproc)
{
    unsigned int hash;
    int couple[2];
    int i;

    couple[0] = src; couple [1] = tgt;

    for (i = 0; i < 2 ; ++i) {
      int hashnum;
      hash = Zoltan_Hash((ZOLTAN_ID_PTR)couple, 2, 2*nnz_loc);
      if (hash < 0) return;

      while((hashtab[hash] != 0) && ((rcvarctab[hashtab[hash]].src != couple[0])
	    || (rcvarctab[hashtab[hash]].tgt != couple[1]))){
	hash++;
	if (hash >= 2*nnz_loc)
	  hash =0;
      }

      if (hashtab[hash] != 0) { /* If already in the hash table, only check weights */
	wgtfct(&rcvarcwgttab[hashtab[hash]*wgt_dim], wgt, wgt_dim);
	return;
      }
      hashnum = hashtab[hash] = *curserrcv;
      (*curserrcv)++;
      rcvarctab[hashnum].src = couple[0];
      rcvarctab[hashnum].tgt = couple[1];
      memcpy (&rcvarcwgttab[hashnum*wgt_dim], wgt, sizeof(float)*wgt_dim);
      rcvarctab[hashnum].wgtoff = hashnum*wgt_dim; /* Usefull to not have to sort wgt array */

      if ((tgt >= vtxdist[myproc]) && (tgt < vtxdist[myproc+1])) {
	couple[0] = tgt; couple [1] = src;
      }
      else { /* Non local vertex */
	sndarctab[*cursersnd].src = tgt;
	sndarctab[*cursersnd].tgt = src;
	sndarctab[*cursersnd].wgtoff = 0;

	memcpy (&sndarcwgttab[(*cursersnd)*wgt_dim], wgt, sizeof(float)*wgt_dim);
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
    int obj_wgt_dim, int * edge_wgt_dim,
    indextype ** vtxdist, indextype **xadj, indextype **adjncy,
    float **ewgts, int **adjproc)
{
  /* Local variables */
  indextype vertsrc;
  int edgenum;
  int nnz_loc;
  int numProc, myproc, currentproc;
  Zoltan_Arcs *sndarctab = NULL, *rcvarctab = NULL;
  int cursersnd, curserrcv;
  int *hashtab;
  int *proctab;
  int rcvsize = 0;
  indextype prevsrc, prevtgt, currentvtx;
  int currentedge;
  int edge;
  WgtFctPtr edgeWgt = &wgtFctAdd;
  float *sndarcwgttab = NULL, *rcvarcwgttab = NULL;
  float baseWgt = 1;
  int msgtag = 8161982;
  int base_wgt_dim = *edge_wgt_dim;
  int new_wgt_dim;

  MPI_Comm comm = zz->Communicator;
  ZOLTAN_COMM_OBJ *comm_plan;

  char *yo = "Zoltan_Symmetrize_Graph";
  int ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);

  numProc = zz->Num_Proc;
  myproc = zz->Proc;

  nnz_loc = (*xadj)[num_obj];
  new_wgt_dim = (base_wgt_dim>0)?base_wgt_dim:1;
  *edge_wgt_dim = new_wgt_dim;
  hashtab = (int*)ZOLTAN_CALLOC(2*nnz_loc, sizeof(int));
  sndarctab = (Zoltan_Arcs*)ZOLTAN_MALLOC(nnz_loc*sizeof(Zoltan_Arcs));
  sndarcwgttab = (float*)ZOLTAN_MALLOC(nnz_loc*sizeof(float)*new_wgt_dim);
  rcvarctab = (Zoltan_Arcs*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(Zoltan_Arcs));
  rcvarcwgttab = (float*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(float)*new_wgt_dim);
  proctab = (int*)ZOLTAN_MALLOC(2*nnz_loc*sizeof(int));

  if ((nnz_loc >0 ) && (hashtab == NULL || sndarctab == NULL || rcvarctab == NULL
			|| proctab == NULL || sndarcwgttab == NULL || rcvarcwgttab == NULL)) {
    ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }


  /* Construct a locally symmetrized graph */
  for (vertsrc = 0, curserrcv=0, cursersnd =0 ; vertsrc < num_obj ; vertsrc ++) { /* For all vertices */
    for (edgenum = (*xadj)[vertsrc] ; edgenum < (*xadj)[vertsrc+1] ; ++edgenum) { /* For all neighbors */
      float * wgt = NULL;

      if (base_wgt_dim != 0)
	wgt = &((*ewgts)[edgenum*new_wgt_dim]);
      else wgt = &baseWgt;

      add_edge(vertsrc + (*vtxdist)[myproc], (*adjncy)[edgenum], wgt, new_wgt_dim, edgeWgt,
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
    rcvarctab = (Zoltan_Arcs*) ZOLTAN_REALLOC(rcvarctab, rcvsize*sizeof(Zoltan_Arcs));
    rcvarcwgttab = (float*) ZOLTAN_REALLOC(rcvarcwgttab, rcvsize*sizeof(float)*new_wgt_dim);
    if (rcvarctab == NULL) {
      Zoltan_Comm_Destroy (&comm_plan);
      ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory (2)");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  Zoltan_Comm_Do (comm_plan, msgtag--, (char*)sndarctab, sizeof(Zoltan_Arcs),
		  (char*)(rcvarctab + curserrcv));
  ZOLTAN_FREE(&sndarctab);

  Zoltan_Comm_Do_Post (comm_plan, msgtag, (char*)sndarcwgttab, sizeof(float)*new_wgt_dim,
		  (char*)(rcvarcwgttab + curserrcv));

  for (edge=curserrcv ; edge < rcvsize ; ++edge) { /* Update wgtindex for received data */
    rcvarctab[edge].wgtoff = edge*new_wgt_dim;
  }

  Zoltan_Comm_Do_Wait (comm_plan, msgtag, (char*)sndarcwgttab, sizeof(float)*new_wgt_dim,
		  (char*)(rcvarcwgttab + curserrcv));

  ZOLTAN_FREE(&sndarcwgttab);
  ZOLTAN_FREE(&proctab);
  Zoltan_Comm_Destroy (&comm_plan);

  if (rcvsize > 0) {
    *adjncy = (indextype*)ZOLTAN_MALLOC(rcvsize*sizeof(indextype));
    *adjproc = (indextype*)ZOLTAN_MALLOC(rcvsize*sizeof(indextype));
    *ewgts = (float*)ZOLTAN_MALLOC(rcvsize*sizeof(float)*new_wgt_dim);

    if ((*adjncy) == NULL || (*adjproc) == NULL ||(*ewgts) == NULL) {
      ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory (3)");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  /* Sort Edges */
  qsort ((void*)rcvarctab, rcvsize, sizeof(Zoltan_Arcs),
	 (int (*)(const void*,const void*))compar_arc);


  /* Now I only need to copy these informations in the edge array and to recompute
     the indirection array */
  prevsrc = (*vtxdist)[myproc] ; prevtgt = -1;
  (*xadj)[0] = 0; currentedge = 0; currentvtx=1;
  currentproc = 0;
  for (edge = 0 ; edge < rcvsize ; ++ edge) {
    Zoltan_Arcs *arc = &(rcvarctab[edge]);

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
      if (prevtgt == arc->tgt) {
	(*edgeWgt)(&((*ewgts)[(currentedge-1)*new_wgt_dim]), &rcvarcwgttab[arc->wgtoff], new_wgt_dim);
	continue;
      }
    }

    prevtgt = arc->tgt;
    (*adjncy)[currentedge] = arc->tgt;
    (*adjproc)[currentedge] = give_proc(arc->tgt, *vtxdist, numProc, &currentproc);
    memcpy(&(*ewgts)[currentedge*new_wgt_dim], &rcvarcwgttab[arc->wgtoff], new_wgt_dim*sizeof(float));

    currentedge++;
  }
  (*xadj)[num_obj] = currentedge;

  /* Realloc arrays that are probably too big to be memory more efficient */
  if (rcvsize > 0) {
    (*adjncy) = (indextype*) ZOLTAN_REALLOC((*adjncy), currentedge*sizeof(indextype));
    (*adjproc) = (int*)  ZOLTAN_REALLOC((*adjproc), currentedge*sizeof(int));
    (*ewgts) = (float*) ZOLTAN_REALLOC((*ewgts), currentedge*sizeof(float)*new_wgt_dim);

    if ((*adjncy) == NULL || (*adjproc) == NULL ||(*ewgts) == NULL) {
      ZOLTAN_PRINT_ERROR (myproc, yo, "Unable to allocate enough memory (4)");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

End:
  ZOLTAN_FREE(&proctab);
  ZOLTAN_FREE(&rcvarctab);
  ZOLTAN_FREE(&rcvarcwgttab);
  ZOLTAN_FREE(&sndarctab);
  ZOLTAN_FREE(&sndarcwgttab);
  ZOLTAN_FREE(&hashtab);


  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}

/**************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
