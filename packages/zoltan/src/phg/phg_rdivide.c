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

#include "phg.h"
#include "phg_distrib.h"

#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}

/* #define SPLIT_PROCESSORS */

static int split_hypergraph (int *pins[2], HGraph*, HGraph*, Partition, int, ZZ*);


/* recursively divides problem into 2 parts until all p found */
int Zoltan_PHG_rdivide(int lo, int hi, Partition final, ZZ *zz, HGraph *hg,
                       PHGPartParams *hgp, int level)
{
  char *yo = "Zoltan_PHG_rdivide";
  int i, j, mid, ierr=ZOLTAN_OK, *pins[2] = {NULL,NULL}, *lpins[2] = {NULL,NULL};
  Partition part=NULL;
  HGraph *left=NULL, *right=NULL;
  PHGComm *hgc = hg->comm;
  float tgpartsize[2]={0.0,0.0};    /* Target partition sizes; dimension is 2 
                                     because we are doing bisection */

  hg->redl = hgp->redl;
  
  /* only one part remaining, record results and exit */
  if (lo == hi) {
    for (i = 0; i < hg->nVtx; ++i)
      final [hg->vmap[i]] = lo;
    return ZOLTAN_OK;
  }

  if (hg->nVtx && !(part = (Partition) ZOLTAN_MALLOC (hg->nVtx * sizeof (int))))
      MEMORY_ERROR;

  /* bipartition current hypergraph with appropriate split ratio */
  mid = (lo+hi)/2;
  tgpartsize[0] = tgpartsize[1] = 0.;
  for (i = lo; i <= mid; i++)  tgpartsize[0] += hgp->part_sizes[i];
  for (i = lo; i <= hi;  i++)  tgpartsize[1] += hgp->part_sizes[i];
  hg->ratio = (double) tgpartsize[0] / (double) tgpartsize[1];
  tgpartsize[0] = hg->ratio;
  tgpartsize[1] = 1. - tgpartsize[0];

  ierr = Zoltan_PHG_Partition (zz, hg, 2, tgpartsize, part, hgp, level);
  if (ierr != ZOLTAN_OK)
      goto End;

  uprintf(hgc, "Rdivide(%d, %d): %.1lf\n", lo, hi, Zoltan_PHG_hcut_size_links(hgc, hg, part, 2));
    
  /* if only two parts total, record results and exit */
  if (lo + 1 == hi)  {
    for (i = 0; i < hg->nVtx; ++i)
      final [hg->vmap[i]] = ((part[i] == 0) ? lo : hi);
    ZOLTAN_FREE (&part);
    return ZOLTAN_OK;
  }

  if (hg->nEdge && (!(pins[0] = (int*) ZOLTAN_CALLOC (2 * hg->nEdge, sizeof(int)))
   || !(lpins[0] = (int*) ZOLTAN_CALLOC (2 * hg->nEdge, sizeof(int)))))
      MEMORY_ERROR;
  if (pins[0] && lpins[0]) {
      pins[1]  = &( pins[0][hg->nEdge]);
      lpins[1] = &(lpins[0][hg->nEdge]);
  }
     
  /* Initial calculation of the local pin distribution  (sigma in UVC's papers)  */
  for (i = 0; i < hg->nEdge; ++i)
      for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j)
        ++(lpins[part[hg->hvertex[j]]][i]);
        
  /* now compute global pin distribution */
  MPI_Allreduce(lpins[0], pins[0], 2*hg->nEdge, MPI_INT, MPI_SUM, hgc->row_comm);
  ZOLTAN_FREE (&lpins[0]);                        /* we don't need lpins */
    
  /* recursively divide in two parts and repartition hypergraph */
  if (mid>lo) { /* only split if we really need it */
      if (!(left = (HGraph*) ZOLTAN_MALLOC (sizeof (HGraph))))
          MEMORY_ERROR;
      Zoltan_HG_HGraph_Init (left);
      
      ierr = split_hypergraph (pins, hg, left, part, 0, zz);
      if (ierr != ZOLTAN_OK) 
          goto End;
  } else {
      for (i = 0; i < hg->nVtx; ++i)
          if (part[i]==0)
              final [hg->vmap[i]] = lo;
  }

  if (hi>mid+1) { /* only split if we need it */
      if (!(right = (HGraph*) ZOLTAN_MALLOC (sizeof (HGraph))))
          MEMORY_ERROR;
      Zoltan_HG_HGraph_Init (right);
      ierr = split_hypergraph (pins, hg, right, part, 1, zz);
  
      ZOLTAN_FREE (&pins[0]); /* we don't need pins */      
      if (ierr != ZOLTAN_OK)
          goto End;
  } else {
      ZOLTAN_FREE (&pins[0]); /* we don't need pins */
      for (i = 0; i < hg->nVtx; ++i)
          if (part[i]==1)
              final [hg->vmap[i]] = hi;
  }

#ifdef SPLIT_PROCESSORS
  if (hgc->nProc>1 && left && right) {
      PHGComm  leftcomm, rightcomm;
      HGraph  newleft, newright;
      int     *leftvmap, *rightvmap, mid;

      Zoltan_HG_HGraph_Init (&newleft);
      Zoltan_HG_HGraph_Init (&newright);

      /* redistribute left and right parts */
      mid = (int)((float) hgc->nProc * (float) left->nPins / (float) hg->nPins);
      Zoltan_PHG_Redistribute(zz, left,
                              0, mid, 
                              &leftcomm, 
                              &newleft,
                              &leftvmap);
      Zoltan_PHG_Redistribute(zz, right,
                              mid+1, hgc->nProc-1, 
                              &rightcomm, 
                              &newright,
                              &rightvmap);
      
      if (hgc->myProc<=mid) {/* I'm on the left part so I should partition newleft */
          ierr = Zoltan_PHG_rdivide (lo, mid, final, zz, &newleft, hgp, level+1);
          Zoltan_HG_HGraph_Free (&newleft);
      } else { /* I'm on the right part so I should partition newright */
          ierr |= Zoltan_PHG_rdivide (mid+1, hi, final, zz, &newright, hgp, level+1);
          Zoltan_HG_HGraph_Free (&newright);
      }
  } else {
#endif
      if (left) 
          ierr = Zoltan_PHG_rdivide (lo, mid, final, zz, left, hgp, level+1);
      if (right) 
          ierr |= Zoltan_PHG_rdivide (mid+1, hi, final, zz, right, hgp, level+1);
#ifdef SPLIT_PROCESSORS      
  }
#endif
  
 End:
  if (left)
      Zoltan_HG_HGraph_Free (left);
  if (right)
      Zoltan_HG_HGraph_Free (right);
  Zoltan_Multifree (__FILE__, __LINE__, 5, &pins[0], &lpins[0], &part, &left, &right);

  return ierr;
}



/* recursively divides problem into 2 parts until all p found */
int Zoltan_PHG_rdivide_NoProcSplit(int lo, int hi, Partition final, ZZ *zz, HGraph *hg,
                       PHGPartParams *hgp, int level)
{
  char *yo = "Zoltan_PHG_rdivide";
  int i, j, mid, ierr=ZOLTAN_OK, *pins[2] = {NULL,NULL}, *lpins[2] = {NULL,NULL};
  Partition part=NULL;
  HGraph *left=NULL, *right=NULL;
  PHGComm *hgc = hg->comm;
  float tgpartsize[2]={0.0,0.0};    /* Target partition sizes; dimension is 2 
                                     because we are doing bisection */

  hg->redl = hgp->redl;
  
  /* only one part remaining, record results and exit */
  if (lo == hi) {
    for (i = 0; i < hg->nVtx; ++i)
      final [hg->vmap[i]] = lo;
    return ZOLTAN_OK;
  }

  if (hg->nVtx && !(part = (Partition) ZOLTAN_MALLOC (hg->nVtx * sizeof (int))))
      MEMORY_ERROR;

  /* bipartition current hypergraph with appropriate split ratio */
  mid = (lo+hi)/2;
  tgpartsize[0] = tgpartsize[1] = 0.;
  for (i = lo; i <= mid; i++)  tgpartsize[0] += hgp->part_sizes[i];
  for (i = lo; i <= hi;  i++)  tgpartsize[1] += hgp->part_sizes[i];
  hg->ratio = (double) tgpartsize[0] / (double) tgpartsize[1];
  tgpartsize[0] = hg->ratio;
  tgpartsize[1] = 1. - tgpartsize[0];

  ierr = Zoltan_PHG_Partition (zz, hg, 2, tgpartsize, part, hgp, level);
  if (ierr != ZOLTAN_OK)
      goto End;

  uprintf(hgc, "Rdivide(%d, %d): %.1lf\n", lo, hi, Zoltan_PHG_hcut_size_links(hgc, hg, part, 2));
    
  /* if only two parts total, record results and exit */
  if (lo + 1 == hi)  {
    for (i = 0; i < hg->nVtx; ++i)
      final [hg->vmap[i]] = ((part[i] == 0) ? lo : hi);
    ZOLTAN_FREE (&part);
    return ZOLTAN_OK;
  }

  if (hg->nEdge && (!(pins[0] = (int*) ZOLTAN_CALLOC (2 * hg->nEdge, sizeof(int)))
   || !(lpins[0] = (int*) ZOLTAN_CALLOC (2 * hg->nEdge, sizeof(int)))))
      MEMORY_ERROR;
  if (pins[0] && lpins[0]) {
      pins[1]  = &( pins[0][hg->nEdge]);
      lpins[1] = &(lpins[0][hg->nEdge]);
  }
     
  /* Initial calculation of the local pin distribution  (sigma in UVC's papers)  */
  for (i = 0; i < hg->nEdge; ++i)
      for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j)
        ++(lpins[part[hg->hvertex[j]]][i]);
        
  /* now compute global pin distribution */
  MPI_Allreduce(lpins[0], pins[0], 2*hg->nEdge, MPI_INT, MPI_SUM, hgc->row_comm);
  ZOLTAN_FREE (&lpins[0]);                        /* we don't need lpins */
    
  /* recursively divide in two parts and repartition hypergraph */
  if (mid>lo) { /* only split if we really need it */
      if (!(left = (HGraph*) ZOLTAN_MALLOC (sizeof (HGraph))))
          MEMORY_ERROR;

      Zoltan_HG_HGraph_Init (left);
      
      ierr = split_hypergraph (pins, hg, left, part, 0, zz);
      if (ierr != ZOLTAN_OK) 
          goto End;

      ierr = Zoltan_PHG_rdivide (lo, mid, final, zz, left, hgp, level+1);
      Zoltan_HG_HGraph_Free (left);
      if (ierr != ZOLTAN_OK) 
          goto End;
  } else {
      for (i = 0; i < hg->nVtx; ++i)
          if (part[i]==0)
              final [hg->vmap[i]] = lo;
  }

  if (hi>mid+1) { /* only split if we need it */
      if (!(right = (HGraph*) ZOLTAN_MALLOC (sizeof (HGraph))))
          MEMORY_ERROR;
      ierr = split_hypergraph (pins, hg, right, part, 1, zz);
  
      ZOLTAN_FREE (&pins[0]); /* we don't need pins */
      
      if (ierr != ZOLTAN_OK)
          goto End;

      ierr = Zoltan_PHG_rdivide (mid+1, hi, final, zz, right, hgp, level+1);
      Zoltan_HG_HGraph_Free (right);
  } else {
      ZOLTAN_FREE (&pins[0]); /* we don't need pins */
      for (i = 0; i < hg->nVtx; ++i)
          if (part[i]==1)
              final [hg->vmap[i]] = hi;
  }
      /* remove alloc'ed structs */
 End:
  Zoltan_Multifree (__FILE__, __LINE__, 5, &pins[0], &lpins[0], &part, &left, &right);

  return ierr;
}


static int split_hypergraph (int *pins[2], HGraph *old, HGraph *new, Partition part,
                             int partid, ZZ *zz)
{
  int *tmap = NULL;  /* temporary array mapping from old HGraph info to new */
  int edge, i, ierr=ZOLTAN_OK;  
  PHGComm *hgc = old->comm;
  char *yo = "split_hypergraph";

  new->comm = old->comm;
  new->info               = old->info;
  new->VtxWeightDim       = old->VtxWeightDim;
  new->EdgeWeightDim      = old->EdgeWeightDim;
  
  /* allocate memory for dynamic arrays in new HGraph and for tmap array */
  if (old->nVtx && (tmap = (int*) ZOLTAN_MALLOC (old->nVtx * sizeof (int)))==NULL)
      MEMORY_ERROR;
  
  /* fill in tmap array, -1 for ignored vertices, otherwise nonnegative int */
  new->nVtx = 0;
  for (i = 0; i < old->nVtx; i++)
      tmap[i] = (part[i] == partid) ? new->nVtx++ : -1; 

  /* save vertex and edge weights if they exist */
  if (old->vwgt && new->VtxWeightDim)
    new->vwgt=(float*)ZOLTAN_MALLOC(new->nVtx*sizeof(float)*new->VtxWeightDim);
  if (new->nVtx && (new->vmap = (int*) ZOLTAN_MALLOC (new->nVtx * sizeof (int)))==NULL)
      MEMORY_ERROR;
  
  for (i = 0; i < old->nVtx; i++) {
      int v=tmap[i];
      if (v!=-1) {
          new->vmap[v] = old->vmap[i];
          if (new->vwgt)
              memcpy(&new->vwgt[v*new->VtxWeightDim], &old->vwgt[i*new->VtxWeightDim],
                     new->VtxWeightDim * sizeof(float));
      }
  }
  
    
  /* fill in hindex and hvertex arrays in new HGraph */
  new->nEdge = 0;
  new->nPins = 0;
  for (edge = 0; edge < old->nEdge; ++edge)
      if (pins[partid][edge] > 1) {
          ++new->nEdge;
          new->nPins += pins[partid][edge];
      }


  /* continue allocating memory for dynamic arrays in new HGraph */
  if (new->nEdge && (new->hindex  = (int*) ZOLTAN_MALLOC ((new->nEdge+1) * sizeof (int)))==NULL)
      MEMORY_ERROR;
  if (new->nPins && (new->hvertex = (int*) ZOLTAN_MALLOC (new->nPins * sizeof (int)))==NULL)
      MEMORY_ERROR;
  if (old->ewgt && new->EdgeWeightDim && new->nEdge)
      if ((new->ewgt=(float*)ZOLTAN_MALLOC(new->nEdge*sizeof(float)*new->EdgeWeightDim))==NULL)
          MEMORY_ERROR;
  
  new->nEdge = 0;
  new->nPins = 0;
  for (edge = 0; edge < old->nEdge; ++edge)
    if (pins[partid][edge] > 1) { /* edge has at least two vertices in partition:
                                        we are skipping size 1 nets */
      new->hindex[new->nEdge] = new->nPins;
      for (i = old->hindex[edge]; i < old->hindex[edge+1]; ++i)
        if (tmap [old->hvertex[i]] >= 0)  {
          new->hvertex[new->nPins] = tmap[old->hvertex[i]];
          new->nPins++;  
        }
        if (new->ewgt)
            memcpy(&new->ewgt[new->nEdge*new->VtxWeightDim], &old->vwgt[edge*new->VtxWeightDim],
                   new->EdgeWeightDim * sizeof(float));
        ++new->nEdge;
    }
  new->hindex[new->nEdge] = new->nPins;

  /* We need to compute dist_x, dist_y */
  if (!(new->dist_x = (int *) ZOLTAN_CALLOC((hgc->nProc_x+1), sizeof(int)))
	 || !(new->dist_y = (int *) ZOLTAN_CALLOC((hgc->nProc_y+1), sizeof(int))))
      MEMORY_ERROR;

  MPI_Scan(&new->nVtx, new->dist_x, 1, MPI_INT, MPI_SUM, hgc->row_comm);
  MPI_Allgather(new->dist_x, 1, MPI_INT, &(new->dist_x[1]), 1, MPI_INT, hgc->row_comm);
  new->dist_x[0] = 0;
  
  MPI_Scan(&new->nEdge, new->dist_y, 1, MPI_INT, MPI_SUM, hgc->col_comm);
  MPI_Allgather(new->dist_y, 1, MPI_INT, &(new->dist_y[1]), 1, MPI_INT, hgc->col_comm);
  new->dist_y[0] = 0;
    
  /* shrink hindex, hvertex arrays to correct size & determine vindex, vedge */
  new->hvertex = (int*)ZOLTAN_REALLOC(new->hvertex, sizeof(int) * new->nPins);
  new->hindex  = (int*)ZOLTAN_REALLOC(new->hindex,  sizeof(int) *(new->nEdge+1));
  
  Zoltan_HG_Create_Mirror (zz, new);
 End:
  ZOLTAN_FREE (&tmap);
  return ierr;
}
