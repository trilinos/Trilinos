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

#include "hypergraph.h"
static float hcut_size_links (ZZ *zz, HGraph *hg, int p, Partition part);



/****************************************************************************/
/* Routine to set function pointers corresponding to input-string options. */
int Zoltan_HG_Set_Part_Options(ZZ *zz, HGPartParams *hgp)
{
  char *yo = "Zoltan_HG_Set_Part_Options";

  if (hgp->bal_tol < 1.0) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid HG_BALANCE_TOLERANCE.");
     return ZOLTAN_FATAL;
     }

  /* Set reduction method. */
  hgp->matching = hgp->matching_opt = NULL;
  hgp->packing  = hgp->packing_opt  = NULL;
  hgp->grouping = hgp->grouping_opt = NULL;

  if (!(Zoltan_HG_Set_Matching_Fn (hgp))
   && !(Zoltan_HG_Set_Packing_Fn  (hgp))
   && !(Zoltan_HG_Set_Grouping_Fn (hgp))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid HG_REDUCTION_METHOD.");
      return ZOLTAN_FATAL;
      }

  /* Set global partitioning method */
  if (!(hgp->global_part = Zoltan_HG_Set_Global_Part_Fn(hgp->global_str))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid HG_GLOBAL_PARTITIONING.");
     return ZOLTAN_FATAL;
     }

  /* Set local refinement method. */
  if (!(hgp->local_ref = Zoltan_HG_Set_Local_Ref_Fn(hgp->local_str))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid HG_LOCAL_REFINEMENT.");
     return ZOLTAN_FATAL;
     }

  return ZOLTAN_OK;
}

/****************************************************************************/



/*  Main partitioning function for hypergraph partitioning. */
int Zoltan_HG_HPart_Lib (
  ZZ *zz,              /* Zoltan data structure */
  HGraph *hg,          /* Input hypergraph to be partitioned */
  int p,               /* Input:  number partitions to be generated */
  Partition part,      /* Output:  partition #s; aligned with vertex arrays. */
  HGPartParams *hgp    /* Input:  parameters for hgraph partitioning. */
)
{
int  i, err = ZOLTAN_OK;
char msg[128];
char *yo = "Zoltan_HG_HPart_Lib";

  /* The partition array has to be allocated prior to this procedure */
  if (!part) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Output partition array is NULL.");
     return ZOLTAN_MEMERR;
     }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
     Zoltan_HG_Plot(zz->Proc, hg->nVtx, p, hg->vindex, hg->vedge, NULL,
      "coarsening plot");

  if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST) {
     printf("START %3d |V|=%6d |E|=%6d |I|=%6d %d/%s/%s-%s p=%d...\n",
      hg->info, hg->nVtx, hg->nEdge, hg->nInput, hgp->redl, hgp->redm_str,
      hgp->global_str, hgp->local_str, p);
     if (zz->Debug_Level > ZOLTAN_DEBUG_LIST) {
        err = Zoltan_HG_Info(zz, hg);
        if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
           return err;
        }
     }

  /* the graph will only be reduced to a size equal to the number of parts */
  if (hgp->redl < p)
      hgp->redl = p;

  /* Something wrong with the part number? */
  if (p <= 0) {
     sprintf(msg, "PART ERROR...p=%d is not a positive number!\n", p);
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
     return ZOLTAN_FATAL;
     }

  /* take care of all special cases first */
  if (p == 1) {            /* everything goes in the one partition */
     for (i =  0; i < hg->nVtx; i++)
        part[i] = 0;
     }
  else if (p >= hg->nVtx) { /* more partitions than vertices, trivial answer */
     for (i = 0; i < hg->nVtx; i++)
        part[i] = i;
     }
  else if (hg->nVtx <= hgp->redl || hg->nEdge == 0
   || (hgp->matching == NULL && hgp->packing == NULL && hgp->grouping == NULL)){
     /* fewer vertices than desired or no edges or no coarsening requested */
     err = Zoltan_HG_Global (zz, hg, p, part, hgp);
     if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
        return err;
     }
  else {        /* normal multilevel situation */
     int *pack = NULL, *LevelMap = NULL, *c_part = NULL, limit;
     HGraph c_hg;

     /* Allocate Packing Array (used for matching, packing & grouping) */
     if (!(pack = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory for pack array.");
        return ZOLTAN_MEMERR;
        }
     for (i = 0; i < hg->nVtx; i++)
        pack[i] = i;

     /* Calculate one of Packing, Grouping or Matching */
     limit = hg->nVtx - hgp->redl;
     if      (hgp->packing)  err= Zoltan_HG_Packing (zz, hg, pack, hgp, &limit);
     else if (hgp->grouping) err= Zoltan_HG_Grouping(zz, hg, pack, hgp, &limit);
     else if (hgp->matching) err= Zoltan_HG_Matching(zz, hg, pack, hgp, &limit);

     if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
        ZOLTAN_FREE((void**) &pack);
        return err;
        }

     /* Allocate and initialize LevelMap */
     if (!(LevelMap = (int*) ZOLTAN_CALLOC (hg->nVtx, sizeof(int)))) {
        ZOLTAN_FREE ((void**) &pack);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory for LevelMap.");
        return ZOLTAN_MEMERR;
        }

     /* Construct coarse hypergraph and LevelMap */
     err = Zoltan_HG_Coarsening(zz, hg, pack, &c_hg, LevelMap);
     if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
        Zoltan_Multifree (__FILE__, __LINE__, 2, &pack, &LevelMap);
        return err;
        }

     /* Check the consistency of the coarsening */
     if (limit != c_hg.nVtx - hgp->redl) {
        sprintf(msg, "limit %d is not %d-%d!\n", limit, c_hg.nVtx, hgp->redl);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
        Zoltan_Multifree (__FILE__, __LINE__, 2, &pack, &LevelMap);
        return ZOLTAN_FATAL;
        }
     else if (c_hg.nVtx < hgp->redl) {
        sprintf(msg, "wanted coarsen to %d vertices, but reached %d vertices\n",
         hgp->redl, c_hg.nVtx);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
        Zoltan_Multifree (__FILE__, __LINE__, 2, &pack, &LevelMap);
        return ZOLTAN_FATAL;
        }

     /* heuristic: stop on diminishing returns */
     if (c_hg.nVtx > (0.9 * hg->nVtx))
        hgp->redl = c_hg.nVtx;

     ZOLTAN_FREE ((void**) &pack);

     /* Allocate Partition of coarse graph */
     if (!(c_part = (int*) ZOLTAN_CALLOC (c_hg.nVtx, sizeof(int)))) {
        ZOLTAN_FREE ((void**) &LevelMap);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
        }

     /* Recursively partition coarse hypergraph */
     err = Zoltan_HG_HPart_Lib (zz, &c_hg, p, c_part, hgp);
     if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
        Zoltan_Multifree (__FILE__, __LINE__, 2, &c_part, &LevelMap);
        return err;
        }

     /* Project coarse partition to fine partition */
     for (i = 0; i < hg->nVtx; i++)
        part[i] = c_part[LevelMap[i]];

     /* Free coarse graph, coarse partition and LevelMap */
     Zoltan_HG_HGraph_Free (&c_hg);
     Zoltan_Multifree (__FILE__, __LINE__, 2, &c_part, &LevelMap);
     }

  /* Locally refine partition */
  err = Zoltan_HG_Local (zz, hg, p, part, hgp);
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
     return err;

  /* print useful information (conditionally) */
  if (zz->Debug_Level > ZOLTAN_DEBUG_LIST) {
     err = Zoltan_HG_HPart_Info (zz, hg, p, part);
     if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
        return err;
     }
  if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST)
    printf("FINAL %3d |V|=%6d |E|=%6d |I|=%6d %d/%s/%s-%s p=%d cutl=%.2f\n",
     hg->info, hg->nVtx, hg->nEdge, hg->nInput, hgp->redl, hgp->redm_str,
     hgp->global_str, hgp->local_str, p, hcut_size_links(zz, hg, p, part));

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
     Zoltan_HG_Plot(zz->Proc, hg->nVtx, p, hg->vindex, hg->vedge, part,
      "partitioned plot");
  return err;
}

/****************************************************************************/



/* Calculates the cutsize of a partition by summing the weight of all edges
   which span more than one part. Time O(|I|). */
float Zoltan_HG_hcut_size_total (HGraph *hg, Partition part)
{
int i, j, hpart;
float cut = 0.0;

  for (i = 0; i < hg->nEdge; i++) {
     hpart = part[hg->hvertex[hg->hindex[i]]];
     for (j = hg->hindex[i] + 1; j < hg->hindex[i+1]
      && part[hg->hvertex[j]] == hpart; j++);
         if (j != hg->hindex[i+1])
            cut += (hg->ewgt ? hg->ewgt[i] : 1.0);
     }
  return cut;
}

/****************************************************************************/



/* Calculates the cutsize of a partition. For each edge it calculates
   the number of parts it spans across. This value minus one is the
   cutsize of this edge and the total cutsize is the sum of the single
   cutsizes. Time O(|I|). */
static float hcut_size_links (ZZ *zz, HGraph *hg, int p, Partition part)
{
int i, j, *parts, nparts;
float cut = 0.0;
char *yo = "hcut_size_links";

  if (!(parts = (int*) ZOLTAN_CALLOC (p,sizeof(int)))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }

  for (i = 0; i < hg->nEdge; i++) {
     nparts = 0;
     for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++) {
        if (parts[part[hg->hvertex[j]]] < i+1)
           nparts++;
        parts[part[hg->hvertex[j]]] = i + 1;
        }
     cut += (nparts-1) * (hg->ewgt ? hg->ewgt[i] : 1.0);
     }
  ZOLTAN_FREE ((void**) &parts);
  return cut;
}

/****************************************************************************/



/* Output procedures to print the minimal and maximal values of an array */
static int hmin_max (ZZ *zz, int P, int *q)
{ int i, values[3];

  if (P > 0) {
     values[0] = INT_MAX;
     values[1] = values[2] = 0;
     for (i = 0; i < P; i++) {
        values[2] += q[i];
        values[0] = MIN(values[0], q[i]);
        values[1] = MAX(values[1], q[i]);
        }
     if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST)
        printf("%9d    %12.2f %9d    %9d\n", values[0], (float)(values[2]) / P,
         values[1],values[2]);
     if (zz->Debug_Level > ZOLTAN_DEBUG_LIST) {
        for (i = 0; i < P; i++)
           printf ("%d ", q[i]);
        printf("\n");
        }
     }
  return values[1];
}



static float hmin_max_float (ZZ *zz, int P, float *q)
{ int i;
  float values[3];

  if (P > 0) {
     values[0] = INT_MAX;
     values[1] = values[2] = 0;
     for (i = 0; i < P; i++) {
        values[2] += q[i];
        values[0] = MIN(values[0], q[i]);
        values[1] = MAX(values[1], q[i]);
        }
     if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST)
        printf("%12.2f %12.2f %12.2f %12.2f\n", values[0], values[2]/P,
         values[1], values[2]);
     if (zz->Debug_Level > ZOLTAN_DEBUG_LIST) {
        for (i = 0; i < P; i++)
            printf ("%.2f ",q[i]);
         printf("\n");
        }
     }
  return values[1];
}



/* Prints important values of the partition on screen */
int Zoltan_HG_HPart_Info (ZZ *zz, HGraph *hg, int p, Partition part)
{
int i, *size, max_size;
char msg[128];
char *yo = "Zoltan_HG_HPart_Info";

  if (zz->Debug_Level < ZOLTAN_DEBUG_LIST)
     return ZOLTAN_OK;

  puts("---------- Partition Information (min/ave/max/tot) ----------------");
  printf ("VERTEX-based:\n");
  if (!(size = (int*) ZOLTAN_CALLOC (p,sizeof(int)))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i = 0; i < hg->nVtx; i++) {
     if (part[i] < 0 || part[i] >= p) {
        sprintf(msg, "PART ERROR...vertex %d has wrong part number %d\n",i,
         part[i]);
        ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
        return ZOLTAN_WARN;
        }
     size[part[i]]++;
     }
  printf (" Size               :"); max_size = hmin_max (zz, p, size);
  printf (" Balance            : %.3f\n", (float) max_size * p / hg->nVtx);
  ZOLTAN_FREE ((void **) &size);

  if (hg->vwgt) {
     float *size_w, max_size_w, tot_w = 0.0;
     if (!(size_w = (float*) ZOLTAN_CALLOC (p, sizeof(float)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
        }
     for (i = 0; i < hg->nVtx; i++) {
        size_w[part[i]] += (hg->vwgt[i]);
        tot_w += hg->vwgt[i];
        }
     printf (" Size w.            :");
     max_size_w = hmin_max_float (zz, p, size_w);
     printf (" Balance w.         : %.3f\n", max_size_w * p / tot_w);
     ZOLTAN_FREE ((void**) &size_w);
     }

  printf ("EDGE-based:\n");
  printf (" Cuts(total/links)  : %.3f %.3f\n",
   Zoltan_HG_hcut_size_total(hg,part),hcut_size_links(zz,hg,p,part));
  printf ("----------------------------------------------------------------\n");

  return ZOLTAN_OK;
}

/****************************************************************************/



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
