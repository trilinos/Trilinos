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

#include "phg.h"
#include <limits.h>

/****************************************************************************/
/* Routine to set function pointers corresponding to input-string options. */
int Zoltan_PHG_Set_Part_Options (ZZ *zz, PHGPartParams *hgp)
{
  char *yo = "Zoltan_PHG_Set_Part_Options";

  if (hgp->bal_tol < 1.0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_BALANCE_TOLERANCE.");
    return ZOLTAN_FATAL;
  }

  /* Set reduction method. */
  hgp->matching = hgp->matching_opt = NULL;
  if (!(Zoltan_PHG_Set_Matching_Fn (hgp))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_REDUCTION_METHOD.");
    return ZOLTAN_FATAL;
  }

  /* Set (serial) coarse partitioning method */
  /* May need parallel partitioning method later if reduction to 1 proc fails */
  if (!(hgp->CoarsePartition 
   = Zoltan_PHG_Set_CoarsePartition_Fn(hgp->coarsepartition_str))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_COARSE_PARTITIONING.");
      return ZOLTAN_FATAL;
  }

  /* Set refinement method. */
  if (!(hgp->Refinement = Zoltan_PHG_Set_Refinement_Fn(hgp->refinement_str))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_REFINEMENT.");
    return ZOLTAN_FATAL;
  }
  return ZOLTAN_OK;
}



/****************************************************************************/
/*  Main partitioning function for hypergraph partitioning. */
int Zoltan_PHG_HPart_Lib (
  ZZ *zz,               /* Zoltan data structure */
  PHGraph *hg,          /* Input hypergraph to be partitioned */
  int p,                /* Input:  number partitions to be generated */
  Partition part,       /* Output:  partition #s; aligned with vertex arrays. */
  PHGPartParams *hgp,   /* Input:  parameters for hgraph partitioning. */
  int level
)
{
  int  i, err = ZOLTAN_OK;
  char msg[128];
  char *yo = "Zoltan_PHG_HPart_Lib";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Check - The partition array must be already allocated */
  if (!part) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Output partition array is NULL.");
    return ZOLTAN_MEMERR;
  }
  
  /* Check - Is something wrong with the part number? */
  if (p <= 0) {
    sprintf(msg, "PART ERROR...p=%d is not a positive number!\n", p);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    return ZOLTAN_FATAL;
  }

  /* Check - the graph must be reduced no smaller than the number of parts */ 
  /* RTHRTH - reduction might be hg->redl / p for parallel */ 
  if (hg->redl < p)
    hg->redl = p;
  
  if (hgp->output_level >= PHG_DEBUG_PLOT)
    Zoltan_PHG_Plot(zz->Proc, hg->nVtx, p, hg->vindex, hg->vedge, NULL,
     "coarsening plot");

  if (hgp->output_level >= PHG_DEBUG_LIST) {
    printf("START %3d |V|=%6d |E|=%6d |I|=%6d %d/%s/%s-%s p=%d...\n",
     hg->info, hg->nVtx, hg->nEdge, hg->nNonZero, hg->redl, hgp->redm_str,
     hgp->coarsepartition_str, hgp->refinement_str, p);
    if (hgp->output_level > PHG_DEBUG_LIST) {
      err = Zoltan_PHG_Info(zz, hg);
      if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
        return err;
    }
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
  else if (hg->nVtx <= hg->redl || hg->nEdge == 0 || hgp->matching == NULL)  {
    /* fewer vertices than desired or no edges or no coarsening requested */
    err = Zoltan_PHG_CoarsePartition (zz, hg, p, part, hgp);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
      return err;
  }
  else {        /* normal multilevel situation */
    int *match = NULL, *LevelMap = NULL, *c_part = NULL, limit;
    PHGraph c_hg;

    /* Allocate and initialize Matching Array */
    if (!(match = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory for Matching array");
      return ZOLTAN_MEMERR;
    }
    for (i = 0; i < hg->nVtx; i++)
      match[i] = i;
          
    /* Calculate matching (packing or grouping) */
    limit = hg->nVtx - hg->redl;
    if (hgp->matching)  {
      err = Zoltan_PHG_Matching (zz, hg, match, hgp);
      if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
        ZOLTAN_FREE ((void**) &match);
        return err;
      }
    }

    /* Allocate and initialize LevelMap */
    if (!(LevelMap = (int*) ZOLTAN_CALLOC (hg->nVtx, sizeof(int)))) {
      ZOLTAN_FREE ((void**) &match);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory for LevelMap.");
      return ZOLTAN_MEMERR;
    }

    /* Construct coarse hypergraph and LevelMap */
    err = Zoltan_PHG_Coarsening (zz, hg, match, &c_hg, LevelMap);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &match, &LevelMap);
      goto End;
      }

    /* Check the consistency of the coarsening */
    /* RTHRTH - needs to be parallel */
    if (limit != c_hg.nVtx - hg->redl) {
      sprintf(msg, "limit %d is not %d-%d!\n", limit, c_hg.nVtx, hg->redl);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      Zoltan_Multifree (__FILE__, __LINE__, 2, &match, &LevelMap);
      err = ZOLTAN_FATAL;
      goto End;
    }
    else if (c_hg.nVtx < hg->redl) {
      sprintf(msg, "wanted coarsen to %d vertices, but reached %d vertices\n",
       hg->redl, c_hg.nVtx);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      Zoltan_Multifree (__FILE__, __LINE__, 2, &match, &LevelMap);
      err = ZOLTAN_FATAL;
      goto End;
    }

    /* heuristic: stop on diminishing returns */
    /* RTHRTH - modify this for parallel version (need global view) */
    if (c_hg.nVtx > 0.9 * hg->nVtx)
      hg->redl = c_hg.redl = c_hg.nVtx;

    ZOLTAN_FREE ((void**) &match);

    /* Allocate Partition of coarse graph */
    if (!(c_part = (int*) ZOLTAN_CALLOC (c_hg.nVtx, sizeof(int)))) {
      ZOLTAN_FREE ((void**) &LevelMap);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      err = ZOLTAN_MEMERR;
      goto End;
    }

    /* Recursively partition coarse hypergraph */
    err = Zoltan_PHG_HPart_Lib (zz, &c_hg, p, c_part, hgp, level);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &c_part, &LevelMap);
      goto End;
    }

    /* Project coarse partition to fine partition */
    for (i = 0; i < hg->nVtx; i++)
      part[i] = c_part[LevelMap[i]];

    /* Free coarse graph, coarse partition and LevelMap */
    Zoltan_PHG_HGraph_Free (&c_hg);
    Zoltan_Multifree (__FILE__, __LINE__, 2, &c_part, &LevelMap);
  }

  /* Locally refine partition */
  err = Zoltan_PHG_Refinement (zz, hg, p, part, hgp);

  /* print useful information (conditionally) */
  if (hgp->output_level > PHG_DEBUG_LIST) {
    err = Zoltan_PHG_HPart_Info (zz, hg, p, part, hgp);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
      goto End;
  }
  if (hgp->output_level >= PHG_DEBUG_LIST)
      printf("FINAL %3d |V|=%6d |E|=%6d |I|=%6d %d/%s/%s-%s p=%d bal=%.2f cutl=%.2f, %.2f\n",
             hg->info, hg->nVtx, hg->nEdge, hg->nNonZero, hg->redl, hgp->redm_str,
             hgp->coarsepartition_str, hgp->refinement_str, p,
             Zoltan_PHG_HPart_balance(zz, hg, p, part),
             Zoltan_PHG_hcut_size_total(&hgp->comm, hg, part, p),
             Zoltan_PHG_hcut_size_links(&hgp->comm, hg, part, p));
  
  if (hgp->output_level >= PHG_DEBUG_PLOT)
    Zoltan_PHG_Plot(zz->Proc, hg->nVtx, p, hg->vindex, hg->vedge, part,
     "partitioned plot");

End:
  ZOLTAN_TRACE_EXIT(zz, yo) ;
  return err;
}



/****************************************************************************/
/* Calculates the cutsize of a partition by summing the weight of all edges
   which span more than one part. Time O(|I|). */
double Zoltan_PHG_hcut_size_total (PHGComm *hgc, PHGraph *hg, Partition part, int p)
{
    int i, j, *netpart, *allparts;    
    double cut = 0.0, totalcut=0.0;
    char *yo = "Zoltan_PHG_hcut_size_total";
    
    if (!(netpart = (int*) ZOLTAN_CALLOC (hg->nEdge, sizeof(int)))) {
        ZOLTAN_PRINT_ERROR(hgc->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }
    
    if (!hgc->myProc_x)
        if (!(allparts = (int*) ZOLTAN_CALLOC (hgc->nProc_x*hg->nEdge, sizeof(int)))) {
            ZOLTAN_PRINT_ERROR(hgc->Proc, yo, "Insufficient memory.");
            return ZOLTAN_MEMERR;
        }

    for (i = 0; i < hg->nEdge; ++i) 
        if (hg->hindex[i]>=hg->hindex[i+1])
            netpart[i] = -1;
        else {
            j = hg->hindex[i];
            netpart[i] = part[hg->hvertex[j]];
            for (++j; j < hg->hindex[i+1]
                     && part[hg->hvertex[j]] == netpart[i]; ++j);
            if (j != hg->hindex[i+1])
                netpart[i] = -2;
    }

    MPI_Gather(netpart, hg->nEdge, MPI_INT, allparts, hg->nEdge, MPI_INT, 0, hgc->row_comm);
    ZOLTAN_FREE ((void**) &netpart);

    if (!hgc->myProc_x) {
        for (i = 0; i < hg->nEdge; ++i) {
            int p=-1;
            for (j = 0; j<hgc->nProc_x; ++j)
                if (allparts[j*hg->nEdge+i]==-2)
                    break;
                else if (allparts[j*hg->nEdge+i]>=0) {
                    if (p==-1)
                        p = allparts[j*hg->nEdge+i];
                    else if (p != allparts[j*hg->nEdge+i])
                        break;
                }            
            if (j<hgc->nProc_x)
                cut += (hg->ewgt ? hg->ewgt[i] : 1.0);
        }
        
        ZOLTAN_FREE ((void**) &allparts);
        MPI_Reduce(&cut, &totalcut, 1, MPI_DOUBLE, MPI_SUM, 0, hgc->col_comm);
    }
    MPI_Bcast(&totalcut, 1, MPI_DOUBLE, 0, hgc->Communicator);
    return totalcut;    
}



/****************************************************************************/
/* Calculates the cutsize of a partition. For each edge it calculates
   the number of parts it spans across. This value minus one is the
   cutsize of this edge and the total cutsize is the sum of the single
   cutsizes. Time O(|I|). */
double Zoltan_PHG_hcut_size_links (PHGComm *hgc, PHGraph *hg, Partition part, int p)
{
    int i, j, *cuts, *rescuts, *parts, nparts;
    double cut = 0.0, totalcut=0.0;
    char *yo = "Zoltan_PHG_hcut_size_links";
    
    if (!(cuts = (int*) ZOLTAN_CALLOC (p*hg->nEdge, sizeof(int)))) {
        ZOLTAN_PRINT_ERROR(hgc->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }   
    if (!hgc->myProc_x)
        if (!(rescuts = (int*) ZOLTAN_CALLOC (p*hg->nEdge, sizeof(int)))) {
            ZOLTAN_PRINT_ERROR(hgc->Proc, yo, "Insufficient memory.");
            return ZOLTAN_MEMERR;
        }
    
    for (i = 0; i < hg->nEdge; ++i) {
        parts = &cuts[i*p];
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j) 
            ++parts[part[hg->hvertex[j]]];
    }
    
    MPI_Reduce(cuts, rescuts, p*hg->nEdge, MPI_INT, MPI_SUM, 0, hgc->row_comm);
    ZOLTAN_FREE ((void**) &cuts);

    if (!hgc->myProc_x) {
        for (i = 0; i < hg->nEdge; ++i) {
            parts = &rescuts[i*p];
            for (j=nparts=0; j<p; ++j)
                if (parts[j])
                    ++nparts;
            cut += (nparts-1) * (hg->ewgt ? hg->ewgt[i] : 1.0);
        }        
        ZOLTAN_FREE ((void**) &rescuts);

        MPI_Reduce(&cut, &totalcut, 1, MPI_DOUBLE, MPI_SUM, 0, hgc->col_comm);
    }
    MPI_Bcast(&totalcut, 1, MPI_DOUBLE, 0, hgc->Communicator);
    return totalcut;
}



/****************************************************************************/
/* Output procedures to print the minimal and maximal values of an array */
static int hmin_max (ZZ *zz, int P, int *q, PHGPartParams *hgp)
{ 
  int i, values[3];

  if (P > 0) {
    values[0] = INT_MAX;
    values[1] = values[2] = 0;
    for (i = 0; i < P; i++) {
      values[2] += q[i];
      values[0] = MIN(values[0], q[i]);
      values[1] = MAX(values[1], q[i]);
    }
    if (hgp->output_level >= PHG_DEBUG_LIST)
      printf("%9d    %12.2f %9d    %9d\n", values[0], (double)(values[2]) / P,
       values[1],values[2]);
    if (hgp->output_level > PHG_DEBUG_LIST) {
      for (i = 0; i < P; i++)
        printf ("%d ", q[i]);
      printf("\n");
    }
  }
  return values[1];
}



static double hmin_max_float (ZZ *zz, int P, double *q, PHGPartParams *hgp)
{ 
  int i;
  double values[3];

  if (P > 0) {
    values[0] = INT_MAX;
    values[1] = values[2] = 0;
    for (i = 0; i < P; i++) {
      values[2] += q[i];
      values[0] = MIN(values[0], q[i]);
      values[1] = MAX(values[1], q[i]);
    }
    if (hgp->output_level >= PHG_DEBUG_LIST)
      printf("%12.2f %12.2f %12.2f %12.2f\n", values[0], values[2]/P,
       values[1], values[2]);
    if (hgp->output_level > PHG_DEBUG_LIST) {
      for (i = 0; i < P; i++)
        printf ("%.2f ",q[i]);
      printf("\n");
    }
  }
  return values[1];
}



/* Prints important values of the partition on screen */
int Zoltan_PHG_HPart_Info (
  ZZ *zz,
  PHGraph *hg,
  int p,
  Partition part,
  PHGPartParams *hgp
)
{
  int i, *size, max_size;
  char msg[128];
  char *yo = "Zoltan_PHG_HPart_Info";

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
  printf (" Size               :");
  max_size = hmin_max (zz, p, size, hgp);
  printf (" Balance            : %.3f\n", (double) max_size * p / hg->nVtx);
  ZOLTAN_FREE ((void **) &size);

  if (hg->vwgt) {
    double *size_w, max_size_w, tot_w = 0.0;
    if (!(size_w = (double*) ZOLTAN_CALLOC (p, sizeof(double)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
    }
    for (i = 0; i < hg->nVtx; i++) {
      size_w[part[i]] += (hg->vwgt[i]);
      tot_w += hg->vwgt[i];
    }
    printf (" Size w.            :");
    max_size_w = hmin_max_float (zz, p, size_w, hgp);
    printf (" Balance w.         : %.3f\n", max_size_w * p / tot_w);
    ZOLTAN_FREE ((void**) &size_w);
  }

  printf ("EDGE-based:\n");
  printf (" Cuts(total/links)  : %.3f %.3f\n",
   Zoltan_PHG_hcut_size_total(&hgp->comm, hg, part, p), Zoltan_PHG_hcut_size_links(&hgp->comm, hg, part, p));
  printf ("----------------------------------------------------------------\n");

  return ZOLTAN_OK;
}



/****************************************************************************/

double Zoltan_PHG_HPart_balance (
  ZZ *zz,
  PHGraph *hg,
  int p,
  Partition part
)
{
  int i;
  char *yo = "Zoltan_PHG_HPart_balance";
  double *size_w, max_size_w=0.0, tot_w = 0.0;

  if (!(size_w = (double*) ZOLTAN_CALLOC (p, sizeof(double)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i = 0; i < hg->nVtx; i++)
    size_w[part[i]] += (hg->vwgt ? hg->vwgt[i] : 1.0);

  for (i = 0; i < p; i++)
    max_size_w = MAX(max_size_w, size_w[i]);

  if (hg->vwgt)
    for (i = 0; i < p; i++)
      tot_w += size_w[i];
  else
    tot_w = (double) (hg->nVtx);

  ZOLTAN_FREE ((void**) &size_w);

  return (max_size_w * p / tot_w);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
