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
  int ierr = ZOLTAN_OK;

  if (hgp->bal_tol < 1.0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_BALANCE_TOLERANCE.");
    return ZOLTAN_FATAL;
  }

  /* Set reduction method. */
  hgp->matching = NULL;
  if (!(Zoltan_PHG_Set_Matching_Fn (hgp))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_REDUCTION_METHOD.");
    return ZOLTAN_FATAL;
  }

  /* Set (serial) coarse partitioning method */
  /* May need parallel partitioning method later if reduction to 1 proc fails */
  
  hgp->CoarsePartition = Zoltan_PHG_Set_CoarsePartition_Fn(hgp, &ierr);
  if (ierr != ZOLTAN_OK) {
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


typedef struct _tagVCycle {
    HGraph             *hg;         /* for finer==NULL, hg and part contains   */
    Partition           part;       /* original hg and part, don't delete them */  
    int                *LevelMap;   /* necessary to uncoarsen */
    int                 LevelCnt;   /* count of external vertices to uncoarsen */
    int                *LevelData;  /* buffer for external vertex information */
    struct _tagVCycle  *finer; 
} VCycle; 

static int allocVCycle(VCycle *v)
{
    int partalloc = 0;
    if (!v->hg)
        return ZOLTAN_OK;
    if (v->hg->nVtx) {
        if (!v->part) {
            if (!(v->part = (int*) ZOLTAN_CALLOC (v->hg->nVtx, sizeof(int))))
                return ZOLTAN_MEMERR;
            partalloc = 1;
        }
        if (!v->LevelMap && 
            !(v->LevelMap = (int*) ZOLTAN_CALLOC (v->hg->nVtx, sizeof(int)))) {
            if (partalloc)
                ZOLTAN_FREE ((void**) &v->part);
            return ZOLTAN_MEMERR;
        }
    }
    return ZOLTAN_OK;
}


static VCycle *newVCycle(HGraph *hg, Partition part, int *levmap, VCycle *finer)
{
    VCycle *vcycle=NULL;
    HGraph *nhg = NULL;
    
    if (!(vcycle = (VCycle*) ZOLTAN_MALLOC (sizeof(VCycle)))) 
        return NULL;
    vcycle->finer = finer;
    if (hg)
        vcycle->hg = hg;
    else {
        if (!(nhg = vcycle->hg = (HGraph *) ZOLTAN_MALLOC (sizeof(HGraph)))) {
            ZOLTAN_FREE ((void**) &vcycle);
            return NULL;
        }
    }
    vcycle->part = part;
    vcycle->LevelMap = levmap;
    if (hg)
        if (allocVCycle(vcycle)) {
            ZOLTAN_FREE ((void**) &vcycle);
            if (nhg)
                ZOLTAN_FREE ((void**) &nhg);
            vcycle = NULL;
        }
    return vcycle;
}

/****************************************************************************/
/*  Main partitioning function for hypergraph partitioning. */
int Zoltan_PHG_Partition (
    ZZ *zz,               /* Zoltan data structure */
    HGraph *hg,           /* Input hypergraph to be partitioned */
    int p,                /* Input:  number partitions to be generated */
    float *part_sizes,    /* Input:  array of length p containing percentages
                             of work to be assigned to each partition */
    Partition Parts,      /* Output:  partition #s; aligned with vtx arrays. */
    PHGPartParams *hgp,   /* Input:  parameters for hgraph partitioning. */
    int level
)
{
    VCycle  *vcycle=NULL, *del=NULL;
    int  i, err = ZOLTAN_OK, prevVcnt=2*hg->dist_x[hg->comm->nProc_x];
    char *yo = "Zoltan_PHG_Partition";
    
    ZOLTAN_TRACE_ENTER(zz, yo);
    
    if (!(vcycle = newVCycle(hg, Parts, NULL, NULL))) {
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "VCycle is NULL.");
        return ZOLTAN_MEMERR;
    }
    
    /****** Coarsening ******/    
    while ((hg->dist_x[hg->comm->nProc_x] > hg->redl)
           && (hg->dist_x[hg->comm->nProc_x] < 0.9 * prevVcnt)
           && (hg->nEdge != 0) && hgp->matching) {
        int *match = NULL;
        VCycle *coarser=NULL;
        
        prevVcnt=hg->dist_x[hg->comm->nProc_x];
        
        if (hgp->output_level >= PHG_DEBUG_LIST) {
            uprintf(hg->comm, "START %3d |V|=%6d |E|=%6d |I|=%6d %d/%s/%s-%s p=%d...\n",
                    hg->info, hg->nVtx, hg->nEdge, hg->nPins, hg->redl, hgp->redm_str,
                    hgp->coarsepartition_str, hgp->refinement_str, p);
            if (hgp->output_level > PHG_DEBUG_LIST) {
                err = Zoltan_HG_Info(zz, hg);
                if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
                    goto End;
            }
        }
        if (hgp->output_level >= PHG_DEBUG_PLOT)
            Zoltan_PHG_Plot(zz->Proc, hg->nVtx, p, hg->vindex, hg->vedge, NULL,
                            "coarsening plot");
        
        
        /* Allocate and initialize Matching Array */
        if (hg->nVtx && !(match = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory for Matching array");
            return ZOLTAN_MEMERR;
        }
        for (i = 0; i < hg->nVtx; i++)
            match[i] = i;
        
        /* Calculate matching (packing or grouping) */
        err = Zoltan_PHG_Matching (zz, hg, match, hgp);
        if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
            ZOLTAN_FREE ((void**) &match);
            goto End;
        }
        
        if (!(coarser = newVCycle(NULL, NULL, NULL, vcycle))) {
            ZOLTAN_FREE ((void**) &match);
            ZOLTAN_PRINT_ERROR (zz->Proc, yo, "coarser is NULL.");
            goto End;
        }
        
        /* Construct coarse hypergraph and LevelMap */
        err = Zoltan_PHG_Coarsening (zz, hg, match, coarser->hg, 
         vcycle->LevelMap, &vcycle->LevelCnt, &vcycle->LevelData);
        if (err != ZOLTAN_OK && err != ZOLTAN_WARN) 
            goto End;
        
        ZOLTAN_FREE ((void**) &match);

        if ((err=allocVCycle(coarser))!= ZOLTAN_OK)
            goto End;
        vcycle = coarser;
        hg = vcycle->hg;
    }


    if (hgp->output_level >= PHG_DEBUG_LIST) {
        uprintf(hg->comm, "START %3d |V|=%6d |E|=%6d |I|=%6d %d/%s/%s-%s p=%d...\n",
                hg->info, hg->nVtx, hg->nEdge, hg->nPins, hg->redl, hgp->redm_str,
                hgp->coarsepartition_str, hgp->refinement_str, p);
        if (hgp->output_level > PHG_DEBUG_LIST) {
            err = Zoltan_HG_Info(zz, hg);
            if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
                goto End;
        }
    }
    if (hgp->output_level >= PHG_DEBUG_PLOT)
        Zoltan_PHG_Plot(zz->Proc, hg->nVtx, p, hg->vindex, hg->vedge, NULL,
                        "coarsening plot");

    
    /****** Coarse Partitioning ******/
    err = Zoltan_PHG_CoarsePartition (zz, hg, p, part_sizes, vcycle->part, hgp);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
        goto End;

    del = vcycle;
    /****** Uncoarsening/Refinement ******/
    while (vcycle) {
        VCycle *finer = vcycle->finer;
        hg = vcycle->hg;

        err = Zoltan_PHG_Refinement (zz, hg, p, vcycle->part, hgp);
        
        if (hgp->output_level >= PHG_DEBUG_LIST)     
            uprintf(hg->comm, "FINAL %3d |V|=%6d |E|=%6d |I|=%6d %d/%s/%s-%s p=%d bal=%.2f cutl=%.2f\n",
                    hg->info, hg->nVtx, hg->nEdge, hg->nPins, hg->redl, hgp->redm_str,
                    hgp->coarsepartition_str, hgp->refinement_str, p,
                    Zoltan_PHG_Compute_Balance(zz, hg, p, vcycle->part),
                    Zoltan_PHG_hcut_size_links(&hgp->comm, hg, vcycle->part, p));

        if (hgp->output_level >= PHG_DEBUG_PLOT)
            Zoltan_PHG_Plot(zz->Proc, hg->nVtx, p, hg->vindex, hg->vedge, Parts,
                            "partitioned plot");
        
        /* Project coarse partition to fine partition */
        if (finer)  { 
           PHGComm *hgc = hg->comm;
           int each_size [hgc->nProc_x], displs[hgc->nProc_x];
           int gno, lpart, size, *ip;
           char *rbuffer;

           for (i = 0; i < finer->hg->nVtx; i++)
             if (finer->LevelMap[i] >= 0)
                finer->part[i] = vcycle->part[finer->LevelMap[i]];
                 
           for (i = 0; i < finer->LevelCnt; i++)  {
              i++;                         /* skip gno entry */
              finer->LevelData[i] = vcycle->part [finer->LevelMap [finer->LevelData[i]]];
           }
                         
           MPI_Allgather (&finer->LevelCnt, 1, MPI_INT, each_size, 1, MPI_INT,
            hgc->row_comm);

           displs[0] = size = 0;
           for (i = 0; i < hgc->nProc_x; i++)
             size += each_size[i];                 /* compute total size of rbuffer */
           for (i = 1; i < hgc->nProc_x; i++)
             displs[i] = displs[i-1] + each_size[i-1];    /* message displacements */
        
           if (!(rbuffer = (char*) ZOLTAN_MALLOC (1 + size * sizeof(int))))    {
             ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
             return ZOLTAN_MEMERR;
           }          
        
           MPI_Allgatherv (finer->LevelData, each_size[hgc->myProc_x], MPI_INT, rbuffer,
             each_size, displs, MPI_INT, hgc->row_comm);
                
           ip = (int*) rbuffer;  
           for (i = 0; i < size; i += 2)   {
              gno   = *ip++;
              lpart = *ip++;
              if (VTX_TO_PROC_X (finer->hg, gno) == hgc->myProc_x)
                 finer->part[VTX_GNO_TO_LNO (finer->hg, gno)] = lpart;
           }          
        }
        vcycle = finer;
    }
    
 End:
    vcycle = del;
    while (vcycle) {
        if (vcycle->finer) {
            Zoltan_HG_HGraph_Free (vcycle->hg);
            Zoltan_Multifree (__FILE__, __LINE__, 2, &vcycle->part, &vcycle->LevelMap);
        } else
            ZOLTAN_FREE(&vcycle->LevelMap);
        
        del = vcycle;
        vcycle = vcycle->finer;
        ZOLTAN_FREE(&del);
    }

    ZOLTAN_TRACE_EXIT(zz, yo) ;
    return err;
}
    


/****************************************************************************/
/* Calculates the cutsize of a partition by summing the weight of all edges
   which span more than one part. Time O(|I|). */
double Zoltan_PHG_hcut_size_total (PHGComm *hgc, HGraph *hg, Partition part, int p)
{
    int i, j, *netpart, *allparts;    
    double cut = 0.0, totalcut=0.0;
    char *yo = "Zoltan_PHG_hcut_size_total";

    if (!(netpart = (int*) ZOLTAN_CALLOC (hg->nEdge, sizeof(int)))) {
        ZOLTAN_PRINT_ERROR(hgc->myProc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }

    if (!hgc->myProc_x)
        if (!(allparts = (int*) ZOLTAN_CALLOC (hgc->nProc_x*hg->nEdge, sizeof(int)))) {
            ZOLTAN_PRINT_ERROR(hgc->myProc, yo, "Insufficient memory.");
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
double Zoltan_PHG_hcut_size_links (PHGComm *hgc, HGraph *hg, Partition part, int p)
{
    int i, j, *cuts=NULL, *rescuts=NULL, *parts, nparts;
    double cut = 0.0, totalcut=0.0;
    char *yo = "Zoltan_PHG_hcut_size_links";
    
    if (hg->nEdge) {
        if (!(cuts = (int*) ZOLTAN_CALLOC (p*hg->nEdge, sizeof(int)))) {
            ZOLTAN_PRINT_ERROR(hgc->myProc, yo, "Insufficient memory.");
            return ZOLTAN_MEMERR;
        }   
        if (!hgc->myProc_x)
            if (!(rescuts = (int*) ZOLTAN_CALLOC (p*hg->nEdge, sizeof(int)))) {
                ZOLTAN_PRINT_ERROR(hgc->myProc, yo, "Insufficient memory.");
                return ZOLTAN_MEMERR;
            }
    
        for (i = 0; i < hg->nEdge; ++i) {
            parts = &cuts[i*p];
            for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j) 
                ++parts[part[hg->hvertex[j]]];
        }
    
        MPI_Reduce(cuts, rescuts, p*hg->nEdge, MPI_INT, MPI_SUM, 0, 
                   hgc->row_comm);
        ZOLTAN_FREE ((void**) &cuts);
    }

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

double Zoltan_PHG_Compute_Balance (
  ZZ *zz,
  HGraph *hg,
  int p,
  Partition part
)
{
  int i;
  char *yo = "Zoltan_PHG_Compute_Balance";
  double *lsize_w, *size_w, max_size_w=0.0, tot_w = 0.0;

  if (!(lsize_w = (double*) ZOLTAN_CALLOC (p, sizeof(double))) 
      ||!(size_w = (double*) ZOLTAN_CALLOC (p, sizeof(double)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
  }
  for (i = 0; i < hg->nVtx; i++)
      lsize_w[part[i]] += (hg->vwgt ? hg->vwgt[i] : 1.0);
  MPI_Allreduce(lsize_w, size_w, p, MPI_DOUBLE, MPI_SUM, hg->comm->row_comm);
  
  for (i = 0; i < p; i++) {
      max_size_w = MAX(max_size_w, size_w[i]);
      tot_w += size_w[i];
  }

  Zoltan_Multifree(__FILE__,__LINE__, 2, &size_w, &lsize_w);

  return (max_size_w * p / tot_w);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
