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

static int split_hypergraph (int *pins[2], PHGraph*, PHGraph*, Partition, int, ZZ*);



/* recursively divides problem into 2 parts until all p found */
int Zoltan_PHG_rdivide (int lo, int hi, Partition final, ZZ *zz, PHGraph *hg,
                        PHGPartParams *hgp, int level)
{
    int i, j, mid, err, *pins[2], *lpins[2];
    Partition part;
    PHGraph *new;
    PHGComm *hgc=hg->comm;
    char *yo = "Zoltan_PHG_rdivide";
    
/*   hg->redl = 2;    */  /* this seems to be the right thing to do, but trying to match previous answers today. */
    hg->redl = hgp->redl;

    uprintf(hgc, "Rdivide(%d, %d)\n", lo, hi);
    /* only one part remaining, record results and exit */
    if (lo == hi) {
        for (i = 0; i < hg->nVtx; i++)
            final [hg->vmap[i]] = lo - 1;
        return ZOLTAN_OK;
    }

    part = (Partition) ZOLTAN_MALLOC (hg->nVtx * sizeof (int));
    if (part == NULL) {
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Unable to allocate memory.");
        return ZOLTAN_MEMERR;
    }

    /* bipartition current hypergraph with appropriate split ratio */
    mid = (lo+hi)/2;
    hg->ratio = (double) (mid-lo+1) / (double) (hi-lo+1);
    /*   hg->redl = 0;  */
    err = Zoltan_PHG_HPart_Lib (zz, hg, 2, part, hgp, level);
    if (err != ZOLTAN_OK) {
        ZOLTAN_FREE (&part);
        return err;
    }
    
    /* if only two parts total, record results and exit */
    if (lo + 1 == hi)  {
        for (i = 0; i < hg->nVtx; i++)
            final [hg->vmap[i]] = ((part[i] == 0) ? (lo-1) : (hi-1));
        ZOLTAN_FREE (&part);
        return ZOLTAN_OK;
    }


    if (!(pins[0]     = (int*) ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int)))
        || !(lpins[0] = (int*) ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int))) ) {
         Zoltan_Multifree(__FILE__,__LINE__, 2, &pins[0], &lpins[0]);
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
    }
    pins[1] = &(pins[0][hg->nEdge]);
    lpins[1] = &(lpins[0][hg->nEdge]);

    /* Initial calculation of the local pin distribution (sigma in UVC's papers)  */
    for (i = 0; i < hg->nEdge; ++i)
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j)
            ++(lpins[part[hg->hvertex[j]]][i]);
    /* now compute global pin distribution */
    MPI_Allreduce(lpins[0], pins[0], 2*hg->nEdge, MPI_INT, MPI_SUM, hgc->row_comm);
    ZOLTAN_FREE (&lpins[0]); /* we don't need lpins */
    
    uprintf(hgc, "before splitting for part 0\n");
    new = (PHGraph*) ZOLTAN_MALLOC (sizeof (PHGraph));
    if (new == NULL)  {
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Unable to allocate memory.");
        return ZOLTAN_MEMERR;
    }
    Zoltan_PHG_HGraph_Init (new);
    
    /* recursively divide in two parts and repartition hypergraph */
    err = split_hypergraph (pins, hg, new, part, 0, zz);
    uprintf(hgc, "after splitting for part 0\n");
    if (err != ZOLTAN_OK)
        return err;
    uprintf(hgc, "before calling RB for 0\n");
    err = Zoltan_PHG_rdivide (lo, mid, final, zz, new, hgp, level+1);
    uprintf(hgc, "after calling RB for 0\n");
    Zoltan_PHG_HGraph_Free (new);
    if (err != ZOLTAN_OK)
        return err;

    uprintf(hgc, "before splitting for part 1\n");
    err = split_hypergraph (pins, hg, new, part, 1, zz);
    ZOLTAN_FREE (&pins[0]); /* we don't need pins */
    uprintf(hgc, "after splitting for part 1\n");
    if (err != ZOLTAN_OK)
        return err;
    uprintf(hgc, "before calling RB for 1\n");
    err = Zoltan_PHG_rdivide (mid+1, hi, final, zz, new, hgp, level+1);
    uprintf(hgc, "after calling RB for 1\n");
    Zoltan_PHG_HGraph_Free (new);
    if (err != ZOLTAN_OK)
      return err;
    
    /* remove alloc'ed structs */
    ZOLTAN_FREE (&part);
    ZOLTAN_FREE (&new);
    
    return ZOLTAN_OK;
}



static int split_hypergraph (int *pins[2], PHGraph *old, PHGraph *new, Partition part,
                             int partid, ZZ *zz)
{
    int *tmap;        /* temporary array mapping from old HGraph info to new */
    int edge, i;                                            /* loop counters */
    char *yo = "split_hypergraph";
    
    /* allocate memory for dynamic arrays in new HGraph and for tmap array */
    new->vmap = (int*) ZOLTAN_MALLOC (old->nVtx * sizeof (int));
    tmap      = (int*) ZOLTAN_MALLOC (old->nVtx * sizeof (int));
    if (new->vmap == NULL || tmap == NULL)  {
        Zoltan_Multifree (__FILE__, __LINE__, 2, &new->vmap, &tmap);
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Unable to allocate memory 1.");
        return ZOLTAN_MEMERR;
    }
    
    /* save vertex and edge weights if they exist */
    if (old->vwgt)
        new->vwgt = (float*) ZOLTAN_MALLOC (old->nVtx  * sizeof(float)
                                            * old->VertexWeightDim);
    if (old->ewgt)
        new->ewgt = (float*) ZOLTAN_MALLOC (old->nEdge * sizeof(float)
                                            * old->EdgeWeightDim);
    
    /* fill in tmap array, -1 for ignored vertices, otherwise nonnegative int */
    new->nVtx = 0;
    for (i = 0; i < old->nVtx; i++)
        if (part[i] == partid)  {
            tmap[i] = new->nVtx;
            new->vmap[new->nVtx] = old->vmap[i];
            if (new->vwgt)
                new->vwgt[new->nVtx] = old->vwgt[i];
            new->nVtx++;
        }
        else
            tmap[i] = -1;

    new->comm = old->comm;
    /* continue allocating memory for dynamic arrays in new HGraph */
    new->vmap    = (int*) ZOLTAN_REALLOC (new->vmap, new->nVtx * sizeof (int));
    new->hindex  = (int*) ZOLTAN_MALLOC ((old->nEdge+1) * sizeof (int));
    new->hvertex = (int*) ZOLTAN_MALLOC (old->nNonZero * sizeof (int));
    if (new->vmap == NULL || new->hindex == NULL || new->hvertex == NULL)  {
        Zoltan_Multifree (__FILE__, __LINE__, 5, &new->vmap, &new->hindex,
                          &new->hvertex, &new->vmap, &tmap);
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Unable to allocate memory 2.");
        return ZOLTAN_MEMERR;
    }
    
    /* fill in hindex and hvertex arrays in new HGraph */
    new->nEdge  = 0;
    new->nNonZero = 0;
    for (edge = 0; edge < old->nEdge; edge++)
        if (pins[partid][edge]) {  /* edge has at least one vertex in partition */
            new->hindex[new->nEdge] = new->nNonZero;
            for (i = old->hindex[edge]; i < old->hindex[edge+1]; i++)
                if (tmap [old->hvertex[i]] >= 0)  {
                    new->hvertex[new->nNonZero] = tmap[old->hvertex[i]];
                    new->nNonZero++;  
                }
            if (new->ewgt)
                new->ewgt[new->nEdge] = old->ewgt[edge];
            new->nEdge++;
        }

    new->hindex[new->nEdge] = new->nNonZero;
    new->info = old->info;
    new->VertexWeightDim = old->VertexWeightDim;
    new->EdgeWeightDim   = old->EdgeWeightDim;
    /* UVCUVC: BUGBUG TODO CHECK we need to compute dist_x, dist_y */
    
    /* shrink hindex, hvertex arrays to correct size & determine vindex, vedge */
    new->hvertex = (int*)ZOLTAN_REALLOC(new->hvertex, sizeof(int) * new->nNonZero);
    new->hindex = (int*)ZOLTAN_REALLOC(new->hindex, sizeof(int) * (new->nEdge+1));
    Zoltan_PHG_Create_Mirror (zz, new);
    
    ZOLTAN_FREE (&tmap);
    return ZOLTAN_OK;
}
