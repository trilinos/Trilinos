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

#include <assert.h>
#include "phg.h"



static ZOLTAN_PHG_MATCHING_FN matching_no;   /* template -- matching */
static ZOLTAN_PHG_MATCHING_FN matching_ipm;  /* inner product matching */


/* static void check_upper_bound_of_matching_weight (Graph*, ZZ*, Matching); */
/* static int graph_connected_components (int, int*, int*, int);             */

static void print_matrix(int*, int*, int, int, int);

/*****************************************************************************/
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams *hgp)
{
    int exist = 1;
    if (!strcasecmp(hgp->redm_str, "no"))       hgp->matching = NULL;
    else if (!strcasecmp(hgp->redm_str, "ipm")) hgp->matching = matching_ipm;
    else {
        exist = 0;
        hgp->matching = NULL;
    }

  return exist;
}


static char *uMe(PHGComm *hgc)
{
    static char msg[1024];

    sprintf(msg, "<%d/%d>: (%d,%d)/[%d,%d] ->", hgc->Proc, hgc->Num_Proc, hgc->myProc_x, hgc->myProc_y, hgc->nProc_x, hgc->nProc_y);
    return msg;
}


/* Removed sim() from serial version at this point */

/*****************************************************************************/

int Zoltan_PHG_Matching (
  ZZ *zz,
  PHGraph *hg,
  Matching match,
  PHGPartParams *hgp)
{
float *old_ewgt = NULL, *new_ewgt = NULL;
int   err;
char  *yo = "Zoltan_PHG_Matching";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hg->vwgt && hgp->ews) {
     if (!(new_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        err = ZOLTAN_MEMERR;
        goto End;
     }
     Zoltan_PHG_Scale_HGraph_Weight (zz, hg, new_ewgt, hgp->ews);
     old_ewgt = hg->ewgt;
     hg->ewgt = new_ewgt;
  }

  /* Do the matching */
  if (hgp->matching) {
     err = hgp->matching (zz, hg, match);
     if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
        goto End;
  }

  /* Removed serial Optimization here */

End:

  /* Restore the old edge weights */
  if (hg->vwgt && hgp->ews)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}



/*****************************************************************************/
/* template for matching, hypergraph version */
static int matching_no (ZZ *zz, PHGraph *hg, Matching match)
{
  return ZOLTAN_OK;
}

/*****************************************************************************/

/*
 * inner product matching
 * Parallelized version of the serial algorithm (see hg_match.c)
 * Based on conversations with Rob Bisseling
 *
 * by Aaron Becker, UIUC, summer 2004
 *
 * For each vertex, we match with the unmatched vertex which has the most
 * hyperedges in common with it.  We first perform horizontal communication
 * and compute partial inner products, then perform vertical communication to
 * compute the complete inner product for each vertex (column).
 *
 * The matching stage is currently unimplemented.
 * 
 * Future improvements:
 *      -perform communication in rounds to reduce memory requirements at the
 *          expense of match quality
 *      -distribute final inner product computation over all processors by
 *          computing product b_ij at processor (i mod nProc_y, t) rather than
 *          doing all computations at (0, t)
 */
static int matching_ipm(ZZ *zz, PHGraph *hg, Matching match)
{
    int   i, j, n, v1, edge;
    int   send_size, err;
    int   *pips, *cips;
    int    cips_size, pips_size;
   
    int   *r_hindex, *r_vindex, *r_hvertex, *r_vedge; /* row-wide hypergraph */
    int   r_hindex_size, r_vindex_size, r_hvertex_size, r_vedge_size;
    char  *yo = "matching_ipm";

    /*
     * note: I assume that the match array comes 
     * in initialized with match[i] = i
     */

    MPI_Comm *comm = &zz->Communicator;
    PHGComm  *hgc = hg->comm;

    printf("Entering matching_ipm on proc %d ",zz->Proc);
    printf("(%d,%d)\n", hgc->myProc_x, hgc->myProc_y);
    
    /* delimit vedge array so processor boundaries are visible in r_vedge */
    hg->vedge[hg->vindex[hg->nVtx] - 1] = 
        -hg->vedge[hg->vindex[hg->nVtx] - 1] - 1;

    /* get row-wide vedge */    
    send_size = hg->vindex[hg->nVtx] * sizeof(int);
    err = Zoltan_PHG_gather_slice(
        hgc->nProc_x, hgc->nProc_y, hgc->myProc_x, hgc->myProc_y, send_size,
        (char*) hg->vedge, &r_vedge_size, (char**) &r_vedge, comm, 1);
    if (err < 0)
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zoltan_PHG_gather_slice failed.");
    r_vedge_size /= sizeof(int);
    
    /* delimit vindex */
    hg->vindex[hg->nVtx] = -hg->vindex[hg->nVtx] - 1;

    /* get row-wide vindex */
    send_size = (hg->nVtx + 1) * sizeof(int);
    err = Zoltan_PHG_gather_slice(
        hgc->nProc_x, hgc->nProc_y, hgc->myProc_x, hgc->myProc_y, send_size,
        (char*) hg->vindex, &r_vindex_size, (char**) &r_vindex, comm, 1);
    if (err < 0)
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zoltan_PHG_gather_slice failed.");
    r_vindex_size /= sizeof(int);


    /* make r_vindex index into r_vedge correctly                    */
    /* also reverse negative entries indicating processor boundaries */
    n = 0;
    for (i = 0; i < r_vindex_size; ++i) {
        while (r_vindex[i] >= 0) {
            r_vindex[i] += n;
            i++;
        }
        r_vindex[i] = -(r_vindex[i] + 1) + n;
        n = r_vindex[i] - n - 1;
        while (r_vedge[n] >= 0)
            n++;
        r_vedge[n] = -(r_vedge[n] + 1);
        n++;
    }

    
    /* build r_hvertex and r_hindex from r_vedge and r_vindex      */
    /* this functionality was copied from Zoltan_PHG_Create_Mirror */

    r_hindex_size = 0;
    for (i = 0; i < r_vedge_size; ++i)
        if (r_vedge[i] > r_hindex_size)
            r_hindex_size = r_vedge[i];
    r_hindex_size++;
    r_hvertex_size = r_vedge_size;
    
    if (!(r_hindex  = (int*) ZOLTAN_MALLOC(r_hindex_size * sizeof(int))) 
     || !(r_hvertex = (int*) ZOLTAN_MALLOC(r_hvertex_size * sizeof(int)))) {
        Zoltan_Multifree(__FILE__, __LINE__, 2, &r_hindex, &r_hvertex);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }
 
    /* count number of thing in the hindex space */
    for (i = 0; i <= r_hindex_size; i++)
        r_hindex[i] = 0;
    for (i = 0; i < r_vindex_size; i++)
        for (j = r_vindex[i];  j < r_vindex[i + 1]; j++)
            (r_hindex[r_vedge[j] + 1])++;

    /* compute partial sum */
    for (i = 1; i < r_hindex_size; i++)
        r_hindex[i] += r_hindex[i-1];

    /* store data in hindex location and increment index */
    for (i = 0; i < r_vindex_size; i++)
        for (j = r_vindex[i]; j < r_vindex[i+1]; j++)
            r_hvertex[(r_hindex[r_vedge[j]])++] = i;

    /* hindex now points past end of sublist. Shift array one position up */
    for (i = r_hindex_size; i > 0; i--)
        r_hindex[i] = r_hindex[i-1];
    r_hindex[0] = 0;

    /* UVCUVC: turned it off tooo much debug output
    print_matrix(r_vindex, r_vedge,   r_vindex_size, r_hindex_size, 0);
    print_matrix(r_hindex, r_hvertex, r_vindex_size, r_hindex_size, 1);
    */
    
    /* compute partial inner products using serial algorithm */
    
    /* 
     * note: pips currently computes partial ips for every vertex in r_vedge,
     * not just vertices belonging to this processor (ie those in hg->vedge).
     *
     * to change this, we need pips_size = hg->nVtx * r_vindex size, and the
     * indexing in the inner product loop must be changed.
     */
    cips = NULL;
    pips_size = r_vindex_size * r_vindex_size;

    printf("**** UVC: %s  r_vindex_size= %d   pips_size= %d \n r_vedge_size= %d",
           uMe(hgc), 
           r_vindex_size, pips_size, r_vedge_size);
    
    if (!(pips = (int*) ZOLTAN_MALLOC(pips_size * sizeof(int)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
        return ZOLTAN_MEMERR;
    }
    
    /* for every vertex on this processor*/
    for (v1 = 0; v1 < r_vindex_size; v1++) {
        
        /* for every hyperedge containing the vertex */
        for (i = r_vindex[v1]; i < r_vindex[v1+1]; i++) {
            edge = r_vedge[i];
                
            /* for every other vertex in the hyperedge */
            for (j = r_hindex[edge]; j < r_hindex[edge+1]; j++) {
                pips[v1 * hg->nVtx + r_hvertex[j]]++;
            }
        }
    }

    /* send results to column root */
    err = Zoltan_PHG_gather_slice_root(
            hgc->nProc_x, hgc->nProc_y, hgc->myProc_x, hgc->myProc_y,
            pips_size * sizeof(int), (char*)pips, 
            &cips_size, (char**)&cips, comm, 0);
    if (err < 0)
        ZOLTAN_PRINT_ERROR(zz->Proc,yo,"Zoltan_PHG_gather_slice_root failed.");
    cips_size /= sizeof(int);
   
    
    if(!hgc->myProc_y) {
        /* column root computes complete inner products, creates matching */
        /* and sends results back                                         */
        
        assert(cips_size == pips_size * hgc->nProc_y);

        /* compute complete inner products */
        for (i = 0; i < pips_size; ++i)
            for (j = 1; j < hgc->nProc_y; ++j)
                cips[i] += cips[j * pips_size];
        
/*****************************************************************************
 * 
 * compute matching
 *
 *  first find all vertices on your processor whose best match is also
 *  on your processor and match them.
 *  
 *  for (i = myProc_start; i < myProc_end; i += num_vertices)
 *      for (j = 0; j < num_vertices; ++j)
 *          if (cips[i + j] < best match so far)
 *              keep track of best match
 *      if (best match on my proc)
 *          modify match array
 *
 *  Send this partial match across the row with Zoltan_PHG_gather_slice
 *  Each processor updates its match array to reflect these matches, and there 
 *  are no conflicts to resolve because all matches are within a single 
 *  processor.
 *
 *  then do a series of matching rounds
 *      each unmatched vertex chooses its preferred remaining vertex as a 
 *      candidate match
 *
 *      candidate = malloc(size of matching)
 *      for (i = 0; i < |matching|; ++i)
 *          candidate[i] = matching[i]
 *      
 *      for (i = myProc_start; i < myProc_end; i+= num_vertices)
 *          for (j = 0; j < num_vertices; ++j)
 *              if (cips[i + j] < best match so far 
 *               && vertex is unmatched in candidate array)
 *                  keep track of best match
 *          modify candidate array
 *      
 *      this candidate matching is sent across the row with 
 *      Zoltan_PHG_gather_slice
 *      store the row-wide candidate array in r_candidate
 *
 *      all columns update their matching array
 *      in case of conflicts, the lower-numbered vertex wins
 *      
 *      for (i = 0; i < nProc_x; ++i)
 *          for (j = 0; j < |matching|; ++j)
 *              if (matching[j] == j)
 *                  matching[j] = r_candidate[j + i * |matching|]
 *  
 *  repeat until no further matching is possible
 *  
 *  send resulting matching array to everyone in column
 *  Zoltan_PHG_gather_slice (need matching call on other processors)
 *  
 ******************************************************************************
 *  
 *  Notes:
 *  Notice that all vertices belonging to the lowest-numbered processor with 
 *  unmatched vertices will always be successfully matched due to the 
 *  tiebreaking mechanism.  Therefore we will need at most nProc_x rounds of 
 *  matching.
 *  
 *  Also, we may want to ignore matches whose inner products are below some
 *  threshold value.  This ignores what might otherwise be poor matches, but
 *  we should be careful use the same rules for ignoring on each processor
 *  and in each round to ensure consistency across the row and to keep nProc_x
 *  as an upper bound on the number of rounds.
 *  
 *  Additionally, we may want a hard limit on the number of rounds we do.
 *  Late matches are unlikely to be of high quality, and might be better
 *  handled farther down the v-cycle.
 *  
 *  This method of tiebreaking is inherently unfair--it favors lower-numbered
 *  vertices.  To ease this problem, we might consider changing the tiebreaking
 *  criterion each round (for example, changing the processor whose vertices
 *  have highest priority).  We need only ensure that the rule is deterministic
 *  and identical on each processor.
 *
 *  The most significant problem to be addressed with this method is the high
 *  memory requirement.  Each processor needs 
 *  number of vertices * number of vertices / processor
 *  space to hold complete inner product information.  We can reduce this
 *  requirement by processing only a subset of vertices at once, but at the
 *  cost of match quality.
 *
 ******************************************************************************
 * 
 *  Rationale:
 *  we cannot send complete information to each processor; the memory 
 *  requirement will not scale
 *
 *  therefore, we must do local decision making with some global rules for 
 *  breaking ties
 *
 *  since matchings within a processor are already favorable, are likely 
 *  common, and will not produce conflicts, do them first
 *
 *  we must now match everything that remains.
 *  this is not an easy problem to do well (see Rob's discussion of stable 
 *  marriage on zoltan-dev)
 *  for now, I recommend the simple-minded strategy of:
 *      try the best match
 *      if this produces a conflict, the lower-numbered vertex wins
 *  this won't produce the highest-quality matching, but it is simple and 
 *  requires no extra communication.
 * 
 *****************************************************************************
 *  Notes by EB 8/6/2004:
 *  The approach above is similar to the LHM and RHM matching schemes.
 *  RHM will pick a random vertex and match it with the "best" vertex,
 *  where best is defined by a similarity function similar to inner products.
 *  LHM is more clever but more expensive, as it computes all similarities
 *  (inner products) and then performs an approximate maximal matching
 *  algorithm (LAM). Note that we should try to leverage earlier
 *  work as much as possible.
 *
 *  As mentioned earlier, computing all inner products is impractical.
 *  Instead, each processor could send a chunk of its data in each 
 *  of several communication rounds; then we do a partial matching using
 *  this partial data after each round. The greedy matching algorithm 
 *  may be used, or possibly LAM.
 *****************************************************************************/
 
    } else
        assert(cips_size == 0);
    
    fflush(NULL);
    MPI_Barrier(*comm);

    Zoltan_Multifree(__FILE__, __LINE__, 6, &pips, 
                                            &cips, 
                                            &r_hindex, 
                                            &r_hvertex,
                                            &r_vindex,
                                            &r_vedge);
    return ZOLTAN_OK;
}


/* 
 * draws the sparse matrix representation of a hypergraph, with vertices as
 * rows and columns as edges.  useful for debugging on small inputs
 */
static void print_matrix(int *index, int *data, int x, int y, int h)
{
    char *matrix = malloc(x * y * sizeof(char));
    int i, j;
 
    fflush(NULL);
 
    for (i = 0; i < y; ++i)
        for (j = 0; j < x; ++j)
            matrix[i * x + j] = '.';

    if (h) {
        for (i = 0; i < y; ++i)
            for (j = index[i]; j < index[i + 1]; ++j)
                matrix[i * x + data[j]] = 'x';
    }
    else {
        for (i = 0; i < x; ++i)
            for (j = index[i]; j < index[i + 1]; ++j)
                matrix[data[j] * x + i] = 'x';
    }
    
    for (i = 0; i < y; ++i) {
        for (j = 0; j < x; ++j)
            printf("%c ", matrix[i * x + j]);
        printf("\n");
    }
    printf("\n");
    fflush(NULL);
    
    free(matrix);
    return;
}


/*****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
