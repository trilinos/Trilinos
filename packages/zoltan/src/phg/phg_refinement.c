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

#include <math.h>
#include <float.h>
#include "phg.h"


    /*

#define _DEBUG
#define _DEBUG2
#define _DEBUG3
    */
    
static ZOLTAN_PHG_REFINEMENT_FN refine_no;
static ZOLTAN_PHG_REFINEMENT_FN refine_fm2;

/****************************************************************************/

ZOLTAN_PHG_REFINEMENT_FN *Zoltan_PHG_Set_Refinement_Fn(char *str)
{
  
  if      (!strcasecmp(str, "fm2"))            return refine_fm2;  
  else if (!strcasecmp(str, "no"))             return refine_no;
  else                                         return NULL;
}



/****************************************************************************/
int Zoltan_PHG_Refinement (ZZ *zz, PHGraph *hg, int p, Partition part,
PHGPartParams *hgp)
{
  return hgp->Refinement(zz, hg, p, part, hgp, hgp->bal_tol);
}



/****************************************************************************/
static int refine_no (
  ZZ *zz,     /* Zoltan data structure */
  PHGraph *hg,
  int p,
  Partition part,
  PHGPartParams *hgp,
  float bal_tol
)
{
  return ZOLTAN_OK;
}



/*****************************************************************************/
/* 2-way Parallel FM refinement. No data movement between processors,
 * just relabeling of vertices. In each FM pass we only move in one
 * direction, from the heavier partition to the lighter one.
 */


static void fm2_move_vertex_oneway(int v, PHGraph *hg, Partition part, float *gain, HEAP *heap,
                                   int *pins[2], int *lpins[2], double *weights, double *lweights,
                                   int *mark, int *adj)
{
    int   pno=part[v], vto=1-pno, adjsz=0, j, i;
    
    mark[v] = 1;  /* mark as moved */
    part[v] = vto;
    weights[pno] -= (hg->vwgt ? hg->vwgt[v] : 1.0);
    weights[vto] += (hg->vwgt ? hg->vwgt[v] : 1.0);
    lweights[pno] -= (hg->vwgt ? hg->vwgt[v] : 1.0);
    lweights[vto] += (hg->vwgt ? hg->vwgt[v] : 1.0);

    for (j = hg->vindex[v]; j < hg->vindex[v+1]; j++) {
        int n = hg->vedge[j];
        float w = hg->ewgt ? hg->ewgt[n] : 1.0;
    
        --pins[pno][n];
        --lpins[pno][n];
        ++pins[vto][n];
        ++lpins[vto][n];

#ifdef _DEBUG
        if (pins[pno][n] < 0)
            errexit("move of %d makes pin[%d][%d]=%d", v, pno, n, pins[pno][n]);
#endif

        if ((pins[0][n] + pins[1][n])==1) /* size 1 net; it is never critical */
            continue;

        if (pins[pno][n]==1) {
            for (i = hg->hindex[n]; i < hg->hindex[n+1]; ++i) {
                int u = hg->hvertex[i]; 
                if (part[u]==pno) {
                    gain[u] += w;
                    if (!mark[u]) {
                        adj[adjsz++] = u;
                        mark[u] = -1;  /* mark neighbors with -1 */
                    }
                }
            }
        }

        if (pins[vto][n]==1) { /* now there is at least one pin here */
            for (i = hg->hindex[n]; i < hg->hindex[n+1]; ++i) {
                int u = hg->hvertex[i];
                if (part[u]==pno) { 
                    gain[u] += w;
                    if (!mark[u]) {
                        adj[adjsz++] = u;
                        mark[u] = -1;  /* mark neighbors with -1 */
                    }
                }
            }
        } 
    }
    
    for (i=0; i<adjsz; i++) {
        int u=adj[i], p=part[u];
        

#ifdef _DEBUG
        if (mark[u]!=-1)
            errexit("hey while moving v=%d mark[%d]=%d", v, u, mark[u]);
        if (part[u] == vto)
            errexit("hey while moving v=%d u=%d is in part %d", v, u, part[u]);
#endif
        mark[u] = 0;
        if (Zoltan_heap_has_elem(&heap[p], u))
            Zoltan_heap_change_value(&heap[p], u, gain[u]);
    }
}


static void fm2_move_vertex_oneway_nonroot(int v, PHGraph *hg, Partition part, int *lpins[2],
                                           double *weights, double *lweights)
{
    int   pno=part[v], vto=1-pno, j;
    
    part[v] = vto;
    weights[pno] -= (hg->vwgt ? hg->vwgt[v] : 1.0);
    weights[vto] += (hg->vwgt ? hg->vwgt[v] : 1.0);
    lweights[pno] -= (hg->vwgt ? hg->vwgt[v] : 1.0);
    lweights[vto] += (hg->vwgt ? hg->vwgt[v] : 1.0);

    for (j = hg->vindex[v]; j < hg->vindex[v+1]; j++) {
        int n = hg->vedge[j];
    
        --lpins[pno][n];
        ++lpins[vto][n];
    }
}


#if 0

typedef int  (*SelectFunc)(HEAP heap[2], double *weights, double *max_weight, double targetw0);

static int fm2_select(HEAP heap[2], double *weights, double *max_weight, double targetw0)
{
    int from;
    /* select a vertex with max gain; if possible */
    if (Zoltan_heap_not_empty(&heap[0]) && Zoltan_heap_not_empty(&heap[1])) {
        if (Zoltan_heap_max_value(&heap[0])==Zoltan_heap_max_value(&heap[1]))
            from = (weights[0] < targetw0) ? 1 : 0;
        else
            from = (Zoltan_heap_max_value(&heap[0])>Zoltan_heap_max_value(&heap[1])) ? 0 : 1;
    } else if (Zoltan_heap_empty(&heap[0])) {
        if (Zoltan_heap_empty(&heap[1])) /* too bad both are empty */
            return -1; /* nothing to select */
        else
            from = 1;
    } else
        from = 0;
    return Zoltan_heap_extract_max(&heap[from]);    
}


static void fm2_move_vertex(int v, PHGraph *hg, Partition part, float *gain, HEAP *heap, int *pins[2], int *lpins[2], double *weights, double *lweights, int *mark, int *adj)
{
    float oldgain=gain[v];
    int   pno=part[v], vto=1-pno, adjsz=0, j, i;
    
    mark[v] = 1;
    part[v] = vto;
    weights[pno] -= (hg->vwgt ? hg->vwgt[v] : 1.0);
    weights[vto] += (hg->vwgt ? hg->vwgt[v] : 1.0);
    lweights[pno] -= (hg->vwgt ? hg->vwgt[v] : 1.0);
    lweights[vto] += (hg->vwgt ? hg->vwgt[v] : 1.0);

    for (j = hg->vindex[v]; j < hg->vindex[v+1]; j++) {
        int n = hg->vedge[j];
        float w = hg->ewgt ? hg->ewgt[n] : 1.0;
    
        --pins[pno][n];
        --lpins[pno][n];
#ifdef _DEBUG2
        if (pins[pno][n] < 0)
            errexit("move of %d makes pin[%d][%d]=%d", v, pno, n, pins[pno][n]);
#endif

        if (!pins[pno][n]) {  /* no pin in source part */ 
            for (i = hg->hindex[n]; i < hg->hindex[n+1]; ++i) {
                int u = hg->hvertex[i]; 
                gain[u] -= w;
                if (!mark[u]) {
                    adj[adjsz++] = u;
                    mark[u] = -1;
		}
	    }
        } else if (pins[pno][n]==1) {
            for (i = hg->hindex[n]; i < hg->hindex[n+1]; ++i) {
                int u = hg->hvertex[i]; 
                if (part[u]==pno) {
                    gain[u] += w;
                    if (!mark[u]) {
                        adj[adjsz++] = u;
                        mark[u] = -1;
                    }
                }
            }
        }

        ++pins[vto][n];
        ++lpins[vto][n];
        if (pins[vto][n]==1) { /* now there is at least one pin here */
            for (i = hg->hindex[n]; i < hg->hindex[n+1]; ++i) {
                int u = hg->hvertex[i];
                if (part[u]==pno) { 
                    gain[u] += w;
                    if (!mark[u]) {
                        adj[adjsz++] = u;
                        mark[u] = -1;
                    }
                }
            }
        } else if (pins[vto][n]==2)
            for (i = hg->hindex[n]; i < hg->hindex[n+1]; ++i) {
                int u = hg->hvertex[i];
                if (part[u]==vto) { 
                    gain[u] -= w;
                    if (!mark[u]) {
                        adj[adjsz++] = u;
                        mark[u] = -1;
                    }
                }
            }
    }

    gain[v] = -oldgain;    
    for (i=0; i<adjsz; i++) {
        int u=adj[i], p=part[u];
        
        mark[u] = 0;
        if (Zoltan_heap_has_elem(&heap[p], u))
            Zoltan_heap_change_value(&heap[p], u, gain[u]);
    }
}
#endif

static int refine_fm2 (
  ZZ *zz,
  PHGraph *hg,
  int p,
  Partition part,
  PHGPartParams *hgp,
  float bal_tol
)
{
    int    i, j,  *pins[2], *lpins[2], *moves=0, *mark=0, *adj=0, round=0, best_cutsizeat, cont;
    double total_weight, ltotal_weight, weights[2], lweights[2], lmax_weight[2];
    double targetw0, ltargetw0, minvw=DBL_MAX;
    double cutsize, best_cutsize, ratio = hg->ratio, best_imbal, best_limbal, imbal, limbal;
    float  *gain=0, *lgain=0;
    HEAP   heap[2];
    char   *yo="local_fm2";
    PHGComm *hgc=&hgp->comm;
    struct {
        int nNonZero;
        int rank;
    } rootin, root;

    /*    SelectFunc select_func = fm2_select;*/
        
#ifdef _DEBUG3    
    FILE *fp;
    char fname[512];

    sprintf(fname, "uvc-phg-debug.%d.txt", zz->Proc);    
    fp = fopen(fname, "w");
    if (!fp) {
        printf ("failed to open debug file '%s'\n", fname);
        return ZOLTAN_FATAL;
    }
    fprintf(fp, "%s\n", uMe(hgc));
    fprintf(fp, "H(%d, %d, %d)\n", hg->nVtx, hg->nEdge, hg->nNonZero);
    fprintf(fp, "p=%d  bal_tol = %.3f\n PartVec:\n", p, bal_tol);    
    for (i = 0; i < hg->nVtx; ++i)
        fprintf(fp, "%d ", part[i]);
    fprintf(fp, "\n\n");
    Zoltan_PHG_Print(zz, hg, fp);
    fclose(fp);
#endif

    if (p != 2) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "p!=2 not allowed for phg_fm2.");
        return ZOLTAN_FATAL;
    }
    
    
    if (hg->nEdge == 0)
        return ZOLTAN_OK;

    /* find the index of the proc in column group with the most #nonzeros; it will be our root
       proc for computing moves since it has better knowedge about global hypergraph */
    rootin.nNonZero = hg->nNonZero; 
    rootin.rank = hgc->myProc_y;
    MPI_Allreduce(&rootin, &root, 1, MPI_2INT, MPI_MAXLOC, hgc->col_comm);
    
    /* Calculate the weights in each partition and total, then maxima */
    weights[0] = weights[1] = 0.0;
    lweights[0] = lweights[1] = 0.0;
    if (hg->vwgt) 
        for (i = 0; i < hg->nVtx; i++) {
            lweights[part[i]] += hg->vwgt[i];
            minvw = (minvw > hg->vwgt[i]) ? hg->vwgt[i] : minvw;
        }
    else {
        minvw = 1.0;
        for (i = 0; i < hg->nVtx; i++)
            lweights[part[i]] += 1.0;
    }

    MPI_Allreduce(lweights, weights, 2, MPI_DOUBLE, MPI_SUM, hgc->row_comm);
    total_weight = weights[0] + weights[1];
    targetw0 = total_weight * ratio;
    ltotal_weight = lweights[0] + lweights[1];
    ltargetw0 = ltotal_weight * ratio;
    lmax_weight[0] = ltotal_weight * bal_tol *      ratio;
    lmax_weight[1] = ltotal_weight * bal_tol * (1 - ratio);

    if (!(pins[0]     = (int*) ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int)))
        || !(lpins[0] = (int*) ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int)))
        || !(moves    = (int*)   ZOLTAN_MALLOC(hg->nVtx * sizeof(int)))
        || !(lgain    = (float*) ZOLTAN_CALLOC(hg->nVtx, sizeof(float))) ) {
         Zoltan_Multifree(__FILE__,__LINE__, 4, &pins[0], &lpins[0], &moves, &lgain);
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
    }
    pins[1] = &(pins[0][hg->nEdge]);
    lpins[1] = &(lpins[0][hg->nEdge]);

    if (hgc->myProc_y==root.rank) { /* only root needs mark, adj, gain and heaps*/
        if (!(mark     = (int*)   ZOLTAN_CALLOC(hg->nVtx, sizeof(int)))
            || !(adj      = (int*)   ZOLTAN_MALLOC(hg->nVtx * sizeof(int)))                
            || !(gain     = (float*) ZOLTAN_CALLOC(hg->nVtx, sizeof(float)))) {
         Zoltan_Multifree(__FILE__,__LINE__, 3, &mark, &adj, &gain);
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
        }        
        Zoltan_heap_init(zz, &heap[0], hg->nVtx);
        Zoltan_heap_init(zz, &heap[1], hg->nVtx);  
    }

    /* Initial calculation of the local pin distribution (sigma in UVC's papers)  */
    for (i = 0; i < hg->nEdge; ++i)
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j)
            ++(lpins[part[hg->hvertex[j]]][i]);

    do {
        int v=1, movecnt=0, neggaincnt=0, from, to;
        int maxneggain = (hgp->fm_max_neg_move < 0) ? hg->nVtx : hgp->fm_max_neg_move;

        /* now compute global pin distribution */
        MPI_Allreduce(lpins[0], pins[0], 2*hg->nEdge, MPI_INT, MPI_SUM, hgc->row_comm);

        /* now we can compute actual cut */
        best_cutsizeat=0;
        cutsize = 0.0;
        for (i=0; i < hg->nEdge; ++i) {
            if (pins[0][i] && pins[1][i])
                cutsize += (hg->ewgt ? hg->ewgt[i] : 1.0);
        }
        MPI_Allreduce(&cutsize, &best_cutsize, 1, MPI_DOUBLE, MPI_SUM, hgc->col_comm);
        cutsize = best_cutsize;
        best_imbal = imbal = fabs(weights[0]-targetw0)/targetw0;
        best_limbal = limbal = fabs(lweights[0]-ltargetw0)/ltargetw0;


        /* decide which way the moves will be in this pass */
        from = (weights[0] < targetw0) ? 1 : 0;
        /* we want to be sure that everybody!!! picks the same source */
        MPI_Bcast(&from, 1, MPI_INT, 0, hgc->Communicator); 
        to = 1-from;
                
#ifdef _DEBUG
    /* Just for debugging */
        best_cutsize = Zoltan_PHG_hcut_size_total(hgc, hg, part, p);
        if (best_cutsize!=cutsize) {
            errexit("%s: Initial cutsize=%.2lf Verify: total=%.2lf\n", uMe(hgc), cutsize,
               best_cutsize);
        }
        if (hgc->myProc_y==root.rank)
            for (i = 0; i< hg->nVtx; ++i)
                if (mark[i])
                    errexit("mark[%d]=%d", i, mark[i]);
        /* debuggging code ends here */
#endif

        /* compute only the gains of the vertices from 'from' part */
        for (i = 0; i < hg->nVtx; ++i) 
            if (part[i]==from) {
                lgain[i] = 0;
                for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++) {
                    int edge = hg->vedge[j];
                    if ((pins[0][edge]+pins[1][edge])>1) { /* if they have at least 2 pins :) */
                        if (pins[part[i]][edge] == 1)
                            lgain[i] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
                        else if (pins[1-part[i]][edge] == 0)
                            lgain[i] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
                    }
                }
            }
        /* now sum up all gains on only root proc */
        MPI_Reduce(lgain, gain, hg->nVtx, MPI_FLOAT, MPI_SUM, root.rank, hgc->col_comm);

        if (hgp->output_level >= PHG_DEBUG_ALL)        
            printf("%s FM Pass %d (%d->%d) Cut=%.2lf W[%5.0lf, %5.0lf] I= %.2lf LW[%5.0lf, %5.0lf] LI= %.2lf\n", uMe(hgc), round, from, to, cutsize, weights[0], weights[1], imbal, lweights[0], lweights[1], limbal);

        if (hgc->myProc_y==root.rank) { /* those are the lucky ones; each proc in column-group
                                 could have compute the same moves concurrently; but for this
                                 version we'll do in in the root procs and broadcast */
            /* Initialize the heaps and fill them with the gain values */
            Zoltan_heap_clear(&heap[from]);  
            for (i = 0; i < hg->nVtx; ++i)
                if (part[i]==from)
                    Zoltan_heap_input(&heap[from], i, gain[i]);
            Zoltan_heap_make(&heap[from]);
            
            while ((v>=0) && (neggaincnt < maxneggain) && ((lweights[to]+minvw) <= lmax_weight[to]) ) {
                if (Zoltan_heap_empty(&heap[from])) /* too bad it is empty */
                    break;
                v = Zoltan_heap_extract_max(&heap[from]);    
                
#ifdef _DEBUG
                if (from != part[v])
                    errexit("hooop from=%d part[%d]=%d", from, v, part[v]);
#endif

                /* Mark vertex we picked from the heap so it is "locked". 
                   For the current strategy, moving only one direction 
                   at a time, the mark information is not critical.
                   Note that the mark array is also used in the move/update 
                   routine so don't remove it! */
                ++mark[v];
                if (lweights[to]+((hg->vwgt)? hg->vwgt[v] : 1.0) > lmax_weight[to]) {
#ifdef _DEBUG2                    
                    printf("%s %4d: %6d (g: %5.1lf), p:%2d [%4.0lf, %4.0lf] NF\n", uMe(hgc), movecnt, v, gain[v], from, weights[0], weights[1]);
#endif
                    /* Negative value in moves array means we have examined 
                       the vertex but couldn't move it. Note offset by one,
                       otherwise zero would be ambiguous. */
                    moves[movecnt++] = -(v+1);
                    continue;
                } 
	
                /* Positive value means we actually moved the vertex. 
                   Offset by one not really necessary here? */
                moves[movecnt] = (v+1);
                ++neggaincnt;
                cutsize -= gain[v];

                fm2_move_vertex_oneway(v, hg, part, gain, heap, pins, lpins, weights, lweights, mark, adj);
                imbal = fabs(weights[0]-targetw0)/targetw0;
                limbal = fabs(lweights[0]-ltargetw0)/ltargetw0;

                /* UVC: note that in the loop we're only using local imbal; hence FM might want to continue
                   to improve local imbalance; but it might be improving the global imbalance at all!
                   Strangely, if we let FM do this; it seems to be finding better local optimums.
                   So I'll leave it as it is.
                */
                
                if ((cutsize<best_cutsize) || (cutsize==best_cutsize && limbal < best_limbal)) {
#ifdef _DEBUG2                    
                    printf("%s %4d: %6d (g: %5.1lf), p:%2d W[%4.0lf, %4.0lf] I:%.2lf LW[%4.0lf, %4.0lf] LI:%.2lf C:%.1lf<-- Best\n", uMe(hgc), movecnt, v, gain[v], from, weights[0], weights[1], imbal, lweights[0], lweights[1], limbal, cutsize); /* after move gain is -oldgain */
#endif
                    
                    best_cutsize = cutsize;
                    best_cutsizeat = movecnt+1;
                    best_limbal = limbal;
                    neggaincnt = 0;
                }
#ifdef _DEBUG2                
                else
                    printf("%s %4d: %6d (g: %5.1lf), p:%2d [%4.0lf, %4.0lf] %.1lf\n", uMe(hgc), movecnt, v, gain[v], from, weights[0], weights[1], cutsize);
#endif
                ++movecnt;
            }

#ifdef _DEBUG
	    if (v<0)
	      uprintf(hgc, "EOLB @ %d there was no vertex to select: v=%d\n", movecnt, v);
	    else if (neggaincnt >= maxneggain) 
	      uprintf(hgc, "EOLB @ %d max neg move reached neggaincnt(%d) >= maxneggain\n", movecnt, neggaincnt, maxneggain);
	    else 
	      uprintf(hgc, "EOLB @ %d balance constraint LW[%.1lf, %.1lf] and MAXW[%.1lf, %.1lf]\n", movecnt, lweights[0], lweights[1], lmax_weight[0], lmax_weight[1]);
#endif
            
            /* roll back the moves without any improvement */
            for (i=movecnt-1; i>=best_cutsizeat; --i) {
                int v = moves[i];
                if (v<0)
                    v = -v-1;
                else {
                    --v;
                    fm2_move_vertex_oneway(v, hg, part, gain, heap, pins, lpins, weights, lweights, mark, adj);
                }
                mark[v] = 0;
            }
            for (i=0; i<best_cutsizeat; ++i){
                int v = (moves[i] < 0 ) ? -moves[i] - 1 : moves[i]-1;
                mark[v] = 0;
            }                
        }
        
        /* now root bcast moves to column procs */
        MPI_Bcast(&best_cutsizeat, 1, MPI_INT, root.rank, hgc->col_comm);
        MPI_Bcast(moves, best_cutsizeat, MPI_INT, root.rank, hgc->col_comm);
        if (hgc->myProc_y!=root.rank) { /* now non-root does move simulation */
            for (i=0; i<best_cutsizeat; ++i) {
                int v = moves[i];
                if (v<0)
                    v = -v-1;
                else {
                    --v;
                    fm2_move_vertex_oneway_nonroot(v, hg, part, lpins, weights, lweights);
                }
            }
        }

#ifdef _DEBUG
        for (i = 0; i < hg->nEdge; ++i) {
            int lp[2];

            lp[0] = lp[1] = 0;
            for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j)
                ++(lp[part[hg->hvertex[j]]]);
            if ((lp[0] != lpins[0][i]) || (lp[1] != lpins[1][i]))
                errexit("for net %d -- lp=[%d, %d] lpins[%d, %d]", i, lp[0], lp[1], lpins[0][i], lpins[1][i]);
        }
#endif

        
        MPI_Allreduce(lweights, weights, 2, MPI_DOUBLE, MPI_SUM, hgc->row_comm);
        cont = 0;
        MPI_Allreduce(&best_cutsizeat, &cont, 1, MPI_INT, MPI_LOR, hgc->row_comm);
#ifdef _DEBUG
    /* Just for debugging */
        best_cutsize = Zoltan_PHG_hcut_size_total(hgc, hg, part, p);
        imbal = fabs(weights[0]-targetw0)/targetw0;
        printf("%s End of Pass %d Comp.Cut=%.2lf RealCut=%.2lf W[%5.0lf, %5.0lf] Imbal=%.2lf\n", uMe(hgc), round, cutsize, best_cutsize, weights[0], weights[1], imbal);
        if (cutsize<best_cutsize) {
            errexit("*** HEY HEY Invalid cut!!!");
        }
        /* debuggging code ends here */
#endif
    } while (cont &&  (++round < hgp->fm_loop_limit));

    if (hgc->myProc_y==root.rank) { /* only root needs mark, adj, gain and heaps*/
        Zoltan_Multifree(__FILE__,__LINE__, 3, &mark, &adj, &gain);
        Zoltan_heap_free(&heap[0]);
        Zoltan_heap_free(&heap[1]);        
    }
    
    Zoltan_Multifree(__FILE__, __LINE__, 4, &pins[0], &lpins[0], &moves, &lgain);
    return ZOLTAN_OK;
}





#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
