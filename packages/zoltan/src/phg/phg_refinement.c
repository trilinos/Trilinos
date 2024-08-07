// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <math.h>
#include <float.h>
#include "phg.h"
#include "zz_heap.h"
#include "zz_const.h"

#define BADBALANCE  2.0
#define HANDLE_ISOLATED_VERTICES    
#define USE_SERIAL_REFINEMENT_ON_ONE_PROC




    
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
  
  if      (!strcasecmp(str, "fm"))             return refine_fm2;  
  else if (!strcasecmp(str, "fm2"))            return refine_fm2;  
  else if (!strcasecmp(str, "no"))             return refine_no;
  else if (!strcasecmp(str, "none"))           return refine_no;
  else                                         return NULL;
}



/****************************************************************************/
int Zoltan_PHG_Refinement (ZZ *zz, HGraph *hg, int p, float *part_sizes, Partition part,
                           PHGPartParams *hgp)
{
    return hgp->Refinement(zz, hg, p, part_sizes, part, hgp, hgp->bal_tol);
}



/****************************************************************************/
static int refine_no (ZZ *zz,     /* Zoltan data structure */
                      HGraph *hg,
                      int p,
                      float *part_sizes,
                      Partition part,
                      PHGPartParams *hgp,
                      float bal_tol
    )
{
  return ZOLTAN_OK;
}



#ifdef USE_SERIAL_REFINEMENT_ON_ONE_PROC


/****************************************************************************/

/* Move a vertex from sour to dest, update cut and gain values. */
/* Note this function is also used by the coarse partitioner. */

int Zoltan_HG_move_vertex (HGraph *hg, int vertex, int sour, int dest,
 int *part, int **cut, double *gain, HEAP *heap)
{
int i, j;
int v, edge;

  gain[vertex] = 0.0;
  part[vertex] = dest;

  for (i = hg->vindex[vertex]; i < hg->vindex[vertex+1]; i++) {
     edge = hg->vedge[i];
     if (cut[sour][edge] == 1) {
        for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
           v = hg->hvertex[j];
           gain[v] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
           if (heap)
              Zoltan_Heap_Change_Value(&heap[part[v]], v, gain[v]);
           }
        }
     else if (cut[sour][edge] == 2) {
        for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
           v = hg->hvertex[j];
           if (part[v] == sour) {
              gain[v] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
              if (heap)
                 Zoltan_Heap_Change_Value(&heap[part[v]], v, gain[v]);
              break;
              }
           }
        }

     if (cut[dest][edge] == 0) {
        for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
           v = hg->hvertex[j];
           gain[v] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
           if (heap)
              Zoltan_Heap_Change_Value(&heap[part[v]], v, gain[v]);
           }
        }
     else if (cut[dest][edge] == 1) {
        for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
           v = hg->hvertex[j];
           if (v != vertex && part[v] == dest) {
              gain[v] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
              if (heap)
                 Zoltan_Heap_Change_Value(&heap[part[v]], v, gain[v]);
              break;
              }
           }
        }
     cut[sour][edge]--;
     cut[dest][edge]++;
     }
  return ZOLTAN_OK;
}

/****************************************************************************/
/* Serial FM 2-way refinement, latest & greatest implementation.            */
/****************************************************************************/
static int serial_fm2 (ZZ *zz,
    HGraph *hg,
    int p,
    float *part_sizes,
    Partition part,
    PHGPartParams *hgp,
    float bal_tol)
{
int    i, j, vertex, edge, *pins[2], *locked = 0, *locked_list = 0, round = 0;
double total_weight, part_weight[2], max_weight[2];
double cutsize_beforepass, best_cutsize, *gain = 0;
HEAP   heap[2];
int    steplimit;
char   *yo="serial_fm2";
int    part_dim = (hg->VtxWeightDim ? hg->VtxWeightDim : 1);
#ifdef HANDLE_ISOLATED_VERTICES    
 int    isocnt=0;
#endif
#ifdef _DEBUG
 double tw0, imbal, cutsize;
#endif

double error, best_error;
int    best_imbalance, imbalance;

  if (p != 2) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "p!=2 not allowed for local_fm2.");
     return ZOLTAN_FATAL;
     }

  if (hg->nEdge == 0)
     return ZOLTAN_OK;

  /* Calculate the weights in each partition and total, then maxima */
  part_weight[0] = 0.0;
  part_weight[1] = 0.0;
  if (hg->vwgt)  {
     for (i = 0; i < hg->nVtx; i++)
        part_weight[part[i]] += hg->vwgt[i*hg->VtxWeightDim];
     total_weight = part_weight[0] + part_weight[1];
  }
  else  {
     total_weight = (double)(hg->nVtx);
     for (i = 0; i < hg->nVtx; i++)
        part_weight[part[i]] += 1.0;
  }
  max_weight[0] = total_weight * bal_tol * part_sizes[0];
  max_weight[1] = total_weight * bal_tol * part_sizes[part_dim];

#ifdef _DEBUG
  tw0 = total_weight * part_sizes[0];
#endif
  
  if (!(pins[0]     = (int*)   ZOLTAN_CALLOC(2*hg->nEdge, sizeof(int)))
   || !(locked      = (int*)   ZOLTAN_CALLOC(hg->nVtx,    sizeof(int)))
   || !(locked_list = (int*)   ZOLTAN_CALLOC(hg->nVtx,    sizeof(int)))
   || !(gain        = (double*)ZOLTAN_CALLOC(hg->nVtx,    sizeof(double))) ) {
     Zoltan_Multifree(__FILE__,__LINE__, 4, &pins[0], &locked, &locked_list,
          &gain);
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
  }
  pins[1] = &(pins[0][hg->nEdge]);

  /* Initial calculation of the pins distribution and gain values */
  for (i = 0; i < hg->nEdge; i++)
     for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
        (pins[part[hg->hvertex[j]]][i])++;
  for (i = 0; i < hg->nVtx; i++)
     for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++) {
        edge = hg->vedge[j];
        if (pins[part[i]][edge] == 1)
           gain[i] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
        else if (pins[1-part[i]][edge] == 0)
           gain[i] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
     }

  /* Initialize the heaps and fill them with the gain values */
  Zoltan_Heap_Init(zz, &heap[0], hg->nVtx);
  Zoltan_Heap_Init(zz, &heap[1], hg->nVtx);  
  for (i = 0; i < hg->nVtx; i++)
      if (!hgp->UseFixedVtx || hg->fixed_part[i]<0) {
#ifdef HANDLE_ISOLATED_VERTICES          
          if (hg->vindex[i+1]==hg->vindex[i]) { /* isolated vertex */
              part_weight[part[i]] -= (hg->vwgt ? hg->vwgt[i*hg->VtxWeightDim] 
                                                : 1.0);
              part[i] = -(part[i]+1); /* remove those vertices from that part*/
              ++isocnt;
          } else
#endif
              Zoltan_Heap_Input(&heap[part[i]], i, gain[i]);
      }
#ifdef _DEBUG
      else {
          int pp = (hg->fixed_part[i] < hg->bisec_split) ? 0 : 1;
          if (part[i]!=pp) 
              errexit("%s: beginning of pass for hg->info=%d vertex %d is fixed at %d bisec_split is %d but its part is %d\n", uMe(hg->comm), hg->info, i, hg->fixed_part[i], hg->bisec_split, part[i]);
          
              
      }
#endif
  Zoltan_Heap_Make(&heap[0]);
  Zoltan_Heap_Make(&heap[1]);

  /* Initialize given partition as best partition */
  best_cutsize = cutsize_beforepass = Zoltan_PHG_Compute_NetCut(hg->comm, hg, part);
  best_error = MAX (part_weight[0]-max_weight[0], part_weight[1]-max_weight[1]);
  best_imbalance = (part_weight[0]>max_weight[0])||(part_weight[1]>max_weight[1]);
  do {
    int step = 0, no_better_steps = 0, number_locked = 0, best_locked = 0;
    int sour, dest;
    double cur_cutsize=best_cutsize;

    round++;
    cutsize_beforepass = best_cutsize;
    if (hgp->output_level > PHG_DEBUG_LIST)
      printf("ROUND %d:\nSTEP VERTEX  PARTS MAX_WGT CHANGE CUTSIZE\n",round);

    steplimit = (hgp->fm_max_neg_move < 0) ? hg->nVtx : hgp->fm_max_neg_move;
    /* steplimit = hg->nVtx/4;  Robsys previous choice */

    while (step < hg->nVtx && no_better_steps < steplimit) {
        step++;
        no_better_steps++;

        if (Zoltan_Heap_Empty(&heap[0]))
           sour = 1;
        else if (Zoltan_Heap_Empty(&heap[1]))
           sour = 0;
        else if (part_weight[0] > max_weight[0])
           sour = 0;
        else if (part_weight[1] > max_weight[1])
           sour = 1;
        else if (Zoltan_Heap_Max_Value(&heap[0])
              >  Zoltan_Heap_Max_Value(&heap[1]))
           sour = 0;
        else
           sour = 1;
        dest = 1-sour;
        vertex = Zoltan_Heap_Extract_Max(&heap[sour]);
        if (vertex<0)
            break;

        locked[vertex] = part[vertex] + 1;
        locked_list[number_locked++] = vertex;
        cur_cutsize -= gain[vertex];        
        
        Zoltan_HG_move_vertex (hg, vertex, sour, dest, part, pins, gain, heap);

#ifdef _DEBUG
        imbal = (tw0==0.0) ? 0.0 : (part_weight[0]-tw0)/tw0;
        uprintf(hg->comm, "%4d: SEQ moving %4d from %d to %d cut=%6.0lf bal=%.3lf\n", step, vertex, sour, dest, cur_cutsize, imbal);
        /* Just for debugging */
        cutsize = Zoltan_PHG_Compute_NetCut(hg->comm, hg, part);
        if (cur_cutsize!=cutsize) {
            errexit("%s: SEQ after move cutsize=%.2lf Verify: total=%.2lf\n", uMe(hg->comm), cur_cutsize,
                    cutsize);
        }
#endif
        
        part_weight[sour] -= (hg->vwgt ? hg->vwgt[vertex*hg->VtxWeightDim] 
                                       : 1.0);
        part_weight[dest] += (hg->vwgt ? hg->vwgt[vertex*hg->VtxWeightDim] 
                                       : 1.0);

        error = MAX (part_weight[0]-max_weight[0],part_weight[1]-max_weight[1]);
        imbalance = (part_weight[0]>max_weight[0])||(part_weight[1]>max_weight[1]);

        if ( ( best_imbalance && (error < best_error))
          || (!imbalance && (cur_cutsize < best_cutsize)))  {
            best_error   = error;
            best_imbalance = imbalance;
            best_locked  = number_locked;
            best_cutsize = cur_cutsize;
            no_better_steps = 0;
        }
        if (hgp->output_level > PHG_DEBUG_LIST+1)
           printf ("%4d %6d %2d->%2d %7.2f %f %f\n", step, vertex, sour, dest,
            error, cur_cutsize - cutsize_beforepass, cur_cutsize);
    }

#ifdef _DEBUG
    uprintf(hg->comm, "SEQ Best CUT=%6.0lf at move %d\n", best_cutsize, best_locked);
#endif
    
    /* rollback */
     while (number_locked != best_locked) {
        vertex = locked_list[--number_locked];
        sour = part[vertex];
        dest = locked[vertex] - 1;

        Zoltan_HG_move_vertex (hg, vertex, sour, dest, part, pins, gain, heap);

        part_weight[sour] -= (hg->vwgt ? hg->vwgt[vertex*hg->VtxWeightDim] 
                                       : 1.0);
        part_weight[dest] += (hg->vwgt ? hg->vwgt[vertex*hg->VtxWeightDim] 
                                       : 1.0);
        Zoltan_Heap_Input(&heap[dest], vertex, gain[vertex]);
        locked[vertex] = 0;
     }

     /* only update data structures if we're going to do another pass */
     if ((best_cutsize < cutsize_beforepass) &&  (round < hgp->fm_loop_limit)) {         
         while (number_locked) {
             vertex = locked_list[--number_locked];
             locked[vertex] = 0;
             Zoltan_Heap_Input(&heap[part[vertex]], vertex, gain[vertex]);
         }
         
         Zoltan_Heap_Make(&(heap[0]));
         Zoltan_Heap_Make(&(heap[1]));
     }
  } while ((best_cutsize < cutsize_beforepass) &&  (round < hgp->fm_loop_limit));

#ifdef HANDLE_ISOLATED_VERTICES
  if (isocnt) {
#ifdef _DEBUG      
      double isoimbalbefore, isoimbal;
#endif
      double targetw0;
      
      targetw0 = total_weight * part_sizes[0];
#ifdef _DEBUG      
      isoimbalbefore = (targetw0==0) ? 0.0 : (part_weight[0] - targetw0)/ targetw0;
#endif
      for (i=0; i < hg->nVtx; ++i)
          if (!hgp->UseFixedVtx || hg->fixed_part[i]<0) {
              if (hg->vindex[i+1]==hg->vindex[i])  { /* go over isolated vertices */
                  int npno = (part_weight[0] <  targetw0) ? 0 : 1;
                  part_weight[npno] += (hg->vwgt ? hg->vwgt[i*hg->VtxWeightDim] 
                                                 : 1.0);                
                  part[i] = npno;
              }
          }
#ifdef _DEBUG      
      isoimbal = (targetw0==0) ? 0.0 : (part_weight[0] - targetw0)/ targetw0;
      uprintf(hg->comm, "SEQ %d isolated vertices, balance before: %.3lf  after: %.3lf\n", isocnt, isoimbalbefore, isoimbal);
#endif
  }
#endif  
  
  /* gain_check (hg, gain, part, pins); */
  Zoltan_Multifree(__FILE__,__LINE__, 4, &pins[0], &locked, &locked_list, &gain);
  Zoltan_Heap_Free(&heap[0]);
  Zoltan_Heap_Free(&heap[1]);

  return ZOLTAN_OK;
}
#endif


/*****************************************************************************/
/* 2-way Parallel FM refinement. No data movement between processors,
 * just relabeling of vertices. In each FM pass we only move in one
 * direction, from the heavier partition to the lighter one. */
/*****************************************************************************/

static void fm2_move_vertex_oneway(int v, HGraph *hg, Partition part, 
                                   float *gain, HEAP *heap,
                                   int *pins[2], int *lpins[2], 
                                   double *weights, double *lweights,
                                   int *mark, int *adj)
{
    int   pno=part[v], vto=1-pno, adjsz=0, j, i;
    
    mark[v] = 1;  /* mark as moved */
    part[v] = vto;
    weights[pno] -= (hg->vwgt ? hg->vwgt[v*hg->VtxWeightDim] : 1.0);
    weights[vto] += (hg->vwgt ? hg->vwgt[v*hg->VtxWeightDim] : 1.0);
    lweights[pno] -= (hg->vwgt ? hg->vwgt[v*hg->VtxWeightDim] : 1.0);
    lweights[vto] += (hg->vwgt ? hg->vwgt[v*hg->VtxWeightDim] : 1.0);

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
        /* UVC: Heap_Change already checks if the element is
           in the heap 
           if (Zoltan_Heap_Has_Elem(&heap[p], u)) */
        Zoltan_Heap_Change_Value(&heap[p], u, gain[u]);
    }
}


static void fm2_move_vertex_oneway_nonroot(int v, HGraph *hg, Partition part, 
                                           int *lpins[2], double *lweights)
{
    int   pno=part[v], vto=1-pno, j;
    
    part[v] = vto;
    lweights[pno] -= (hg->vwgt ? hg->vwgt[v*hg->VtxWeightDim] : 1.0);
    lweights[vto] += (hg->vwgt ? hg->vwgt[v*hg->VtxWeightDim] : 1.0);

    for (j = hg->vindex[v]; j < hg->vindex[v+1]; j++) {
        int n = hg->vedge[j];
    
        --lpins[pno][n];
        ++lpins[vto][n];
    }
}



static int refine_fm2 (ZZ *zz,
                       HGraph *hg,
                       int p,
                       float *part_sizes,
                       Partition part,
                       PHGPartParams *hgp,
                       float bal_tol
    )
{
    int    i, j, ierr=ZOLTAN_OK, *pins[2]={NULL,NULL}, *lpins[2]={NULL,NULL};
    int    *moves=NULL, *mark=NULL, *adj=NULL, passcnt=0;
    float  *gain=NULL, *lgain=NULL;
    int    best_cutsizeat, cont, successivefails=0;
    double total_weight, weights[2], total_lweight, lweights[2], lwadjust[2],
        max_weight[2], lmax_weight[2], avail[2], gavail[2];
    int availcnt[2], gavailcnt[2];
    double targetw0, ltargetw0, minvw=DBL_MAX;
    double cutsize, best_cutsize, 
        best_limbal, imbal, limbal;
    HEAP   heap[2];
    char   *yo="refine_fm2";
    int    part_dim = (hg->VtxWeightDim ? hg->VtxWeightDim : 1);
#ifdef HANDLE_ISOLATED_VERTICES    
    int    isocnt=hg->nVtx; /* only root uses isocnt, isolated vertices
                               are kept at the end of moves array */
    int    *deg=NULL, *ldeg=NULL;
#if 0
    double best_imbal;
#endif
#endif
    PHGComm *hgc=hg->comm;
    int rootRank;
    
    struct phg_timer_indices *timer = Zoltan_PHG_LB_Data_timers(zz);
    int do_timing = (hgp->use_timers > 2);
    int detail_timing = (hgp->use_timers > 3);

    ZOLTAN_TRACE_ENTER(zz, yo);

    if (p != 2) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "p!=2 not allowed for refine_fm2.");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_FATAL;
    }

    /* return only if globally there is no edge or vertex */
    if (!hg->dist_y[hgc->nProc_y] || hg->dist_x[hgc->nProc_x] == 0) {
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_OK;
    }


#ifdef USE_SERIAL_REFINEMENT_ON_ONE_PROC
    if (hgc->nProc==1){ /* only one proc? use serial code */
        ZOLTAN_TRACE_EXIT(zz, yo);
        return serial_fm2 (zz, hg, p, part_sizes, part, hgp, bal_tol);
    }
#endif

    if (do_timing) { 
        if (timer->rfrefine < 0) 
            timer->rfrefine = Zoltan_Timer_Init(zz->ZTime, 1, "Ref_P_Total");
        ZOLTAN_TIMER_START(zz->ZTime, timer->rfrefine, hgc->Communicator);
    }
    if (detail_timing) {
        if (timer->rfpins < 0) 
            timer->rfpins = Zoltan_Timer_Init(zz->ZTime, 0, "Ref_P_Pins");
        if (timer->rfiso < 0) 
            timer->rfiso = Zoltan_Timer_Init(zz->ZTime, 0, "Ref_P_IsolatedVert");
        if (timer->rfgain < 0) 
            timer->rfgain = Zoltan_Timer_Init(zz->ZTime, 0, "Ref_P_Gain");
        if (timer->rfheap < 0) 
            timer->rfheap = Zoltan_Timer_Init(zz->ZTime, 0, "Ref_P_Heap");
        if (timer->rfpass < 0) 
            timer->rfpass = Zoltan_Timer_Init(zz->ZTime, 0, "Ref_P_Pass");
        if (timer->rfroll < 0) 
            timer->rfroll = Zoltan_Timer_Init(zz->ZTime, 0, "Ref_P_Roll");
        if (timer->rfnonroot < 0) 
            timer->rfnonroot = Zoltan_Timer_Init(zz->ZTime, 0, "Ref_P_NonRoot");
    }
    
    
    /* find the index of the proc in column group with 
       the most #nonzeros; it will be our root
       proc for computing moves since it has better 
       knowedge about global hypergraph.
       We ignore returned #pins (i) in root */
    Zoltan_PHG_Find_Root(hg->nPins, hgc->myProc_y, hgc->col_comm, 
                         &i, &rootRank);
    
    /* Calculate the weights in each partition and total, then maxima */
    weights[0] = weights[1] = 0.0;
    lweights[0] = lweights[1] = 0.0;
    if (hg->vwgt) 
        for (i = 0; i < hg->nVtx; i++) {
            lweights[part[i]] += hg->vwgt[i*hg->VtxWeightDim];
            minvw = (minvw > hg->vwgt[i*hg->VtxWeightDim]) 
                  ? hg->vwgt[i*hg->VtxWeightDim] 
                  : minvw;
        }
    else {
        minvw = 1.0;
        for (i = 0; i < hg->nVtx; i++)
            lweights[part[i]] += 1.0;
    }

    MPI_Allreduce(lweights, weights, 2, MPI_DOUBLE, MPI_SUM, hgc->row_comm);
    total_weight = weights[0] + weights[1];
    targetw0 = total_weight * part_sizes[0]; /* global target weight for part 0 */

    max_weight[0] = total_weight * bal_tol * part_sizes[0];
    max_weight[1] = total_weight * bal_tol * part_sizes[part_dim]; /* should be (1 - part_sizes[0]) */


    if (weights[0]==0.0) {
        ltargetw0 = targetw0 / hgc->nProc_x;
        lmax_weight[0] = max_weight[0] / hgc->nProc_x;
    } else {
        lmax_weight[0] = (weights[0]==0.0) ? 0.0 : lweights[0] +
            (max_weight[0] - weights[0]) * ( lweights[0] / weights[0] );
        ltargetw0 = targetw0 * ( lweights[0] / weights[0] ); /* local target weight */
    }
    if (weights[1]==0.0)
        lmax_weight[1] = max_weight[1] / hgc->nProc_x;
    else
        lmax_weight[1] = (weights[1]==0.0) ? 0.0 : lweights[1] +
            (max_weight[1] - weights[1]) * ( lweights[1] / weights[1] );

    total_lweight = lweights[0]+lweights[1];
    
    avail[0] = MAX(0.0, lmax_weight[0]-total_lweight);
    avail[1] = MAX(0.0, lmax_weight[1]-total_lweight);
    availcnt[0] = (avail[0] == 0) ? 1 : 0;
    availcnt[1] = (avail[1] == 0) ? 1 : 0; 
    MPI_Allreduce(avail, gavail, 2, MPI_DOUBLE, MPI_SUM, hgc->row_comm);
    MPI_Allreduce(availcnt, gavailcnt, 2, MPI_INT, MPI_SUM, hgc->row_comm);

#ifdef _DEBUG
    if (gavailcnt[0] || gavailcnt[1])
        uprintf(hgc, "before adjustment, LMW[%.1lf, %.1lf]\n", lmax_weight[0], lmax_weight[1]);
#endif

    if (gavailcnt[0]) 
        lmax_weight[0] += gavail[0] / (double) gavailcnt[0];
    
    if (gavailcnt[1]) 
        lmax_weight[1] += gavail[1] / (double) gavailcnt[1];
    
    /* Our strategy is to stay close to the current local weight balance.
       We do not need the same local balance on each proc, as long as
       we achieve approximate global balance.                            */

#ifdef _DEBUG
    imbal = (targetw0==0.0) ? 0.0 : fabs(weights[0]-targetw0)/targetw0;
    limbal = (ltargetw0==0.0) ? 0.0 : fabs(lweights[0]-ltargetw0)/ltargetw0;
    uprintf(hgc, "H(%d, %d, %d), FM2: W[%.1lf, %.1lf] MW:[%.1lf, %.1lf] I=%.3lf  LW[%.1lf, %.1lf] LMW[%.1lf, %.1lf] LI=%.3lf\n", hg->nVtx, hg->nEdge, hg->nPins, weights[0], weights[1], max_weight[0], max_weight[1], imbal, lweights[0], lweights[1], lmax_weight[0], lmax_weight[1], limbal);
#endif

    
    if ((hg->nEdge && (!(pins[0]    = (int*) ZOLTAN_MALLOC(2 * hg->nEdge * sizeof(int)))
                      || !(lpins[0] = (int*) ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int))))) ||
        (hg->nVtx && (!(moves   = (int*)   ZOLTAN_MALLOC(hg->nVtx * sizeof(int)))
                     || !(lgain = (float*) ZOLTAN_MALLOC(hg->nVtx * sizeof(float))))))
        MEMORY_ERROR;

    if (hg->nEdge) {
        pins[1] = &(pins[0][hg->nEdge]);
        lpins[1] = &(lpins[0][hg->nEdge]);
    }

    if (hgc->myProc_y==rootRank) { /* only root needs mark, adj, gain and heaps*/
        Zoltan_Heap_Init(zz, &heap[0], hg->nVtx);
        Zoltan_Heap_Init(zz, &heap[1], hg->nVtx);  
        if (hg->nVtx &&
            (!(mark     = (int*)   ZOLTAN_CALLOC(hg->nVtx, sizeof(int)))
             || !(adj   = (int*)   ZOLTAN_MALLOC(hg->nVtx * sizeof(int)))   
             || !(gain  = (float*) ZOLTAN_MALLOC(hg->nVtx * sizeof(float)))))
            MEMORY_ERROR;
    }

    /* Initial calculation of the local pin distribution (sigma in UVC's papers)  */
    if (detail_timing)         
        ZOLTAN_TIMER_START(zz->ZTime, timer->rfpins, hgc->Communicator);                        
    for (i = 0; i < hg->nEdge; ++i)
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j){
            ++(lpins[part[hg->hvertex[j]]][i]);
        }
    if (detail_timing)         
        ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfpins, hgc->Communicator);                    
    

#ifdef HANDLE_ISOLATED_VERTICES        
    /* first compute vertex degree to find any isolated vertices
       we use lgain and gain, as ldeg, deg.*/
    if (hg->nVtx) {
        if (detail_timing)         
            ZOLTAN_TIMER_START(zz->ZTime, timer->rfiso, hgc->Communicator);        
        ldeg = (int *) lgain;
        deg = (int *) gain; /* null for non-root but that is fine */
        for (i = 0; i < hg->nVtx; ++i)
            ldeg[i] = hg->vindex[i+1] - hg->vindex[i];
        MPI_Reduce(ldeg, deg, hg->nVtx, MPI_INT, MPI_SUM, rootRank,
                   hg->comm->col_comm);

        if (hgc->myProc_y==rootRank) { /* root marks isolated vertices */
            for (i=0; i<hg->nVtx; ++i)
                if (!hgp->UseFixedVtx || hg->fixed_part[i]<0) {
                    if (!deg[i]) {
                        moves[--isocnt] = i;
                        part[i] = -(part[i]+1); /* remove those vertices from that part*/
                    }
                }
        }   
        if (detail_timing)         
            ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfiso, hgc->Communicator);        

    }
#endif
    
    do {
        int v=1, movecnt=0, neggaincnt=0, from, to;
        int maxneggain = (hgp->fm_max_neg_move < 0) ? hg->nVtx : hgp->fm_max_neg_move;
        int notfeasible=(weights[0]>max_weight[0]) || (weights[1]>max_weight[1]);
    
        /* now compute global pin distribution */
        if (hg->nEdge) {
            if (detail_timing)         
                ZOLTAN_TIMER_START(zz->ZTime, timer->rfpins, hgc->Communicator);                    
            MPI_Allreduce(lpins[0], pins[0], 2*hg->nEdge, MPI_INT, MPI_SUM, 
                          hgc->row_comm);
            if (detail_timing)         
                ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfpins, hgc->Communicator);                    
        }

        /* now we can compute actual cut */
        best_cutsizeat=0;
        cutsize = 0.0;
        for (i=0; i < hg->nEdge; ++i) {
            if (pins[0][i] && pins[1][i])
                cutsize += (hg->ewgt ? hg->ewgt[i] : 1.0);
        }
        MPI_Allreduce(&cutsize, &best_cutsize, 1, MPI_DOUBLE, MPI_SUM, hgc->col_comm);
        cutsize = best_cutsize;

        imbal = (targetw0==0.0) ? 0.0 : fabs(weights[0]-targetw0)/targetw0;        
        best_limbal = limbal = (ltargetw0==0.0) ? 0.0
            : fabs(lweights[0]-ltargetw0)/ltargetw0;

        /* UVCUVC: it looks like instead of moving always from overloaded
           part, alternating the 'from' part gives better results.
           Hence if the imbal is not really bad (2x worse) we use that approach  */
        if (imbal > BADBALANCE*(bal_tol-1.0) ) /* decide which way the moves will be in this pass */
            from = (weights[0] < targetw0) ? 1 : 0;
        else 
            from = passcnt % 2; 
        /* we want to be sure that everybody!!! picks the same source */
        MPI_Bcast(&from, 1, MPI_INT, 0, hgc->Communicator); 

        to = 1-from;
        
#ifdef _DEBUG
        /* Just for debugging */
        best_cutsize = Zoltan_PHG_Compute_NetCut(hgc, hg, part);
        if (best_cutsize!=cutsize) {
            errexit("%s: Initial cutsize=%.2lf Verify: total=%.2lf\n", uMe(hgc), cutsize,
                    best_cutsize);
        }
        if (hgc->myProc_y==rootRank)
            for (i = 0; i< hg->nVtx; ++i)
                if (mark[i])
                    errexit("mark[%d]=%d", i, mark[i]);
        /* debuggging code ends here */
#endif

        /* compute only the gains of the vertices from 'from' part */
        if (detail_timing)         
            ZOLTAN_TIMER_START(zz->ZTime, timer->rfgain, hgc->Communicator);                    
        
        for (i = 0; i < hg->nVtx; ++i) {
            lgain[i] = 0.0;
            if ((part[i]==from) && (!hgp->UseFixedVtx || hg->fixed_part[i]<0))
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
        if (hg->nVtx)
            MPI_Reduce(lgain, gain, hg->nVtx, MPI_FLOAT, MPI_SUM, rootRank, 
                       hgc->col_comm);
        if (detail_timing)         
            ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfgain, hgc->Communicator);                    
        

        if (hgp->output_level >= PHG_DEBUG_ALL) {
            imbal = (targetw0==0.0) ? 0.0 : fabs(weights[0]-targetw0)/targetw0;
            printf("%s FM Pass %d (%d->%d) Cut=%.2f W[%5.0f, %5.0f] I= %.2f LW[%5.0f, %5.0f] LI= %.2f\n", uMe(hgc), passcnt, from, to, cutsize, weights[0], weights[1], imbal, lweights[0], lweights[1], limbal);
        }

        if (hgc->myProc_y==rootRank) {
            /* those are the lucky ones; each proc in column-group
               could have compute the same moves concurrently; but for this
               version we'll do it in the root procs and broadcast */

#ifdef HANDLE_ISOLATED_VERTICES
            if (detail_timing)         
                ZOLTAN_TIMER_START(zz->ZTime, timer->rfiso, hgc->Communicator);                    
            lwadjust[0] = lwadjust[1] = 0.0;
            for (i=isocnt; i < hg->nVtx; ++i) { /* go over isolated vertices */
                int   u=moves[i], pno=-part[u]-1;
                float w=(hg->vwgt ? hg->vwgt[u*hg->VtxWeightDim] : 1.0);

                if (pno<0 || pno>1)
                    errexit("heeeey pno=%d", pno);
                /* let's remove it from its part */
                lwadjust[pno] -= w;                
            }
            lweights[0] += lwadjust[0];
            lweights[1] += lwadjust[1];
            if (detail_timing)         
                ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfiso, hgc->Communicator);                    
#endif

            if (detail_timing)         
                ZOLTAN_TIMER_START(zz->ZTime, timer->rfheap, hgc->Communicator);                    
            
            /* Initialize the heaps and fill them with the gain values */
            Zoltan_Heap_Clear(&heap[from]);  
            for (i = 0; i < hg->nVtx; ++i)
                if ((part[i]==from) && (!hgp->UseFixedVtx || hg->fixed_part[i]<0))
                    Zoltan_Heap_Input(&heap[from], i, gain[i]);
            Zoltan_Heap_Make(&heap[from]);
            if (detail_timing) {
                ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfheap, hgc->Communicator);
                ZOLTAN_TIMER_START(zz->ZTime, timer->rfpass, hgc->Communicator);
            }

            while ((neggaincnt < maxneggain) && ((lweights[to]+minvw) <= lmax_weight[to]) ) {
                if (Zoltan_Heap_Empty(&heap[from])) { /* too bad it is empty */
                    v = -1;
                    break;
                }
                
                v = Zoltan_Heap_Extract_Max(&heap[from]);    
                
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
                if (lweights[to]+((hg->vwgt)?hg->vwgt[v*hg->VtxWeightDim]:1.0) > lmax_weight[to]) {
#ifdef _DEBUG2                    
                    printf("%s %4d: %6d (g: %5.1lf), p:%2d [%4.0lf, %4.0lf] NF\n", uMe(hgc), movecnt, v, gain[v], from, weights[0], weights[1]);
#endif
                    /* Negative value in moves array means we have examined 
                       the vertex but couldn't move it. Note offset by one,
                       otherwise zero would be ambiguous. */
                    moves[movecnt++] = -(v+1);
                    continue;
                } 

                    
                moves[movecnt] = v;
                ++neggaincnt;
                cutsize -= gain[v];

                fm2_move_vertex_oneway(v, hg, part, gain, heap, pins, lpins, weights, lweights, mark, adj);
                imbal = (targetw0==0.0) ? 0.0
                    : fabs(weights[0]-targetw0)/targetw0;
                limbal = (ltargetw0==0.0) ? 0.0
                    : fabs(lweights[0]-ltargetw0)/ltargetw0;

                if (notfeasible || (cutsize<best_cutsize) ||
                                   (cutsize==best_cutsize && limbal < best_limbal)) {
#ifdef _DEBUG2                    
                    printf("%s %4d: %6d (g: %5.1lf), p:%2d W[%4.0lf, %4.0lf] I:%.2lf LW[%4.0lf, %4.0lf] LI:%.2lf C:%.1lf<-- Best\n", uMe(hgc), movecnt, v, gain[v], from, weights[0], weights[1], imbal, lweights[0], lweights[1], limbal, cutsize); /* after move gain is -oldgain */
#endif
                    notfeasible = weights[from]>max_weight[from];
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

            
            if (detail_timing) {
                ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfpass, hgc->Communicator);
                ZOLTAN_TIMER_START(zz->ZTime, timer->rfroll, hgc->Communicator);
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
                int vv = moves[i];
                if (vv<0)
                    vv = -vv-1;
                else /* we don't need to roll pins, or weights etc; rolling local ones suffices */
                    fm2_move_vertex_oneway_nonroot(vv, hg, part, lpins, lweights);
                mark[vv] = 0;
            }
            for (i=0; i<best_cutsizeat; ++i){
                int vv = (moves[i] < 0 ) ? -moves[i] - 1 : moves[i];
                mark[vv] = 0;
            }
            if (detail_timing) 
                ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfroll, hgc->Communicator);            
        }

        if (detail_timing) 
            ZOLTAN_TIMER_START(zz->ZTime, timer->rfnonroot, hgc->Communicator);            
        
        /* now root bcast moves to column procs */
        MPI_Bcast(&best_cutsizeat, 1, MPI_INT, rootRank, hgc->col_comm);
        MPI_Bcast(moves, best_cutsizeat, MPI_INT, rootRank, hgc->col_comm);
        if (hgc->myProc_y!=rootRank) { /* now non-root does move simulation */
            for (i=0; i<best_cutsizeat; ++i) {
                int vv = moves[i];
                if (vv>=0)
                    fm2_move_vertex_oneway_nonroot(vv, hg, part, lpins, lweights);
            }
        }
        if (detail_timing) 
            ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfnonroot, hgc->Communicator);            

        
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


#ifdef HANDLE_ISOLATED_VERTICES
        if (detail_timing)         
            ZOLTAN_TIMER_START(zz->ZTime, timer->rfiso, hgc->Communicator);        
        
#if 0
        MPI_Allreduce(lweights, weights, 2, MPI_DOUBLE, MPI_SUM, hgc->row_comm);        
        best_imbal = (targetw0==0.0) ? 0.0 : fabs(weights[0]-targetw0)/targetw0;
        if (hgc->myProc_y==rootRank)             
            uprintf(hgc, "BEFORE ISOLATED VERTEX HANDLING WE *THINK* GLOBAL IMBALANCE is %.3lf\n", best_imbal);
#endif
        
        if (hgc->myProc_y==rootRank) {
            best_limbal = (ltargetw0==0.0) ? 0.0
                : fabs(lweights[0]-ltargetw0)/ltargetw0;
            
            for (i=isocnt; i < hg->nVtx; ++i) { /* go over isolated vertices */
                int u = moves[i], npno;
                float w=(hg->vwgt ? hg->vwgt[u*hg->VtxWeightDim] : 1.0);

                npno = (lweights[0] < ltargetw0) ? 0 : 1;
                lweights[npno] += w;
                lwadjust[npno] += w;
                part[u] = -(npno+1); /* move to npno (might be same as pno;
                                        so it may not be a real move */
            }
            limbal = (ltargetw0==0.0) ? 0.0
                : fabs(lweights[0]-ltargetw0)/ltargetw0;
#if 0           
            uprintf(hgc, "before binpacking of %d isolated vertices balance was: %.3lf now: %.3lf\n", hg->nVtx-isocnt, best_limbal, limbal);
#endif
        }

        MPI_Bcast(lwadjust, 2, MPI_DOUBLE, rootRank, hgc->col_comm);
        if (hgc->myProc_y!=rootRank) {
            lweights[0] += lwadjust[0];
            lweights[1] += lwadjust[1];
        }
        if (detail_timing)         
            ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfiso, hgc->Communicator);                
#endif        
        
        MPI_Allreduce(lweights, weights, 2, MPI_DOUBLE, MPI_SUM, hgc->row_comm);
#if 0       
        best_imbal = (targetw0==0.0) ? 0.0 : fabs(weights[0]-targetw0)/targetw0;
        if (hgc->myProc_y==rootRank)             
            uprintf(hgc, "NEW GLOBAL IMBALANCE is %.3lf\n", best_imbal);
#endif
        
        if (weights[0]==0.0) 
            ltargetw0 = lmax_weight[0] = 0.0;
        else {
            lmax_weight[0] = lweights[0] +
                (max_weight[0] - weights[0]) * ( lweights[0] / weights[0] );
            ltargetw0 = targetw0 * ( lweights[0] / weights[0] ); /* local target weight */
        }
        lmax_weight[1] = (weights[1]==0.0) ? 0.0 : lweights[1] +
            (max_weight[1] - weights[1]) * ( lweights[1] / weights[1] );
        
        cont = 0;
        MPI_Allreduce(&best_cutsizeat, &cont, 1, MPI_INT, MPI_LOR, hgc->row_comm);

        /* since we're only moving in one direction; make sure two successive
           pass didn't produce any improvement before terminating */
        if (!cont)
            ++successivefails; 
        else
            successivefails = 0; 
#ifdef _DEBUG
        /* Just for debugging */
        best_cutsize = Zoltan_PHG_Compute_NetCut(hgc, hg, part);
        imbal = (targetw0 == 0.0) ? 0.0 : fabs(weights[0]-targetw0)/targetw0;
        printf("%s End of Pass %d Comp.Cut=%.2lf RealCut=%.2lf W[%5.0lf, %5.0lf] Imbal=%.2lf\n", uMe(hgc), passcnt, cutsize, best_cutsize, weights[0], weights[1], imbal);
        /* debuggging code ends here */
#endif
    } while (successivefails<2 &&  (++passcnt < hgp->fm_loop_limit));


#ifdef HANDLE_ISOLATED_VERTICES
    if (detail_timing)         
        ZOLTAN_TIMER_START(zz->ZTime, timer->rfiso, hgc->Communicator);            
    /* now root sneds the final part no's of isolated vertices; if any */
    MPI_Bcast(&isocnt, 1, MPI_INT, rootRank, hgc->col_comm);
    if (isocnt<hg->nVtx) {
        deg = (int *) lgain; /* we'll use for part no's of isolated vertices */
        if (hgc->myProc_y==rootRank) 
            for (i=isocnt; i < hg->nVtx; ++i) { /* go over isolated vertices */
                int u = moves[i];
                deg[i] = part[u] = -part[u]-1; 
            }
            
        MPI_Bcast(&moves[isocnt], hg->nVtx-isocnt, MPI_INT, rootRank, hgc->col_comm);
        MPI_Bcast(&deg[isocnt], hg->nVtx-isocnt, MPI_INT, rootRank, hgc->col_comm);
        if (hgc->myProc_y!=rootRank) 
            for (i=isocnt; i < hg->nVtx; ++i)  /* go over isolated vertices */
                part[moves[i]] = deg[i];
    }
    if (detail_timing)         
        ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfiso, hgc->Communicator);            
#endif
 End:    

    if (hgc->myProc_y==rootRank) { /* only root needs mark, adj, gain and heaps*/        
        Zoltan_Multifree(__FILE__,__LINE__, 3, &mark, &adj, &gain);
        Zoltan_Heap_Free(&heap[0]);
        Zoltan_Heap_Free(&heap[1]);        
    }
    
    Zoltan_Multifree(__FILE__, __LINE__, 4, &pins[0], &lpins[0], &moves, &lgain);

    if (do_timing) 
        ZOLTAN_TIMER_STOP(zz->ZTime, timer->rfrefine, hgc->Communicator);
    
    
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ierr;
}





#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
