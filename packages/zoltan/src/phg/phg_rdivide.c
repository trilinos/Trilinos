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


#include "phg.h"
#include "phg_tree.h"
#include "phg_distrib.h"
#include "zz_util_const.h"

/*
#define _DEBUG1
#define _DEBUG2
#define _DEBUG3
*/

/* if you want to disable processor split, set PHG_PROC_SPLIT to 0 */
#define PHG_PROC_SPLIT 1


static int split_hypergraph(int *pins[2], HGraph*, HGraph*, PHGPartParams*,
                            Partition, int,
                            ZZ*, double *, double *, int connectivitycut);


static int rdivide_and_prepsend(int, int, Partition, ZZ *, HGraph *,
                                PHGPartParams *, int, int *, int *, int *,
                                int *, int *, int, int);
static float balanceTol(PHGPartParams *hgp, int part_dim, int pno, float *ratios,
                        float tot, float pw);

/* Recursively divides both the problem and the processes (if enabled)
   into 2 parts until all parts are found */
int Zoltan_PHG_rdivide(
  int lo,                 /* Lowest partition number to be found */
  int hi,                 /* Highest partition number to be found */
  Partition final,        /* Input:  initial partition assignments for vtxs;
                             Output:  computed partition assignments. */
  ZZ *zz, 
  HGraph *hg,
  PHGPartParams *hgp, 
  int level,
  int father
)
{
  char *yo = "Zoltan_PHG_rdivide";
  int i, j, mid, ierr=ZOLTAN_OK, leftend, rightstart; 
  int *pins[2] = {NULL,NULL}, *lpins[2] = {NULL,NULL};
  Partition part=NULL;
  HGraph *left=NULL, *right=NULL;
  int    *proclist=NULL, *sendbuf=NULL, *recvbuf=NULL, nsend, msg_tag=7777;
  PHGComm *hgc = hg->comm;
  int nVtx = hg->nVtx;
  ZOLTAN_GNO_TYPE gnVtx = hg->dist_x[hgc->nProc_x];
  double leftw=0.0, rightw=0.0;
  float  bal_tol = hgp->bal_tol;
  float  *bisec_part_sizes=NULL;   /* Target partition sizes; dimension is 2*hg->VtxWeightDim  
                                        because we are doing bisection */
  int part_dim = hg->VtxWeightDim ? hg->VtxWeightDim : 1;
  struct phg_timer_indices *timer = Zoltan_PHG_LB_Data_timers(zz);
  int do_timing = (hgp->use_timers > 1);
  int detail_timing = (hgp->use_timers > 3);

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (hgp->keep_tree) Zoltan_PHG_Tree_Set(zz, father, lo, hi);

  if (!gnVtx) { /* UVC: no vertex; no need for recursion!? */
      if (level>0)
          Zoltan_HG_HGraph_Free(hg);

      return ierr;
  }
  
  if (do_timing) { 
    if (timer->rdrdivide < 0) 
      timer->rdrdivide = Zoltan_Timer_Init(zz->ZTime, 1, "Rdivide");
    ZOLTAN_TIMER_START(zz->ZTime, timer->rdrdivide, hgc->Communicator);
  }
  if (detail_timing) {
    if (timer->rdbefore < 0) 
      timer->rdbefore = Zoltan_Timer_Init(zz->ZTime, 0, "Rdivide_BefPart");
    ZOLTAN_TIMER_START(zz->ZTime, timer->rdbefore, hgc->Communicator);
    if (timer->rdafter < 0) 
      timer->rdafter = Zoltan_Timer_Init(zz->ZTime, 0, "Rdivide_AftPart");
    if (timer->rdsplit < 0) 
      timer->rdsplit = Zoltan_Timer_Init(zz->ZTime, 0, "Rdivide_Split");
    if (timer->rdredist < 0) 
      timer->rdredist = Zoltan_Timer_Init(zz->ZTime, 0, "Rdivide_Redist");
    if (timer->rdsend < 0) 
      timer->rdsend = Zoltan_Timer_Init(zz->ZTime, 0, "Rdivide_Send");
    if (timer->rdwait < 0) 
      timer->rdwait = Zoltan_Timer_Init(zz->ZTime, 0, "Rdivide_LoadBalWait");
  }
  hg->redl = hgp->redl;
  
  /* only one part remaining, record results and exit */
  if (lo == hi) {
    for (i = 0; i < hg->nVtx; ++i)
      final[hg->vmap[i]] = lo;
    goto End;
  }

  if (hg->nVtx && !(part = (Partition) ZOLTAN_MALLOC (hg->nVtx * sizeof (int))))
      MEMORY_ERROR;
  for (i = 0; i < hg->nVtx; i++)
    part[i] = final[i];

  /* bipartition current hypergraph with appropriate split ratio */
  if (!(bisec_part_sizes = (float *) ZOLTAN_MALLOC(2*part_dim*sizeof(float))))
      MEMORY_ERROR;
  mid = (lo+hi)/2;
  for (i=0; i<2*part_dim; ++i) bisec_part_sizes[i] = 0.;
  for (i = lo; i <= mid; ++i)
      for (j=0; j<part_dim; ++j)
          bisec_part_sizes[j] += hgp->part_sizes[i*part_dim+j];
  for (i = lo; i <= hi;  ++i)
      for (j=0; j<part_dim; ++j)
          bisec_part_sizes[part_dim+j] += hgp->part_sizes[i*part_dim+j];
  for (j=0; j<part_dim; ++j) {
      bisec_part_sizes[j] = (float) ((double) bisec_part_sizes[j] / (double) bisec_part_sizes[part_dim+j]);
      bisec_part_sizes[part_dim+j] = 1. - bisec_part_sizes[j];
  }


/*  uprintf(hgc, "RDIVIDE: [%d, %d] mid=%d  target sizes: [%.3lf, %.3lf]\n", lo, hi, mid, bisec_part_sizes[0], bisec_part_sizes[1]);  */
  
  if (hgp->bal_tol_adjustment>1.0) {
      float q = (float) ceil(log((double)1+hi-lo) / log(2.0));
      /* uprintf(hgc, "for k=%d q=%.1f\n", 1+hi-lo, q);*/
      hgp->bal_tol = MIN(bal_tol,
                         1.0 + hgp->bal_tol_adjustment*(bal_tol-1.0)/q);
  } else 
      hgp->bal_tol = (hi==lo+1) ? bal_tol 
                                : 1.0 + hgp->bal_tol_adjustment*(bal_tol-1.0);

  if (do_timing)  /* Don't include partitioning time in rdivide */
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdrdivide, hgc->Communicator);
  if (detail_timing)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdbefore, hgc->Communicator);

  /*uprintf(hgc, "OLD MAxImbal: %.3f   New MaxImbal: %.3f\n", bal_tol, hgp->bal_tol);*/
  if (hgp->UseFixedVtx || hgp->UsePrefPart)
      hg->bisec_split = mid+1;

  ierr = Zoltan_PHG_Partition (zz, hg, 2, bisec_part_sizes, part, hgp);

  if (do_timing)  /* Restart rdivide timer */
    ZOLTAN_TIMER_START(zz->ZTime, timer->rdrdivide, hgc->Communicator);
  if (detail_timing)
    ZOLTAN_TIMER_START(zz->ZTime, timer->rdafter, hgc->Communicator);

  hgp->bal_tol = bal_tol;
  if (ierr != ZOLTAN_OK)
      goto End;

  if (hgp->output_level)
    uprintf(hgc, "Rdivide(%d, %d): %.1lf\n", lo, hi, 
                 Zoltan_PHG_Compute_ConCut(hgc, hg, part, 2, &ierr));
    
  /* if only two parts total, record results and exit */
  if (lo + 1 == hi)  {
    if (hgp->keep_tree) {
      Zoltan_PHG_Tree_Set(zz, 2*father, lo, lo);
      Zoltan_PHG_Tree_Set(zz, 2*father+1, hi, hi);
    }
    for (i = 0; i < hg->nVtx; ++i)
      final[hg->vmap[i]] = ((part[i] == 0) ? lo : hi);
    ZOLTAN_FREE (&part);
    if (detail_timing) 
      ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdafter, hgc->Communicator);
    goto End;
  }

  if (hg->nEdge && (!(pins[0] = (int*) ZOLTAN_CALLOC(2*hg->nEdge, sizeof(int)))
   || !(lpins[0] = (int*) ZOLTAN_CALLOC (2 * hg->nEdge, sizeof(int)))))
      MEMORY_ERROR;
  if (pins[0] && lpins[0]) {
      pins[1]  = &( pins[0][hg->nEdge]);
      lpins[1] = &(lpins[0][hg->nEdge]);
  }
     
  /* Initial calculation of the local pin distribution (sigma in UVC's papers)*/
  for (i = 0; i < hg->nEdge; ++i)
      for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j)
        ++(lpins[part[hg->hvertex[j]]][i]);
        
  /* now compute global pin distribution */
  if (hg->nEdge)
      MPI_Allreduce(lpins[0], pins[0], 2*hg->nEdge, MPI_INT, MPI_SUM, 
                    hgc->row_comm);
  ZOLTAN_FREE (&lpins[0]);   /* we don't need lpins anymore */
    
  if (detail_timing) {
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdafter, hgc->Communicator);
    ZOLTAN_TIMER_START(zz->ZTime, timer->rdsplit, hgc->Communicator);
  }
  
  /* recursively divide in two parts and repartition hypergraph */
  if (mid>lo) { /* only split if we really need it */
      if (!(left = (HGraph*) ZOLTAN_MALLOC (sizeof (HGraph))))
          MEMORY_ERROR;
      
      ierr = split_hypergraph (pins, hg, left, hgp, part, 0, zz, 
                               &leftw, &rightw, hgp->connectivity_cut);
      if (ierr != ZOLTAN_OK) 
          goto End;
      if (!left->dist_x[hgc->nProc_x]) { /* left is empty */
          Zoltan_HG_HGraph_Free (left);
          ZOLTAN_FREE(&left);
      }          
  } else {
      for (i = 0; i < hg->nVtx; ++i)
          if (part[i]==0)
              final[hg->vmap[i]] = lo;
      /* No recursion for the tree */
      if (hgp->keep_tree) Zoltan_PHG_Tree_Set(zz, 2*father, lo, lo);
  }

  if (hi>mid+1) { /* only split if we need it */
      if (!(right = (HGraph*) ZOLTAN_MALLOC (sizeof (HGraph))))
          MEMORY_ERROR;
      ierr = split_hypergraph (pins, hg, right, hgp, part, 1, zz, 
                               &rightw, &leftw, hgp->connectivity_cut);
  
      if (ierr != ZOLTAN_OK)
          goto End;
      if (!right->dist_x[hgc->nProc_x]) { /* right is empty */
          Zoltan_HG_HGraph_Free (right);
          ZOLTAN_FREE(&right);
      }
  } else {
      for (i = 0; i < hg->nVtx; ++i)
          if (part[i]==1)
              final[hg->vmap[i]] = hi;
      /* No recursion for the tree */
      if (hgp->keep_tree) Zoltan_PHG_Tree_Set(zz, 2*father+1, hi, hi);
  }

  
#ifdef _DEBUG1
  for (i=0; i<hg->nVtx; ++i)
      if (part[i]<0 || part[i]>1)
          errexit("During recursive bisection for [%d, %d] invalid partvec[%d]=%d\n", lo, hi, i, part[i]);
  if (left || right) {
      double ltotw=0, totw=0.0, imbal, targetw0;
      for (i=0; i<hg->nVtx; ++i)
          ltotw += hg->vwgt[i*hg->VtxWeightDim];
      MPI_Allreduce(&ltotw, &totw, 1, MPI_DOUBLE, MPI_SUM, hgc->row_comm);
      targetw0=totw*bisec_part_sizes[0];
      imbal = (targetw0==0.0) ? 0.0 : fabs(leftw-targetw0)/targetw0;      
      uprintf(hgc, "Total Vertex Weight=%.2lf   LeftW=%.2lf   RightW=%.2lf  tps=(%.2lf, %.2lf)  ib=%.2lf\n", totw, leftw, rightw, bisec_part_sizes[0], bisec_part_sizes[1], imbal);
      if ((totw - (leftw+rightw)) > 0.000001)
          errexit("During recursive bisection for [%d, %d] totw (%.2lf) !=  Left(%.2lf) + Right(%.2lf)\n", lo, hi, totw, leftw, rightw);
  }
#endif

  if (level>0) { /* some output functions in phg depends on hg being exist for
                    top level hypergraph hence we only free the rest after we split
                    because in theory we shouldn't need them */      
      Zoltan_HG_HGraph_Free(hg);
      hg = NULL;
  }  
  Zoltan_Multifree (__FILE__, __LINE__, 2, &pins[0], &part);

  
  if (detail_timing) 
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdsplit, hgc->Communicator);


#if PHG_PROC_SPLIT 
  if ((hgc->nProc>1) && left && right) {
      /* redistribute left and right parts */
      leftend = (int)((float) (hgc->nProc-1) 
                    * (float) left->dist_x[hgc->nProc_x] 
                    / (float) gnVtx);
      if ((leftend+1)>=hgc->nProc) /* just to be sure :) */
          leftend =  hgc->nProc-1;    
      
      if (hgp->nProc_x_req != 1 && hgp->nProc_y_req != 1)  { /* Want 2D decomp */
          if ((leftend+1) > SMALL_PRIME && Zoltan_PHG_isPrime(leftend+1))
              --leftend; /* if it was prime just use one less #procs
                            (since it should be bigger than 7 it is safe to decrement)  */
          rightstart = leftend + 1;
          if ((hgc->nProc-rightstart) > SMALL_PRIME && Zoltan_PHG_isPrime(hgc->nProc-rightstart))
          ++rightstart; /* #procs for right was prime so reduce the number of procs
                           by increasing offset by one */
      } else
          rightstart = leftend + 1;
      
#ifdef _DEBUG1
      if (leftend<0 || (leftend+1)>=hgc->nProc ||
          rightstart<0 || rightstart>=hgc->nProc)
          errexit("hey hey Proc Number range is [0, %d] leftend=%d rightstart=%d"
                  "for left #pins=%zd nPins=%zd", hgc->nProc-1, leftend, rightstart,
                  left->dist_x[hgc->nProc_x] , gnVtx);
      uprintf(hgc, "before redistribute for left nProc=%d leftend=%d  rightstart=%d  ---------------\n",
              hgc->nProc, leftend, rightstart);
#endif
  } else {
      rightstart = leftend = -1;
  }
#else
  rightstart = leftend = -1;  
#endif

      /* if we want proc_split and there are more than one procs avail, and
         the left and right exist and they have enough vertex to distribute */
  if (rightstart != -1 && 
      (left->dist_x[hgc->nProc_x]>2*(leftend+1)) &&
      (right->dist_x[hgc->nProc_x]>2*(hgc->nProc-rightstart))) {
      PHGComm  leftcomm, rightcomm;
      HGraph  newleft, newright;
      int *leftvmap=NULL, *rightvmap=NULL, 
          *leftdest=NULL, *rightdest=NULL;
      ZOLTAN_COMM_OBJ *plan=NULL;    

      if (detail_timing) 
        ZOLTAN_TIMER_START(zz->ZTime, timer->rdredist, hgc->Communicator);
      
#ifdef _DEBUG1
      uprintf(hgc, "before redistribute for left leftend=%d ---------------\n",
                    leftend);
#endif
      Zoltan_PHG_Redistribute(zz, hgp, left, 0, leftend, &leftcomm, 
                              &newleft, &leftvmap, &leftdest);
      if (hgp->output_level >= PHG_DEBUG_LIST)     
          uprintf(hgc, "Left: H(%d, %d, %d) ----> H(%d, %d, %d) Weights=(%.2lf, %.2lf)\n",
                  left->nVtx, left->nEdge, left->nPins,
                  newleft.nVtx, newleft.nEdge, newleft.nPins, leftw, rightw);
      Zoltan_HG_HGraph_Free (left);
      
#ifdef _DEBUG1
      uprintf(hgc, "before redistribute for right ++++++++++++++++++++++\n");
#endif
      Zoltan_PHG_Redistribute(zz, hgp, right, rightstart, hgc->nProc-1,
                              &rightcomm, &newright, &rightvmap, &rightdest);
      if (hgp->output_level >= PHG_DEBUG_LIST)     
          uprintf(hgc, "Right: H(%d, %d, %d) ----> H(%d, %d, %d) Weights=(%.2lf, %.2lf)\n",
                  right->nVtx, right->nEdge, right->nPins, newright.nVtx, 
                  newright.nEdge, newright.nPins, leftw, rightw);
      Zoltan_HG_HGraph_Free (right);
      
      if (detail_timing) 
          ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdredist, hgc->Communicator);
      
      nsend = MAX(newleft.nVtx, newright.nVtx);
      part = (Partition) ZOLTAN_MALLOC (nsend * sizeof (int));
      proclist = (int *) ZOLTAN_MALLOC (nsend * sizeof (int));
      sendbuf =  (int *) ZOLTAN_MALLOC (nsend * 2 * sizeof (int));
      recvbuf =  (int *) ZOLTAN_MALLOC (nVtx * 2 * sizeof (int));
      if ((nsend && (!proclist || !sendbuf || !part)) ||
          (nVtx && !recvbuf))
          MEMORY_ERROR;
      if (hgc->myProc<=leftend) {
          float save_bal_tol=hgp->bal_tol;

          /* I'm on the left part so I should partition newleft */
          hgp->bal_tol = balanceTol(hgp, part_dim, 0, bisec_part_sizes,
                                    leftw+rightw, leftw);
          if (hgp->output_level >= PHG_DEBUG_LIST)     
              uprintf(hgc, "Left: H(%d, %d, %d) OldI: %.2lf NewI: %.2lf\n",
                      newleft.nVtx, newleft.nEdge, newleft.nPins,
                      save_bal_tol, hgp->bal_tol);          
          
          ierr = rdivide_and_prepsend (lo, mid, part, zz, &newleft, hgp, 
                                       level+1, proclist, sendbuf, 
                                       leftdest, leftvmap, &nsend,
                                       (timer ? timer->rdrdivide : -1),
				       2*father);

          hgp->bal_tol = save_bal_tol;
          Zoltan_HG_HGraph_Free (&newright); /* free dist_x and dist_y
                                                allocated in Redistribute*/
          if (detail_timing) {
              ZOLTAN_TIMER_START(zz->ZTime, timer->rdwait, hgc->Communicator);
          }          
      } else if (hgc->myProc>=rightstart) {
          float save_bal_tol=hgp->bal_tol;

          /* I'm on the right part so I should partition newright */
          hgp->bal_tol = balanceTol(hgp, part_dim, 1, bisec_part_sizes,
                                    leftw+rightw, rightw);
          if (hgp->output_level >= PHG_DEBUG_LIST)     
              uprintf(hgc, "Right: H(%d, %d, %d) OldI: %.2lf NewI: %.2lf\n",
                      newright.nVtx, newright.nEdge, newright.nPins,
                      save_bal_tol, hgp->bal_tol);

          ierr |= rdivide_and_prepsend (mid+1, hi, part, zz, &newright, hgp, 
                                        level+1, proclist, sendbuf, 
                                        rightdest, rightvmap, &nsend,
                                        (timer ? timer->rdrdivide : -1),
					2*father+1);

          hgp->bal_tol = save_bal_tol;          
          Zoltan_HG_HGraph_Free (&newleft); /* free dist_x and dist_y
                                               allocated in Redistribute*/
          if (detail_timing) {
              ZOLTAN_TIMER_START(zz->ZTime, timer->rdwait, hgc->Communicator);
          }          
      } else {
          nsend = 0; 
          Zoltan_HG_HGraph_Free (&newleft);  /* free dist_x and dist_y */
          Zoltan_HG_HGraph_Free (&newright); /* allocated in Redistribute*/
      }

      if (detail_timing) {
          MPI_Barrier(hgc->Communicator);
          ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdwait, hgc->Communicator);
      }

      if (detail_timing) {
          ZOLTAN_TIMER_START(zz->ZTime, timer->rdsend, hgc->Communicator);
      }
      --msg_tag;
      ierr |= Zoltan_Comm_Create(&plan, nsend, proclist, hgc->Communicator, msg_tag, &i);

#ifdef _DEBUG1
      if (!hgc->myProc_y) {
          if (i!=nVtx) 
              errexit("(%d,%d) I should be receiving nVtx(%d) part info but received %d", hgc->myProc_x, hgc->myProc_y, nVtx, i);
      } else {
          if (i)
              errexit("I'm not in the first row; why I'm receiving %d vertices?", i); 
      }
#endif
      
      --msg_tag;
      Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, 2*sizeof(int), (char *) recvbuf);

      MPI_Bcast(recvbuf, nVtx*2, MPI_INT, 0, hgc->col_comm);

      for (i=0; i<nVtx; ++i) {
#ifdef _DEBUG1
          int p=recvbuf[i*2+1];
          int v=recvbuf[i*2];

          if (v<0 || v>nVtx)
              errexit("sanity check failed for v=%d nVtx=%d\n", v, nVtx);
          if (p<lo || p>hi)
              errexit("sanity check failed for v=%d p=%d lo=%d hi=%d\n", v, p, lo, hi);
#endif
          final[recvbuf[i*2]] = recvbuf[i*2+1];
      }
      
      Zoltan_Comm_Destroy(&plan);
      if (detail_timing) {
        ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdsend, hgc->Communicator);
      }
  } else {
      if (left) {
          float save_bal_tol=hgp->bal_tol;
          
          hgp->bal_tol = balanceTol(hgp, part_dim, 0, bisec_part_sizes,
                                    leftw+rightw, leftw);
          if (hgp->output_level >= PHG_DEBUG_LIST)     
              uprintf(hgc, "Left: H(%d, %d, %d) OldI: %.2lf NewI: %.2lf\n",
                      left->nVtx, left->nEdge, left->nPins,
                      save_bal_tol, hgp->bal_tol);          
          
          if (do_timing)  /* Stop timer before recursion */
              ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdrdivide,
                                hgc->Communicator);

          ierr = Zoltan_PHG_rdivide(lo, mid, final, zz, left, hgp, level+1, 2*father);
          /* rdivide call will free "left" */

          if (do_timing)  /* Restart timer after recursion */
              ZOLTAN_TIMER_START(zz->ZTime, timer->rdrdivide, 
                                hgc->Communicator);

          hgp->bal_tol = save_bal_tol;                    
      }
      if (right) {
          float save_bal_tol=hgp->bal_tol;
          
          hgp->bal_tol = balanceTol(hgp, part_dim, 1, bisec_part_sizes,
                                    leftw+rightw, rightw);
          if (hgp->output_level >= PHG_DEBUG_LIST)     
              uprintf(hgc, "Right: H(%d, %d, %d) OldI: %.2lf NewI: %.2lf\n",
                      right->nVtx, right->nEdge, right->nPins,
                      save_bal_tol, hgp->bal_tol);

          
          if (do_timing)  /* Stop timer before recursion */
              ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdrdivide,
                                hgc->Communicator);
          
          ierr |= Zoltan_PHG_rdivide(mid+1, hi, final, zz, right, hgp, level+1, 2*father+1);
          /* rdivide call will free "right" */
          
          if (do_timing)  /* Restart timer after recursion */
              ZOLTAN_TIMER_START(zz->ZTime, timer->rdrdivide, 
                                hgc->Communicator);

          hgp->bal_tol = save_bal_tol;          
      }
  }
  

End:
  if (level>0)
      Zoltan_HG_HGraph_Free(hg);
    
  Zoltan_Multifree (__FILE__, __LINE__, 9, &pins[0], &lpins[0], &part, 
                    &left, &right, &proclist, &sendbuf, &recvbuf, 
                    &bisec_part_sizes);

  if (do_timing) 
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->rdrdivide, hgc->Communicator);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}


static int rdivide_and_prepsend(int lo, int hi, Partition final, ZZ *zz,
                                HGraph *hg,
                                PHGPartParams *hgp, int level,
                                int *proclist, int *sendbuf, int *dest,
                                int *vmap, int *nsend, int timer_rdivide, 
				int father)
{
    int      ierr=ZOLTAN_OK, i, nVtx=hg->nVtx;
    PHGComm  *hgc=hg->comm;
    int do_timing = (timer_rdivide > -1);
#ifdef _DEBUG1
    Zoltan_HG_Check(zz, hg);
#endif

    if (do_timing)  /* Stop timer before recursion */
        ZOLTAN_TIMER_STOP(zz->ZTime, timer_rdivide, hgc->Communicator);

    ierr = Zoltan_PHG_rdivide (lo, hi, final, zz, hg, hgp, level, father);
    /* rdivide will free the content of "hg" */
    
    if (do_timing) /* Restart rdivide timer */
        ZOLTAN_TIMER_START(zz->ZTime, timer_rdivide, hgc->Communicator);

    *nsend = 0;
    if (!hgc->myProc_y) { /* only first row sends the part vector */
        for (i=0; i<nVtx; ++i) {
            proclist[*nsend] = dest[i];
            sendbuf[(*nsend)*2] = vmap[i];
            sendbuf[(*nsend)*2+1] = final[i];
            ++(*nsend);
        }
    }
    
    ZOLTAN_FREE(&vmap);
    ZOLTAN_FREE(&dest);
#ifdef _DEBUG1
    if (hgc->col_comm==MPI_COMM_NULL || hgc->row_comm==MPI_COMM_NULL || hgc->Communicator==MPI_COMM_NULL)
              errexit("hey comm is NULL com=%x col=%x row=%x", hgc->Communicator, hgc->col_comm, hgc->row_comm);
#endif
    MPI_Comm_free(&hgc->col_comm);
    MPI_Comm_free(&hgc->row_comm);
    MPI_Comm_free(&hgc->Communicator);
    return ierr;
}




static int split_hypergraph (int *pins[2], HGraph *ohg, HGraph *nhg, 
                             PHGPartParams *hgp, Partition part, int partid, 
                             ZZ *zz, double *splitpw, double *otherpw, 
                             int connectivitycut)
{
  int *tmap = NULL;  /* temporary array mapping from old HGraph info to new */
  int edge, i, j, ierr=ZOLTAN_OK;  
  PHGComm *hgc = ohg->comm;
  char *yo = "split_hypergraph";
  double pw[2], tpw[2];
  ZOLTAN_GNO_TYPE tmp_gno;
  MPI_Datatype zoltan_gno_mpi_type;

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

  pw[0] = pw[1] = 0; /* 0 is the part being splitted, 1 is the other part(s) */
  Zoltan_HG_HGraph_Init (nhg);
  nhg->comm = ohg->comm;
  nhg->info               = ohg->info;
  nhg->nDim               = ohg->nDim;
  nhg->VtxWeightDim       = ohg->VtxWeightDim;
  nhg->EdgeWeightDim      = ohg->EdgeWeightDim;
  
  /* allocate memory for dynamic arrays in new HGraph and for tmap array */
  if (ohg->nVtx && (tmap = (int*) ZOLTAN_MALLOC (ohg->nVtx * sizeof (int)))==NULL)
      MEMORY_ERROR;
  
  /* fill in tmap array, -1 for ignored vertices, otherwise nonnegative int */
  nhg->nVtx = 0;
  for (i = 0; i < ohg->nVtx; i++)
      tmap[i] = (part[i] == partid) ? nhg->nVtx++ : -1; 

  /* save vertex weights, edge weights, and coordinates if they exist */
  if (nhg->nVtx && ohg->vwgt && nhg->VtxWeightDim &&  /* KDDKDD DO WE NEED CHECK OF ohg->vwgt HERE? */
      !(nhg->vwgt=(float*)ZOLTAN_MALLOC(nhg->nVtx*sizeof(float)*nhg->VtxWeightDim)))
      MEMORY_ERROR;
  if (nhg->nVtx && nhg->nDim &&
      !(nhg->coor=(double *)ZOLTAN_MALLOC(nhg->nVtx * sizeof(double) * nhg->nDim)))
      MEMORY_ERROR;
  if (nhg->nVtx && hgp->UseFixedVtx &&
      !(nhg->fixed_part = (int*)ZOLTAN_MALLOC(nhg->nVtx*sizeof(int))))
      MEMORY_ERROR;
  if (nhg->nVtx && hgp->UsePrefPart &&
      !(nhg->pref_part = (int*)ZOLTAN_MALLOC(nhg->nVtx*sizeof(int))))
      MEMORY_ERROR;  
  if (nhg->nVtx && (nhg->vmap = (int*) ZOLTAN_MALLOC (nhg->nVtx * sizeof (int)))==NULL)
      MEMORY_ERROR;
  
  for (i = 0; i < ohg->nVtx; i++) {
      int v=tmap[i];
      if (v!=-1) {
          nhg->vmap[v] = ohg->vmap[i];
          if (hgp->UseFixedVtx)
              nhg->fixed_part[v] = ohg->fixed_part[i];
          if (hgp->UsePrefPart) {
              if (ohg->pref_part[i]<0) /* was free */
                  nhg->pref_part[v] = -1; /* will be free */
              else { /* if it had a prefered part but moved to
                        other side of bisection; set it free
                        otherwise keep the preferred part */
                  nhg->pref_part[v] = (partid!=((ohg->pref_part[i] < ohg->bisec_split) ? 0 : 1) ) ?
                      -1 : ohg->pref_part[i];
              }              
          }
          if (nhg->VtxWeightDim) {
              /* UVC: TODO CHECK we're only using 1st weight! Right now this will be used
                 to compute balance ratio! Check this code when multiconstraint is added!
              */              
              pw[0] += ohg->vwgt[i*nhg->VtxWeightDim];
              memcpy(&nhg->vwgt[v*nhg->VtxWeightDim], &ohg->vwgt[i*nhg->VtxWeightDim], 
                     nhg->VtxWeightDim * sizeof(float)); 
          } else
              pw[0] += 1.0;
	  if (nhg->nDim) 
	    for (j = 0; j < nhg->nDim; j++)
	      nhg->coor[v*nhg->nDim + j] = ohg->coor[i*nhg->nDim + j];
      } else {
          pw[1] += (nhg->VtxWeightDim) ? ohg->vwgt[i*nhg->VtxWeightDim] : 1.0;
      }
      
  }
  MPI_Allreduce(pw, tpw, 2, MPI_DOUBLE, MPI_SUM, hgc->row_comm);
  *splitpw = tpw[0];  *otherpw = tpw[1];
  
  /* fill in hindex and hvertex arrays in new HGraph */
  nhg->nEdge = 0;
  nhg->nPins = 0;
  for (edge = 0; edge < ohg->nEdge; ++edge)
      if (pins[partid][edge] > 1) {
          ++nhg->nEdge;
          nhg->nPins += pins[partid][edge];
      }

  /* continue allocating memory for dynamic arrays in new HGraph */
  if (nhg->nEdge && !(nhg->hindex  = (int*) ZOLTAN_MALLOC ((nhg->nEdge+1) * sizeof (int))))
      MEMORY_ERROR;
  if (nhg->nPins && !(nhg->hvertex = (int*) ZOLTAN_MALLOC (nhg->nPins * sizeof (int))))
      MEMORY_ERROR;
  if (ohg->ewgt && nhg->EdgeWeightDim && nhg->nEdge &&
      !(nhg->ewgt=(float*)ZOLTAN_MALLOC(nhg->nEdge*sizeof(float)*nhg->EdgeWeightDim)))
      MEMORY_ERROR;
  
  nhg->nEdge = 0;
  nhg->nPins = 0;
  for (edge = 0; edge < ohg->nEdge; ++edge)
    if ((pins[partid][edge] > 1)         /* edge has at least two vertices in partition:
                                            we are skipping size 1 nets */
        && (connectivitycut              /* if connectivity-1 metric; split the edge */
            || !pins[1-partid][edge])) { /* if cut-net metric only take if
                                            it is not a cut net */
        nhg->hindex[nhg->nEdge] = nhg->nPins;
        for (i = ohg->hindex[edge]; i < ohg->hindex[edge+1]; ++i)
          if (tmap [ohg->hvertex[i]] >= 0)  {
            nhg->hvertex[nhg->nPins] = tmap[ohg->hvertex[i]];
            nhg->nPins++;  
          }
        if (nhg->ewgt)
            memcpy(&nhg->ewgt[nhg->nEdge*nhg->EdgeWeightDim], &ohg->ewgt[edge*nhg->EdgeWeightDim],
                   nhg->EdgeWeightDim * sizeof(float));
        ++nhg->nEdge;
    }
  if (nhg->hindex)
      nhg->hindex[nhg->nEdge] = nhg->nPins;

  /* We need to compute dist_x, dist_y */
  if (!(nhg->dist_x = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC((hgc->nProc_x+1), sizeof(ZOLTAN_GNO_TYPE)))
	 || !(nhg->dist_y = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC((hgc->nProc_y+1), sizeof(ZOLTAN_GNO_TYPE))))
      MEMORY_ERROR;

  tmp_gno = (ZOLTAN_GNO_TYPE)nhg->nVtx;

  MPI_Scan(&tmp_gno, nhg->dist_x, 1, zoltan_gno_mpi_type, MPI_SUM, hgc->row_comm);
  MPI_Allgather(nhg->dist_x, 1, zoltan_gno_mpi_type, &(nhg->dist_x[1]), 1, zoltan_gno_mpi_type, hgc->row_comm);
  nhg->dist_x[0] = 0;
  
  tmp_gno = (ZOLTAN_GNO_TYPE)nhg->nEdge;

  MPI_Scan(&tmp_gno, nhg->dist_y, 1, zoltan_gno_mpi_type, MPI_SUM, hgc->col_comm);
  MPI_Allgather(nhg->dist_y, 1, zoltan_gno_mpi_type, &(nhg->dist_y[1]), 1, zoltan_gno_mpi_type, hgc->col_comm);
  nhg->dist_y[0] = 0;
    
  ierr = Zoltan_HG_Create_Mirror (zz, nhg);
 End:
  ZOLTAN_FREE (&tmap);
  return ierr;
}


#ifdef __INTEL_COMPILER
#pragma optimize("",off)
#endif
/* KDD 3/3/16:  Release builds with Intel 15.0.2 compiler were failing nightly testing
 * for hg_ibm03/zdrive.inp.phg.ipm.fixed07.  I traced the differences in solution to the
 * optimization of this function.  All the input arguments were identical, yet the function
 * returned different values for release builds than for debug builds.  The debug builds'
 * results matched the results we get with every other compiler (release or debug), so I am
 * disabling optimization of this function only to restore passing tests with this compiler. 
 */
/* UVCUVC: CHECK currently only uses 1st weight */
static float balanceTol(PHGPartParams *hgp, int part_dim, int pno, float *ratios,
                        float tot, float pw)
{
    float ntol=(pw==0.0) ? 0.0 : (tot*hgp->bal_tol*ratios[part_dim*pno])/pw;


/*    printf("%s: TW=%.1lf pw=%.1lf (%.3lf) old_tol=%.2f  part_s=(%.3f, %.3f) and new tol=%.2f\n", (pno==0) ? "LEFT" : "RIGHT", tot, pw, pw/tot, hgp->bal_tol, ratios[0], ratios[1], ntol);*/
    return ntol;
}
#ifdef __INTEL_COMPILER
#pragma optimize("",on)
#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
