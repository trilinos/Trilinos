/*
 * Copyright 1998, Regents of the University of Minnesota
 *
 * tstadpt.c
 * 
 * This file contains code for testing teh adaptive partitioning routines
 *
 * Started 5/19/97
 * George
 *
 * $Id$
 *
 */

#include <parmetis.h>


/***********************************************************************************
* This function is the testing routine for the adaptive multilevel partitioning code.
* It computes a partition from scratch, it then moves the graph and changes some
* of the vertex weights and then call the adaptive code.
************************************************************************************/
void TestAdaptiveMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, 
                idxtype *part, int *options, int afactor, MPI_Comm comm)
{
  int i, j, k, min, max;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;
  int CoarsenTo;
  GraphType *mgraph;
  idxtype *rpart;
  int numflag = 0, wgtflag = 3, edgecut;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);


  SetUpCtrl(&ctrl, npes, options, comm);
  ctrl.CoarsenTo = amin(vtxdist[npes]-1, 25*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, NULL, adjncy, NULL, 0);

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  if (npes%2 == 0)
    FldGlobal_Partition(&ctrl, graph, &wspace, 0);
  else
    Global_Partition(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, rprintf(&ctrl, "Final Cut: %6d, \tBalance: %6.3f [%d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*graph->gnvtxs), 
          graph->gpwgts[idxamax(npes, graph->gpwgts)], idxsum(npes, graph->gpwgts), graph->gnvtxs));

  mgraph = MoveGraph(&ctrl, graph, &wspace);

  FreeInitialGraphAndRemap(graph, 0);
  FreeWSpace(&wspace);

  AdaptGraph2(mgraph, afactor, comm); 

  rpart = idxmalloc(mgraph->nvtxs, "rpart");

  rprintf(&ctrl, "Testing ParMETIS_RepartLDiffusion...\n");
  ParMETIS_RepartLDiffusion(mgraph->vtxdist, mgraph->xadj, mgraph->adjncy, mgraph->vwgt, 
     mgraph->adjwgt, &wgtflag, &numflag, options, &edgecut, rpart, &comm);

  rprintf(&ctrl, "Testing ParMETIS_RepartGDiffusion...\n");
  ParMETIS_RepartGDiffusion(mgraph->vtxdist, mgraph->xadj, mgraph->adjncy, mgraph->vwgt, 
     mgraph->adjwgt, &wgtflag, &numflag, options, &edgecut, rpart, &comm);

  rprintf(&ctrl, "Testing ParMETIS_RepartRemap...\n");
  ParMETIS_RepartRemap(mgraph->vtxdist, mgraph->xadj, mgraph->adjncy, mgraph->vwgt, 
     mgraph->adjwgt, &wgtflag, &numflag, options, &edgecut, rpart, &comm);

  rprintf(&ctrl, "Testing ParMETIS_RepartMLRemap...\n");
  ParMETIS_RepartMLRemap(mgraph->vtxdist, mgraph->xadj, mgraph->adjncy, mgraph->vwgt, 
     mgraph->adjwgt, &wgtflag, &numflag, options, &edgecut, rpart, &comm);
}



/*************************************************************************
* This function implements a simple graph adaption strategy.
**************************************************************************/
void AdaptGraph(GraphType *graph, int afactor, MPI_Comm comm)
{
  int i, j, k, nvtxs, nadapt, firstvtx, lastvtx;
  int npes, mype, mypwgt, max, min, sum;
  idxtype *vwgt, *xadj, *adjncy, *adjwgt, *perm;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  srand(mype*afactor);
  srand48(mype*afactor);

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  if (graph->adjwgt == NULL)
    adjwgt = graph->adjwgt = idxsmalloc(graph->nedges, 1, "AdaptGraph: adjwgt");
  else
    adjwgt = graph->adjwgt;
  vwgt = graph->vwgt;

  firstvtx = graph->vtxdist[mype];
  lastvtx = graph->vtxdist[mype+1];

  perm = idxmalloc(nvtxs, "AdaptGraph: perm");
  FastRandomPermute(nvtxs, perm, 1);

  nadapt = RandomInRange(nvtxs);

  for (i=0; i<nadapt; i++)
    vwgt[perm[i]] = afactor*vwgt[perm[i]];

/*
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (k >= firstvtx && k < lastvtx) {
	adjwgt[j] = (int)pow(1.0*(amin(vwgt[i],vwgt[k-firstvtx])), .6667);
        if (adjwgt[j] == 0)
          adjwgt[j] = 1;
      }
    }
  }
*/

  mypwgt = idxsum(nvtxs, vwgt);

  MPI_Allreduce((void *)&mypwgt, (void *)&max, 1, MPI_INT, MPI_MAX, comm);
  MPI_Allreduce((void *)&mypwgt, (void *)&min, 1, MPI_INT, MPI_MIN, comm);
  MPI_Allreduce((void *)&mypwgt, (void *)&sum, 1, MPI_INT, MPI_SUM, comm);

  if (mype == 0)
    printf("Initial Load Imbalance: %5.4f, [%5d %5d %5d]\n", (1.0*max*npes)/(1.0*sum), min, max, sum);

  free(perm);
}


/*************************************************************************
* This function implements a simple graph adaption strategy.
**************************************************************************/
void AdaptGraph2(GraphType *graph, int afactor, MPI_Comm comm)
{
  int i, j, k, nvtxs, nadapt, firstvtx, lastvtx;
  int npes, mype, mypwgt, max, min, sum;
  idxtype *vwgt, *xadj, *adjncy, *adjwgt, *perm;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  srand(mype*afactor);
  srand48(mype*afactor);

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  if (graph->adjwgt == NULL)
    adjwgt = graph->adjwgt = idxsmalloc(graph->nedges, 1, "AdaptGraph: adjwgt");
  else
    adjwgt = graph->adjwgt;
  vwgt = graph->vwgt;

  firstvtx = graph->vtxdist[mype];
  lastvtx = graph->vtxdist[mype+1];


/*  if (RandomInRange(npes+1) < .05*npes) { */ 
  if (RandomInRange(npes+1) < 2) { 
    printf("[%d] is adapting\n", mype);
    for (i=0; i<nvtxs; i++)
      vwgt[i] = afactor*vwgt[i];
  }

  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (k >= firstvtx && k < lastvtx) {
	adjwgt[j] = (int)pow(1.0*(amin(vwgt[i],vwgt[k-firstvtx])), .6667);
        if (adjwgt[j] == 0)
          adjwgt[j] = 1;
      }
    }
  }
      
  mypwgt = idxsum(nvtxs, vwgt);

  MPI_Allreduce((void *)&mypwgt, (void *)&max, 1, MPI_INT, MPI_MAX, comm);
  MPI_Allreduce((void *)&mypwgt, (void *)&min, 1, MPI_INT, MPI_MIN, comm);
  MPI_Allreduce((void *)&mypwgt, (void *)&sum, 1, MPI_INT, MPI_SUM, comm);

  if (mype == 0)
    printf("Initial Load Imbalance: %5.4f, [%5d %5d %5d]\n", (1.0*max*npes)/(1.0*sum), min, max, sum);

}

