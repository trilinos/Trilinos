/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ptest.c
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
void TestParMetis(char *filename, MPI_Comm comm)
{
  int i, j, k, l, nparts, npes, mype, realcut;
  int opt1, opt2;
  GraphType graph, mgraph;
  idxtype *part, *mpart, *order, *sizes;
  int numflag=0, wgtflag=0, options[5], edgecut, ndims;
  float *xyz;


  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  ndims = 2;

  ReadTestGraph(&graph, filename, comm);
  xyz = ReadTestCoordinates(&graph, filename, 2, comm);

  part = idxmalloc(graph.nvtxs, "TestParMetis: part");


  /*======================================================================
  / ParMETIS_PartKway
  /=======================================================================*/
  options[0] = 1;
  options[3] = 3;
  wgtflag = 0;
  numflag = 0;
  for (nparts=2*npes; nparts>=npes/2; nparts = nparts/2) {
    for (opt1=1; opt1<=2; opt1++) {
      for (opt2=0; opt2<=300; opt2+=100) {
        options[1] = opt1;
        options[2] = opt2;

        if (mype == 0)
          printf("\nTesting ParMETIS_PartKway with options[1-3] = {%2d %2d %2d}, Nparts: %d\n", options[1], options[2], options[3], nparts);

        ParMETIS_PartKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, &wgtflag,
                 &numflag, &nparts, options, &edgecut, part, &comm);

        realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
        if (mype == 0) {
          if (realcut == edgecut)
            printf("ParMETIS_PartKway reported a cut of %d which is correct!\n", edgecut);
          else
            printf("ParMETIS_PartKway reported a cut of %d which is incorrect (realcut = %d)!\n", edgecut, realcut);
        }
      }
    }
  }


  /*======================================================================
  / ParMETIS_PartGeomKway 
  /=======================================================================*/
  options[0] = 1;
  options[3] = 3;
  wgtflag = 0;
  numflag = 0;
  for (nparts=2*npes; nparts>=npes/2; nparts = nparts/2) {
    for (opt1=1; opt1<=2; opt1++) {
      for (opt2=0; opt2<=300; opt2+=100) {
        options[1] = opt1;
        options[2] = opt2;

        if (mype == 0)
          printf("\nTesting ParMETIS_PartGeomKway with options[1-3] = {%2d %2d %2d}, Nparts: %d\n", options[1], options[2], options[3], nparts);

        ParMETIS_PartGeomKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, &wgtflag,
                 &numflag, &ndims, xyz, &nparts, options, &edgecut, part, &comm);

        realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
        if (mype == 0) {
          if (realcut == edgecut)
            printf("ParMETIS_PartGeomKway reported a cut of %d which is correct!\n", edgecut);
          else
            printf("ParMETIS_PartGeomKway reported a cut of %d which is incorrect (realcut = %d)!\n", edgecut, realcut);
        }
      }
    }
  }



  /*======================================================================
  / ParMETIS_PartGeom 
  /=======================================================================*/
  wgtflag = 0;
  numflag = 0;
  if (mype == 0)
    printf("\nTesting ParMETIS_PartGeom\n");

  ParMETIS_PartGeom(graph.vtxdist, &ndims, xyz, part, &comm);

  realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
  if (mype == 0) 
    printf("ParMETIS_PartGeom reported a cut of %d\n", realcut);


  /*======================================================================
  / ParMETIS_RefineKway 
  /=======================================================================*/
  options[0] = 1;
  options[3] = 3;
  wgtflag = 0;
  numflag = 0;
  if (mype == 0)
    printf("\nTesting ParMETIS_RefineKway with default options\n");

  ParMETIS_RefineKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, &wgtflag,
           &numflag, options, &edgecut, part, &comm);

  realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
  if (mype == 0) {
    if (realcut == edgecut)
      printf("ParMETIS_RefineKway reported a cut of %d which is correct!\n", edgecut);
    else
      printf("ParMETIS_RefineKway reported a cut of %d which is incorrect (realcut = %d)!\n", edgecut, realcut);
  }

  MALLOC_CHECK(NULL);

  /* Compute a good partition and move the graph. Do so quietly! */
  options[0] = 0;
  ParMETIS_PartKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, &wgtflag,
           &numflag, &npes, options, &edgecut, part, &comm);
  TestMoveGraph(&graph, &mgraph, part, comm);
  mpart = idxmalloc(mgraph.nvtxs, "TestParMetis: mpart");

  MALLOC_CHECK(NULL);

  /* Test #2 */
  options[0] = 1;
  options[3] = 3;

  if (mype == 0)
    printf("\nTesting ParMETIS_RefineKway with defualt options (after move)\n");

  ParMETIS_RefineKway(mgraph.vtxdist, mgraph.xadj, mgraph.adjncy, NULL, NULL, &wgtflag,
           &numflag, options, &edgecut, mpart, &comm);

  MALLOC_CHECK(NULL);

  realcut = ComputeRealCut2(graph.vtxdist, mgraph.vtxdist, part, mpart, filename, comm);
  if (mype == 0) {
    if (realcut == edgecut)
      printf("ParMETIS_RefineKway reported a cut of %d which is correct!\n", edgecut);
    else
      printf("ParMETIS_RefineKway reported a cut of %d which is incorrect (realcut = %d)!\n", edgecut, realcut);
  }


  /*======================================================================
  / Adaptive Routines
  /=======================================================================*/
  mgraph.vwgt = idxsmalloc(mgraph.nvtxs, 1, "TestParKMetis: mgraph.vwgt");
  AdaptGraph(&mgraph, 4, comm); 

  options[0] = 1;
  options[3] = 3;
  wgtflag = 2;
  numflag = 0;

  if (mype == 0)
    printf("\nTesting ParMETIS_RepartLDiffusion with defualt options\n");

  ParMETIS_RepartLDiffusion(mgraph.vtxdist, mgraph.xadj, mgraph.adjncy, mgraph.vwgt, NULL, 
           &wgtflag, &numflag, options, &edgecut, mpart, &comm);

  realcut = ComputeRealCut2(graph.vtxdist, mgraph.vtxdist, part, mpart, filename, comm);
  if (mype == 0) {
    if (realcut == edgecut)
      printf("ParMETIS_RepartLDiffusion reported a cut of %d which is correct!\n", edgecut);
    else
      printf("ParMETIS_RepartLDiffusion reported a cut of %d which is incorrect (realcut = %d)!\n", edgecut, realcut);
  }

  if (mype == 0)
    printf("\nTesting ParMETIS_RepartGDiffusion with defualt options\n");

  ParMETIS_RepartGDiffusion(mgraph.vtxdist, mgraph.xadj, mgraph.adjncy, mgraph.vwgt, NULL, 
           &wgtflag, &numflag, options, &edgecut, mpart, &comm);

  realcut = ComputeRealCut2(graph.vtxdist, mgraph.vtxdist, part, mpart, filename, comm);
  if (mype == 0) {
    if (realcut == edgecut)
      printf("ParMETIS_RepartGDiffusion reported a cut of %d which is correct!\n", edgecut);
    else
      printf("ParMETIS_RepartGDiffusion reported a cut of %d which is incorrect (realcut = %d)!\n", edgecut, realcut);
  }

  if (mype == 0)
    printf("\nTesting ParMETIS_RepartRemap with defualt options\n");

  ParMETIS_RepartRemap(mgraph.vtxdist, mgraph.xadj, mgraph.adjncy, mgraph.vwgt, NULL, 
           &wgtflag, &numflag, options, &edgecut, mpart, &comm);

  realcut = ComputeRealCut2(graph.vtxdist, mgraph.vtxdist, part, mpart, filename, comm);
  if (mype == 0) {
    if (realcut == edgecut)
      printf("ParMETIS_RepartRemap reported a cut of %d which is correct!\n", edgecut);
    else
      printf("ParMETIS_RepartRemap reported a cut of %d which is incorrect (realcut = %d)!\n", edgecut, realcut);
  }

  if (mype == 0)
    printf("\nTesting ParMETIS_RepartMLRemap with defualt options\n");

  ParMETIS_RepartMLRemap(mgraph.vtxdist, mgraph.xadj, mgraph.adjncy, mgraph.vwgt, NULL, 
           &wgtflag, &numflag, options, &edgecut, mpart, &comm);

  realcut = ComputeRealCut2(graph.vtxdist, mgraph.vtxdist, part, mpart, filename, comm);
  if (mype == 0) {
    if (realcut == edgecut)
      printf("ParMETIS_RepartMLRemap reported a cut of %d which is correct!\n", edgecut);
    else
      printf("ParMETIS_RepartMLRemap reported a cut of %d which is incorrect (realcut = %d)!\n", edgecut, realcut);
  }





  /*======================================================================
  / ParMETIS_NodeND 
  /=======================================================================*/
  sizes = idxmalloc(2*npes, "TestParMetis: sizes");
  order = idxmalloc(graph.nvtxs, "TestParMetis: sizes");

  options[0] = 1;
  options[3] = 3;
  numflag = 0;

  for (opt2=1; opt2<=2; opt2++) {
    options[1] = opt2;

    if (mype == 0)
      printf("\nTesting ParMETIS_NodeND with options[1,3] = {%d %d}\n", options[1], options[3]);

    ParMETIS_NodeND(graph.vtxdist, graph.xadj, graph.adjncy, &numflag, options,
             order, sizes, &comm);
  }


  GKfree(&part, &mpart, &order, &sizes, LTERM);

}


/******************************************************************************
* This function takes a partition vector that is distributed and reads in
* the original graph and computes the edgecut
*******************************************************************************/
int ComputeRealCut(idxtype *vtxdist, idxtype *part, char *filename, MPI_Comm comm)
{
  int i, j, k, nvtxs, nedges, mype, npes, cut;
  idxtype *xadj, *adjncy, *gpart, scratch;
  MPI_Status status;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (mype != 0) {
    MPI_Send((void *)part, vtxdist[mype+1]-vtxdist[mype], IDX_DATATYPE, 0, 1, comm);
  }
  else {  /* Processor 0 does all the rest */
    gpart = idxmalloc(vtxdist[npes], "ComputeRealCut: gpart");
    idxcopy(vtxdist[1], part, gpart);

    for (i=1; i<npes; i++) 
      MPI_Recv((void *)(gpart+vtxdist[i]), vtxdist[i+1]-vtxdist[i], IDX_DATATYPE, i, 1, comm, &status);

    ReadMetisGraph(filename, &nvtxs, &xadj, &adjncy);

    /* OK, now compute the cut */
    for (cut=0, i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (gpart[i] != gpart[adjncy[j]])
          cut++;
      }
    }
    cut = cut/2;

    GKfree(&gpart, &xadj, &adjncy, LTERM);

    return cut;
  }
}


/******************************************************************************
* This function takes a partition vector that is distributed and reads in
* the original graph and computes the edgecut
*******************************************************************************/
int ComputeRealCut2(idxtype *vtxdist, idxtype *mvtxdist, idxtype *part, idxtype *mpart, char *filename, MPI_Comm comm)
{
  int i, j, k, nvtxs, nedges, mype, npes, cut;
  idxtype *xadj, *adjncy, *gpart, *gmpart, *perm, *sizes, scratch;
  MPI_Status status;


  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (mype != 0) {
    MPI_Send((void *)part, vtxdist[mype+1]-vtxdist[mype], IDX_DATATYPE, 0, 1, comm);
    MPI_Send((void *)mpart, mvtxdist[mype+1]-mvtxdist[mype], IDX_DATATYPE, 0, 1, comm);
  }
  else {  /* Processor 0 does all the rest */
    gpart = idxmalloc(vtxdist[npes], "ComputeRealCut: gpart");
    idxcopy(vtxdist[1], part, gpart);
    gmpart = idxmalloc(mvtxdist[npes], "ComputeRealCut: gmpart");
    idxcopy(mvtxdist[1], mpart, gmpart);

    for (i=1; i<npes; i++) {
      MPI_Recv((void *)(gpart+vtxdist[i]), vtxdist[i+1]-vtxdist[i], IDX_DATATYPE, i, 1, comm, &status);
      MPI_Recv((void *)(gmpart+mvtxdist[i]), mvtxdist[i+1]-mvtxdist[i], IDX_DATATYPE, i, 1, comm, &status);
    }

    /* OK, now go and reconstruct the permutation to go from the graph to mgraph */
    perm = idxmalloc(vtxdist[npes], "ComputeRealCut: perm");
    sizes = idxsmalloc(npes+1, 0, "ComputeRealCut: sizes");

    for (i=0; i<vtxdist[npes]; i++)
      sizes[gpart[i]]++;
    MAKECSR(i, npes, sizes);
    for (i=0; i<vtxdist[npes]; i++)
      perm[i] = sizes[gpart[i]]++;

    /* Ok, now read the graph from the file */
    ReadMetisGraph(filename, &nvtxs, &xadj, &adjncy);

    /* OK, now compute the cut */
    for (cut=0, i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (gmpart[perm[i]] != gmpart[perm[adjncy[j]]])
          cut++;
      }
    }
    cut = cut/2;

    GKfree(&gpart, &gmpart, &perm, &sizes, &xadj, &adjncy, LTERM);

    return cut;
  }
}



/******************************************************************************
* This function takes a graph and its partition vector and creates a new
* graph corresponding to the one after the movement
*******************************************************************************/
void TestMoveGraph(GraphType *ograph, GraphType *omgraph, idxtype *part, MPI_Comm comm)
{
  int i, j, k;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;
  int options[5] = {0, 0, 1, 0, 0};

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  SetUpCtrl(&ctrl, npes, options, comm); 
  ctrl.CoarsenTo = 1;  /* Needed by SetUpGraph, otherwise we can FP errors */
  graph = SetUpGraph(&ctrl, ograph->vtxdist, ograph->xadj, NULL, ograph->adjncy, NULL, 0);
  PreAllocateMemory(&ctrl, graph, &wspace);

  SetUp(&ctrl, graph, &wspace);
  graph->where = part;
  mgraph = MoveGraph(&ctrl, graph, &wspace);

  omgraph->gnvtxs = mgraph->gnvtxs;
  omgraph->nvtxs = mgraph->nvtxs;
  omgraph->nedges = mgraph->nedges;
  omgraph->vtxdist = mgraph->vtxdist;
  omgraph->xadj = mgraph->xadj;
  omgraph->adjncy = mgraph->adjncy;
  mgraph->vtxdist = NULL;
  mgraph->xadj = NULL;
  mgraph->adjncy = NULL;
  FreeGraph(mgraph);

  graph->where = NULL;
  FreeInitialGraphAndRemap(graph, 0);
  FreeWSpace(&wspace);
}  


