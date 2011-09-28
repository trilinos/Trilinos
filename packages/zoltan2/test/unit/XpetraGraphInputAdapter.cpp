// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Testing *petra graph input adapters.
//  TODO
//  Find a few small mtx files for graphs.
//  for each graph:
//     Read it in with EpetraExt
//     Create an Epetra::CrsGraph
//     Instaniate a Zoltan2::GraphInputAdapter
//     Test all queries, and verify.
//     Test copy and assignment.         TODO
//     Test support of vertex and edge weights.
//     Test support of vertex coordinates.
//     Local IDs should be optional (TODO). Test without providing local IDs.
//     Ditto for Tpetra.
//     Ditto for direct Xpetra.
//
//   We don't support changing the graph in an adapter once
//   it is set up.
//
// Xpetra does not support const input, TODO - ask for const support
//
//  TODO this is a pretty boring test.  All the vertex and
//    and edge weights are 1.  We should get weights from
//    the mtx file and check against them.
//
#
#include <Zoltan2_config.h>
#include <string>
#include <vector>
#include <iostream>
#include <Zoltan2_TpetraCrsGraphInput.hpp>
#include <Zoltan2_EpetraCrsGraphInput.hpp>
#include <Zoltan2_Util.hpp>
#include <EpetraExt_CrsMatrixIn.h>
#include <Teuchos_ArrayRCP.hpp>

// temporarily
#include <Zoltan2_EpetraCrsMatrixInput.hpp>

// For Epetra tests
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>

// For Tpetra tests
#include <Teuchos_DefaultMpiComm.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>

using namespace std;
using namespace Teuchos;

int z2_rank;
int z2_globalToken;

#define TEST_FAIL_AND_EXIT(comm, f, s){ \
  comm.SumAll(&f, &z2_globalToken, 1); \
  if (z2_globalToken){ \
    if (z2_rank == 0) { \
     std::cerr << s << std::endl; \
     std::cout << "FAIL" << std::endl; \
    } \
    exit(1); \
  } \
}

#define CATCH_EXCEPTION(s)    \
  catch(std::exception &e){    \
    std::cerr << z2_rank << ") " << s << ", " << e.what() << std::endl;    \
    fail=1;    \
  }

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm c = MPI_COMM_WORLD;
  RCP<MpiComm<int> > comm = Zoltan2::getTeuchosMpiComm<int>(c);
  Epetra_MpiComm ecomm(MPI_COMM_WORLD);
#else
  RCP<SerialComm<int> > comm = rcp(new SerialComm<int>);
  Epetra_SerialComm ecomm;
#endif

  z2_rank = comm->getRank();

  std::vector<std::string> mtxFiles;
  
  mtxFiles.push_back("../data/simple.mtx");
  mtxFiles.push_back("../data/cage10.mtx");
  mtxFiles.push_back("../data/diag500_4.mtx");

  // Tests of EpetraCrsGraphInput, TpetraCrsGraphInput adapter.

  for (unsigned int fileNum=0; fileNum < mtxFiles.size(); fileNum++){

    Epetra_CrsMatrix *M=NULL;

    int fail = EpetraExt::MatrixMarketFileToCrsMatrix(
      mtxFiles[fileNum].c_str(), ecomm, M, 0, 0);

    TEST_FAIL_AND_EXIT(ecomm, fail, "reading mtx file");

    const Epetra_CrsGraph &g = M->Graph();
    const Epetra_BlockMap &map = g.RowMap();
    int lidMin = map.MinLID();
    int lidMax = map.MaxLID();
    int numVtx = g.NumMyRows();
    int numNZ = g.NumMyEntries();
    int base = map.IndexBase();
 
    //////////////////////////////////////////////
    // Create an EpetraCrsGraphInput adapter.
    //////////////////////////////////////////////

    RCP<Epetra_CrsGraph> graph = rcp(const_cast<Epetra_CrsGraph *>(&g));
    graph.release();  // we're not responsible for deleting graph

    Zoltan2::EpetraCrsGraphInput<float> *adapterE = NULL;

    try {
      adapterE = new Zoltan2::EpetraCrsGraphInput<float>(graph);
    }
    CATCH_EXCEPTION("adapterE constructor");

    // Just want to see if this will compile

    Zoltan2::EpetraCrsMatrixInput *adapterM = NULL;
    try {
      adapterM = new Zoltan2::EpetraCrsMatrixInput(Tpetra::rcp<Epetra_CrsMatrix>(M));
    }
    CATCH_EXCEPTION("adapterM constructor");

    //////////////////////////////////////////////
    // Create a TpetraCrsGraphInput adapter 
    // representing the same graph.
    //////////////////////////////////////////////

    RCP<Tpetra::Map<int, int> > tmap = rcp(
      new Tpetra::Map<int, int>(map.NumGlobalElements(), numVtx, base, comm));

    ArrayRCP<size_t> entriesPerRow(numVtx);
    for (int i=0; i < numVtx; i++){
      entriesPerRow[i] = g.NumMyIndices(i + base);
    }

    ArrayRCP<int> nonZeroIds(g.MaxNumIndices());

    RCP<Tpetra::CrsGraph<int, int> > tpetraVersion =
      rcp(new Tpetra::CrsGraph<int, int>(tmap, entriesPerRow));

    const Epetra_BlockMap &colMap = g.ColMap();
    for (int i=0; i < numVtx; i++){
      int numEdges;
      int *edgeIds;
      g.ExtractMyRowView(i+base, numEdges, edgeIds);
      for (int j=0; j < numEdges; j++){
        nonZeroIds[j] = colMap.GID(edgeIds[j]);
      }
      tpetraVersion->insertGlobalIndices(map.GID(i+base), 
        nonZeroIds.view(0, numEdges));
    }

    tpetraVersion->fillComplete();

    //Zoltan2::TpetraCrsGraphInput<float, int, int> adapterT(tpetraVersion);

    //////////////////////////////////////////////
    // Now test both adapters 
    //////////////////////////////////////////////

    size_t lidBase;
    bool consecLids = adapterE->haveConsecutiveLocalIds(lidBase);
    if (!consecLids || (lidBase != lidMin))
      fail = 1;
    
    TEST_FAIL_AND_EXIT(ecomm, fail, "lid base vs lid min");

    if (adapterE->inputAdapterName() != std::string("EpetraCrsGraph"))
      fail = 1;

    TEST_FAIL_AND_EXIT(ecomm, fail, "inputAdapterName");

    if (adapterE->getLocalNumVertices() != numVtx)
      fail = 1;

    TEST_FAIL_AND_EXIT(ecomm, fail, "getLocalNumVertices");

    // create some vertex weights and coordinates
    std::vector<int> lidList(numVtx);
    std::vector<float> vwgtList(numVtx, 1.0);
    std::vector<float> xyzList(numVtx*3);
    for (int i=lidMin,j=0,k=0; i <= lidMax;  i++){
      lidList[j++] = i;
      xyzList[k++] = i;
      xyzList[k++] = 2*i;
      xyzList[k++] = 3*i;
    }

    try{
      adapterE->setVertexWeights(lidList, vwgtList);
    }
    CATCH_EXCEPTION("setting vertex weights");

    TEST_FAIL_AND_EXIT(ecomm, fail, "setVertexWeights");

    try{
      adapterE->setVertexCoordinates(lidList, xyzList);
    }
    CATCH_EXCEPTION("setting vertex coords");

    TEST_FAIL_AND_EXIT(ecomm, fail, "setVertexCoordinates");

    if (adapterE->getLocalNumEdges() != numNZ)
      fail = 1;

    TEST_FAIL_AND_EXIT(ecomm, fail, "getLocalNumEdges");

    // create some edge weights
    std::vector<int> numNbors(numVtx);
    std::vector<int> nborGID(numNZ);
    std::vector<float> ewgtList(numNZ, 1.0);
    for (int i=lidMin,j=0, k=0; i <= lidMax;  i++){
      int numIndices;
      int *nzList;
      fail = g.ExtractMyRowView(i, numIndices, nzList);
      if (fail)
        break;

      numNbors[j++] = numIndices;
      for(int n=0; n < numIndices; n++)
        nborGID[k++] = nzList[n];
    }

    TEST_FAIL_AND_EXIT(ecomm, fail, "g.ExtractGlobalRowView");

    try{
      adapterE->setEdgeWeights(lidList, numNbors, nborGID, ewgtList);
    }
    CATCH_EXCEPTION("setting edge weights");

    TEST_FAIL_AND_EXIT(ecomm, fail, "setEdgeWeights");

    // Test that get methods are correct

    if (adapterE->getVertexWeightDim() != 1)
      fail = 1;

    TEST_FAIL_AND_EXIT(ecomm, fail, "getVertexWeightDim");

    if (adapterE->getEdgeWeightDim() != 1)
      fail = 1;

    TEST_FAIL_AND_EXIT(ecomm, fail, "getEdgeWeightDim");

    if (adapterE->getCoordinateDim() != 3)
      fail = 1;

    TEST_FAIL_AND_EXIT(ecomm, fail, "getCoordinateDim");

    std::vector<int> adapterEGids;
    std::vector<int> adapterELids;
    std::vector<float> adapterECoords;
    std::vector<float> adapterEVwgts;

    try{
      adapterE->getVertexListCopy(adapterEGids, adapterELids,
        adapterECoords, adapterEVwgts);
    }
    CATCH_EXCEPTION("getting vertex copy");

    TEST_FAIL_AND_EXIT(ecomm, fail, "getVertexListCopy");

    if ((adapterEGids.size() != numVtx) ||
        (adapterECoords.size() != numVtx*3) ||
        (adapterEVwgts.size() != numVtx) ){
      fail = 1;
    }
    else{
      for (int i=0, k=0; i < numVtx; i++, k+=3){
        int lid = i + lidBase;
        if ((adapterEGids[i] != map.GID(lid)) ||
            (adapterECoords[k] != lid) ||
            (adapterECoords[k+1] != 2*lid) ||
            (adapterECoords[k+2] != 3*lid) ||
            (adapterEVwgts[i] != 1.0) ){
          fail=1;
          break;
        }
      }
    }

    TEST_FAIL_AND_EXIT(ecomm, fail, "getVertexListCopy results");

    const int *gidView=NULL, *lidView=NULL;
    const float *wgtView=NULL, *coordView=NULL;
    int nv=0;

    try{
      nv = adapterE->getVertexListView(gidView, lidView, coordView, wgtView);
    }
    CATCH_EXCEPTION("getting vertex view");

    TEST_FAIL_AND_EXIT(ecomm, fail, "getVertexListView");

    if (nv != numVtx) 
      fail=1;
    else{
      for (int i=0, k=0; i < numVtx; i++, k+=3){
        int lid = i + lidBase;
        if ((gidView[i] != map.GID(lid)) ||
            (coordView[k] != lid) ||
            (coordView[k+1] != 2*lid) ||
            (coordView[k+2] != 3*lid) ||
            (wgtView[i] != 1.0) ){
          fail=1;
          break;
        }
      }
    }

    TEST_FAIL_AND_EXIT(ecomm, fail, "getVertexListView results");

    for (int lid=lidMin, i=0, k=0; lid <= lidMax; lid++, i++){
      unsigned nnbors = numNbors[i];
      std::vector<int> id;
      std::vector<float> wgt;

      try{
        adapterE->getVertexEdgeCopy(map.GID(lid), lid, id, wgt);
      }
      CATCH_EXCEPTION("getting edge copy");

      if (fail)
        break;
     

      if ((id.size() != nnbors) || (wgt.size() != nnbors))
        fail=1;
      else{
        for (int j=0; j < nnbors; j++, k++){
          if ((id[j] != nborGID[k]) || 
              (wgt[j] != 1.0) ){
            fail = 1;
            break;
          }
        }
      }
    }

    TEST_FAIL_AND_EXIT(ecomm, fail, "getVertexEdgeCopy/View results");
  }  // Next graph

  if (!z2_rank)
    std::cout << "PASS" << std::endl;

  MPI_Finalize();
  return 0;
}
