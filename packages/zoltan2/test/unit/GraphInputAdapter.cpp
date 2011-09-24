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
//     Test copy and assignment.
//     Test support of vertex and edge weights.
//     Test support of vertex coordinates.
//     Local IDs should be optional (TODO). Test without providing local IDs.
//     Ditto for Tpetra.
//     Ditto for direct Xpetra.
//
//   We don't support changing the graph in an adapter once
//   it is set up.
//
#include <Zoltan2_config.h>
#include <string>
#include <vector>
#include <iostream>
#include <Zoltan2_TpetraCrsGraphInput.hpp>
#include <Zoltan2_EpetraCrsGraphInput.hpp>
#include <Zoltan2_Util.hpp>
#include <EpetraExt_CrsMatrixIn.h>

// For Epetra tests
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>

// For Tpetra tests
#include <Teuchos_DefaultMpiComm.hpp>

using namespace std;
using namespace Teuchos;

int z2_rank;
int z2_globalToken;

#define TEST_FAIL_AND_RETURN(comm, f, s){ \
  comm.SumAll(&f, &z2_globalToken, 1); \
  if (z2_globalToken){ \
    if (z2_rank == 0) { \
     std::cerr << s << std::endl; \
     std::cout << "FAIL" << std::endl; \
    } \
    return 1; \
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

  for (unsigned int fileNum=0; fileNum < mtxFiles.size(); fileNum++){

    Epetra_CrsMatrix *M=NULL;

    int fail = EpetraExt::MatrixMarketFileToCrsMatrix(
      mtxFiles[fileNum].c_str(), ecomm, M, 0, 0);

    TEST_FAIL_AND_RETURN(ecomm, fail, "reading mtx file");

    const Epetra_CrsGraph &g = M->Graph();
    const Epetra_BlockMap &map = g.RowMap();
    int lidMin = map.MinLID();
    int lidMax = map.MaxLID();
    int numVtx = g.NumMyRows();
    int numNZ = g.NumMyEntries();

    // Xpetra does not support const input, because it is also
    // used to change Epetra and Tpetra objects.  TODO - ask for const support

    RCP<Epetra_CrsGraph> graph = rcp(const_cast<Epetra_CrsGraph *>(&g));
    Zoltan2::EpetraCrsGraphInput<float> adapter;

    try {
      adapter.setGraph(graph);
    }
    CATCH_EXCEPTION("adapter constructor");
   

    if (adapter.inputAdapterName() != std::string("EpetraCrsGraph"))
      fail = 1;

    TEST_FAIL_AND_RETURN(ecomm, fail, "inputAdapterName");

    if (adapter.getLocalNumVertices() != numVtx)
      fail = 1;

    TEST_FAIL_AND_RETURN(ecomm, fail, "getLocalNumVertices");

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
      adapter.setVertexWeights(lidList, vwgtList);
    }
    CATCH_EXCEPTION("setting vertex weights");

    TEST_FAIL_AND_RETURN(ecomm, fail, "setVertexWeights");

    try{
      adapter.setVertexCoordinates(lidList, xyzList);
    }
    CATCH_EXCEPTION("setting vertex coords");

    TEST_FAIL_AND_RETURN(ecomm, fail, "setVertexCoordinates");
     
    if (adapter.getLocalNumEdges() != numNZ);
      fail = 1;

    TEST_FAIL_AND_RETURN(ecomm, fail, "getLocalNumEdges");

    // create some edge weights
    std::vector<int> numNbors(numVtx);
    std::vector<int> nborGID(numNZ);
    std::vector<float> ewgtList(numNZ, 1.0);
    for (int i=lidMin,j=0, k=0; i <= lidMax;  i++){
      int numIndices;
      int *nzList;
      fail = g.ExtractGlobalRowView(map.GID(i), numIndices, nzList);
      if (fail)
        break;

      numNbors[j++] = numIndices;
      for(int n=0; n < numIndices; n++)
        nborGID[k++] = nzList[n];
    }

    TEST_FAIL_AND_RETURN(ecomm, fail, "g.ExtractGlobalRowView");

    try{
      adapter.setEdgeWeights(lidList, numNbors, nborGID, ewgtList);
    }
    CATCH_EXCEPTION("setting edge weights");

    TEST_FAIL_AND_RETURN(ecomm, fail, "setEdgeWeights");

    // Test that get methods are correct

    if (adapter.getVertexWeightDim() != 1)
      fail = 1;

    TEST_FAIL_AND_RETURN(ecomm, fail, "getVertexWeightDim");

    if (adapter.getEdgeWeightDim() != 1)
      fail = 1;

    TEST_FAIL_AND_RETURN(ecomm, fail, "getEdgeWeightDim");

    if (adapter.getCoordinateDim() != 3)
      fail = 1;

    TEST_FAIL_AND_RETURN(ecomm, fail, "getCoordinateDim");

    std::vector<int> adapterGids;
    std::vector<int> adapterLids;
    std::vector<float> adapterCoords;
    std::vector<float> adapterVwgts;

    try{
      adapter.getVertexListCopy(adapterGids, adapterLids,
        adapterCoords, adapterVwgts);
    }
    CATCH_EXCEPTION("getting vertex copy");

    TEST_FAIL_AND_RETURN(ecomm, fail, "getVertexListCopy");

    if ((adapterGids.size() != numVtx) ||
        (adapterLids.size() != numVtx) ||
        (adapterCoords.size() != numVtx*3) ||
        (adapterVwgts.size() != numVtx) ){
      fail = 1;
    }
    else{
      for (int i=0, k=0; i < numVtx; i++, k+=3){
        int lid = adapterLids[i];
        if ((lid < lidMin) || (lid > lidMax)){
          fail = 1;
          break;
        }
        if ((adapterGids[i] != map.GID(lid)) ||
            (adapterCoords[k] != lid) ||
            (adapterCoords[k+1] != 2*lid) ||
            (adapterCoords[k+2] != 3*lid) ||
            (adapterVwgts[i] != 1.0) ){
          fail=1;
          break;
        }
      }
    }

    TEST_FAIL_AND_RETURN(ecomm, fail, "getVertexListCopy results");

    const int *gidView=NULL, *lidView=NULL;
    const float *wgtView=NULL, *coordView=NULL;
    int nv=0;

    try{
      nv = adapter.getVertexListView(gidView, lidView, wgtView, coordView);
    }
    CATCH_EXCEPTION("getting vertex view");

    TEST_FAIL_AND_RETURN(ecomm, fail, "getVertexListView");

    if (nv != numVtx) 
      fail=1;
    else{
      for (int i=0, k=0; i < numVtx; i++, k+=3){
        int lid = lidView[i];
        if ((lid < lidMin) || (lid > lidMax)){
          fail = 1;
          break;
        }
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

    TEST_FAIL_AND_RETURN(ecomm, fail, "getVertexListView results");

    for (int lid=lidMin, i=0, k=0; lid <= lidMax; lid++, i++){
      int num=-1;
      unsigned nnbors = numNbors[i];
      std::vector<int> id;
      std::vector<float> wgt;

      try{
        adapter.getVertexEdgeCopy(map.GID(lid), lid, id, wgt);
      }
      CATCH_EXCEPTION("getting edge copy");

      TEST_FAIL_AND_RETURN(ecomm, fail, "getVertexEdgeCopy");

      try{
        num = adapter.getVertexEdgeView(map.GID(lid), lid, gidView, wgtView);
      }
      CATCH_EXCEPTION("getting edge view");

      TEST_FAIL_AND_RETURN(ecomm, fail, "getVertexEdgeView");

      if ((id.size() != nnbors) || (num != nnbors) || (wgt.size() != nnbors))
        fail=1;
      else{
        for (int j=0; j < nnbors; j++, k++){
          if ((id[j] != nborGID[k]) || 
              (wgt[j] != 1.0) ||
              (gidView[j] != nborGID[k]) ||
              (wgtView[j] != 1.0) ){
            fail = 1;
            break;
          }
        }
      }
      TEST_FAIL_AND_RETURN(ecomm, fail, "getVertexEdgeCopy/View results");
    }
  }  // Next graph

  if (!z2_rank)
    std::cout << "PASS" << std::endl;

  return 0;
}
