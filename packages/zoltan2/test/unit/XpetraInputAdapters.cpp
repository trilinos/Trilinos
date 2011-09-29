// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Basic testing xpetra graph input adapters.  

#include <iostream>
#include <vector>
#include <string>
#include <Zoltan2_config.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>

#include <EpetraExt_CrsMatrixIn.h>

#include <Teuchos_DefaultComm.hpp>

#include <Epetra_CrsGraph.h>
#include <Zoltan2_EpetraCrsGraphInput.hpp>

#include <Tpetra_CrsGraph.hpp>
#include <Zoltan2_TpetraCrsGraphInput.hpp>

#include <Epetra_CrsMatrix.h>
#include <Zoltan2_EpetraCrsMatrixInput.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Zoltan2_TpetraCrsMatrixInput.hpp>

#include <Xpetra_CrsMatrix.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>

#include <Xpetra_CrsGraph.hpp>
#include <Zoltan2_XpetraCrsGraphInput.hpp>

#include <Tpetra_Map.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::rcp_const_cast;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

int z2_localToken, z2_globalToken, z2_rank;

#ifdef HAVE_MPI
#define TEST_FAIL_AND_EXIT(comm, ok, s){ \
  z2_localToken = ( (ok) ? 0 : 1);       \
  comm.SumAll(&z2_localToken, &z2_globalToken, 1); \
  if (z2_globalToken){ \
    if (z2_rank == 0) { \
     std::cerr << s << std::endl; \
     std::cout << "FAIL" << std::endl; \
    } \
    MPI_Finalize(); \
    exit(1); \
  } \
}
#else
#define TEST_FAIL_AND_EXIT(comm, ok, s){ \
  z2_localToken = ( (ok) ? 0 : 1);       \
  if (z2_localToken){ \
   std::cerr << s << std::endl; \
   std::cout << "FAIL" << std::endl; \
  } \
  exit(1); \
}
#endif

#define CATCH_EXCEPTION(s)    \
  catch(std::exception &e){    \
    std::cerr << z2_rank << ") " << s << ", " << e.what() << std::endl;    \
    fail=1;    \
  }

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm ecomm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm ecomm;
#endif

  RCP<const Teuchos::Comm<int> > tcomm =
     Teuchos::DefaultComm<int>::getComm();

  z2_rank = tcomm->getRank();

  std::vector<std::string> mtxFiles;
  
  mtxFiles.push_back("../data/simple.mtx");
  mtxFiles.push_back("../data/cage10.mtx");

  // To use this matrix we would need to pass a domain map
  // to FillComplete.  So we skip it for now.  TODO
  // mtxFiles.push_back("../data/diag500_4.mtx");

  for (unsigned int fileNum=0; fileNum < mtxFiles.size(); fileNum++){

    Epetra_CrsMatrix *M=NULL;

    int fail = EpetraExt::MatrixMarketFileToCrsMatrix(
      mtxFiles[fileNum].c_str(), ecomm, M, 0, 0);

    TEST_FAIL_AND_EXIT(ecomm, fail==0, "reading mtx file");

    // Epetra matrix input

    Zoltan2::EpetraCrsMatrixInput emi(rcp(M));

    // Epetra graph input

    const Epetra_CrsGraph &g = M->Graph();
    RCP<Epetra_CrsGraph> grcp =
      rcp_const_cast<Epetra_CrsGraph>(rcpFromRef(g));
    Zoltan2::EpetraCrsGraphInput<float> egi(grcp); 

    // Tpetra graph input

    const Epetra_Map &map = M->RowMap();
    int gNumVtx= map.NumGlobalElements();
    int numVtx = map.NumMyElements();
    int base = map.IndexBase();

    ArrayRCP<int> myGIDs(numVtx);
    for (int i=0; i < numVtx; i++)
      myGIDs[i] = map.GID(i+base);

    ArrayRCP<size_t> entriesPerRow(numVtx);
    for (int i=0; i < numVtx; i++){
      entriesPerRow[i] = g.NumMyIndices(i+base);
    }

    RCP<Tpetra::Map<int, int> > tmap = rcp(
      new Tpetra::Map<int, int>(gNumVtx, myGIDs.view(0,numVtx), base, tcomm));

    RCP<Tpetra::CrsGraph<int, int> > tpetraGraph =
      rcp(new Tpetra::CrsGraph<int, int>(
        tmap, entriesPerRow, Tpetra::StaticProfile));

    Array<int> nonZeroIds(g.MaxNumIndices());
    const Epetra_BlockMap &colMap = g.ColMap();
    for (int i=0; i < numVtx; i++){
      int numEdges;
      int *edgeIds;
      g.ExtractMyRowView(i+base, numEdges, edgeIds);
      for (int j=0; j < numEdges; j++){
        nonZeroIds[j] = colMap.GID(edgeIds[j]);
      }
      tpetraGraph->insertGlobalIndices(map.GID(i+base),
        nonZeroIds.view(0, numEdges));
    }

    tpetraGraph->fillComplete();

    Zoltan2::TpetraCrsGraphInput<int, int> tgi(tpetraGraph);  

    // Tpetra matrix input

    RCP<Tpetra::CrsMatrix<float, int, int> > tpetraMatrix =
     rcp(new Tpetra::CrsMatrix<float, int, int>(tpetraGraph));

    tpetraMatrix->setAllToScalar(1.0);
    tpetraMatrix->fillComplete();

    Zoltan2::TpetraCrsMatrixInput<float, int, int> tmi(tpetraMatrix);

    ///////////////////////////////////////////////////////
    // Test graph adapters using graph input interface
    ///////////////////////////////////////////////////////

    size_t num = egi.getLocalNumVertices();
    TEST_FAIL_AND_EXIT(ecomm, num==numVtx, "egi.getLocalNumVertices");

    num = tgi.getLocalNumVertices();
    TEST_FAIL_AND_EXIT(ecomm, num==numVtx, "tgi.getLocalNumVertices");

    num = egi.getGlobalNumVertices();
    TEST_FAIL_AND_EXIT(ecomm, num==gNumVtx, "egi.getGlobalNumVertices");

    num = tgi.getGlobalNumVertices();
    TEST_FAIL_AND_EXIT(ecomm, num==gNumVtx, "tgi.getGlobalNumVertices");

    num = egi.getLocalNumEdges();
    TEST_FAIL_AND_EXIT(ecomm, num==g.NumMyEntries(), "egi.getLocalNumEdges");

    num = tgi.getLocalNumEdges();
    TEST_FAIL_AND_EXIT(ecomm, num==g.NumMyEntries(), "tgi.getLocalNumEdges");

    num = egi.getGlobalNumEdges();
    TEST_FAIL_AND_EXIT(ecomm, num==g.NumGlobalEntries(), "egi.getGlobalNumEdges");

    num = tgi.getGlobalNumEdges();
    TEST_FAIL_AND_EXIT(ecomm, num==g.NumGlobalEntries(), "tgi.getGlobalNumEdges");

    // TODO the rest of the methods

    ///////////////////////////////////////////////////////
    // Test matrix adapters using matrix input interface
    ///////////////////////////////////////////////////////

    num = emi.getLocalNumRows();
    TEST_FAIL_AND_EXIT(ecomm, num==numVtx, "emi.getLocalNumRows");

    num = tmi.getLocalNumRows();
    TEST_FAIL_AND_EXIT(ecomm, num==numVtx, "tmi.getLocalNumRows");

    num = emi.getGlobalNumRows();
    TEST_FAIL_AND_EXIT(ecomm, num==gNumVtx, "emi.getGlobalNumRows");

    num = tmi.getGlobalNumRows();
    TEST_FAIL_AND_EXIT(ecomm, num==gNumVtx, "tmi.getGlobalNumRows");

    num = emi.getLocalNumColumns();
    TEST_FAIL_AND_EXIT(ecomm, num==g.NumMyCols(), "emi.getLocalNumColumns");

    num = tmi.getLocalNumColumns();
    TEST_FAIL_AND_EXIT(ecomm, num==g.NumMyCols(), "tmi.getLocalNumColumns");

    num = emi.getGlobalNumColumns();
    TEST_FAIL_AND_EXIT(ecomm, num==g.NumGlobalCols(), "emi.getGlobalNumColumns");

    num = tmi.getGlobalNumColumns();
    TEST_FAIL_AND_EXIT(ecomm, num==g.NumGlobalCols(), "tmi.getGlobalNumColumns");

    // TODO the rest of the methods

    if (!z2_rank){
      std::cout << "Processed " << mtxFiles[fileNum] << ": ";
      std::cout << gNumVtx << " vertices, ";
      std::cout << g.NumGlobalEntries() << " edges." << std::endl;
    }
    //delete M;
  }

  if (!z2_rank)
    std::cout << "PASS" << std::endl;

  MPI_Finalize();
  return 0;
}
