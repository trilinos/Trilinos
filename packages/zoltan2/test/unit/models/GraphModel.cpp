// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Testing of GraphModel built from Xpetra matrix input adapters.
//
//    TODO test GraphModel for a matrix that is not Xpetra, that
//         that global IDs that are not Teuchos::Ordinals.
//    TODO test for input with gids that are not consecutive, but
//              ask graph model to map them to consecutive
//    TODO test Epetra inputs
//    TODO this test does not require MPI, but uses TestAdapters
//               which does.  Modify TestAdapters and all unit
//               tests to run both Serial and MPI.

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>

#include <UserInputForTests.hpp>

#include <string>
#include <vector>
#include <iostream>

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayView.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::ArrayView;

using std::string;
using std::vector;

template<typename LNO, typename VertexT, typename EdgeT>
  void printGraph(LNO nrows, const VertexT *v, const EdgeT *e, 
    const int *owner, const LNO *idx, const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  if (rank == 0)
    if (owner)
      std::cout << "Global graph:" << std::endl;
    else
      std::cout << "Local graph:" << std::endl;

  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << "Rank " << p << std::endl;
      if (owner){
        for (LNO i=0; i < nrows; i++){
          std::cout << "  Vtx " << *v++ << ": ";
          for (LNO j=idx[i]; j < idx[i+1]; j++)
            std::cout << *e++ << " (" << *owner++ << ") ";
          std::cout << std::endl;
        }
        std::cout.flush();
      }
      else{
        for (LNO i=0; i < nrows; i++){
          std::cout << "  Vtx " << i << ": ";
          for (LNO j=idx[i]; j < idx[i+1]; j++)
            std::cout << *e++ << " ";
          std::cout << std::endl;
        }
        std::cout.flush();
      }
    }
    comm->barrier();
  }
  comm->barrier();
}

template <typename Scalar, typename LNO, typename GNO>
  void checkModel(RCP<const Tpetra::CrsMatrix<Scalar, LNO, GNO> > &M,
    const RCP<const Comm<int> > &comm,
    bool consecutiveIdsRequested, bool noSelfEdges)
{
  int fail=0;
  int rank = comm->getRank();
  RCP<const Zoltan2::Environment> default_env = 
    Zoltan2::getDefaultEnvironment();

  if (rank==0){
    std::cout << "     Request consecutive IDs " << consecutiveIdsRequested;
    std::cout << ", remove self edges " << noSelfEdges << std::endl;;
  }

  typedef Zoltan2::default_node_t node_t;
  typedef Tpetra::CrsMatrix<Scalar, LNO, GNO, node_t> tcrsMatrix_t;
  typedef Zoltan2::MatrixInput<tcrsMatrix_t> base_adapter_t;
  typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> adapter_t;

  adapter_t tmi(M);

  Zoltan2::GraphModel<base_adapter_t> *model = NULL;
  const base_adapter_t *baseTmi = &tmi;

  try{
    model = new Zoltan2::GraphModel<base_adapter_t>(baseTmi, default_env, 
      comm, consecutiveIdsRequested, noSelfEdges);
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "Creating xpetra graph model", 1)

  // Test the GraphModel interface

  LNO nLocalRows = M->getNodeNumRows();
  LNO nLocalNonZeros = M->getNodeNumEntries();
  GNO nGlobalRows =  M->getGlobalNumRows();
  GNO nGlobalNonZeros = M->getGlobalNumEntries();

  if (model->getLocalNumVertices() != nLocalRows)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumVertices", 1)

  if (model->getGlobalNumVertices() != nGlobalRows)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getGlobalNumVertices", 1)

  if (noSelfEdges){
    if (model->getGlobalNumEdges() >  nGlobalNonZeros)
      fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "getGlobalNumEdges", 1)

    if (model->getLocalNumEdges() > nLocalNonZeros)
      fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumEdges", 1)
  }
  else{
    if (model->getGlobalNumEdges() !=  nGlobalNonZeros)
      fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "getGlobalNumEdges", 1)

    if (model->getLocalNumEdges() != nLocalNonZeros)
      fail = 1;
    TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalNumEdges", 1)
  }

  ArrayView<const GNO> vertexGids;
  ArrayView<const Scalar> coords;  // not implemented yet
  ArrayView<const Scalar> wgts;    // not implemented yet

  try{
    model->getVertexList(vertexGids, coords, wgts);
  }
  catch (std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexList", 1)

  if (vertexGids.size() != nLocalRows)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, !fail, "getVertexList", 1)
  
  ArrayView<const GNO> edgeGids;
  ArrayView<const int> procIds;
  ArrayView<const LNO> offsets;
  size_t numEdges=0;

  try{
    numEdges = model->getEdgeList(edgeGids, procIds, offsets, wgts);
  }
  catch(std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getEdgeList", 1)

  LNO numLocalEdges = 0;
  for (LNO i=0; i < numEdges; i++){
    if (procIds[i] == rank)
      numLocalEdges++;
  }

  if (numEdges != model->getLocalNumEdges())
    fail = 1;

  TEST_FAIL_AND_EXIT(*comm, !fail, "getEdgeList size", 1)

  if (nGlobalRows < 200){
    printGraph<LNO, GNO, GNO>(nLocalRows, vertexGids.getRawPtr(), 
      edgeGids.getRawPtr(), procIds.getRawPtr(), offsets.getRawPtr(), comm);
  }
  else{
    if (rank==0) 
      std::cout << "    " << nGlobalRows << " total rows" << std::endl;
  }

  // Get graph restricted to this process

  ArrayView<const LNO> localEdges;
  ArrayView<const LNO> localOffsets;
  size_t numLocalNeighbors=0;

  try{
    numLocalNeighbors= model->getLocalEdgeList(localEdges, localOffsets, wgts);
  }
  catch(std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalEdgeList", 1)

  if (numLocalEdges != numLocalNeighbors)
    fail = 1;

  TEST_FAIL_AND_EXIT(*comm, !fail, "getLocalEdgeList size", 1)

  if (nGlobalRows < 200){
    printGraph<LNO, GNO, LNO>(nLocalRows, vertexGids.getRawPtr(), 
      localEdges.getRawPtr(), NULL, localOffsets.getRawPtr(), comm);
  }

  // Get graph restricted to this process


  delete model;

  if (rank==0) std::cout << "    OK" << std::endl;
}

template <typename Scalar, typename LNO, typename GNO>
  void testGraphModel(string fname, GNO xdim, GNO ydim, GNO zdim,
    const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  typedef Zoltan2::default_node_t node_t;
  typedef Tpetra::CrsMatrix<Scalar, LNO, GNO, node_t> tcrsMatrix_t;

  // Input generator
  UserInputForTests<Scalar,LNO,GNO> *input;

  if (fname.size() > 0)
    input = new UserInputForTests<Scalar,LNO,GNO>(fname, comm);
  else
    input = new UserInputForTests<Scalar,LNO,GNO>(xdim,ydim,zdim,comm);

  RCP<tcrsMatrix_t> M = input->getTpetraCrsMatrix();

  bool consecutiveIds=true;
  bool noSelfEdges=true;

  // Row Ids of test input are already consecutive

  if (rank == 0)
    std::cout << "  Matrix row IDs are globally consecutive." << std::endl;

  RCP<const tcrsMatrix_t> Mconsec = rcp_const_cast<const tcrsMatrix_t>(M);

  checkModel<Scalar,LNO,GNO>(Mconsec, comm, !consecutiveIds, !noSelfEdges);
  checkModel<Scalar,LNO,GNO>(Mconsec, comm, !consecutiveIds, noSelfEdges);
  checkModel<Scalar,LNO,GNO>(Mconsec, comm, consecutiveIds, noSelfEdges);

  // Do a round robin migration so that global IDs are not consecutive.

  Array<GNO> myNewRows;
  for (int i=rank; i < Mconsec->getGlobalNumRows(); i+=nprocs)
    myNewRows.push_back(i);

  RCP<const tcrsMatrix_t> Mnonconsec = 
    Zoltan2::XpetraTraits<tcrsMatrix_t>::doMigration(
      Mconsec, myNewRows.size(), myNewRows.getRawPtr());

  if (rank == 0)
    std::cout << "  Matrix row IDs are not globally consecutive." << std::endl;

  checkModel<Scalar,LNO,GNO>(Mnonconsec, comm, !consecutiveIds, !noSelfEdges);
  checkModel<Scalar,LNO,GNO>(Mnonconsec, comm, !consecutiveIds, noSelfEdges);
  checkModel<Scalar,LNO,GNO>(Mnonconsec, comm, consecutiveIds, noSelfEdges);

  delete input;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  string nullString;
  vector<string> mtxFiles;
  
  mtxFiles.push_back("../data/simple.mtx");
  mtxFiles.push_back("../data/cage10.mtx");

  for (int fileNum=0; fileNum < mtxFiles.size(); fileNum++){
    if (rank==0)
      std::cout << mtxFiles[fileNum] << ", float, int, int:" << std::endl;
    testGraphModel<float, int, int>(mtxFiles[fileNum], 0, 0, 0, comm);
  }

  if (rank==0)
    std::cout << "4x4x4 mesh generated matrix, float, int, long:" << std::endl;
  testGraphModel<float, int, long>(nullString, 4, 4, 4, comm);

  if (rank==0)
    std::cout << "PASS" << std::endl;

  return 0;
}

