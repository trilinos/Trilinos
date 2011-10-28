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

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <TestAdapters.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_IdentifierTraits.hpp>

#ifdef HAVE_MPI
#include <Teuchos_MPISession.hpp>
#endif

using namespace std;
using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::OrdinalTraits;
using Teuchos::ScalarTraits;
using Teuchos::ArrayView;

#define COMMENT(s) {if (rank==0) {std::cout << s << std::endl;}}

template <typename Scalar, typename LNO, typename GNO, typename Node>
  void testGraphModel(std::string fname, GNO xdim, GNO ydim, GNO zdim,
    const RCP<const Comm<int> > &comm,
    bool verbose, bool consecutiveIds)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0;

  // A default environment 
  RCP<const Zoltan2::Environment> default_env = 
    Teuchos::rcp(new Zoltan2::Environment);

  // User data
  typedef Tpetra::CrsMatrix<Scalar, LNO, GNO> tcrsMatrix_t;

  // Zoltan2 user data input adapter
  typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> adapter_t;

  // Test adapter generator
  TestAdapters<Scalar,LNO,GNO,LNO,GNO,Node> *input;

  if (fname.size() > 0)
    input = new TestAdapters<Scalar,LNO,GNO,LNO,GNO,Node>(fname, comm);
  else
    input = new TestAdapters<Scalar,LNO,GNO,LNO,GNO,Node>(xdim,ydim,zdim,comm);

  // TODO return by reference
  RCP<adapter_t> tmi = input->getTpetraCrsMatrixInputAdapter();

  // Question: Are the matrix global IDs consecutive (locally)?

  const GNO *rowIds, *colIds;
  const LNO *off, *rowLocalIds;
  size_t n = tmi->getRowListView(rowIds, rowLocalIds, off, colIds);
  bool inARow = Zoltan2::IdentifierTraits<GNO>::areConsecutive(rowIds, n);

  if (inARow)
    COMMENT("Matrix row IDs are locally consecutive")
  else
    COMMENT("Matrix row IDs are locally not consecutive")

  // Get original matrix for debugging.
  RCP<tcrsMatrix_t > M = input->getTpetraMatrix();
  LNO nLocalRows = M->getNodeNumRows();
  LNO nLocalNonZeros = M->getNodeNumEntries();
  GNO nGlobalRows =  M->getGlobalNumRows();
  GNO nGlobalNonZeros = M->getGlobalNumEntries();

  // Create a graph model with this input
  Zoltan2::GraphModel<adapter_t> *model = NULL;

  try{
    model = new Zoltan2::GraphModel<adapter_t>(tmi, comm, default_env);
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "Creating xpetra graph model", 1)

  // Test the GraphModel interface

  if (model->getLocalNumVertices() != nLocalRows)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getLocalNumVertices", 1)

  if (model->getGlobalNumVertices() != nGlobalRows)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getGlobalNumVertices", 1)

  if (model->getGlobalNumEdges() !=  nGlobalNonZeros)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getGlobalNumEdges", 1)

  if (model->getLocalNumEdges() != nLocalNonZeros)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getLocalNumEdges", 1)

  if (verbose && (nLocalRows > 100)){
    if (rank==0)
       std::cout << "Note: Turned off verbose for big graph" << std::endl;
    verbose=false;
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
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getVertexList", 1)

  if (vertexGids.size() != nLocalRows)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getVertexList", 1)
  
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
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getEdgeList", 1)

  if (numEdges != nLocalNonZeros)
    fail = 1;

  TEST_FAIL_AND_EXIT(*comm, fail==0, "getEdgeList size", 1)

  if (verbose){
    const GNO *v = vertexGids.getRawPtr();
    const GNO *e = edgeGids.getRawPtr();
    const int *owner = procIds.getRawPtr();
    const LNO *idx = offsets.getRawPtr();
    comm->barrier();
    for (int p=0; p < nprocs; p++){
      if (p == rank){
        std::cout << "Rank " << p << std::endl;
        for (LNO i=0; i < nLocalRows; i++){
          std::cout << "  Vtx " << *v++ << ": ";
          for (LNO j=idx[i]; j < idx[i+1]; j++)
            std::cout << *e++ << " (" << *owner++ << ") ";
          std::cout << std::endl;
        }
        std::cout.flush();
      }
      comm->barrier();
    }
    comm->barrier();
  }

  // Test getting the edge list for one vertex

  for (LNO i=0; (fail==0) && (i < nLocalRows); i++){
     LNO vertexLid =  i;
     GNO vertexGid =  vertexGids[i];
     ArrayView<const GNO> edges;
     ArrayView<const int> procs;
     ArrayView<const Scalar> wgts;

     size_t nedges = model->getVertexGlobalEdge(vertexGid, edges, procs, wgts);

     LNO idx0 = offsets[i];
     LNO idx1 = offsets[i+1];

     if (nedges != idx1 - idx0)
       fail = 1;
     else{
       for (LNO j=idx0,k=0; !fail && (j < idx1); j++,k++){
         if ((edges[k] != edgeGids[j]) || (procs[k] != procIds[j]))
           fail = 2;
       }
     }

     if (!fail){
       nedges = model->getVertexLocalEdge(vertexLid, edges, procs, wgts);

       if (nedges != idx1 - idx0)
         fail = 3;
       else{
         for (LNO j=idx0,k=0; !fail && (j < idx1); j++,k++){
           if ((edges[k] != edgeGids[j]) || (procs[k] != procIds[j]))
             fail = 4;
         }
       }
     }
  }

  delete model;
  delete input;
  ostringstream oss;
  oss << "Fail code " << fail << " getting vertex edges";
  TEST_FAIL_AND_EXIT(*comm, fail==0, oss.str(), 1)
}

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  Teuchos::GlobalMPISession session(&argc, &argv);
#endif
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  std::string nullString;
  std::vector<std::string> mtxFiles;
  
  mtxFiles.push_back("../data/simple.mtx");
  mtxFiles.push_back("../data/cage10.mtx");

  bool verbose = true;
  bool wishConsecutiveIds = true;

  for (unsigned int fileNum=0; fileNum < mtxFiles.size(); fileNum++){

    COMMENT(mtxFiles[fileNum]+" double, int, int, !wishConsecutiveIds");
    testGraphModel<double, int, int, Zoltan2::default_node_t>(
      mtxFiles[fileNum], 0,0,0,comm, verbose, !wishConsecutiveIds);

    COMMENT(mtxFiles[fileNum]+" float, int, long, !wishConsecutiveIds");
    testGraphModel<float, int, long, Zoltan2::default_node_t>(
      mtxFiles[fileNum], 0,0,0,comm,  !verbose, !wishConsecutiveIds);

    COMMENT(mtxFiles[fileNum]+" float, int, int, wishConsecutiveIds");
    testGraphModel<float, int, long, Zoltan2::default_node_t>(
      mtxFiles[fileNum], 0,0,0,comm,  !verbose, wishConsecutiveIds);
  }

  COMMENT("5x5x5 mesh, double, int, int, !wishConsecutiveIds");
  testGraphModel<double, int, int, Zoltan2::default_node_t>(
    nullString, 5, 5, 5, comm, verbose, !wishConsecutiveIds);

  COMMENT("5x5x5 mesh, double, int, int, wishConsecutiveIds");
  testGraphModel<double, int, int, Zoltan2::default_node_t>(
    nullString, 5, 5, 5, comm, verbose, wishConsecutiveIds);

  COMMENT("PASS");

  return 0;
}
