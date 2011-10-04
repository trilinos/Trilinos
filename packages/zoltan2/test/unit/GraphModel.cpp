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

#include <string>
#include <TestAdapters.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_GraphModel.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::OrdinalTraits;
using Teuchos::ScalarTraits;
using Teuchos::ArrayView;

//
// TODO - Add a MatrixInput adapter for CSR arrays, and use GIDs
//   that are not Ordinals - like std::pair<int,int>
//#include <Zoltan2_IdentifierTraits.hpp>
//using Zoltan2::IdentifierTraits;

template <typename Scalar, typename LNO, typename GNO>
  void checkGraph(
    Zoltan2::GraphModel<Scalar, LNO, GNO, 
      Zoltan2::MatrixInput<Scalar, LNO, GNO> > &graph, 
    RCP<Tpetra::CrsMatrix<Scalar, LNO, GNO> > M, std::string errMsg)
{
  RCP<Comm<int> > &comm = M->getComm();
  RCP<const Zoltan2::Environment> default_env = 
    Teuchos::rcp(new Zoltan2::Environment);

  // Values to check against
  size_t nRows = M->getNodeNumRows();
  size_t ngRows = M->getGlobalNumRows();
  size_t nEntries = M->getNodeNumEntries();
  size_t ngEntries = M->getGlobalNumEntries();
  size_t nColumns = M->getNodeNumCols();
  size_t ngColumns = M->getGlobalNumCols();

  int fail = 0;

  if (graph.getLocalNumVertices() != nRows)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getLocalNumVertices"+errMsg, 1)

  if (graph.getGlobalNumVertices() != ngRows)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getGlobalNumVertices"+errMsg, 1)

  if (graph.getGlobalNumEdges() != ngEntries)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getGlobalNumEdges"+errMsg, 1)

  if (graph.getLocalNumEdges() != nEntries)
    fail = 1;
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getLocalNumEdges"+errMsg, 1)

  ArrayView<GNO> vertexGids;
  ArrayView<Scalar> coords;     // not implemented
  ArrayView<Scalar> wgts;     // not implemented

  try{
    graph.getVertexList(vertexGids, coords, wgts);
  }
  catch (std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getVertexList"+errMsg, 1)

  size_t *edgeSize = new size_t [nRows];
  GNO *edgeId = new GNO [nEntries];
  int *procOwner = new int [nEntries];

  ArrayView<GNO> edgeGids;
  ArrayView<int> procIds;
  ArrayView<Scalar> wgts;           // not implemented

  for (LNO lid=0; lid < vertexGids.size(); lid++){
    edgeSize[lid] = graph.getVertexLocalEdge(lid, edgesGids, procIds, wgts);
    edgeId[lid] = edgeGids.getPtr();
    procOwner[lid] = procIds.getPtr();
  }

  for (LNO lid=0; lid < vertexGids.size(); lid++){
    edgeSize[lid] = 
      graph.getVertexGlobalEdge(vertexGids[lid], edgesGids, procIds, wgts);
    edgeId[lid] = edgeGids.getPtr();
    procOwner[lid] = procIds.getPtr();
  }
  // TODO check that graph returned represents M
}

template <typename Scalar, typename LNO, typename GNO>
  void testGraphFromXpetraMatrix(std::string fname, RCP<const Comm<int> > comm)
{
  int rank = comm->getRank;
  bool includeEpetra = false;
  std::string intName(OrdinalTraits<int>::name());
  std::string ScalarName(ScalarTraits<Scalar>::name());
  std::string LNOName(OrdinalTraits<LNO>::name());
  std::string GNOName(OrdinalTraits<GNO>::name());

  // A default environment 

  RCP<const Zoltan2::Environment> default_env = 
    Teuchos::rcp(new Zoltan2::Environment);

  if ((LNOName == intName) && (GNOName == intName))
    includeEpetra = true;

  // Create some matrix input adapters

  typedef Zoltan2::EpetraCrsMatrixInput epetraMatrix_t;
  typedef Zoltan2::TpetraCrsMatrixInput<Scalar, LNO, GNO> tpetraMatrix_t;
  typedef Zoltan2::XpetraCrsMatrixInput<Scalar, LNO, GNO> xpetraMatrix_t;

  RCP<epetraMatrix_t> emi;
  RCP<tpetraMatrix_t> tmi;
  RCP<xpetraMatrix_t> xmi;

  TestAdapters<Scalar,LNO,GNO> input(fname);

#ifdef HAVE_MPI
  input.setMpiCommunicator(MPI_COMM_WORLD);
#endif

  if (includeEpetra){
    emi = input.getEpetraCrsMatrixInputAdapter();
  }

  tmi = input.getTpetraCrsMatrixInputAdapter();

  xmi = input.getXpetraCrsMatrixInputAdapter();

  // Get original matrix for debugging.
  RCP<Tpetra::CrsMatrix<Scalar, LNO, GNO> > M = input.getMatrix();

  // Create a graph model with each and test it.

  int fail = 0;
  if (includeEpetra){

    // GraphModel built with EpetraCrsMatrixInput

    typedef Zoltan2::GraphModel<Scalar, LNO, GNO, epetraMatrix_t> epetraMatrixGraphModel_t;

    epetraMatrixGraphModel_t *graph = NULL;

    try{
      graph = new epetraMatrixGraphModel_t(emi, comm, default_env);
    }
    catch (std::exception &e){
      std::cerr << rank << ") Error " << e.what() << std::endl;
      fail = 1;
    }
    TEST_FAIL_AND_EXIT(*comm, fail==0, "Creating epetra graph model", 1)

    checkGraph(*graph, M, " (epetra matrix input)");

    delete graph;
  }

  // GraphModel built with TpetraCrsMatrixInput

  typedef Zoltan2::GraphModel<Scalar, LNO, GNO, tpetraMatrix_t> tpetraMatrixGraphModel_t;

  tpetraMatrixGraphModel_t *tgraph=NULL;
  try{
    tgraph = new tpetraMatrixGraphModel_t(tmi, comm, default_env);
  }
  catch (std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "Creating tpetra graph model", 1)

  checkGraph(*tgraph, M, " (tpetra matrix input)");

  delete tgraph;

  // GraphModel built with XpetraCrsMatrixInput

  typedef Zoltan2::GraphModel<Scalar, LNO, GNO, xpetraMatrix_t> xpetraMatrixGraphModel_t;

  xpetraMatrixGraphModel_t *xgraph=NULL;
  try{
    xgraph = new xpetraMatrixGraphModel_t(xmi, comm, default_env);
  }
  catch (std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "Creating xpetra graph model", 1)

  checkGraph(*xgraph, M, " (xpetra matrix input)");

  delete xgraph;

  if (!rank){
    std::cout << "Processed " << fname << ": ";
    std::cout << "Scalar = " << ScalarName << std::endl;
    std::cout << "LNO = " << LNOName << std::endl;
    std::cout << "GNO = " << GNOName << std::endl;
    std::cout << std::endl;
  }
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();

  std::vector<std::string> mtxFiles;
  
  mtxFiles.push_back("../data/simple.mtx");
  mtxFiles.push_back("../data/cage10.mtx");

  // To use this matrix we would need to pass a domain map
  // to FillComplete.  So we skip it for now.  TODO
  // mtxFiles.push_back("../data/diag500_4.mtx");

  for (unsigned int fileNum=0; fileNum < mtxFiles.size(); fileNum++){
    testGraphFromXpetraMatrix<double, int, int>(mtxFiles[fileNum], comm);
    testGraphFromXpetraMatrix<double, int, long>(mtxFiles[fileNum], comm);
  }

  return 0;
}
