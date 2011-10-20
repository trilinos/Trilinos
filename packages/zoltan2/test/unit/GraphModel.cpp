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

#include <string>
#include <vector>
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


template <typename Model, typename Matrix>
  void checkGraph(Model &graph, RCP<Matrix> M, std::string errMsg)
{
  typedef typename Model::scalar_t Scalar;
  typedef typename Model::lno_t LNO;
  typedef typename Model::gno_t GNO;

  const RCP<const Comm<int> > &comm = M->getComm();
  int rank = comm->getRank();
  RCP<const Zoltan2::Environment> default_env = 
    Teuchos::rcp(new Zoltan2::Environment);

  // Values to check against
  size_t nRows = M->getNodeNumRows();
  size_t ngRows = M->getGlobalNumRows();
  size_t nEntries = M->getNodeNumEntries();
  size_t ngEntries = M->getGlobalNumEntries();

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

  ArrayView<const GNO> vertexGids;
  ArrayView<const Scalar> coords;     // not implemented
  ArrayView<const Scalar> wgts;     // not implemented

  try{
    graph.getVertexList(vertexGids, coords, wgts);
  }
  catch (std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "getVertexList"+errMsg, 1)

  ArrayView<const GNO> edgeGids;
  ArrayView<const int> procIds;

#if 0
   // TODO won't compile
  for (LNO lid=0; lid < vertexGids.size(); lid++){
    size_t edgeSize = 
       graph.getVertexLocalEdge(lid, edgeGids, procIds, wgts);
    //size_t edgeSize = graph.getVertexLocalEdge(lid, edgeGids, procIds, wgts);
    // TODO check that these make sense
  }
#endif

#if 0
  // graph.getVertexGlobalEdge not implemented for now
  for (LNO lid=0; lid < vertexGids.size(); lid++){
    edgeSize[lid] = 
      graph.getVertexGlobalEdge(vertexGids[lid], edgesGids, procIds, wgts);
    edgeId[lid] = edgeGids.getRawPtr();
    procOwner[lid] = procIds.getRawPtr();
  }
#endif
  // TODO check that graph returned represents M
  // TODO test graph.getEdgeList()
}

template <typename Scalar, typename LNO, typename GNO, typename Node>
  void testGraphInputFromMatrix(std::string fname, 
  const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();

  // If LNO and GNO are int we can test an Epetra input adapter.
  bool includeEpetra = false;
  std::string intName(OrdinalTraits<int>::name());
  std::string ScalarName(ScalarTraits<Scalar>::name());
  std::string LNOName(OrdinalTraits<LNO>::name());
  std::string GNOName(OrdinalTraits<GNO>::name());

  // A default environment 

  RCP<const Zoltan2::Environment> default_env = 
    Teuchos::rcp(new Zoltan2::Environment);

  // User input data types

  typedef Tpetra::CrsMatrix<Scalar, LNO, GNO> tcrsMatrix_t;
  typedef Epetra_CrsMatrix                    ecrsMatrix_t;
  typedef Xpetra::CrsMatrix<Scalar, LNO, GNO> xcrsMatrix_t;

  // Zoltan2 user data input adapters

  typedef Zoltan2::XpetraCrsMatrixInput<ecrsMatrix_t> EpetraCrsMatrixInput;
  typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> TpetraCrsMatrixInput;
  typedef Zoltan2::XpetraCrsMatrixInput<xcrsMatrix_t> XpetraCrsMatrixInput;

  // Test adapters

  RCP<TpetraCrsMatrixInput> tmi;
  RCP<EpetraCrsMatrixInput> emi;
  RCP<XpetraCrsMatrixInput> xmi;

  TestAdapters<Scalar,LNO,GNO,LNO,GNO,Node> input(fname, comm);

  if (includeEpetra){
    emi = input.getEpetraCrsMatrixInputAdapter();
  }

  tmi = input.getTpetraCrsMatrixInputAdapter();

  xmi = input.getXpetraCrsMatrixInputAdapter();

  // Get original matrix for debugging.
  RCP<tcrsMatrix_t > M = input.getTpetraMatrix();

  // Create graph models with different user objects and test them.

  typedef Zoltan2::XpetraCrsMatrixInput<ecrsMatrix_t> EpetraCrsMatrixUpcast;
  typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> TpetraCrsMatrixUpcast;

  typedef Zoltan2::GraphModel<XpetraCrsMatrixInput> xGraphModel_t;
  typedef Zoltan2::GraphModel<EpetraCrsMatrixUpcast> eGraphModel_t;
  typedef Zoltan2::GraphModel<TpetraCrsMatrixUpcast> tGraphModel_t;

  int fail = 0;
  if (includeEpetra){
    eGraphModel_t *graph = NULL;
    try{
      RCP<EpetraCrsMatrixUpcast> upcastEmi = 
        Teuchos::rcp_implicit_cast<EpetraCrsMatrixUpcast>(emi);
      graph = new eGraphModel_t(upcastEmi, comm, default_env);
    }
    catch (std::exception &e){
      std::cerr << rank << ") Error " << e.what() << std::endl;
      fail = 1;
    }
    TEST_FAIL_AND_EXIT(*comm, fail==0, "Creating epetra graph model", 1)

    checkGraph<eGraphModel_t, tcrsMatrix_t>(*graph, M, " (epetra matrix input)");

    delete graph;
  }

  // GraphModel built with TpetraCrsMatrixInput

  tGraphModel_t *tgraph=NULL;
  try{
    RCP<TpetraCrsMatrixUpcast> upcastTmi = 
      Teuchos::rcp_implicit_cast<TpetraCrsMatrixUpcast>(tmi);
    tgraph = new tGraphModel_t(upcastTmi, comm, default_env);
  }
  catch (std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "Creating tpetra graph model", 1)

  checkGraph<tGraphModel_t, tcrsMatrix_t>(*tgraph, M, " (tpetra matrix input)");

  delete tgraph;

  // GraphModel built with XpetraCrsMatrixInput

  xGraphModel_t *xgraph=NULL;
  try{
    xgraph = new xGraphModel_t(xmi, comm, default_env);
  }
  catch (std::exception &e){
    std::cerr << rank << ") Error " << e.what() << std::endl;
    fail = 1;
  }
  TEST_FAIL_AND_EXIT(*comm, fail==0, "Creating xpetra graph model", 1)

  checkGraph<xGraphModel_t, tcrsMatrix_t>(*xgraph, M, " (xpetra matrix input)");

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
  const RCP<const Comm<int> > &comm = DefaultComm<int>::getComm();

  std::vector<std::string> mtxFiles;
  
  mtxFiles.push_back("../data/simple.mtx");
  mtxFiles.push_back("../data/cage10.mtx");

  // To use this matrix we would need to pass a domain map
  // to FillComplete.  So we skip it for now.  TODO
  // mtxFiles.push_back("../data/diag500_4.mtx");

  for (unsigned int fileNum=0; fileNum < mtxFiles.size(); fileNum++){

    testGraphInputFromMatrix<
      double, int, int, Zoltan2::default_node_t>(mtxFiles[fileNum], comm);

    testGraphInputFromMatrix<
      double, int, long, Zoltan2::default_node_t>(mtxFiles[fileNum], comm);
  }

  if (comm->getRank() == 0)
    std::cout << "PASS" << std::endl;

  return 0;
}
