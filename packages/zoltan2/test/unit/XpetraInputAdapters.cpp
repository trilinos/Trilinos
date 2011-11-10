// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Basic testing of xpetra and tpetra matrix and graph input adapters
//
// TODO Test with epetra matrix and graphs.

#include <TestAdapters.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <string>

using namespace std;
using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::OrdinalTraits;
using Teuchos::ScalarTraits;

template <typename GraphAdapter, typename Graph>
  void testGraphAdapter(RCP<GraphAdapter > input, RCP<const Graph> g, 
    RCP<const Comm<int> > comm)
{
  typedef typename GraphAdapter::lid_t LID;
  typedef typename GraphAdapter::gid_t GID;
  typedef typename GraphAdapter::lno_t LNO;
  typedef typename GraphAdapter::gno_t GNO;
  typedef typename GraphAdapter::node_t Node;

  size_t nRows = g->getNodeNumRows();
  size_t ngRows = g->getGlobalNumRows();
  size_t nEntries = g->getNodeNumEntries();
  size_t ngEntries = g->getGlobalNumEntries();

  size_t num = input->getLocalNumVertices();
  TEST_FAIL_AND_EXIT(*comm, num==nRows, "getLocalNumVertices", 1);
  num = input->getGlobalNumVertices();
  TEST_FAIL_AND_EXIT(*comm, num==ngRows, "getGlobalNumVertices", 1);
  num = input->getLocalNumEdges();
  TEST_FAIL_AND_EXIT(*comm, num==nEntries, "getLocalNumEdges", 1);
  num = input->getGlobalNumEdges();
  TEST_FAIL_AND_EXIT(*comm, num==ngEntries, "getGlobalNumEdges", 1);
}

template <typename MatrixAdapter, typename Matrix>
  void testMatrixAdapter(RCP<MatrixAdapter > input, RCP<Matrix> m, RCP<const Comm<int> >&comm)
{
  typedef typename MatrixAdapter::scalar_t Scalar;
  typedef typename MatrixAdapter::lid_t LID;
  typedef typename MatrixAdapter::gid_t GID;
  typedef typename MatrixAdapter::lno_t LNO;
  typedef typename MatrixAdapter::gno_t GNO;
  typedef typename MatrixAdapter::node_t Node;

  size_t nRows = m->getNodeNumRows();
  size_t ngRows = m->getGlobalNumRows();
  size_t nColumns = m->getNodeNumCols();
  size_t ngColumns = m->getGlobalNumCols();

  size_t num = input->getLocalNumRows();
  TEST_FAIL_AND_EXIT(*comm, num==nRows, "getLocalNumRows", 1);
  num = input->getGlobalNumRows();
  TEST_FAIL_AND_EXIT(*comm, num==ngRows, "getGlobalNumRows", 1);
  num = input->getLocalNumColumns();
  TEST_FAIL_AND_EXIT(*comm, num==nColumns, "getLocalNumColumns", 1);
  num = input->getGlobalNumColumns();
  TEST_FAIL_AND_EXIT(*comm, num==ngColumns, "getGlobalNumColumns", 1);
}

template <typename Scalar, typename LNO, typename GNO>
  void testInputAdapters(GNO xdim, GNO ydim, GNO zdim,
    std::string &fname, RCP<const Comm<int> > &comm, int rank)
{
  bool include64BitIds = false;

  std::string ScalarName(ScalarTraits<Scalar>::name());
  std::string LNOName(OrdinalTraits<LNO>::name());
  std::string GNOName(OrdinalTraits<GNO>::name());

//
//  64-bit global IDs need to wait for Tpetra problem to be resolved.
//  if (sizeof(GNO) >= 8)
//    include64BitIds = true;
//

  // User input data types

  typedef Tpetra::CrsMatrix<Scalar, LNO, GNO> tcrsMatrix_t;
  typedef Epetra_CrsMatrix                    ecrsMatrix_t;
  typedef Xpetra::CrsMatrix<Scalar, LNO, GNO> xcrsMatrix_t;
  typedef Tpetra::CrsGraph<LNO, GNO>          tcrsGraph_t;
  typedef Epetra_CrsGraph                     ecrsGraph_t;
  typedef Xpetra::CrsGraph<LNO, GNO>          xcrsGraph_t;

  // Zoltan2 user data input adapters

  typedef Zoltan2::XpetraCrsMatrixInput<ecrsMatrix_t> EpetraCrsMatrixInput;
  typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> TpetraCrsMatrixInput;
  typedef Zoltan2::XpetraCrsMatrixInput<xcrsMatrix_t> XpetraCrsMatrixInput;
  typedef Zoltan2::XpetraCrsGraphInput<ecrsGraph_t> EpetraCrsGraphInput;
  typedef Zoltan2::XpetraCrsGraphInput<tcrsGraph_t> TpetraCrsGraphInput;
  typedef Zoltan2::XpetraCrsGraphInput<xcrsGraph_t> XpetraCrsGraphInput;

  // Test adapters

  RCP<TpetraCrsMatrixInput> tmi;
  RCP<EpetraCrsMatrixInput> emi;
  RCP<XpetraCrsMatrixInput> xmi;
  RCP<TpetraCrsGraphInput> tgi;
  RCP<EpetraCrsGraphInput> egi;
  RCP<XpetraCrsGraphInput> xgi;

  TestAdapters<Scalar,LNO,GNO> *input = NULL;

  if (fname.size() > 0)
    input = new TestAdapters<Scalar,LNO,GNO>(fname, comm);
  else
    input = new TestAdapters<Scalar,LNO,GNO>(xdim, ydim, zdim, comm);

  // Get the original matrix in order to compare answers to it.
  RCP<tcrsMatrix_t> M = input->getTpetraMatrix();
  RCP<const tcrsGraph_t> G = M->getCrsGraph();

  cout << "Testing with Tpetra::CrsMatrix" << endl;
  tmi = input->getTpetraCrsMatrixInputAdapter();
  testMatrixAdapter<TpetraCrsMatrixInput, tcrsMatrix_t>(tmi, M, comm);

  cout << "Testing with Tpetra::CrsGraph" << endl;
  tgi = input->getTpetraCrsGraphInputAdapter();
  testGraphAdapter<TpetraCrsGraphInput, tcrsGraph_t>(tgi, G, comm);

  RCP<xcrsMatrix_t> xM = input->getXpetraMatrix();
  RCP<const xcrsGraph_t> xG = xM->getCrsGraph();

  cout << "Testing with Xpetra::TpetraCrsMatrix" << endl;
  xmi = input->getXpetraCrsMatrixInputAdapter();
  testMatrixAdapter<XpetraCrsMatrixInput, xcrsMatrix_t>(xmi, xM, comm);

  cout << "Testing with Xpetra::TpetraCrsGraph" << endl;
  xgi = input->getXpetraCrsGraphInputAdapter();
  testGraphAdapter<XpetraCrsGraphInput, xcrsGraph_t>(xgi, xG, comm);

  if (include64BitIds){
    // This matrix has ids that use the high order 4 bytes of
    // of the 8 byte ID.  true: delete the original M after making a new one.
    //RCP<tpetraM_t> tmi64 = input->getTpetraCrsMatrixInputAdapter64(true);
    //testMatrixAdapter<TpetraCrsMatrixInput, tcrsMatrix_t>(tmi64, M, comm);
  }

  delete input;

  if (!rank){
    if (fname.size())
      std::cout << "Processed " << fname << ": ";
    else{
      std::cout << "Created input adapter for matrix/graph of size ";
      std::cout << xdim << " x " << ydim << " x " << zdim << std::endl;
    }
    std::cout << "Scalar = " << ScalarName;
    std::cout << ", LNO = " << LNOName;
    std::cout << ", GNO = " << GNOName << std::endl;
    if (include64BitIds)
      std::cout << "Including a Tpetra::CrsMatrix with 64-bit Ids." << std::endl;
    std::cout << std::endl;
  }
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = session.getRank();

  // Test input adapters created from Matrix Market files.

  std::vector<std::string> mtxFiles;
  
  mtxFiles.push_back("../data/simple.mtx");
  mtxFiles.push_back("../data/cage10.mtx");

  for (unsigned int fileNum=0; fileNum < mtxFiles.size(); fileNum++){
    testInputAdapters<double, int, int>(0,0,0,mtxFiles[fileNum], comm, rank);
    testInputAdapters<double, int, long>(0,0,0,mtxFiles[fileNum],comm, rank);
  }

  // Test input adapters created the Muelu Gallery.

  std::string nullString;

  testInputAdapters<double, int, long>(10, 20, 10, nullString, comm, rank);

  if (rank == 0)
    std::cout << "PASS" << std::endl;

  return 0;
}
