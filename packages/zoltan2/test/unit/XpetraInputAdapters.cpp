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
//

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

template <typename Scalar, typename LNO, typename GNO>
  void testInputAdapters(
    std::string fname, RCP<const Comm<int> > comm, int rank)
{
  bool includeEpetra = false;
  std::string intName(OrdinalTraits<int>::name());
  std::string ScalarName(ScalarTraits<Scalar>::name());
  std::string LNOName(OrdinalTraits<LNO>::name());
  std::string GNOName(OrdinalTraits<GNO>::name());

  if ((LNOName == intName) && (GNOName == intName))
    includeEpetra = true;

  RCP<Zoltan2::EpetraCrsMatrixInput> emi;
  RCP<Zoltan2::EpetraCrsGraphInput> egi;

  // Create some input adapters for xpetra objects

  TestAdapters<Scalar,LNO,GNO> input(fname);

#ifdef HAVE_MPI
  input.setMpiCommunicator(MPI_COMM_WORLD);
#endif

  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LNO, GNO> > M = input.getMatrix();
  size_t nRows = M->getNodeNumRows();
  size_t ngRows = M->getGlobalNumRows();
  size_t nEntries = M->getNodeNumEntries();
  size_t ngEntries = M->getGlobalNumEntries();
  size_t nColumns = M->getNodeNumCols();
  size_t ngColumns = M->getGlobalNumCols();

  if (includeEpetra){
    emi = input.getEpetraCrsMatrixInputAdapter();
    egi = input.getEpetraCrsGraphInputAdapter();
  }

  RCP<Zoltan2::TpetraCrsMatrixInput<Scalar, LNO, GNO> > tmi = 
    input.getTpetraCrsMatrixInputAdapter();

  RCP<Zoltan2::TpetraCrsGraphInput<LNO, GNO> > tgi = 
    input.getTpetraCrsGraphInputAdapter();

  RCP<Zoltan2::XpetraCrsMatrixInput<Scalar, LNO, GNO> > xmi = 
    input.getXpetraCrsMatrixInputAdapter();

  RCP<Zoltan2::XpetraCrsGraphInput<LNO, GNO> > xgi = 
    input.getXpetraCrsGraphInputAdapter();


  ///////////////////////////////////////////////////////
  // Test graph adapters using graph input interface
  ///////////////////////////////////////////////////////

  size_t num;

  if (includeEpetra){
    num = egi->getLocalNumVertices();
    TEST_FAIL_AND_EXIT(*comm, num==nRows, "egi->getLocalNumVertices", 1);
    num = egi->getGlobalNumVertices();
    TEST_FAIL_AND_EXIT(*comm, num==ngRows, "egi->getGlobalNumVertices", 1);
    num = egi->getLocalNumEdges();
    TEST_FAIL_AND_EXIT(*comm, num==nEntries, "egi->getLocalNumEdges", 1);
    num = egi->getGlobalNumEdges();
    TEST_FAIL_AND_EXIT(*comm, num==ngEntries, "egi->getGlobalNumEdges", 1);
  }

  num = tgi->getLocalNumVertices();
  TEST_FAIL_AND_EXIT(*comm, num==nRows, "tgi->getLocalNumVertices", 1);
  num = tgi->getGlobalNumVertices();
  TEST_FAIL_AND_EXIT(*comm, num==ngRows, "tgi->getGlobalNumVertices", 1);
  num = tgi->getLocalNumEdges();
  TEST_FAIL_AND_EXIT(*comm, num==nEntries, "tgi->getLocalNumEdges", 1);
  num = tgi->getGlobalNumEdges();
  TEST_FAIL_AND_EXIT(*comm, num==ngEntries, "tgi->getGlobalNumEdges", 1);

  // TODO the rest of the methods

  ///////////////////////////////////////////////////////
  // Test matrix adapters using matrix input interface
  ///////////////////////////////////////////////////////

  if (includeEpetra){
    num = emi->getLocalNumRows();
    TEST_FAIL_AND_EXIT(*comm, num==nRows, "emi->getLocalNumRows", 1);
    num = emi->getGlobalNumRows();
    TEST_FAIL_AND_EXIT(*comm, num==ngRows, "emi->getGlobalNumRows", 1);
    num = emi->getLocalNumColumns();
    TEST_FAIL_AND_EXIT(*comm, num==nColumns, "emi->getLocalNumColumns", 1);
    num = emi->getGlobalNumColumns();
    TEST_FAIL_AND_EXIT(*comm, num==ngColumns, "emi->getGlobalNumColumns", 1);
  }

  num = tmi->getLocalNumRows();
  TEST_FAIL_AND_EXIT(*comm, num==nRows, "tmi->getLocalNumRows", 1);
  num = tmi->getGlobalNumRows();
  TEST_FAIL_AND_EXIT(*comm, num==ngRows, "tmi->getGlobalNumRows", 1);
  num = tmi->getLocalNumColumns();
  TEST_FAIL_AND_EXIT(*comm, num==nColumns, "tmi->getLocalNumColumns", 1);
  num = tmi->getGlobalNumColumns();
  TEST_FAIL_AND_EXIT(*comm, num==ngColumns, "tmi->getGlobalNumColumns", 1);

  // TODO the rest of the methods

  // TODO test the Xpetra adapters too

  if (!rank){
    std::cout << "Processed " << fname << ": ";
    std::cout << ngRows << " vertices, ";
    std::cout << ngEntries << " edges." << std::endl;
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
  int rank = session.getRank();

  std::vector<std::string> mtxFiles;
  
  mtxFiles.push_back("../data/simple.mtx");
  mtxFiles.push_back("../data/cage10.mtx");

  // To use this matrix we would need to pass a domain map
  // to FillComplete.  So we skip it for now.  TODO
  // mtxFiles.push_back("../data/diag500_4.mtx");

  for (unsigned int fileNum=0; fileNum < mtxFiles.size(); fileNum++){
    testInputAdapters<double, int, int>(mtxFiles[fileNum], comm, rank);
    testInputAdapters<double, int, long>(mtxFiles[fileNum], comm, rank);
  }


  return 0;
}
