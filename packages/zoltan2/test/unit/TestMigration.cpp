// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Test migration of a Tpetra::CrsMatrix with through the input adapter.
//
// For now testing that it compiles.  TODO - test that it's correct. 

#include <TestAdapters.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::OrdinalTraits;
using Teuchos::ScalarTraits;

#define Scalar float
#define LNO    int
#define GNO    long

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = session.getRank();
  int nprocs = session.getNProc();

  TestAdapters<Scalar,LNO,GNO> *input = 
    new TestAdapters<Scalar,LNO,GNO>(10, 10, 10, comm);

  typedef Zoltan2::default_node_t node_t;
  typedef Tpetra::CrsMatrix<Scalar, LNO, GNO, node_t> tmatrix_t;
  typedef Zoltan2::XpetraCrsMatrixInput<tmatrix_t> adapter_t;

  RCP<adapter_t> adapter = input->getTpetraCrsMatrixInputAdapter();

#if 0
LNO testval=5;
Array<GNO> vals(5, 0);
Array<LNO> lnovals(5, 0);
adapter->simpleFunc(testval, &vals[2], lnovals.getRawPtr(),
     &lnovals[3]);
#endif

  RCP<tmatrix_t> M = input->getTpetraMatrix();

  // reassign matrix rows in a round robin fashion

  size_t nrows = M->getGlobalNumRows();

  const GNO *rowIds, *colIds;
  const LNO *localIds, *offsets;

  // Since local ids are consecutive beginning with 0, adapter returns
  // with localIds=NULL.

  size_t n = adapter->getRowListView(rowIds, localIds, offsets, colIds);

  int *newpartList = new int [n];
  int *localIdList = new int [n];

  for (LNO i=0; i < n; i++){
    newpartList[i] = i % nprocs;
    localIdList[i] = i;
  }

  tmatrix_t *N=NULL;

  try{
    N = adapter->applyPartitioningSolution(*M, n, nprocs, rowIds, 
      localIdList, newpartList);
  }
  catch (std::exception &e){
    std::cerr << "broken" << e.what() << std::endl;
  }

  // TODO - verify that the N == M

  if (rank == 0)
    std::cout << "PASS" << std::endl;

  return 0;
}
