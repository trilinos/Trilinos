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

#define Scalar float
#define LNO    int
#define GNO    long
#define Node   Zoltan2::default_node_t

int main(int argc, char *argv[])
{
  typedef Tpetra::Map<LNO, GNO, Node> map_t;
  typedef Tpetra::CrsMatrix<Scalar, LNO, GNO, Node> matrix_t;
  typedef Zoltan2::XpetraCrsMatrixInput<matrix_t> adapter_t;
  
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = session.getRank();
  int nprocs = session.getNProc();

  // create test input adapters using Muelu gallery.

  TestAdapters<Scalar,LNO,GNO> *input = 
    new TestAdapters<Scalar,LNO,GNO>(10, 10, 10, comm);

  RCP<adapter_t> adapter = input->getTpetraCrsMatrixInputAdapter();

  // The matrix created by the Muelu gallery
  RCP<matrix_t> M = input->getTpetraMatrix();

  // reassign matrix rows in a round robin fashion
  const RCP<const map_t> &rowMap = M->getRowMap();
  size_t localrows = M->getNodeNumRows();

  ArrayView<const GNO> gids = rowMap->getNodeElementList();
  Array<LNO> newpart(localrows);
  Array<LNO> localIds(localrows);

  for (LNO i=0; i < localrows; i++){
    newpart[i] = gids[i] % nprocs;
    localIds[i] = i;
  }

  matrix_t *newM = NULL;

  adapter->applyPartitioningSolution( *M, newM,
    *comm, LNO(localrows), nprocs, 
    gids.getRawPtr(), localIds.getRawPtr(), newpart.getRawPtr());
  
  // TODO - test that the import worked and newM is the same
  //   matrix as old M

  delete newM;

  delete input;

  if (rank == 0)
    std::cout << "PASS" << std::endl;

  return 0;
}
