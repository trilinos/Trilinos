// Requires MPI
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Xpetra_TpetraCrsGraph.hpp>

using namespace Teuchos;

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  Teuchos::RCP<Teuchos::OpaqueWrapper<MPI_Comm> >handle =
    Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD);
  Teuchos::RCP<Teuchos::MpiComm<int> > comm(new Teuchos::MpiComm<int>(handle));

  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int numVtx = 100;
  int globalNumVtx = 100*nprocs;
  int gidOffset = rank * numVtx;
  int base=0;

  RCP<Tpetra::Map<int, int> > tmap = rcp(
    new Tpetra::Map<int, int>(globalNumVtx, numVtx, base, comm));

  ArrayRCP<size_t> entriesPerRow(numVtx);

  for (int i=0; i < numVtx; i++){
    int numEdges=0;
    for (int j=gidOffset + i-2; j <= gidOffset + i+2; j++)
      if ( (j >= 0) && (j < globalNumVtx))
        numEdges++; 
    entriesPerRow[i] = numEdges;
  }

  RCP<Tpetra::CrsGraph<int, int> > tpetraVersion =
    rcp(new Tpetra::CrsGraph<int, int>(tmap, entriesPerRow));

  ArrayRCP<int> nonZeroIds(5);

  for (int i=0; i < numVtx; i++){
    int numEdges=0;
    for (int j=gidOffset + i-2; j <= gidOffset + i+2; j++)
      if ( (j >= 0) && (j < globalNumVtx))
        nonZeroIds[numEdges++] = j;

    tpetraVersion->insertGlobalIndices(i + gidOffset, 
      nonZeroIds.view(0, numEdges));
  }

  tpetraVersion->fillComplete();

  Xpetra::TpetraCrsGraph<int, int> xgraph(tpetraVersion); 
}
