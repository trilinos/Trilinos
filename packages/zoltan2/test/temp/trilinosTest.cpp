#include <Zoltan2_TestHelpers.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <iostream>

int main(int argc, char *argv[])
{
Teuchos::GlobalMPISession session(&argc, &argv);
RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
int nprocs = comm->getSize();
int rank = comm->getRank();

  gno_t baseId = 10000000;

  lno_t numLocalIds = 1000;

  gno_t *ids = new gno_t [numLocalIds];
  gno_t firstId = baseId + (rank * numLocalIds);
  for (lno_t i=0; i < numLocalIds; i++){
    ids[i] = firstId++;
  }

  Tpetra::Map<lno_t, gno_t> map1(
    numLocalIds*nprocs,
    ArrayView<gno_t>(ids, numLocalIds),
    baseId,
    comm);

  Tpetra::Map<lno_t, gno_t> map2(
    numLocalIds*nprocs,
    ArrayView<gno_t>(ids, numLocalIds),
    0,
    comm);

}
