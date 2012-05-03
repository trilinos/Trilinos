#include <Zoltan2_TestHelpers.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <iostream>

#ifdef SHOW_ZOLTAN2_LINUX_MEMORY
extern "C"{
static char *z2_meminfo=NULL;
extern void Zoltan_get_linux_meminfo(char *msg, char **result);
}
#endif

int main(int argc, char *argv[])
{
Teuchos::GlobalMPISession session(&argc, &argv);
RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
int nprocs = comm->getSize();
int rank = comm->getRank();

#ifdef SHOW_ZOLTAN2_LINUX_MEMORY
  if (rank==0){
    Zoltan_get_linux_meminfo("Starting ", &z2_meminfo);
    if (z2_meminfo){
      std::cout << z2_meminfo << std::endl;
      free(z2_meminfo);
      z2_meminfo=NULL;
    }
  }
#endif

  gno_t baseId = 10000000;

  lno_t numLocalIds = 1000;

  gno_t *ids = new gno_t [numLocalIds];
  gno_t firstId = baseId + (rank * numLocalIds);
  for (lno_t i=0; i < numLocalIds; i++){
    ids[i] = firstId++;
  }

#ifdef SHOW_ZOLTAN2_LINUX_MEMORY
  if (rank==0){
    Zoltan_get_linux_meminfo("After creating global id list ", &z2_meminfo);
    if (z2_meminfo){
      std::cout << z2_meminfo << std::endl;
      free(z2_meminfo);
      z2_meminfo=NULL;
    }
  }
#endif

  Tpetra::Map<lno_t, gno_t> map1(
    numLocalIds*nprocs,
    ArrayView<gno_t>(ids, numLocalIds),
    baseId,
    comm);

#ifdef SHOW_ZOLTAN2_LINUX_MEMORY
  if (rank==0){
    Zoltan_get_linux_meminfo("After creating map with base 10000000 ", &z2_meminfo);
    if (z2_meminfo){
      std::cout << z2_meminfo << std::endl;
      free(z2_meminfo);
      z2_meminfo=NULL;
    }
  }
#endif


  

  Tpetra::Map<lno_t, gno_t> map2(
    numLocalIds*nprocs,
    ArrayView<gno_t>(ids, numLocalIds),
    0,
    comm);

#ifdef SHOW_ZOLTAN2_LINUX_MEMORY
  if (rank==0){
    Zoltan_get_linux_meminfo("After creating map with base 0 ", &z2_meminfo);
    if (z2_meminfo){
      std::cout << z2_meminfo << std::endl;
      free(z2_meminfo);
      z2_meminfo=NULL;
    }
  }
#endif


}
