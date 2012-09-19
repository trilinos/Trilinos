// Test to determine whether Teuchos::Comm's createSubcommunicator leaks memory.
// Since I don't see calls to MPI_Comm_free anywhere in Teuchos::Comm, I 
// suspect memory could be leaked when the subcommunicator falls out of scope.
// The RCPs are satisfied, but any memory allocated behind the scenes in MPI
// will be lost since MPI_Comm_free is not called.
// Running with purify confirms this result.

#include <stdio.h>
#include <mpi.h>
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"


int main(int narg, char **arg)
{
  Teuchos::GlobalMPISession mpiSession(&narg,&arg);

  Teuchos::RCP<const Teuchos::Comm<int> >
    comm = Teuchos::DefaultComm<int>::getComm();
  int me = comm->getRank();
  int np = comm->getSize();

  if (me == 0) 
    printf("Usage:  Zoltan2_teuchosSubcommTest.exe [#_of_subcomms_to_create]\n" 
           "        default number is 4\n");

  int niter = 4;
  if (narg > 1) niter = atoi(arg[1]);

  int *ids = NULL;
  ids = new int[np/2+1];
  ids[0] = me;
  for (int i = 1; i < np/2+1; i++) {
    ids[i] = (i != me ? i : 0);
  }
  Teuchos::ArrayView<const int> list(ids, np/2+1);
  for (int i = 0; i < niter; i++) {
    Teuchos::RCP<const Teuchos::Comm<int> > a 
                                            = comm->createSubcommunicator(list);
    printf("weak: %d  strong: %d total: %d\n",
            a.weak_count(), a.strong_count(), a.total_count());
  }
  delete [] ids;
  if (me == 0)
    printf("\nPASS\n");

  return 0;
}
