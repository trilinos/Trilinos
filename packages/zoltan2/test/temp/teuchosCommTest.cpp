#include <mpi.h>
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"

typedef int Ordinal;

int main(int narg, char **arg)
{
  Teuchos::GlobalMPISession mpiSession(&narg,&arg);

  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();

  int niter = 4000;
  if (narg > 1) niter = atoi(arg[1]);

  int me = comm->getRank();
  double tstart, tend;

  int iin = me, iout;
  double din = me * 2., dout;

  tstart = MPI_Wtime();
  for (int i = 0; i < niter; i++) {
    reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &iin, &iout);
    reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &din, &dout);
  }
  tend = MPI_Wtime();
  if (me == 0)
    std::cout << "reduceAll time using Teuchos::Comm = " 
              << tend - tstart << std::endl;
    
  tstart = MPI_Wtime();
  for (int i = 0; i < niter; i++) {
    MPI_Allreduce(&iin, &iout, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&din, &dout, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  tend = MPI_Wtime();
  if (me == 0)
    std::cout << "Allreduce time using MPI_Allreduce = " 
              << tend - tstart << std::endl;

  if (me == 0)
    std::cout << "PASS" << std::endl;

  return 0;
}
