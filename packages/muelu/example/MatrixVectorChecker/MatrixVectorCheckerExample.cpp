#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_ParameterList.hpp>
#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>

#include "MatrixVectorChecker.hpp"

using namespace Teuchos;

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", 10 * Comm.NumProc());
  galeriList.set("ny", 10);

  galeriList.set("mx", Comm.NumProc());
  galeriList.set("my", 1);

  RCP<Epetra_Map> map = rcp( Galeri::CreateMap("Cartesian2D", Comm, galeriList) );
  RCP<Epetra_CrsMatrix> matrix = rcp( Galeri::CreateCrsMatrix("Laplace2D", map.get(), galeriList) );

  MatrixVectorChecker(*matrix);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);

}
