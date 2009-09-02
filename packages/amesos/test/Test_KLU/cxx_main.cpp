#include "Amesos_ConfigDefs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Amesos_Klu.h"
#include "Amesos_TestRowMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"

using namespace Galeri;

//============ //
// main driver //
//============ //

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::ParameterList GaleriList;
  GaleriList.set("n", 5);

  Epetra_Map* Map = CreateMap("Random", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix = CreateCrsMatrix("Minij", Map, GaleriList);
 
  int NumVectors = 2;
  Amesos_TestRowMatrix A(Matrix);
  Epetra_MultiVector x(*Map,NumVectors);
  Epetra_MultiVector x_exact(*Map,NumVectors);
  Epetra_MultiVector b(*Map,NumVectors);
  x_exact.Random();
  A.Multiply(false,x_exact,b);

  // =========== //
  // AMESOS PART //
  // =========== //

  Epetra_LinearProblem Problem(&A, &x, &b);
  Amesos_Klu Solver(Problem);

  AMESOS_CHK_ERR(Solver.SymbolicFactorization());
  AMESOS_CHK_ERR(Solver.NumericFactorization());
  AMESOS_CHK_ERR(Solver.Solve());

  double norm = ComputeNorm(Matrix, &x_exact, &b);
  if (Comm.MyPID() == 0)
    std::cout << "norm = " << norm << std::endl;

  if (norm > 1e-5)
    exit(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
