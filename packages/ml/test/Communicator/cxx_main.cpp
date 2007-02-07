/* 
 * Goal of this test:
 * - verify that ML is not using MPI functions outside the given communicator.
 *
 * This test will:
 * - Create a matrix on processor 0 only
 * - run ML on processor 0 only
 *
 * \date 27-Oct-05
 *
 * \author Marzio Sala, ETHZ/COLAB
 *
 */

#include <iostream>
#include <math.h>
#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_MPI) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_GALERI)

// epetra objects
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
// required to build the example matrix
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
// required by the linear system solver
#include "AztecOO.h"
// required by ML
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Galeri;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  Epetra_MpiComm SerialComm(MPI_COMM_SELF);

  if (Comm.MyPID() == 0) {

    ParameterList GaleriList;
    GaleriList.set("nx", 10);
    GaleriList.set("ny", 10);
    GaleriList.set("mx", 1);
    GaleriList.set("my", 1);
    
    Epetra_Map* Map = CreateMap("Cartesian2D", SerialComm, GaleriList);
    Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace2D", Map, GaleriList);

    Epetra_Vector LHS(A->OperatorDomainMap());
    Epetra_Vector RHS(A->OperatorRangeMap());

    LHS.PutScalar(0.0);
    RHS.Random();

    Epetra_LinearProblem problem(A, &LHS, &RHS);

    AztecOO solver(problem);

    ParameterList MLList;

    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("ML output", 0);

    ML_Epetra::MultiLevelPreconditioner* MLPrec = 
      new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

    solver.SetPrecOperator(MLPrec);
    solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
    solver.SetAztecOption(AZ_output, 32);
    solver.Iterate(500, 1e-12);

    // destroy the preconditioner
    delete MLPrec;

    delete A;
    delete Map;
  }

  MPI_Finalize();
  return(0);

} // main driver 

#else

#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  puts("This test requires:");
  puts("--enable-epetra");
  puts("--enable-aztecoo");
  puts("--enable-mpi");
  puts("--enable-ifpack");
  puts("--enable-teuchos");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // not to break tests
  return(EXIT_SUCCESS);
}

#endif
