/* 
 * Goal of this test:
 * - compare the two converters from ML_Operator to Epetra_RowMatrix.
 *   ML offers two ways to convert from ML_Operator to an Epetra_RowMatrix
 *   derived class. The first way is ML_Operator2EpetraCrsMatrix, which
 *   creates a new Epetra_CrsMatrix, and fills it with the elements of the
 *   input ML_Operator. This is an expensive conversion.
 *   The other conversion is given by class ML_Epetra::RowMatrix, that defines
 *   suitable wraps from the ML data format, to the Epetra data format.
 *
 * This test will:
 * - create an Aztec matrix;
 * - convert this matrix into ML_Operator;
 * - create an Epetra_CrsMatrix;
 * - create an ML_Epetra::RowMatrix;
 * - multiply those two matrices by a random vector, and compare
 *   the results.
 *
 * \date 29-Aug-04
 *
 * \author Marzio Sala, SNL 9214
 *
 */

#include <iostream>
#include <math.h>
#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_MPI) && defined(HAVE_ML_IFPACK)

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
#include "Trilinos_Util_CrsMatrixGallery.h"
// required by the linear system solver
#include "AztecOO.h"
// required by ML
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Trilinos_Util;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  Epetra_MpiComm SerialComm(MPI_COMM_SELF);

  if (Comm.MyPID() == 0) {

    CrsMatrixGallery Gallery("laplace_2d", SerialComm);
    Gallery.Set("problem_size", 100);
    Epetra_RowMatrix* A = Gallery.GetMatrix();
    //Epetra_CrsMatrix* CrsA = dynamic_cast<Epetra_CrsMatrix*>(A);

    Epetra_Vector LHS(A->OperatorDomainMap());
    Epetra_Vector RHS(A->OperatorRangeMap());

    LHS.PutScalar(0.0);
    RHS.Random();

    Epetra_LinearProblem problem(A, &LHS, &RHS);

    AztecOO solver(problem);

    ParameterList MLList;

    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("output", 0);

    ML_Epetra::MultiLevelPreconditioner* MLPrec = 
      new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

    solver.SetPrecOperator(MLPrec);
    solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
    solver.SetAztecOption(AZ_output, 32);
    solver.Iterate(500, 1e-12);

    // destroy the preconditioner
    delete MLPrec;

  }

  MPI_Finalize();
  return(0);

} // main driver 

#else

int main(int argc, char *argv[])
{
  puts("This test requires:");
  puts("--enable-epetra");
  puts("--enable-aztecoo");
  puts("--enable-mpi");
  puts("--enable-ifpack");
  puts("--enable-teuchos");

  // not to break tests
  return(EXIT_SUCCESS);
}

#endif
