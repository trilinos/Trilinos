#include <iostream>
#include <math.h>
#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_MPI) && defined(HAVE_ML_IFPACK)

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
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "AztecOO.h"
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_Ifpack_ML.h"

using namespace Teuchos;
using namespace Trilinos_Util;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 100);
  Epetra_RowMatrix* A = Gallery.GetMatrix();
  //Epetra_CrsMatrix* CrsA = dynamic_cast<Epetra_CrsMatrix*>(A);

  Epetra_Vector LHS(A->OperatorDomainMap());
  Epetra_Vector RHS(A->OperatorRangeMap());

  LHS.PutScalar(0.0);
  RHS.Random();

  Epetra_LinearProblem problem(A, &LHS, &RHS);

  AztecOO solver(problem);

  Teuchos::ParameterList List;
  ML_Epetra::SetDefaults("SA", List);
  List.set("output", 0);
  List.set("schwarz: combine mode", Add);
  List.set("cycle applications", 10);

  int Overlap = 2;
  Ifpack_AdditiveSchwarz<ML_Epetra::Ifpack_ML> Prec(A, Overlap);
  Prec.SetParameters(List);
  ML_CHK_ERR(Prec.Initialize());
  ML_CHK_ERR(Prec.Compute());

  solver.SetPrecOperator(&Prec);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

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
