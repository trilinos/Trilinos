#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"


#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include "Teuchos_ParameterList.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"

#include "Trilinos_Util_CrsMatrixGallery.h"

void PrintLine() 
{
  cout << endl;
  for( int i=0 ; i<80 ; ++i )
    cout << "=";
  cout << endl;
  cout << endl;
  
  return;
}

int TestMultiLevelPreconditioner(char ProblemType[],
				 Teuchos::ParameterList & MLList,
				 Epetra_LinearProblem & Problem, double & TotalErrorResidual,
				 double & TotalErrorExactSol)
{
  
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();
  Epetra_RowMatrix* A = Problem.GetMatrix();
  
  // ======================================== //
  // create a rhs corresponding to lhs or 1's //
  // ======================================== //
  
  lhs->PutScalar(1.0);
  A->Multiply(false,*lhs,*rhs);

  lhs->PutScalar(0.0);
  
  Epetra_Time Time(A->Comm());
  
  // =================== //
  // call ML and AztecOO //
  // =================== //
  
  AztecOO solver(Problem);
  
  MLList.set("output", 0);

  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);
  
  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.SetAztecOption(AZ_kspace, 160);
  
  solver.Iterate(1550, 1e-12);
  
  delete MLPrec;
  
  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //
  
  double d = 0.0, d_tot = 0.0;
  
  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - 1.0) * ((*lhs)[0][i] - 1.0);
  
  A->Comm().SumAll(&d,&d_tot,1);
  
  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A->Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);
  
  string msg = ProblemType;
  
  if (A->Comm().MyPID() == 0) {
    cout << msg << "......Using " << A->Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }
  
  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;
  
  return( solver.NumIters() );
  
}

using namespace Trilinos_Util;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // initialize the random number generator

  srandom((unsigned int)1);

  // ===================== //
  // create linear problem //
  // ===================== //

  CrsMatrixGallery Gallery("laplace_3d", Comm);
  Gallery.Set("problem_size", 27000);

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();
  Teuchos::ParameterList MLList;
  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;

  // ====================== //
  // default options for SA //
  // ====================== //

  PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Gauss-Seidel");
  TestMultiLevelPreconditioner("SA", MLList, *Problem, 
                               TotalErrorResidual, TotalErrorExactSol);

  // ============================== //
  // default options for SA, Jacobi //
  // ============================== //

  PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "Jacobi");

  TestMultiLevelPreconditioner("SA", MLList, *Problem, TotalErrorResidual, 
                               TotalErrorExactSol);

  // =========================== //
  // default options for SA, MLS //
  // =========================== //

  PrintLine();

  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type", "MLS");

  TestMultiLevelPreconditioner("SA", MLList, *Problem, 
                               TotalErrorResidual, TotalErrorExactSol);

  // ===================== //
  // print out total error //
  // ===================== //

  if (Comm.MyPID() == 0) {
    cout << endl;
    cout << "......Total error for residual        = " << TotalErrorResidual << endl;
    cout << "......Total error for exact solution  = " << TotalErrorExactSol << endl;
    cout << endl;
  }

  if (TotalErrorResidual > 1e-8) {
    cerr << "Error: `MultiLevelPrecoditioner_Sym.exe' failed!" << endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (Comm.MyPID() == 0)
    cerr << "`MultiLevelPrecoditioner_Sym.exe' passed!" << endl;

  return (EXIT_SUCCESS);

}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
  // still need to deal with MPI, some architecture don't like
  // an exit(0) without MPI_Finalize()
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
