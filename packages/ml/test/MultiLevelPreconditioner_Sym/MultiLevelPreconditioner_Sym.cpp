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

static bool verbose = false;
static bool ml_verbose = false;

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
  
  Epetra_MultiVector   * lhs     = Problem.GetLHS();
  Epetra_MultiVector   * rhs     = Problem.GetRHS();
  Epetra_RowMatrix     * A       = Problem.GetMatrix();
  
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
  
  if (ml_verbose)
    MLList.set("output", 10);
  else
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
  
  if( A->Comm().MyPID() == 0 ) {
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

  if (Comm.MyPID() == 0)
    verbose = true;

  if (argc >= 2)
    if (strcmp(argv[1], "-v") == 0) 
      ml_verbose = true;

  // initialize the random number generator

  srandom((unsigned int)1);

  int NumProcs = Comm.NumProc();
  
  // ===================== //
  // create linear problem //
  // ===================== //
  
  CrsMatrixGallery Gallery("laplace_3d", Comm);
  Gallery.Set("problem_size", 27000);
  
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  int TotalFailed = 0;
  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;

  int iters;
  
  // ====================== //
  // default options for SA //
  // ====================== //

  if (true) {

    if (verbose) PrintLine();

    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("SA",MLList);
    iters = TestMultiLevelPreconditioner("SA", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    if (verbose) {
#ifdef HAVE_ML_AMESOS
      // expected iterations
      switch( NumProcs ) {
      case 1:
        if( iters != 16 ) {
          ++TotalFailed;
          cout << endl << "### TEST FAILED : expecting 16 iterations, got "
            << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      case 4:
        if( iters != 16 ) {
          ++TotalFailed;
          cout << endl << "### TEST FAILED : expecting 16 iterations, got "
            << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      }
#else
      cout << "### Cannot check the number of iterations" << endl;
      cout << "### (no `enable-amesos')" << endl;
      cout << "### Checking linear system residual only (at the bottom of the file)" << endl;
#endif
    }
  }
  
  // ============================== //
  // default options for SA, Jacobi //
  // ============================== //

  if (true) {

    if (verbose) PrintLine();

    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("smoother: type", "Jacobi");

    iters = TestMultiLevelPreconditioner("SA", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    if (verbose) {
#ifdef HAVE_ML_AMESOS
      // expected iterations
      switch( NumProcs ) {
      case 1:
        if( iters != 17 ) {
          ++TotalFailed;
          cout << endl << "### TEST FAILED : expecting 17 iterations, got "
            << iters << endl << endl;
        }
        else
          cout << "### TEST PASSED" << endl;
        break;
      case 4:
        if( iters != 17 ) {
          ++TotalFailed;
          cout << endl << "### TEST FAILED : expecting 17 iterations, got "
            << iters << endl << endl;
        }
        else
          cout << "### TEST PASSED" << endl;
        break;
      }
#else
      cout << "### Cannot check the number of iterations" << endl;
      cout << "### (no `enable-amesos')" << endl;
      cout << "### Checking linear system residual only (at the bottom of the file)" << endl;
#endif
    }
  }

  // =========================== //
  // default options for SA, MLS //
  // =========================== //

  if (true) {

    if (verbose) PrintLine();

    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("smoother: type", "MLS");

    iters = TestMultiLevelPreconditioner("SA", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    if (verbose) {
#ifdef HAVE_ML_AMESOS
      // expected iterations
      switch( NumProcs ) {
      case 1:
        if( iters != 14 ) {
          ++TotalFailed;
            cout << endl << "### TEST FAILED : expecting 14 iterations, got "
                 << iters << endl << endl;
        }
        else
          cout << "### TEST PASSED" << endl;
        break;
      case 4:
        if( iters != 14 ) {
          ++TotalFailed;
            cout << endl << "### TEST FAILED : expecting 14 iterations, got "
                 << iters << endl << endl;
        }
        else
          cout << "### TEST PASSED" << endl;
        break;
      }
#else
      cout << "### Cannot check the number of iterations" << endl;
      cout << "### (no `enable-amesos')" << endl;
      cout << "### Checking linear system residual only (at the bottom of the file)" << endl;
#endif
    }
  }

  // ===================== //
  // print out total error //
  // ===================== //

  if( Comm.MyPID() == 0 ) {
    cout << endl;
    cout << "......Total error for residual        = " << TotalErrorResidual << endl;
    cout << "......Total error for exact solution  = " << TotalErrorExactSol << endl;
    cout << "......Total # of failed tests         = " << TotalFailed << endl;
    cout << endl;
  }

  if (TotalErrorResidual > 1e-8) exit(EXIT_FAILURE);

#ifdef HAVE_ML_AMESOS
  if (TotalFailed) exit(EXIT_FAILURE);
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (verbose)
    cout << "### ALL TESTS PASSED" << endl;

  return( EXIT_SUCCESS );

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
