#include "Amesos_ConfigDefs.h"
#ifdef HAVE_AMESOS_TRIUTILS

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Amesos.h"
#include "Teuchos_ParameterList.hpp"
#include "Trilinos_Util_CrsMatrixGallery.h"

void TestAmesos(char ProblemType[],
		Teuchos::ParameterList & AmesosList,
		bool UseTranspose,
		Epetra_CrsMatrix* A,
		Epetra_MultiVector* lhs,
		Epetra_MultiVector* rhs,
		double & TotalErrorResidual,
		double & TotalErrorExactSol)
{

  // create an empty linear problem
  Epetra_LinearProblem Problem;

  Amesos A_Factory;

  Amesos_BaseSolver * Solver = A_Factory.Create(ProblemType, Problem);

  assert (Solver != 0);

  // Both sentences should work
  //    Solver->SetUseTranspose(UseTranspose);
  AmesosList.set("UseTranspose",UseTranspose);

  Solver->SetParameters(AmesosList);

  // create a rhs corresponding to lhs or 1's
  lhs->PutScalar(1.0);
  A->Multiply(UseTranspose,*lhs,*rhs);

  Epetra_Time Time(A->Comm());
  Epetra_Time StartTime(A->Comm());

  // set the problem here
  Problem.SetOperator(A);

  Solver->SymbolicFactorization();
  double TimeForSymbolicFactorization = Time.ElapsedTime();

  Time.ResetStartTime();
  Solver->NumericFactorization();
  double TimeForNumericFactorization = Time.ElapsedTime();

  // set sol and rhs here
  Problem.SetLHS(lhs);
  Problem.SetRHS(rhs);

  Time.ResetStartTime();
  Solver->Solve();
  double TimeForSolve = Time.ElapsedTime();

  // compute difference between exact solution and Amesos
  double d = 0.0, d_tot = 0.0;

  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    for (int j = 0 ; j < lhs->NumVectors() ; ++j)
      d += ((*lhs)[j][i] - 1.0) * ((*lhs)[j][i] - 1.0);

  A->Comm().SumAll(&d,&d_tot,1);

  // compute ||Ax - b||
  double* Norm;
  Norm = new double[rhs->NumVectors()];

  Epetra_MultiVector Ax(*rhs);
  A->Multiply(UseTranspose, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(Norm);

  string msg = ProblemType;

  if( A->Comm().MyPID() == 0 ) {
    cout << msg << "......Using " << A->Comm().NumProc() << " processes, UseTranspose = " << UseTranspose << endl;
    for (int j = 0 ; j < rhs->NumVectors() ; ++j)
      cout << msg << "...... eq " << j 
	   << ", ||A x - b||_2 = " << Norm[j] << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Time for Symbolic Factorization = " << TimeForSymbolicFactorization << endl;
    cout << msg << "......Time for Numeric Factorization  = " << TimeForNumericFactorization << endl;
    cout << msg << "......Time for Solve                  = " << TimeForSolve << endl;
    cout << msg << "......Total Time = " << StartTime.ElapsedTime() << endl;

    if( Norm[0] > 1e-9 ) {
      cerr << endl << msg << " WARNING : TEST FAILED!" << endl << endl;
    }
  }

  TotalErrorExactSol += sqrt(d_tot);
  for (int j = 0 ; j < rhs->NumVectors() ; ++j)
    TotalErrorResidual += Norm[j];

  delete Solver;
  delete [] Norm;

  return;
  
}

using namespace Trilinos_Util;

#include <vector>

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // =========================================== //
  // non-symmetric case (for all except DSCPACK) //
  // =========================================== //

  CrsMatrixGallery Gallery("recirc_2d", Comm);
  Gallery.Set("problem_size", 900);
  Gallery.Set("num_vectors", 5);

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();
  Epetra_RowMatrix* RowA = Problem->GetMatrix();
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(RowA);

  Epetra_MultiVector* LHS = Problem->GetLHS();
  Epetra_MultiVector* RHS = Problem->GetRHS();

  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;
  Amesos Factory;  
  
  vector<string> SolverType;
  SolverType.push_back("Amesos_Lapack");
  SolverType.push_back("Amesos_Klu");
  SolverType.push_back("Amesos_Umfpack");
  SolverType.push_back("Amesos_Superlu");
  SolverType.push_back("Amesos_Superludist");
  SolverType.push_back("Amesos_Mumps");

  for (int i = 0 ; i < SolverType.size() ; ++i) {
    string Solver = SolverType[i];
    if (Factory.Query((char*)Solver.c_str())) {
      {
	// solve with matrix
	Teuchos::ParameterList AmesosList;
	TestAmesos((char*)Solver.c_str(), AmesosList, false, A, LHS, RHS,
		 TotalErrorResidual, TotalErrorExactSol );
      }
      {
	// solve with matrix
	if (Solver != "Amesos_Superludist") {// still not implementes
	  Teuchos::ParameterList AmesosList;
	  TestAmesos((char*)Solver.c_str(), AmesosList, false, A, LHS, RHS,
		     TotalErrorResidual, TotalErrorExactSol );
	}
      }
    } else
      if (Comm.MyPID() == 0) {
	cerr << endl;
	cerr << "WARNING: SOLVER `" << Solver << "' NOT TESTED" << endl;
	cerr << endl;
      }
  }
   
  // =============================== //
  // symmetric problem (for DSCPACK) //
  // =============================== //

  CrsMatrixGallery SymGallery("laplace_3d", Comm);
  SymGallery.Set("problem_size", 8000);
  SymGallery.Set("num_vectors", 5);

  Epetra_LinearProblem* SymProblem = SymGallery.GetLinearProblem();
  Epetra_RowMatrix* SymRowA = SymProblem->GetMatrix();
  Epetra_CrsMatrix* SymA = dynamic_cast<Epetra_CrsMatrix*>(SymRowA);

  Epetra_MultiVector* SymLHS = SymProblem->GetLHS();
  Epetra_MultiVector* SymRHS = SymProblem->GetRHS();

  if (Factory.Query("Amesos_Dscpack")) {
    
    Teuchos::ParameterList AmesosList;
    TestAmesos("Amesos_Dscpack", AmesosList, false, SymA, SymLHS, SymRHS,
	       TotalErrorResidual, TotalErrorExactSol );
  } else
    if (Comm.MyPID() == 0) {
      cerr << endl;
      cerr << "WARNING: SOLVER `Amesos_Dscpack' NOT TESTED" << endl;
      cerr << endl;
    }

  // print out total error
  
  if( Comm.MyPID() == 0 ) {
    cout << endl;
    cout << "......Total error for residual = " << TotalErrorResidual << endl;
    cout << "......Total error for exact solution  = " << TotalErrorExactSol << endl;
    cout << endl;
 }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (TotalErrorResidual < 1e-9) 
    return( EXIT_SUCCESS );
  else
    return( EXIT_FAILURE );

}

#else

// Triutils is not available. Sorry, we have to give up.

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#else
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  puts("Please configure AMESOS with --enable-triutils");
  puts("to run this example");
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(0);
}

#endif


