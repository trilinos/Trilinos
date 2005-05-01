#include "Amesos_ConfigDefs.h"
#ifdef HAVE_AMESOS_TRIUTILS

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Amesos.h"
#include "Amesos_TestRowMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include <vector>
using namespace Trilinos_Util;

// ====================================================================== 
// this function tests two things:
// - the return code from Amesos;
// - the true residual.
// The return value is `true' if both tests are passed, `false' otherwise.
//
// \author Marzio Sala, SNL 9214.
//
// \date Last updated on 21-Apr-04.
// ====================================================================== 

bool TestAmesos(bool quiet, char ProblemType[], Teuchos::ParameterList& AmesosList,
                bool UseTranspose, Epetra_RowMatrix* A, Epetra_MultiVector* lhs,
                Epetra_MultiVector* rhs)
{
  Epetra_LinearProblem Problem;
  Amesos A_Factory;

  Amesos_BaseSolver * Solver = A_Factory.Create(ProblemType, Problem);
  assert (Solver != 0);

  // Both sentences should work
  Solver->SetUseTranspose(UseTranspose);
  //AmesosList.set("UseTranspose",UseTranspose);

  Solver->SetParameters(AmesosList);

  // create a rhs corresponding to lhs or 1's
  lhs->PutScalar(1.0);
  A->Multiply(UseTranspose,*lhs,*rhs);
  lhs->PutScalar(0.0);

  Problem.SetOperator(A);


  if (Solver->SymbolicFactorization()) return(false);
  if (Solver->NumericFactorization()) return(false);

  // set sol and rhs here
  Problem.SetLHS(lhs);
  Problem.SetRHS(rhs);

  if (Solver->Solve()) return(false);

  // compute difference between exact solution and Amesos
  double d = 0.0, d_tot = 0.0;

  for (int i=0 ; i<lhs->Map().NumMyElements() ; ++i)
    for (int j = 0 ; j < lhs->NumVectors() ; ++j)
      d += ((*lhs)[j][i] - 1.0) * ((*lhs)[j][i] - 1.0);

  A->Comm().SumAll(&d,&d_tot,1);

  // compute ||Ax - b||
  vector<double> Norm(rhs->NumVectors());


  Epetra_MultiVector Ax(*rhs);
  A->Multiply(UseTranspose, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm[0]);

  string msg = ProblemType;

  if (!quiet && !A->Comm().MyPID()) 
  {
    cout << endl;
    cout << msg << " : Using " << A->Comm().NumProc() << " processes, UseTranspose = " << UseTranspose << endl;
    for (int j = 0 ; j < rhs->NumVectors() ; ++j)
      cout << msg << " : eq " << j 
	   << ", ||A x - b||_2 = " << Norm[j] << endl;
    cout << msg << " : ||x_exact - x||_2 = " << sqrt(d_tot) << endl;
  }

  if (Norm[0] > 1e-9)  
  {
    cerr << endl << msg << " WARNING : TEST FAILED!" << endl;
    return(false);
  }

  delete Solver;

  return(true);
}

// =========== //
// main driver //
// =========== //
//
int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool quiet = false ;
  if ( argc > 1 && argv[1][0] == '-' &&  argv[1][1] == 'q' ) quiet = true ;

  CrsMatrixGallery Gallery("recirc_2d", Comm);
  Gallery.Set("problem_size", 64);
  Gallery.Set("num_vectors", 2);

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();
  Epetra_RowMatrix* RowA = Problem->GetMatrix();
  // I wanna use this class to be usre that the solver cannot cast
  // to Epetra_CrsMatrix or Epetra_VbrMatrix.
  Amesos_TestRowMatrix A(RowA);

  Epetra_MultiVector* LHS = Problem->GetLHS();
  Epetra_MultiVector* RHS = Problem->GetRHS();

  Amesos Factory;  
  
  vector<string> SolverType;
  SolverType.push_back("Amesos_Lapack");
  SolverType.push_back("Amesos_Klu");
  SolverType.push_back("Amesos_Umfpack");
  SolverType.push_back("Amesos_Superlu");
  SolverType.push_back("Amesos_Superludist");  // Does not seem to work
  SolverType.push_back("Amesos_Mumps");
  SolverType.push_back("Amesos_Scalapack");
  // NOTE: DSCPACK does not support RowMatrix's.

  bool res;

  // If a given test fails, than the code stops, bue to the assert()
  // statement.
  for (unsigned int i = 0 ; i < SolverType.size() ; ++i) 
  {
    string Solver = SolverType[i];

    if (Factory.Query((char*)Solver.c_str())) 
    {
      if (1) {
	// solve with matrix
	Teuchos::ParameterList AmesosList;
	res = TestAmesos(quiet,(char*)Solver.c_str(), AmesosList, false, 
                         &A, LHS, RHS);
        assert (res == true);
      }
      if (1) {
	// solve transpose with matrix
	if (Solver != "Amesos_Superludist") {// still not implementes
	  Teuchos::ParameterList AmesosList;
	  res  = TestAmesos(quiet,(char*)Solver.c_str(), AmesosList, true, 
                            &A, LHS, RHS);
          assert (res == true);
	}
      }
    } 
    else
      if (!Comm.MyPID()) 
      {
	cerr << endl;
	cerr << "WARNING: SOLVER `" << Solver << "' NOT TESTED" << endl;
	cerr << endl;
      }
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(0); 
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


