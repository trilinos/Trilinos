// This example compares all the available Amesos solvers
// for the solution of the same linear system. 
//
// The example can be run in serial and in parallel.

#include "Amesos_ConfigDefs.h"
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
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Trilinos_Util;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = (Comm.MyPID() == 0);
  double TotalResidual = 0.0;

  // NOTE: matrix must be symmetric if DSCPACK is used.
  //       Please refer to the Trilinos tutorial (chapter on triutils)
  //       for an overview of available matrices.
  char* ProblemType = "laplace_3d";
  int NumGlobalRows = 1000; // must be a perfect cube
  
  CrsMatrixGallery Gallery(ProblemType, Comm);
  Gallery.Set("problem_size",NumGlobalRows);
  Gallery.Set("num_vectors",1);

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();
  assert (Problem != 0);
  // retrive pointers to solution (lhs), right-hand side (rhs)
  // and matrix itself (A)
  Epetra_MultiVector* lhs = Problem->GetLHS();
  Epetra_MultiVector* rhs = Problem->GetRHS();
  Epetra_RowMatrix* A = Problem->GetMatrix();
    
  // use this list to set up parameters, now it is required
  // to explicitly compute the residual after solution
  Teuchos::ParameterList List;
  List.set("ComputeTruResidual",true);

  vector<string> SolverType;
  SolverType.push_back("Amesos_Lapack");
  SolverType.push_back("Amesos_Klu");
  SolverType.push_back("Amesos_Umfpack");
  SolverType.push_back("Amesos_Superlu");
  SolverType.push_back("Amesos_Superludist");
  SolverType.push_back("Amesos_Mumps");
  SolverType.push_back("Amesos_Dscpack");

  Epetra_Time Time(Comm);
  
  // this is the Amesos factory object that will create 
  // a specific Amesos solver.
  Amesos Factory;

  // Cycle over all solvers.
  // Only installed solvers will be tested.
  for (int i = 0 ; i < SolverType.size() ; ++i) {

    // Check whether the solver is available or not
    if (Factory.Query(SolverType[i])) {
 
      // 1.- set exact solution (constant vector)
      lhs->PutScalar(1.0);
 
      // 2.- create corresponding rhs
      A->Multiply(false,*lhs,*rhs);
 
      // 3.- randomize solution vector
      lhs->Random();
 
      // 4.- create the amesos solver object
      Amesos_BaseSolver* Solver = Factory.Create(SolverType[i], *Problem);
      assert (Solver != 0);

      // 5.- factorize and solve
      
      Time.ResetStartTime();
      AMESOS_CHK_ERR(Solver->SymbolicFactorization());
      if (verbose) 
	cout << endl
	     << "Solver " << SolverType[i] 
	     << ", symbolic factorization time = " 
	     << Time.ElapsedTime() << endl;

      AMESOS_CHK_ERR(Solver->NumericFactorization());
      if (verbose) 
	cout << "Solver " << SolverType[i] 
	     << ", numeric factorization time = " 
	     << Time.ElapsedTime() << endl;

      AMESOS_CHK_ERR(Solver->Solve());
      if (verbose) 
	cout << "Solver " << SolverType[i] 
	     << ", solve time = " 
	     << Time.ElapsedTime() << endl;
  
      // 6.- compute difference between exact solution and Amesos one
      //     (there are other ways of doing this in Epetra, but let's
      //     keep it simple)
      double d = 0.0, d_tot = 0.0;
      for (int j = 0 ; j< lhs->Map().NumMyElements() ; ++j)
	d += ((*lhs)[0][j] - 1.0) * ((*lhs)[0][j] - 1.0);

      Comm.SumAll(&d,&d_tot,1);
      if (verbose)
	cout << "Solver " << SolverType[i] << ", ||x - x_exact||_2 = " 
	     << sqrt(d_tot) << endl;

      // 7.- delete the object
      delete Solver;

      TotalResidual += d_tot;
    }
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (TotalResidual < 1e-9) 
    exit(EXIT_SUCCESS);
  else
    exit(EXIT_FAILURE);

}
