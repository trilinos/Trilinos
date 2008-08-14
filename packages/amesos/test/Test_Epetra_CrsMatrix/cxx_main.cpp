#include "Amesos_ConfigDefs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
// I have no idea why we have to include Epetra_Map.h here (and not 
// elsewhere, such as Test_EpteraRowMatrix/cxx_main.cpp) 
// but this was added as a workaround to bug #1869   
#include "Epetra_Map.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Amesos.h"
#include "Teuchos_ParameterList.hpp"

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

bool TestAmesos(char ProblemType[], Teuchos::ParameterList& AmesosList,
                bool UseTranspose, Epetra_RowMatrix* A, Epetra_MultiVector* lhs,
                Epetra_MultiVector* rhs)
{
  Epetra_LinearProblem Problem;
  Amesos A_Factory;

  Amesos_BaseSolver * Solver = A_Factory.Create(ProblemType, Problem);
  assert (Solver != 0);

  // Both sentences should work
  Solver->SetUseTranspose(UseTranspose);

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
  std::vector<double> Norm(rhs->NumVectors());

  Epetra_MultiVector Ax(*rhs);
  A->Multiply(UseTranspose, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm[0]);

  std::string msg = ProblemType;

  if (!A->Comm().MyPID()) 
  {
    std::cout << std::endl;
    std::cout << msg << " : Using " << A->Comm().NumProc() << " processes, UseTranspose = " << UseTranspose << std::endl;
    for (int j = 0 ; j < rhs->NumVectors() ; ++j)
      std::cout << msg << " : eq " << j 
	   << ", ||A x - b||_2 = " << Norm[j] << std::endl;
    std::cout << msg << " : ||x_exact - x||_2 = " << sqrt(d_tot) << std::endl;
  }

  if (Norm[0] > 1e-9)  
  {
    std::cerr << std::endl << msg << " WARNING : TEST FAILED!" << std::endl;
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

  int NumGlobalRows = 100 * Comm.NumProc();

  Epetra_Map Map(NumGlobalRows, 0, Comm);

  Epetra_CrsMatrix Matrix(Copy, Map, 0);

  int NumMyRows = Map.NumMyElements();
  int* MyGlobalElements = Map.MyGlobalElements();

  int Indices[3];
  double Values[3];
  int NumEntries;

  for (int LRID = 0 ; LRID < NumMyRows ; ++LRID)
  {
    int GRID = MyGlobalElements[LRID];

    Indices[0] = GRID;
    Values[0] = 2.0;
    NumEntries = 1;

    if (GRID != 0)
    {
      Indices[NumEntries] = GRID - 1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if (GRID != NumGlobalRows - 1)
    {
      Indices[NumEntries] = GRID + 1;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }

    Matrix.InsertGlobalValues(GRID, NumEntries, Values, Indices);
  }

  Matrix.FillComplete();

  Epetra_MultiVector LHS(Map, 3);
  Epetra_MultiVector RHS(Map, 3);

  Amesos Factory;  
  
  std::vector<std::string> SolverType;
  SolverType.push_back("Amesos_Lapack");
  SolverType.push_back("Amesos_Klu");
  SolverType.push_back("Amesos_Umfpack");
  //  SolverType.push_back("Amesos_Pardiso");   bug #1994
  SolverType.push_back("Amesos_Taucs");
  SolverType.push_back("Amesos_Superlu");
  SolverType.push_back("Amesos_Superludist");
  SolverType.push_back("Amesos_Mumps");
  SolverType.push_back("Amesos_Dscpack");
  SolverType.push_back("Amesos_Scalapack");

  bool res;

  // If a given test fails, than the code stops, due to the assert()
  // statement.
  for (unsigned int i = 0 ; i < SolverType.size() ; ++i) 
  {
    std::string Solver = SolverType[i];

    if (Factory.Query((char*)Solver.c_str())) 
    {
      if (1) {
	// solve with matrix
	Teuchos::ParameterList AmesosList;
	res = TestAmesos((char*)Solver.c_str(), AmesosList, false, 
                         &Matrix, &LHS, &RHS);
        assert (res == true);
      }
      if (1) {
	// solve transpose with matrix
	if (Solver != "Amesos_Superludist") {// still not implementes
	  Teuchos::ParameterList AmesosList;
	  res  = TestAmesos((char*)Solver.c_str(), AmesosList, true, 
                            &Matrix, &LHS, &RHS);
          assert (res == true);
	}
      }
    } 
    else
      if (!Comm.MyPID()) 
      {
	std::cerr << std::endl;
	std::cerr << "WARNING: SOLVER `" << Solver << "' NOT TESTED" << std::endl;
	std::cerr << std::endl;
      }
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(0); 
}
