#include "Amesos_ConfigDefs.h"

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
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"

using namespace Galeri;

bool quiet = false ;

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

  if (!quiet && !A->Comm().MyPID()) 
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

void driver(Epetra_Comm& Comm, const bool IsSymmetric, const bool UseTranspose, 
            std::vector<std::string>& SolverType)
{
  std::string ProblemType;
  if (IsSymmetric)
    ProblemType = "Laplace2D";
  else
    ProblemType = "Recirc2D";

  Teuchos::ParameterList GaleriList;
  GaleriList.set("nx", 8);
  GaleriList.set("ny", 8 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix = CreateCrsMatrix(ProblemType, Map, GaleriList);
  Epetra_MultiVector LHS(*Map, 3); LHS.PutScalar(0.0);
  Epetra_MultiVector RHS(*Map, 3); RHS.Random();

  Amesos_TestRowMatrix A(Matrix);

  Amesos Factory;  
  
  bool res;

  // If a given test fails, than the code stops, bue to the assert()
  // statement.
  for (unsigned int i = 0 ; i < SolverType.size() ; ++i) 
  {
    std::string Solver = SolverType[i];

    if (Factory.Query((char*)Solver.c_str())) 
    {
      // solve with matrix
      Teuchos::ParameterList AmesosList;
      AmesosList.set("Redistribute",true);
      res = TestAmesos((char*)Solver.c_str(), AmesosList, false, 
                       &A, &LHS, &RHS);
      assert (res == true);
      if (UseTranspose) {
        // solve transpose with matrix
        Teuchos::ParameterList AmesosList;
        res  = TestAmesos((char*)Solver.c_str(), AmesosList, true, 
                          &A, &LHS, &RHS);
        assert (res == true);
      }
    } 
    else
      if (!quiet && !Comm.MyPID()) 
      {
        std::cerr << std::endl;
        std::cerr << "WARNING: SOLVER `" << Solver << "' NOT TESTED" << std::endl;
        std::cerr << std::endl;
      }
  }

  delete Matrix;
  delete Map;
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

  if ( argc > 1 && argv[1][0] == '-' &&  argv[1][1] == 'q' ) quiet = true ;

  // NOTE: DSCPACK does not support RowMatrix's.
  
  if (true)
  {
    // non-symmetric matrix, test A and A^T
    std::vector<std::string> SolverType;
    SolverType.push_back("Amesos_Lapack");
    SolverType.push_back("Amesos_Klu");
    SolverType.push_back("Amesos_Umfpack");
    SolverType.push_back("Amesos_Superlu");
//    SolverType.push_back("Amesos_Mumps");       Bug #1896 
    SolverType.push_back("Amesos_Scalapack");
    driver(Comm, false, true, SolverType);
  }
 
  if (true)
  {
    // non-symmetric matrix, test only A
    std::vector<std::string> SolverType;
    //  
    //    SolverType.push_back("Amesos_Pardiso");    bug #1994
    SolverType.push_back("Amesos_Superludist");
    driver(Comm, false, false, SolverType);
  }

  // I have some problems with Taucs with LAM
  if (true)
  {
    // symmetric
    std::vector<std::string> SolverType;
    SolverType.push_back("Amesos_Taucs");
    driver(Comm, true, false, SolverType);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(0); 
}
