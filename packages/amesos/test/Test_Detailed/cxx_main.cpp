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
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_TestRowMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include <vector>
using namespace Trilinos_Util;

//=============================================================================
bool CheckError(const string SolverType,
		const string Descriptor,
		const Epetra_RowMatrix& A,
		const Epetra_MultiVector& x,
		const Epetra_MultiVector& b,
		const Epetra_MultiVector& x_exact)
{
  if (A.Comm().MyPID() == 0)
    cout << endl << "*** " << SolverType << " : " 
         << Descriptor << " ***" << endl;
  vector<double> Norm;
  int NumVectors = x.NumVectors();
  Norm.resize(NumVectors);
  Epetra_MultiVector Ax(x);
  A.Multiply(false,x,Ax);
  Ax.Update(1.0,b,-1.0);
  Ax.Norm2(&Norm[0]);
  bool TestPassed = false;
  double TotalNorm = 0.0;
  for (int i = 0 ; i < NumVectors ; ++i) {
    TotalNorm += Norm[i];
  }
  if (A.Comm().MyPID() == 0)
    cout << "||Ax - b||  = " << TotalNorm << endl;
  if (TotalNorm < 1e-5 )
    TestPassed = true;
  else
    TestPassed = false;

  Ax.Update (1.0,x,-1.0,x_exact,0.0);
  Ax.Norm2(&Norm[0]);
  for (int i = 0 ; i < NumVectors ; ++i) {
    TotalNorm += Norm[i];
  }
  if (A.Comm().MyPID() == 0)
    cout << "||x - x_exact||  = " << TotalNorm << endl;
  if (TotalNorm < 1e-5 )
    TestPassed = true;
  else
    TestPassed = false;

  return(TestPassed);
}

//=============================================================================
bool Test(char* SolverType,
	  Epetra_RowMatrix& A, Epetra_MultiVector& x_A,
	  Epetra_MultiVector& b_A, Epetra_MultiVector& x_exactA,
	  Epetra_RowMatrix& B, Epetra_MultiVector& x_B,
	  Epetra_MultiVector& b_B, Epetra_MultiVector& x_exactB,
	  Epetra_RowMatrix& C, Epetra_MultiVector& x_C,
	  Epetra_MultiVector& b_C, Epetra_MultiVector& x_exactC)
{
  bool TestPassed = true;
  Epetra_LinearProblem ProblemA;
  Epetra_LinearProblem ProblemB;
  Epetra_LinearProblem ProblemC;
  Amesos Factory;
  Amesos_BaseSolver* Solver;

  // Test simple usage:
  // - set problem (empty)
  // - set matrix and vector
  // - call Solve() directly
  {
    x_A.PutScalar(0.0);
    ProblemA.SetOperator((Epetra_RowMatrix*)0);
    ProblemA.SetLHS((Epetra_MultiVector*)0);
    ProblemA.SetRHS((Epetra_MultiVector*)0);

    Solver = Factory.Create(SolverType,ProblemA);
  
    ProblemA.SetOperator(&A);
    ProblemA.SetLHS(&x_A);
    ProblemA.SetRHS(&b_A);

    AMESOS_CHK_ERR(Solver->Solve());

    TestPassed = TestPassed && 
      CheckError(SolverType, "Solve() only", A,x_A,b_A,x_exactA);
  }

  // Test almost simple usage:
  // - set problem (empty)
  // - set matrix and vector
  // - call NumericFactorization()
  // - call Solve() directly
  {
    x_A.PutScalar(0.0);
    ProblemA.SetOperator((Epetra_RowMatrix*)0);
    ProblemA.SetLHS((Epetra_MultiVector*)0);
    ProblemA.SetRHS((Epetra_MultiVector*)0);

    Solver = Factory.Create(SolverType,ProblemA);
  
    ProblemA.SetOperator(&A);
    AMESOS_CHK_ERR(Solver->NumericFactorization());
    ProblemA.SetLHS(&x_A);
    ProblemA.SetRHS(&b_A);
    AMESOS_CHK_ERR(Solver->Solve());

    TestPassed = TestPassed && 
      CheckError(SolverType, "NumFact() + Solve()", A,x_A,b_A,x_exactA);
  }
  // Test normal usage:
  // - set problem (empty)
  // - set matrix
  // - call SymbolicFactorization() several times
  // - call NumericFactorization() several times
  // - set vectors
  // - call Solve() several times
  // Repeated calls should *not* break the object. Data are supposed
  // to be re-computed as needed.
  {
    x_A.PutScalar(0.0);
    ProblemA.SetOperator((Epetra_RowMatrix*)0);
    ProblemA.SetLHS((Epetra_MultiVector*)0);
    ProblemA.SetRHS((Epetra_MultiVector*)0);

    Solver = Factory.Create(SolverType,ProblemA);
  
    ProblemA.SetOperator(&A);
    AMESOS_CHK_ERR(Solver->SymbolicFactorization());
    AMESOS_CHK_ERR(Solver->SymbolicFactorization());
    AMESOS_CHK_ERR(Solver->SymbolicFactorization());
    AMESOS_CHK_ERR(Solver->NumericFactorization());
    AMESOS_CHK_ERR(Solver->NumericFactorization());
    AMESOS_CHK_ERR(Solver->NumericFactorization());
    AMESOS_CHK_ERR(Solver->NumericFactorization());
    ProblemA.SetLHS(&x_A);
    ProblemA.SetRHS(&b_A);

    AMESOS_CHK_ERR(Solver->Solve());
    AMESOS_CHK_ERR(Solver->Solve());
    AMESOS_CHK_ERR(Solver->Solve());
    AMESOS_CHK_ERR(Solver->Solve());
    AMESOS_CHK_ERR(Solver->Solve());
    AMESOS_CHK_ERR(Solver->Solve());

    TestPassed = TestPassed && 
      CheckError(SolverType, "SymFact() + NumFact() + Solve()", 
		 A,x_A,b_A,x_exactA);
  }

  // Test normal usage:
  // - set problem (empty)
  // - set matrix (as A)
  // - call SymbolicFactorization()
  // - set pointer to B
  // - call NumericFactorization()
  // - set vectors for B
  // - call Solve()
  {
    x_A.PutScalar(0.0);
    ProblemA.SetOperator((Epetra_RowMatrix*)0);
    ProblemA.SetLHS((Epetra_MultiVector*)0);
    ProblemA.SetRHS((Epetra_MultiVector*)0);

    Solver = Factory.Create(SolverType,ProblemA);
 
    ProblemA.SetOperator(&A);
    AMESOS_CHK_ERR(Solver->SymbolicFactorization());
    ProblemA.SetOperator(&B);
    AMESOS_CHK_ERR(Solver->NumericFactorization());
    ProblemA.SetLHS(&x_B);
    ProblemA.SetRHS(&b_B);

    AMESOS_CHK_ERR(Solver->Solve());

    TestPassed = TestPassed && 
      CheckError(SolverType, "Set A, solve B", B,x_B,b_B,x_exactB);
  }

  // Construct Solver with filled ProblemA.
  // Then, stick pointers for different vectors.
  // a different map as well.
  {
    
    x_C.PutScalar(0.0);
    ProblemA.SetOperator(&A);
    ProblemA.SetLHS(&x_A);
    ProblemA.SetRHS(&b_A);

    Solver = Factory.Create(SolverType,ProblemA);

    ProblemA.SetOperator(&C);
    AMESOS_CHK_ERR(Solver->SymbolicFactorization());
    AMESOS_CHK_ERR(Solver->NumericFactorization());

    ProblemA.SetLHS(&x_C);
    ProblemA.SetRHS(&b_C);
    AMESOS_CHK_ERR(Solver->Solve());

    TestPassed = TestPassed && 
      CheckError(SolverType, "Set A, Solve C", C,x_C,b_C,x_exactC);
  }

  // Construct Solver with filled ProblemA, call Solve().
  // Then, replace the pointers for matrix and vectors,
  // and call Solve() again.
  {
    x_C.PutScalar(0.0);
    ProblemA.SetOperator(&A);
    ProblemA.SetLHS(&x_A);
    ProblemA.SetRHS(&b_A);

    Solver = Factory.Create(SolverType,ProblemA);

    AMESOS_CHK_ERR(Solver->Solve());

    ProblemA.SetOperator(&C);
    AMESOS_CHK_ERR(Solver->SymbolicFactorization());
    AMESOS_CHK_ERR(Solver->NumericFactorization());

    ProblemA.SetLHS(&x_C);
    ProblemA.SetRHS(&b_C);
    AMESOS_CHK_ERR(Solver->Solve());

    TestPassed = TestPassed && 
      CheckError(SolverType, "Solve A + Solve C", C,x_C,b_C,x_exactC);
  }

  return(TestPassed);
}

//=============================================================================
int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Creation of data
  // A and B refer to two different matrices, with the *same*
  // structure and *different* values. Clearly, the map is the same.
  //
  // C refers to a completely different matrix
 
  int Size_AB = 900; // must be square
  int Size_C = 30; // must be square
  int NumVectors_AB = 7;
  int NumVectors_C = 13;
  
  CrsMatrixGallery GalleryA("recirc_2d", Comm);
  GalleryA.Set("problem_size", Size_AB);
  GalleryA.Set("map_type", "interlaced");
  GalleryA.Set("num_vectors", NumVectors_AB);

  CrsMatrixGallery GalleryB("laplace_2d", Comm);
  GalleryB.Set("problem_size", Size_AB);
  GalleryB.Set("map_type", "interlaced");
  GalleryB.Set("num_vectors", NumVectors_AB);

  CrsMatrixGallery GalleryC("minij", Comm);
  GalleryC.Set("problem_size", Size_C);
  GalleryC.Set("num_vectors", NumVectors_C);

  Epetra_RowMatrix* RowA = GalleryA.GetLinearProblem()->GetMatrix();
  Epetra_RowMatrix* RowB = GalleryB.GetLinearProblem()->GetMatrix();
  Epetra_RowMatrix* RowC = GalleryC.GetLinearProblem()->GetMatrix();

  Amesos_TestRowMatrix A(RowA);
  Amesos_TestRowMatrix B(RowB);
  Amesos_TestRowMatrix C(RowC);

  Epetra_MultiVector x_A(A.OperatorDomainMap(),NumVectors_AB);
  Epetra_MultiVector x_exactA(A.OperatorDomainMap(),NumVectors_AB);
  Epetra_MultiVector b_A(A.OperatorRangeMap(),NumVectors_AB);
  x_exactA.Random();
  A.Multiply(false,x_exactA,b_A);

  Epetra_MultiVector x_B(B.OperatorDomainMap(),NumVectors_AB);
  Epetra_MultiVector x_exactB(B.OperatorDomainMap(),NumVectors_AB);
  Epetra_MultiVector b_B(B.OperatorRangeMap(),NumVectors_AB);
  x_exactB.Random();
  B.Multiply(false,x_exactB,b_B);

  Epetra_MultiVector x_C(C.OperatorDomainMap(),NumVectors_C);
  Epetra_MultiVector x_exactC(C.OperatorDomainMap(),NumVectors_C);
  Epetra_MultiVector b_C(C.OperatorRangeMap(),NumVectors_C);
  x_exactC.Random();
  C.Multiply(false,x_exactC,b_C);

  vector<string> SolverType;
  SolverType.push_back("Amesos_Lapack");
  SolverType.push_back("Amesos_Klu");
  SolverType.push_back("Amesos_Umfpack");
  SolverType.push_back("Amesos_Superlu");
  SolverType.push_back("Amesos_Superludist");
  SolverType.push_back("Amesos_Mumps");
//  SolverType.push_back("Amesos_Dscpack");

  Amesos Factory;
  bool TestPassed = true;

  for (int i = 0 ; i < SolverType.size() ; ++i) {
    string Solver = SolverType[i];
    if (Factory.Query((char*)Solver.c_str())) {
      bool ok = Test((char*)Solver.c_str(),
		     A, x_A, b_A, x_exactA,
		     B, x_B, b_B, x_exactB,
		     C, x_C, b_C, x_exactC);
      TestPassed = TestPassed && ok;
    }
    else
      cout << "Solver " << Solver << " not available" << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (TestPassed) {
    if (Comm.MyPID() == 0)
      cout << endl << "TEST PASSED" << endl << endl;
    exit(EXIT_SUCCESS);
  }
  else {
    if (Comm.MyPID() == 0)
      cout << endl << "TEST FAILED" << endl << endl;
    exit(EXIT_FAILURE);
  }

}

#else

// Triutils is not available. Sorry, we have to give up.

#include <stdlib.h>
#ifdef HAVE_MPI
#include "mpi.h"
#else
#endif
#include <stdio.h>

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


