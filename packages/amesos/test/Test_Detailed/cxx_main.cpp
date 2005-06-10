#include "Amesos_ConfigDefs.h"
#ifdef HAVE_AMESOS_TRIUTILS

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_TestRowMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include <vector>

#ifdef HAVE_VALGRIND_H
#include <valgrind/valgrind.h>
#define HAVE_VALGRIND
#else
#ifdef HAVE_VALGRIND_VALGRIND_H
#include <valgrind/valgrind.h>
#define HAVE_VALGRIND
#endif 
#endif 

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
  if (true) {
    x_A.PutScalar(0.0);
    ProblemA.SetOperator((Epetra_RowMatrix*)0);
    ProblemA.SetLHS((Epetra_MultiVector*)0);
    ProblemA.SetRHS((Epetra_MultiVector*)0);

    Solver = Factory.Create(SolverType,ProblemA);
  
    ProblemA.SetOperator(&A);
    ProblemA.SetLHS(&x_A);
    ProblemA.SetRHS(&b_A);

    string ST = SolverType ;    
    if (! ( ST == "Amesos_Superludist" ) ) {    // Kludge see bug #1141 and bug #1138
      AMESOS_CHK_ERR(Solver->Solve());

      TestPassed = TestPassed && 
	CheckError(SolverType, "Solve() only", A,x_A,b_A,x_exactA);
    }
    delete Solver; 
  }

  // Test almost simple usage:
  // - set problem (empty)
  // - set matrix and vector
  // - call NumericFactorization()
  // - call Solve() directly
  
  if (true) 
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

    delete Solver; 
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
  
  if (true) 
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

    delete Solver; 
  }

  // Test normal usage:
  // - set problem (empty)
  // - set matrix (as A)
  // - call SymbolicFactorization()
  // - set pointer to B
  // - call NumericFactorization()
  // - set vectors for B
  // - call Solve()
  
  if (true)
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

    delete Solver; 
  }

  // Construct Solver with filled ProblemA.
  // Then, stick pointers for different vectors.
  // a different map as well.
  
  if (true)
  {
    x_C.PutScalar(0.0);
    ProblemA.SetOperator(&A);
    ProblemA.SetLHS(&x_A);
    ProblemA.SetRHS(&b_A);

    Solver = Factory.Create(SolverType,ProblemA);

    ProblemA.SetOperator(&C);

    string ST = SolverType ; 
    if (! ( ST == "Amesos_Superludist" ) ) { // Kludge see bug #1141 and bug #1138

      AMESOS_CHK_ERR(Solver->SymbolicFactorization());
      AMESOS_CHK_ERR(Solver->NumericFactorization());
      
      ProblemA.SetLHS(&x_C);
      ProblemA.SetRHS(&b_C);
      AMESOS_CHK_ERR(Solver->Solve());
      
      TestPassed = TestPassed && 
	CheckError(SolverType, "Set A, Solve C", C,x_C,b_C,x_exactC);
    }
    delete Solver; 
  }

  // Construct Solver with filled ProblemA, call Solve().
  // Then, replace the pointers for matrix and vectors,
  // and call Solve() again.
  
  if (true)
  {
    x_C.PutScalar(0.0);
    ProblemA.SetOperator(&A);
    ProblemA.SetLHS(&x_A);
    ProblemA.SetRHS(&b_A);

    Solver = Factory.Create(SolverType,ProblemA);

    string ST = SolverType ; 
    if (! ( ST == "Amesos_Superludist" ) ) {   // bug #1141 and bug #1138
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
    delete Solver; 
  }

  return(TestPassed);
}

//=============================================================================

int SubMain( Epetra_Comm &Comm ) {
  // Creation of data
  // A and B refer to two different matrices, with the *same*
  // structure and *different* values. Clearly, the map is the same.
  //
  // C refers to a completely different matrix
 
  int Size_AB = 900; // must be square
  int Size_C = 30; // must be square
  int NumVectors_AB = 7;
  int NumVectors_C = 13;
  
#ifdef HAVE_VALGRIND 
  if ( RUNNING_ON_VALGRIND ) {
   Size_AB = 36; // must be square
   Size_C = 6; 
   NumVectors_AB = 2;
   NumVectors_C = 3;
  }
#endif

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
  //  cout << " A= " ; 
  //  dynamic_cast<Epetra_CrsMatrix *>(RowA)->Print( cout ) ; 
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
  SolverType.push_back("Amesos_Klu");
  SolverType.push_back("Amesos_Lapack");
  SolverType.push_back("Amesos_Umfpack");
  SolverType.push_back("Amesos_Superlu");
  SolverType.push_back("Amesos_Superludist");
  SolverType.push_back("Amesos_Mumps");
  // NOTE: DSCPACK does not support Epetra_RowMatrix's

  Amesos Factory;
  bool TestPassed = true;

  for (unsigned int i = 0 ; i < SolverType.size() ; ++i) {
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

  if (TestPassed) {
    if (Comm.MyPID() == 0)
      cout << endl << "TEST PASSED" << endl << endl;
    return(EXIT_SUCCESS);
  }
  else {
    if (Comm.MyPID() == 0)
      cout << endl << "TEST FAILED" << endl << endl;
    return(EXIT_FAILURE);
  }

}

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

#if 0
  if ( Comm.MyPID() == 0 ) {
    cout << "Enter a char to continue" ;
    char any;
    cin >> any ; 
  }
  Comm.Barrier();
#endif

  int retval = SubMain( Comm ) ;   // 
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return retval ; 
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


