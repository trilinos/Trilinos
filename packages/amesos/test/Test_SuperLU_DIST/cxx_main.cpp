#include "Amesos_ConfigDefs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Amesos_Superludist.h"
#include "Amesos_TestRowMatrix.h"
#include "Teuchos_ParameterList.hpp"

//=============================================================================
bool CheckError(bool verbose,
		const Epetra_RowMatrix& A,
		const Epetra_MultiVector& x,
		const Epetra_MultiVector& b,
		const Epetra_MultiVector& x_exact)
{
  std::vector<double> Norm;
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
  if (verbose && A.Comm().MyPID() == 0)
    std::cout << "||Ax - b||  = " << TotalNorm << std::endl;
  if (TotalNorm < 1e-5 )
    TestPassed = true;
  else
    TestPassed = false;

  Ax.Update (1.0,x,-1.0,x_exact,0.0);
  Ax.Norm2(&Norm[0]);
  for (int i = 0 ; i < NumVectors ; ++i) {
    TotalNorm += Norm[i];
  }
  if (verbose && A.Comm().MyPID() == 0)
    std::cout << "||Ax - b||  = " << TotalNorm << std::endl;
  if (TotalNorm < 1e-5 )
    TestPassed = true;
  else
    TestPassed = false;

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

  int NumGlobalElements = 1000;   // kludge  see bug #1142 - when set to 7, and Min_jGlobal is set to zero, 
                               // Test_SuperLU_DIST.exe dies during Amesos_Superludist destruction
  int NumVectors = 7;
  
  bool verbose = true ;  
  if ( argc > 1 && argv[1][0] == '-' &&  argv[1][1] == 'q' ) verbose = false ;

  // =================== //
  // create a random map //
  // =================== //
 
  int* part = new int[NumGlobalElements];

  if (Comm.MyPID() == 0) {
    Epetra_Util Util;

    for( int i=0 ; i<NumGlobalElements ; ++i ) {
      unsigned int r = Util.RandomInt();	
      part[i] = r%(Comm.NumProc());
    }
  }

  Comm.Broadcast(part,NumGlobalElements,0);

  // count the elements assigned to this proc
  int NumMyElements = 0;
  for (int i = 0 ; i < NumGlobalElements ; ++i) {
    if (part[i] == Comm.MyPID()) 
      NumMyElements++;
  }

  // get the loc2global list
  int* MyGlobalElements = new int[NumMyElements];
  int count = 0;
  for (int i = 0 ; i < NumGlobalElements ; ++i) {
    if (part[i] == Comm.MyPID() ) 
      MyGlobalElements[count++] = i;
  }

  Epetra_Map Map(NumGlobalElements,NumMyElements,MyGlobalElements,
		 0,Comm);

  delete [] part;

  // ===================== //
  // Create a dense matrix //
  // ===================== //
 
  Epetra_CrsMatrix Matrix(Copy,Map,NumGlobalElements);

  int* Indices = new int[NumGlobalElements];
  double* Values = new double[NumGlobalElements];

  for (int i = 0 ; i < NumGlobalElements ; ++i) 
    Indices[i] = i;

  for (int i = 0 ; i < NumMyElements ; ++i) 
  {
    int iGlobal = MyGlobalElements[i];
    const int MakeNotDense = 1;  // kludge  see bug #1142 - set to zero to demonstrate bug #1142 on atlantis
    int Min_jGlobal = std::min(i,MakeNotDense );
    for (int jGlobal = Min_jGlobal ; jGlobal < NumGlobalElements ; ++jGlobal) {
      if (iGlobal == jGlobal) 
	Values[jGlobal-Min_jGlobal] = 1.0 * (NumGlobalElements + 1 ) *
	  (NumGlobalElements + 1);
      else if (iGlobal > jGlobal)
	Values[jGlobal-Min_jGlobal] = -1.0*(jGlobal+1);
      else
	Values[jGlobal-Min_jGlobal] = 1.0*(iGlobal+1);
    }
    Matrix.InsertGlobalValues(MyGlobalElements[i],
                              NumGlobalElements-MakeNotDense, Values, &Indices[Min_jGlobal]);

  }

  Matrix.FillComplete();

  delete[] MyGlobalElements;
  delete[] Indices;
  delete[] Values;
 
  // ======================== //
  // other data for this test //
  // ======================== //

  Amesos_TestRowMatrix A(&Matrix);
  Epetra_MultiVector x(Map,NumVectors);
  Epetra_MultiVector x_exact(Map,NumVectors);
  Epetra_MultiVector b(Map,NumVectors);
  x_exact.Random();
  A.Multiply(false,x_exact,b);

  // =========== //
  // AMESOS PART //
  // =========== //

  Epetra_LinearProblem Problem(&Matrix, &x, &b);
  Amesos_Superludist* Solver = new Amesos_Superludist(Problem);

  Teuchos::ParameterList List;
  List.set("MaxProcs", Comm.NumProc());

  Solver->SetParameters(List);

  AMESOS_CHK_ERR(Solver->SymbolicFactorization());
  AMESOS_CHK_ERR(Solver->NumericFactorization());
  AMESOS_CHK_ERR(Solver->Solve());

  if (!CheckError(verbose,A,x,b,x_exact))
    AMESOS_CHK_ERR(-1);

  delete Solver;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
}
