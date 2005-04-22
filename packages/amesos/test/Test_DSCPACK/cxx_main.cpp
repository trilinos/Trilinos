//
//  This test does not exercise multiple processes on DSCPACK.
//  DSCPACK uses only one process for dense matrices.
//
//
#include "Amesos_ConfigDefs.h"
#ifdef HAVE_AMESOS_DSCPACK
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Amesos_Dscpack.h"
#include "Amesos_TestRowMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include <vector>

//=============================================================================
bool CheckError(const Epetra_RowMatrix& A,
		const Epetra_MultiVector& x,
		const Epetra_MultiVector& b,
		const Epetra_MultiVector& x_exact)
{
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
    cout << "||Ax - b||  = " << TotalNorm << endl;
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

  int NumGlobalElements = 1000;
  int NumVectors = 7;
  
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

  // =========================== //
  // Create a tridiagonal matrix //
  // =========================== //
 
  double Values[3];
  // Right now, we put zeros only in the matrix.
  Values[0] = 0.0;
  Values[1] = 0.0;
  Values[2] = 0.0;
  int Indices[3];
  int NumEntries;

  // At this point we simply set the nonzero structure of A.
  // Actual values will be inserted later (now all zeros)
  for (int i = 0; i < NumMyElements; i++)
  {
    if (MyGlobalElements[i] == 0)
    {
      Indices[0] = 0;
      Indices[1] = 1;
      NumEntries = 2;
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1)
    {
      Indices[0] = NumGlobalElements-1;
      Indices[1] = NumGlobalElements-2;
      NumEntries = 2;
    }
    else
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i];
      Indices[2] = MyGlobalElements[i]+1;
      NumEntries = 3;
    }

    AMESOS_CHK_ERR(Matrix.InsertGlobalValues(MyGlobalElements[i],
                                             NumEntries, Values, Indices));
  }

  // Finish up.
  Matrix.FillComplete();

  // ======================== //
  // other data for this test //
  // ======================== //

//  Amesos_TestRowMatrix A(&Matrix);
  Epetra_MultiVector x(Map,NumVectors);
  Epetra_MultiVector x_exact(Map,NumVectors);
  Epetra_MultiVector b(Map,NumVectors);
  x_exact.Random();
  Matrix.Multiply(false,x_exact,b);

  // =========== //
  // AMESOS PART //
  // =========== //

  Epetra_LinearProblem Problem(&Matrix, &x, &b);
  Amesos_Dscpack* Solver = new Amesos_Dscpack(Problem);

  Teuchos::ParameterList ParamList;
  ParamList.set("MaxProcs", Comm.NumProc());
  ParamList.set("PrintStatus", true);

  EPETRA_CHK_ERR(Solver->SetParameters(ParamList)); 

  AMESOS_CHK_ERR(Solver->SymbolicFactorization());
  AMESOS_CHK_ERR(Solver->NumericFactorization());
  AMESOS_CHK_ERR(Solver->Solve());

  if (!CheckError(Matrix,x,b,x_exact))
    AMESOS_CHK_ERR(-1);

  delete Solver;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#else

// DSCPACK is not available. Sorry, we have to give up.

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

  puts("Please configure Amesos with:");
  puts("--enable-amesos-dscpack");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
}
#endif
