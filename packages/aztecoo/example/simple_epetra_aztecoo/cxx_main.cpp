#include <cstdio>
#include <cstdlib>
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "AztecOO.h"

int main(int argc, char *argv[])
{

  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
  cout << Comm <<endl;

  int NumMyElements = 1000;

  // Construct a Map that puts same number of equations on each processor
  Epetra_Map Map(-1, NumMyElements, 0, Comm);
  int NumGlobalElements = Map.NumGlobalElements();

  // Create a Epetra_Matrix
  Epetra_CrsMatrix A(Copy, Map, 3);
  
  // Add  rows one-at-a-time
  double negOne = -1.0;
  double posTwo = 2.0;
  for (int i=0; i<NumMyElements; i++) {
    int GlobalRow = A.GRID(i); int RowLess1 = GlobalRow - 1; int RowPlus1 = GlobalRow + 1;

    if (RowLess1!=-1) A.InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess1);
    if (RowPlus1!=NumGlobalElements) A.InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus1);
    A.InsertGlobalValues(GlobalRow, 1, &posTwo, &GlobalRow);
  }
  
  // Finish up
  A.TransformToLocal();

  // Create x and b vectors
  Epetra_Vector x(Map);
  Epetra_Vector b(Map);

  b.Random();

  // Create Linear Problem
  Epetra_LinearProblem problem(&A, &x, &b);
  // Create AztecOO instance
  AztecOO solver(problem);

  solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  solver.Iterate(1000, 1.0E-8);

  cout << "Solver performed " << solver.NumIters() << " iterations." << endl
       << "Norm of true residual = " << solver.TrueResidual() << endl;

  MPI_Finalize() ;
  return 0;
}

