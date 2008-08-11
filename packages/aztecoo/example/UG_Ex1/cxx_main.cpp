#include "AztecOO_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "AztecOO.h"

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif
  cout << " Testing on " << Comm <<endl;

  int NumMyElements = 100;
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
  A.FillComplete();

  // Create x and b vectors
  Epetra_Vector x(Map);
  Epetra_Vector b(Map);
  b.Random(); // Fill b with random values

  // Create Linear Problem
  Epetra_LinearProblem problem(&A, &x, &b);
  // Create AztecOO instance
  AztecOO solver(problem);
  
  cout << "Starting iterations for: " << endl;
  solver.SetLabel("AztecOO User Guide Example 1");
  cout << solver.GetLabel() << endl;

  solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  solver.Iterate(1000, 1.0E-4);

  cout << "Solver performed " << solver.NumIters() << " iterations." << endl
       << "Norm of true residual = " << solver.TrueResidual() << endl;
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return 0;
}
