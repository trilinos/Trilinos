#include <stdio.h>
#include <stdlib.h>
#ifdef UG_EX1_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
// prototype
double power_method(const Epetra_CrsMatrix& A);
int main(int argc, char *argv[]) {
#ifdef UG_EX1_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  cout << Comm << endl; // Print out process information
  // Get the number of global equations from the command line
  if (argc!=2) { 
    cout << "Usage: " << argv[0] << " number_of_equations" << endl;
    exit(1);
   }
  int NumGlobalElements = atoi(argv[1]);
  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Epetra_Map Map(NumGlobalElements, 0, Comm);
  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map.NumMyElements();
  // Create an Epetra_CrsMatrix
  Epetra_CrsMatrix A(Copy, Map, 1);
 // Add  rows one-at-a-time
  double two = 2.0;
  for (int i=0; i<NumMyElements; i++) {
    int index = Map.GID(i);
    // Put in the diagonal entry
    A.InsertGlobalValues(index, 1, &two, &index);  
  }
  double four = 4.0;
  int index = 0;
  A.ReplaceGlobalValues(index, 1, &four, &index);  
  // Finish up
 A.TransformToLocal();
  // Iterate
  double lambda = power_method(A);
  cout << "Estimate of Dominant Eigenvalue = " << lambda << endl;		
#ifdef UG_EX1_MPI
  MPI_Finalize() ;
#endif
return 0;
}
