#include <stdio.h>
#include <stdlib.h>
#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Version.h"
                                            
// prototype
double power_method(const Epetra_CrsMatrix& A);
int main(int argc, char *argv[]) {
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  if (Comm.MyPID()==0)
    cout << Epetra_Version() << endl << endl;

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
  Epetra_SerialDenseVector DiagVal(1); 
  DiagVal[0] = 2.0; // We set all diagonal values to 2
  Epetra_IntSerialDenseVector ColIndices(1);
  for (int i=0; i<NumMyElements; i++) {
    int RowIndex = Map.GID(i);
    ColIndices[0] = RowIndex;
    // Put in the diagonal entry
    A.InsertGlobalValues(RowIndex, DiagVal.Length(), 
			 DiagVal.Values(), ColIndices.Values());  
  }
  if (Map.MyGID(0)) { // Change the first global diagonal value to 4.0
    DiagVal[0] = 4.0;
    int RowIndex = 0;
    ColIndices[0] = RowIndex;
  A.ReplaceGlobalValues(RowIndex, DiagVal.Length(), 
			DiagVal.Values(), ColIndices.Values());
  }
  // Finish up
 A.TransformToLocal();
  // Iterate
  double lambda = power_method(A);
  if (Comm.MyPID()==0) 
    cout << endl << "Estimate of Dominant Eigenvalue = " << lambda << endl;		
#ifdef UG_EX1_MPI
  MPI_Finalize() ;
#endif
return 0;
}
