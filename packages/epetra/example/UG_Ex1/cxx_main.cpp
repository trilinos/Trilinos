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

// Simple Power method algorithm
double power_method(const Epetra_CrsMatrix& A) {  
  // variable needed for iteration
  double lambda = 0.0;
  int niters = A.RowMap().NumGlobalElements()*10;
  double tolerance = 1.0e-10;
  // Create vectors
  Epetra_Vector q(A.RowMap());
  Epetra_Vector z(A.RowMap());
  Epetra_Vector resid(A.RowMap());
  // Fill z with random Numbers
  z.Random();
  // variable needed for iteration
  double normz;
  double residual = 1.0 + tolerance;
  int iter = 0;
  while (iter < niters && residual > tolerance) {
    z.Norm2(&normz); // Compute 2-norm of z
    q.Scale(1.0/normz, z);
    A.Multiply(false, q, z); // Compute z = A*q
    q.Dot(z, &lambda); // Approximate maximum eigenvalue
    if (iter%10==0 || iter+1==niters) {
      resid.Update(1.0, z, -lambda, q, 0.0); // Compute A*q - lambda*q
      resid.Norm2(&residual);
      cout << "Iter = " << iter << "  Lambda = " << lambda 
	   << "  Two-norm of A*q - lambda*q = " 
	   << residual << endl;
    } 
    iter++;
  }
  return(lambda);
}
