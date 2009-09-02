#ifdef EPETRA_MPI
#include <mpi.h>
#endif
#include "Epetra_C_wrappers.h"

int main(int argc, char *argv[]) {
  int i;
  EPETRA_OBJECT_PTR Comm, Map, X, Y;
  int MyPID, NumProc;
  int NumGlobalElements;
  int NumMyElements;

#ifdef EPETRA_MPI
  /* Initialize MPI */
  MPI_Init(&argc,&argv);
  Comm = epetra_mpicomm_create2( MPI_COMM_WORLD );
#else
  Comm = epetra_serialcomm_create();
#endif

   MyPID = epetra_comm_mypid(Comm);
   NumProc = epetra_comm_numproc(Comm);


   /* Construct a Map that puts 2 elements on each PE */

   NumGlobalElements = 2*NumProc;
   Map = epetra_map_create1(NumGlobalElements, 0, Comm);
  
 
  X = epetra_vector_create1(Map);
  Y = epetra_vector_create1(Map);

  epetra_vector_random(X);
  epetra_vector_random(Y);
  printf("Contents of X vector\n");
  epetra_vector_print(X);


  printf("Contents of Y vector\n");
  epetra_vector_print(Y);

  /* Add X and Y (need to pass Y twice for now, since this is the only update 
     interface wrapped by C at this time) */
  epetra_vector_update(X, 1.0, Y, 0.0, Y, 1.0);

  printf("Sum of X and Y vectors\n");
  epetra_vector_print(X);

  epetra_vector_destroy(X);
  epetra_vector_destroy(Y);
  epetra_map_destroy(Map);
  epetra_comm_destroy(Comm);
  
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return 0 ;
}
