#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

  // Total number of elements in the vector
  int NumElements = 10;

  // Construct a Map with NumElements and index base of 0
  Epetra_Map Map(NumElements, 0, Comm);

  // Create x as an Epetra_vector
  Epetra_Vector x(Map);
 
  // get the local size of the vector
  int MyLength = x.MyLength();

  /* First way to define the vector:   */
  /* use the [] operator on the object */

  for (int i = 0 ; i < MyLength ; ++i) x[i] = 1.0 * i;

  // need a double pointer because this works with multi-vectors
  double * pointer;
  
  x.ExtractView( &pointer );

  for (int i = 0; i < MyLength; ++i)
    cout << "on proc " << Comm.MyPID() << ", x["
	 << i << "] = " << pointer[i] << endl;
  
  // now modify the values
  for (int i = 0; i < MyLength; ++i)
    pointer[i] *= 10;

  // this affects the object x
  cout << x;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
} /* main */
