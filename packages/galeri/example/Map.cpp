#include "Galeri_Maps.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"

using namespace Galeri;

int main(int argv, char* argc[])
{
#ifdef HAVE_MPI
  MPI_Init(&argv, &argc);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = (Comm.MyPID() == 0);

  // Creates an Epetra_Map corresponding to a 2D Cartesian grid
  // on the unit square. For parallel runs, the nodes are divided into
  // strips, so that the total number of subdomains is Comm.NumProc() x 1.
  
  // Pointer to the object to be created
  Epetra_Map* Map = 0; 
  // Type of the object
  string MapType = "Cartesian2D";
  // Container for parameters
  Teuchos::ParameterList GaleriList;
  GaleriList.set("nx", 2 * Comm.NumProc()); 
  GaleriList.set("ny", 2);
  GaleriList.set("mx", Comm.NumProc());
  GaleriList.set("my", 1);

  try
  {
    // Creation of the map
    Map = CreateMap("Cartesian2D", Comm, GaleriList);

    // print out the map
    cout << *Map;

    // To created object must be free'd using delete
    delete Map;
  }
  catch (Exception& rhs)
  {
    if (Comm.MyPID() == 0)
      rhs.Print();
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
