#include "Galeri_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <vector>
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Galeri_HDF5.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"

using namespace Galeri;

// Showing the usage of HDF5 I/O.
// This example can be run with any number of processors.
//
// \author Marzio Sala, D-INFK/ETHZ.
//
// \date Last modified on 09-Mar-06.

int main (int argc, char **argv)
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // create a map, two vectors, a matrix using Galeri
  Teuchos::ParameterList GaleriList;
  int nx = 10 * Comm.NumProc();
  int ny = 10 * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix = CreateCrsMatrix("Biharmonic2D", Map, GaleriList);
  Epetra_Vector x(Matrix->RowMatrixRowMap());
  Epetra_Vector b(Matrix->RowMatrixRowMap());
  x.Random();
  b.Random();

  // Creates the HDF5 file manager
  Galeri::HDF5 HDF5(Comm);

  HDF5.Create("myfile.h5");

  // =========================== //
  // P A R T   I:  W R I T I N G //
  // =========================== //
  
  // We first write the map, whose name contains the number of processors
  HDF5.Write("map-" + toString(Comm.NumProc()), *Map);
  // Then we write the matrix...
  HDF5.Write("matrix", *Matrix);
  // ... and the x, b vectors
  HDF5.Write("x", x);
  HDF5.Write("b", b);
  // Now we write the parameter list we have used to create the matrix
  HDF5.Write("GaleriList", GaleriList);
  // We can also write integers/doubles/arrays. All these quantities are
  // supposed to have the same value on all processors.
  int IntFlag = 12;
  HDF5.Write("IntFlag", IntFlag);

  double DoubleFlag = 13.0;
  HDF5.Write("IntFlag", IntFlag);

  vector<int> IntArray(3); 
  IntArray[0] = 0, IntArray[1] = 1; IntArray[2] = 2;
  HDF5.Write("IntArray", H5T_NATIVE_INT, 3, &IntArray[0]);

  vector<double> DoubleArray(3);
  DoubleArray[0] = 0, DoubleArray[1] = 1; DoubleArray[2] = 2;
  HDF5.Write("DoubleArray", H5T_NATIVE_DOUBLE, 3, &DoubleArray[0]);
  
  // To analyze the content of the file, one can use 
  // "h5dump filename.h5" or "h5dump filename.h5 -H" o

  // ============================ //
  // P A R T   II:  R E A D I N G //
  // ============================ //
  
  Epetra_Map* NewMap = 0;
  Epetra_CrsMatrix* NewMatrix = 0;

  // Check if the map is there (in this case it is). If it is, read the
  // matrix using the map, otherwise read the matrix with a linear map.
  if (HDF5.IsContained("map-" + toString(Comm.NumProc())))
  {
    HDF5.Read("map-" + toString(Comm.NumProc()), NewMap);
    HDF5.Read("matrix", *NewMap, *NewMap, NewMatrix);
  }
  else
    HDF5.Read("matrix", *NewMap, *NewMap, NewMatrix);

  // read the number of nonzeros from file, compare them with the
  // actual number of NewMatrix
  int NewNumGlobalNonzeros;
  HDF5.Read("matrix", "NumGlobalNonzeros", NewNumGlobalNonzeros);

  cout << Matrix->NumGlobalNonzeros() << " vs. " << NewNumGlobalNonzeros << endl;

  // We finally close the file. Better to close it before calling
  // MPI_Finalize() to avoid MPI-related errors, since Close() might call MPI
  // functions.
  HDF5.Close();

  // delete memory
  delete Matrix;
  delete Map;
  if (NewMap) delete NewMap;
  if (NewMatrix) delete NewMatrix;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
