#include "Galeri_ConfigDefs.h"
#include "hdf5.h"
#include <iostream>
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <vector>
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Galeri_HDF5.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
using namespace Galeri;

int main (int argc, char **argv)
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MySize = 1 + 2 * Comm.MyPID();
  int GlobalSize;
  Comm.SumAll(&MySize, &GlobalSize, 1);

  // Initialize data buffer 
  vector<int> data(MySize);
  for (int i=0; i < MySize; i++) {
    data[i] = Comm.MyPID() + 10;
  }

  Galeri::HDF5 HDF5(Comm);

  Teuchos::ParameterList GaleriList;
  int nx = 5 * Comm.NumProc();
  int ny = 5 * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());


  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix   = CreateCrsMatrix("Biharmonic2D", Map, GaleriList);
  Epetra_Vector x(Matrix->RowMatrixRowMap());
  Epetra_MultiVector x2(Matrix->RowMatrixRowMap(), 2);
  x.Random();
  x2.Random();

  HDF5.Open("myfile.h5");       // opens the file
  HDF5.Write("map", *Map);        // writes the map 
  HDF5.Write("matrix", *Matrix);  // writes the matrix...
  HDF5.Write("x", x);             // and a vector
  HDF5.Write("x2", x2);             // and a vector
  HDF5.Write("nx", nx);
  HDF5.Write("GaleriList", GaleriList);
  Teuchos::ParameterList NewList;

  HDF5.Read("GaleriList", NewList);
  cout << NewList;

  Epetra_MultiVector* x3;
  Epetra_Map* NewMap;
  Epetra_CrsMatrix* NewMatrix;
  HDF5.Read("map", NewMap);
  HDF5.Read("x2", *NewMap, x3);
  HDF5.Read("matrix", *NewMap, NewMatrix);
  cout << *NewMatrix;

  HDF5.Close();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
