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
  int nx = 10 * Comm.NumProc();
  int ny = 10 * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map;
  Epetra_CrsMatrix* Matrix;
  Epetra_Vector* x = 0,* readb = 0,* readxexact = 0;
  Galeri::ReadHB(ProblemType.c_str(), Comm, Map, Matrix, x, readb, readxexact);

#if 0
  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix   = CreateCrsMatrix("Biharmonic2D", Map, GaleriList);
  Epetra_Vector x(Matrix->RowMatrixRowMap());
  x.Random();
  Teuchos::ParameterList NewList;
#endif

  HDF5.Create("myfile.h5");       // opens the file
  HDF5.Write("map-" + toString(Comm.NumProc()), *Map);        // writes the map 
  HDF5.Write("matrix", *Matrix);  // writes the matrix...
  HDF5.Write("x", x);             // and a vector
  HDF5.Write("GaleriList", GaleriList);

#if 0
  cout << HDF5.IsDataSet("nx") << endl;
  cout << HDF5.IsDataSet("n") << endl;

  vector<double> array(10);
  for (int i = 0; i < 10; ++i)
    array[i] = 1000.0;
  HDF5.Write("/array", H5T_NATIVE_DOUBLE, 10, &array[0]);

  HDF5.Read("GaleriList", NewList);
  cout << NewList;

  Epetra_Map* NewMap;
  if (HDF5.IsDataSet("map-" + toString(Comm.NumProc())))
  {
    HDF5.Read("map-" + toString(Comm.NumProc()), NewMap);
    cout << *NewMap;
  }
  Epetra_CrsMatrix* NewMatrix;
  Epetra_Vector* NewX;
  HDF5.Read("x", NewX);
 // cout << *NewX;
  //HDF5.Read("x", *NewMap, NewX);
  HDF5.Read("matrix", NewMatrix);
//  cout << *NewMatrix;
#endif

  HDF5.Close();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
