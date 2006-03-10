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
#include "Epetra_MultiVector.h"
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

  try {

    // create a map, two vectors, a matrix using Galeri
    Teuchos::ParameterList GaleriList;
    int nx = 2;
    int ny = 2 * Comm.NumProc();
    GaleriList.set("nx", nx);
    GaleriList.set("ny", ny);
    GaleriList.set("mx", 1);
    GaleriList.set("my", Comm.NumProc());

    Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
    Epetra_CrsMatrix* Matrix = CreateCrsMatrix("Biharmonic2D", Map, GaleriList);
    Epetra_MultiVector x(Matrix->RowMatrixRowMap(), 2);
    Epetra_MultiVector b(Matrix->RowMatrixRowMap(), 2);
    x.Random();
    b.Random();

    // This is the HDF5 file manager
    Galeri::HDF5 HDF5(Comm);

    // creates a new file. To open an existing file, use Open("myfile.h5")
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
    // we can associate basic data types with a given group, for example:
    HDF5.Write("matrix", "integration order", 1);
    HDF5.Write("matrix", "numerical drop", 0.1);
    // or put them in a new group
    HDF5.Write("my parameters", "latitude", 12);
    HDF5.Write("my parameters", "longitude", 67);

    // Now we write the parameter list we have used to create the matrix
    HDF5.Write("GaleriList", GaleriList);
    // We can also write integers/doubles/arrays. All these quantities are
    // supposed to have the same value on all processors.
    vector<int> iarray(3); 
    iarray[0] = 0, iarray[1] = 1; iarray[2] = 2;
    HDF5.Write("my parameters", "int array", H5T_NATIVE_INT, 3, &iarray[0]);

    vector<double> darray(3); 
    darray[0] = 0.1, darray[1] = 1.1; darray[2] = 2.1;
    HDF5.Write("my parameters", "double array", H5T_NATIVE_DOUBLE, 3, &darray[0]);

    // To analyze the content of the file, one can use 
    // "h5dump filename.h5" or "h5dump filename.h5 -H" o

    // ============================ //
    // P A R T   II:  R E A D I N G //
    // ============================ //

    Epetra_Map* NewMap = 0;
    Epetra_CrsMatrix* NewMatrix = 0;
    Epetra_MultiVector* NewX,* NewB;

    // Check if the map is there (in this case it is). If it is, read the
    // matrix using the map, otherwise read the matrix with a linear map.
    if (HDF5.IsContained("map-" + toString(Comm.NumProc())))
    {
      HDF5.Read("map-" + toString(Comm.NumProc()), NewMap);
      HDF5.Read("matrix", *NewMap, *NewMap, NewMatrix);
    }
    else
    {
      int NumGlobalRows;
      HDF5.Read("matrix", "NumGlobalRows", NumGlobalRows);
      NewMap = new Epetra_Map(NumGlobalRows, 0, Comm);
      HDF5.Read("matrix", *NewMap, *NewMap, NewMatrix);
    }

    // read the number of nonzeros from file, compare them with the
    // actual number of NewMatrix
    int NewNumGlobalNonzeros;
    HDF5.Read("matrix", "NumGlobalNonzeros", NewNumGlobalNonzeros);

    if (Comm.MyPID() == 0)
      cout << Matrix->NumGlobalNonzeros() << " vs. " << NewNumGlobalNonzeros << endl;

    // Now read the MultiVector's
    HDF5.Read("x", NewX);
    HDF5.Read("b", NewB);

    // the int/double values, and the int/double arrays
    int new_latitude, new_longitude;
    vector<int>    new_iarray(3);
    vector<double> new_darray(3);
    HDF5.Read("my parameters", "latitude", new_latitude);
    HDF5.Read("my parameters", "longitude", new_longitude);
    HDF5.Read("my parameters", "int array", H5T_NATIVE_INT, 3, &new_iarray[0]);
    HDF5.Read("my parameters", "double array", H5T_NATIVE_DOUBLE, 3, &new_darray[0]);

    if (Comm.MyPID() == 0)
    {
      for (int i = 0; i < 3; ++i)
        cout << "iarray[" << i << "] = " << iarray[i] << " vs. " << new_iarray[i] << endl;
      for (int i = 0; i < 3; ++i)
        cout << "darray[" << i << "] = " << darray[i] << " vs. " << new_darray[i] << endl;
    }

    // We finally close the file. Better to close it before calling
    // MPI_Finalize() to avoid MPI-related errors, since Close() might call MPI
    // functions.
    HDF5.Close();

    // delete memory
    delete Matrix;
    delete Map;
    if (NewMap) delete NewMap;
    if (NewMatrix) delete NewMatrix;
    if (NewX) delete NewX;
    if (NewB) delete NewB;

  }
  catch(Exception& rhs) 
  {
    rhs.Print();
  }
  catch (...) {
    cerr << "Caught generic exception" << endl;
  }


#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
