#include "EpetraExt_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <vector>
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_HDF5.h"
#include "EpetraExt_Utils.h"
#include "EpetraExt_Exception.h"

// Read a matrix from file "matlab.h5". The matrix has been generated
// as sparse matrix in MATLAB, then save to file using MATLAB's commands;
// see the Doxygen documentation of class EpetraExt::HDF5 for the MATLAB
// commands used.
//
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

  try 
  {
    // This is the HDF5 file manager
    EpetraExt::HDF5 HDF5(Comm);

    // creates a new file. To open an existing file, use Open("myfile.h5")
    // This file contains:
    // - a sparse (diagonal) matrix, whose group name is "speye"
    // - a multivector, whose group name is "x"
    // - a map for 2-processor run, whose group name is "map-2"

    HDF5.Open("matlab.h5");

    if (Comm.MyPID() == 0)
      cout << endl;
      cout << "*) Reading Epetra_CrsMatrix from HDF5 file matlab.h5..." << endl;
      cout << endl;

    // first query for matrix properties:
    int NumGlobalRows, NumGlobalCols, NumGlobalNonzeros;
    int NumGlobalDiagonals, MaxNumEntries;
    double NormOne, NormInf;

    HDF5.ReadCrsMatrixProperties("speye", NumGlobalRows, NumGlobalCols,
                                 NumGlobalNonzeros, NumGlobalDiagonals, 
                                 MaxNumEntries, NormOne, NormInf);

    if (Comm.MyPID() == 0)
    {
      cout << "Matrix information as given by ReadCrsMatrixProperties()";
      cout << endl << endl;
      cout << "NumGlobalRows = " << NumGlobalRows << endl;
      cout << "NumGlobalCols = " << NumGlobalCols << endl;
      cout << "NumGlobalNonzeros = " << NumGlobalNonzeros << endl;
      cout << "NumGlobalDiagonals = " << NumGlobalDiagonals << endl;
      cout << "MaxNumEntries = " << MaxNumEntries << endl;
      cout << "NormOne = " << NormOne << endl;
      cout << "NormInf = " << NormInf << endl;
    }

    // the reading the actual matrix, with a linear map, since no map 
    // has been specified.
    Epetra_CrsMatrix* Matrix = 0;
    HDF5.Read("speye", Matrix);

    cout << *Matrix;

    if (Comm.MyPID() == 0)
    {
      cout << endl;
      cout << "*) Reading Epetra_MultiVector from HDF5 file matlab.h5..." << endl;
      cout << endl;
    }

    Epetra_MultiVector* x;
    HDF5.Read("x", x);
    cout << *x;

    if (Comm.NumProc() == 2)
    {
      if (Comm.MyPID() == 0)
      {
        cout << endl;
        cout << "*) Reading Epetra_Map from HDF5 file matlab.h5..." << endl;
        cout << endl;
      }

      Epetra_Map* Map;
      HDF5.Read("map-2", Map);
      cout << *Map;
    }

    // We finally close the file. Better to close it before calling
    // MPI_Finalize() to avoid MPI-related errors, since Close() might call MPI
    // functions.
    HDF5.Close();

    // delete memory
    if (Matrix) delete Matrix;
  }
  catch(EpetraExt::Exception& rhs) 
  {
    rhs.Print();
  }
  catch (...) 
  {
    cerr << "Caught generic exception" << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
