#include "Galeri_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <vector>
#include "Teuchos_CommandLineProcessor.hpp"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Galeri_HDF5.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_CrsMatrixIn.h"

// converts file from EpetraExt format into HDF5 format.
//
// \author Marzio Sala, D-INFK/ETHZ
//
// \date Last updated on 09-Mar-06.

using namespace Galeri;

int main (int argc, char **argv)
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Creating an empty command line processor looks like:
  Teuchos::CommandLineProcessor CLP;

  string MapFileName    = "not-set";
  string XFileName      = "not-set";
  string MatrixFileName = "not-set";
  string HDF5FileName   = "myfile.f5";
  string MapHDF5Name    = "map";
  string XHDF5Name      = "X";
  string MatrixHDF5Name = "matrix";


  CLP.setOption("in-map",    &MapFileName,    "map file name");
  CLP.setOption("in-matrix", &MatrixFileName, "matrix file name");
  CLP.setOption("in-vector", &XFileName,      "vector file name");
  CLP.setOption("output",    &HDF5FileName,   "name of HDF5 file");
  CLP.setOption("out-map",    &MapHDF5Name,    "map name in HDF5 file");
  CLP.setOption("out-matrix", &MatrixHDF5Name, "matrix name in HDF5 file");
  CLP.setOption("out-vector", &XHDF5Name,      "vector name in HDF5 file");


  CLP.throwExceptions(false);
  CLP.parse(argc,argv);

  Epetra_Map* Map = 0;
  Epetra_CrsMatrix* Matrix = 0;
  Epetra_MultiVector* X = 0;

  if (MapFileName != "not-set")
  {
    if (Comm.MyPID() == 0)
      cout << "Reading map from " << MapFileName << endl;

    EpetraExt::MatrixMarketFileToMap(MapFileName.c_str(), Comm, Map);
  }
  else
  {
    cerr << "You need to specify a map, sorry" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }

  if (XFileName != "not-set")
  {
    if (Comm.MyPID() == 0)
      cout << "Reading vector from " << XFileName << endl;

    EpetraExt::MatrixMarketFileToMultiVector(XFileName.c_str(), *Map, X);
  }

  if (MatrixFileName != "not-set")
  {
    if (Comm.MyPID() == 0)
      cout << "Reading matrix from " << MatrixFileName << endl;

    EpetraExt::MatrixMarketFileToCrsMatrix(MatrixFileName.c_str(), *Map, Matrix);
  }

  // ================================= //
  // Open HDF5 file and append data in //
  // ================================= //
  
  Galeri::HDF5 HDF5(Comm);

  HDF5.Create(HDF5FileName);

  if (Map)
    HDF5.Write("map-" + toString(Comm.NumProc()), *Map);
  if (Matrix)
    HDF5.Write("matrix", *Matrix);
  if (X)
    HDF5.Write("X", *(*X)(0));
  HDF5.Close();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
