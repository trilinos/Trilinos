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
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
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
  Epetra_BlockMap* BlockMap;
  Epetra_CrsMatrix* Matrix;
  Epetra_Vector* x = 0,* readb = 0,* readxexact = 0;
//  Galeri::ReadHB(ProblemType.c_str(), Comm, Map, Matrix, x, readb, readxexact);
  Epetra_Time Time(Comm);
  Time.ResetStartTime();
  EpetraExt::MatrixMarketFileToBlockMap("/Users/marzio/matrices/PREMO/Falcon-step2/Falcon_ss_1_map", Comm, BlockMap);
  EpetraExt::BlockMapToMatrixMarketFile("map.mm", *BlockMap);

  Map = new Epetra_Map(-1, BlockMap->NumMyElements(), BlockMap->MyGlobalElements(), 0, Comm);

  cout << "TIME FOR EPETRAEXT MAP = " << Time.ElapsedTime() << endl;

#if 0
  Time.ResetStartTime();
  EpetraExt::MatrixMarketFileToCrsMatrix("/Users/marzio/matrices/PREMO/Falcon-step2/Falcon_ss_4_matrix", *Map, Matrix);
  cout << "TIME FOR EPETRAEXT READING MATRIX = " << Time.ElapsedTime() << endl;

  Time.ResetStartTime();
  Matrix->FillComplete();
#endif

#if 0
  EpetraExt::RowMatrixToMatrixMarketFile("matrix.mm", *Matrix);
  cout << "TIME FOR EPETRAEXT MATRIX = " << Time.ElapsedTime() << endl;
#endif

#if 0
  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* Matrix   = CreateCrsMatrix("Biharmonic2D", Map, GaleriList);
  Epetra_Vector x(Matrix->RowMatrixRowMap());
  x.Random();
  Teuchos::ParameterList NewList;
#endif


  Epetra_Map* NewMap;

  Time.ResetStartTime();
  HDF5.Open("myfile.h5");       // opens the file
  //HDF5.Write("map-" + toString(Comm.NumProc()), *Map);        // writes the map 
  //HDF5.Read("map-" + toString(Comm.NumProc()), NewMap);
  //cout << "HDF TIME FOR MAP = " << Time.ElapsedTime() << endl;
  //Time.ResetStartTime();

  //HDF5.Write("matrix", *Matrix);  // writes the matrix...
  //cout << "TIME FOR WRITING " << Time.ElapsedTime() << endl;
  //Time.ResetStartTime();
  HDF5.Read("matrix", Matrix);
  cout << "TIME FOR READ " << Time.ElapsedTime() << endl;

  //HDF5.Write("x", x);             // and a vector
  //HDF5.Write("GaleriList", GaleriList);

  cout << HDF5.IsContained("nx") << endl;
  cout << HDF5.IsContained("n") << endl;

#if 0
  vector<int> array(10);
  for (int i = 0; i < 10; ++i)
    array[i] = 1000;
  HDF5.Write("/array", H5T_NATIVE_INT, 10, &array[0]);

  HDF5.Read("GaleriList", NewList);
  cout << NewList;

  if (HDF5.IsContained("map-" + toString(Comm.NumProc())))
  {
    HDF5.Read("map-" + toString(Comm.NumProc()), NewMap);
    cout << *NewMap;
  }
  Epetra_CrsMatrix* NewMatrix;
  Epetra_Vector* NewX;
  HDF5.Read("x", NewX);
 // cout << *NewX;
  HDF5.Read("x", *NewMap, NewX);
  HDF5.Read("matrix", NewMatrix);
//  cout << *NewMatrix;
#endif
  cout << "=================================================" << endl;
  HDF5.Close();
  cout << "TIME FOR EPETRAEXT = " << Time.ElapsedTime() << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
