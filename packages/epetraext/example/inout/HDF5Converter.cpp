/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "EpetraExt_ConfigDefs.h"
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
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "EpetraExt_Exception.h"
#include "EpetraExt_Utils.h"
#include "EpetraExt_HDF5.h"
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

int main (int argc, char **argv)
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  if (Comm.MyPID() == 0)
  {
    cout << "Converter from MatrixMarket files to HDF5 files" << endl;
    cout << "For notes on the usage, execute" << endl;
    cout << "  ./HDF5Converter.exe --help" << endl;
    cout << endl;
  }

  // Creating an empty command line processor looks like:
  Teuchos::CommandLineProcessor CLP;

  string MapFileName    = "not-set";
  string XFileName      = "not-set";
  string BFileName      = "not-set";
  string MatrixFileName = "not-set";
  string HDF5FileName   = "myfile.f5";
  string MapHDF5Name    = "map";
  string XHDF5Name      = "X";
  string BHDF5Name      = "B";
  string MatrixHDF5Name = "matrix";

  CLP.setOption("in-map",    &MapFileName,    "map file name");
  CLP.setOption("in-matrix", &MatrixFileName, "matrix file name");
  CLP.setOption("in-x",      &XFileName,      "x vector file name");
  CLP.setOption("in-b",      &BFileName,      "b vector file name");
  CLP.setOption("output",    &HDF5FileName,   "name of HDF5 file");
  CLP.setOption("out-map",    &MapHDF5Name,    "map name in HDF5 file");
  CLP.setOption("out-matrix", &MatrixHDF5Name, "matrix name in HDF5 file");
  CLP.setOption("out-x",      &XHDF5Name,      "x vector name in HDF5 file");
  CLP.setOption("out-b",      &BHDF5Name,      "b vector name in HDF5 file");

  CLP.throwExceptions(false);
  CLP.parse(argc,argv);

  Epetra_Map* Map = 0;
  Epetra_CrsMatrix* Matrix = 0;
  Epetra_MultiVector* X = 0;
  Epetra_MultiVector* B = 0;

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

  if (BFileName != "not-set")
  {
    if (Comm.MyPID() == 0)
      cout << "Reading vector from " << BFileName << endl;

    EpetraExt::MatrixMarketFileToMultiVector(BFileName.c_str(), *Map, B);
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
  
  EpetraExt::HDF5 HDF5(Comm);

  HDF5.Create(HDF5FileName);

  if (Map)
    HDF5.Write(MapHDF5Name + EpetraExt::toString(Comm.NumProc()), *Map);
  if (Matrix)
    HDF5.Write(MatrixHDF5Name, *Matrix);
  if (X)
    HDF5.Write(XHDF5Name, *X);
  if (B)
    HDF5.Write(BHDF5Name, *B);
  HDF5.Close();

  if (Map) delete Map;
  if (Matrix) delete Matrix;
  if (X) delete X;
  if (B) delete B;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
