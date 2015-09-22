/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#include "Ifpack_ConfigDefs.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Ifpack_OverlappingRowMatrix.h"
#include "Ifpack_LocalFilter.h"
#include "Ifpack_Utils.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  if (Comm.NumProc() == 1)
  {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    cout << "Test `TestOverlappingRowMatrix.exe' passed!" << endl;
    exit(EXIT_SUCCESS);
  }

  Teuchos::ParameterList GaleriList;
  int nx = 100;
  GaleriList.set("n", nx * nx);
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx);
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap64("Linear", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( Galeri::CreateCrsMatrix("Laplace2D", &*Map, GaleriList) );

  int OverlapLevel = 5;
  Epetra_Time Time(Comm);

  // ======================================== //
  // Build the overlapping matrix using class //
  // Ifpack_OverlappingRowMatrix.             //
  // ======================================== //

  Time.ResetStartTime();
  Ifpack_OverlappingRowMatrix B(A,OverlapLevel);
  if (Comm.MyPID() == 0)
    cout << "Time to create B = " << Time.ElapsedTime() << endl;

  long long NumGlobalRowsB = B.NumGlobalRows64();
  long long NumGlobalNonzerosB = B.NumGlobalNonzeros64();

  Epetra_Vector X(A->RowMatrixRowMap());
  Epetra_Vector Y(A->RowMatrixRowMap());
  for (int i = 0 ; i < A->NumMyRows() ; ++i)
    X[i] = 1.0* A->RowMatrixRowMap().GID64(i);
  Y.PutScalar(0.0);

  Epetra_Vector ExtX_B(B.RowMatrixRowMap());
  Epetra_Vector ExtY_B(B.RowMatrixRowMap());
  ExtY_B.PutScalar(0.0);

  IFPACK_CHK_ERR(B.ImportMultiVector(X,ExtX_B));
  IFPACK_CHK_ERR(B.Multiply(false,ExtX_B,ExtY_B));
  IFPACK_CHK_ERR(B.ExportMultiVector(ExtY_B,Y,Add));

  double Norm_B;
  Y.Norm2(&Norm_B);
  if (Comm.MyPID() == 0)
    cout << "Norm of Y using B = " << Norm_B << endl;

  // ================================================== //
  //Build the overlapping matrix as an Epetra_CrsMatrix //
  // ================================================== //

  Time.ResetStartTime();
  Epetra_CrsMatrix& C =
    *(Ifpack_CreateOverlappingCrsMatrix(&*A,OverlapLevel));
  if (Comm.MyPID() == 0)
    cout << "Time to create C = " << Time.ElapsedTime() << endl;

  // simple checks on global quantities
  long long NumGlobalRowsC = C.NumGlobalRows64();
  long long NumGlobalNonzerosC = C.NumGlobalNonzeros64();
  if (NumGlobalRowsB != NumGlobalRowsC) {
    std::ostringstream os;
    os << "NumGlobalRowsB = " << NumGlobalRowsB
       << " != NumGlobalRowsC = " << NumGlobalRowsC << ".";
    throw std::logic_error (os.str ());
  }
  if (NumGlobalNonzerosB != NumGlobalNonzerosC) {
    std::ostringstream os;
    os << "NumGlobalNonzerosB = " << NumGlobalNonzerosB
       << " != NumGlobalNonzerosC = " << NumGlobalNonzerosC << ".";
    throw std::logic_error (os.str ());
  }

  Epetra_Vector ExtX_C(C.RowMatrixRowMap());
  Epetra_Vector ExtY_C(C.RowMatrixRowMap());
  ExtY_C.PutScalar(0.0);
  Y.PutScalar(0.0);

  IFPACK_CHK_ERR(C.Multiply(false,X,Y));

  double Norm_C;
  Y.Norm2(&Norm_C);
  if (Comm.MyPID() == 0)
    cout << "Norm of Y using C = " << Norm_C << endl;

  if (IFPACK_ABS(Norm_B - Norm_C) > 1e-5)
    IFPACK_CHK_ERR(-1);

  // ======================= //
  // now localize the matrix //
  // ======================= //

  Ifpack_LocalFilter D(Teuchos::rcp(&B, false));

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (Comm.MyPID() == 0)
    cout << "Test `TestOverlappingRowMatrix.exe' passed!" << endl;

  return(EXIT_SUCCESS);
}
