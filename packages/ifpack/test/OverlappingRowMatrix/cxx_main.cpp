// @HEADER
// ***********************************************************************
// 
//                IFPACK
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Ifpack_ConfigDefs.h"
// FIXME: this test has to be re-enabled
#if defined(HAVE_IFPACK_AZTECOO) && defined(HAVE_IFPACK_AMESOS) && defined(HAVE_IFPACK_TEUCHOS) 
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
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack_OverlappingRowMatrix.h"
#include "Ifpack_LocalFilter.h"
#include "Ifpack_Utils.h"

using namespace Trilinos_Util;

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  if (Comm.NumProc() == 1)
    exit(EXIT_SUCCESS);

  // size of the global matrix. 
  const int NumPoints = 10000;

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type", "linear");

  Epetra_RowMatrix& A = *(Gallery.GetMatrix());
  int OverlapLevel = 5;
  Epetra_Time Time(Comm);

  // ======================================== //
  // Build the overlapping matrix using class //
  // Ifpack_OverlappingRowMatrix.             //
  // ======================================== //
 
  Time.ResetStartTime();
  Ifpack_OverlappingRowMatrix B(&A,OverlapLevel);
  if (Comm.MyPID() == 0)
    cout << "Time to create B = " << Time.ElapsedTime() << endl;

  int NumGlobalRowsB = B.NumGlobalRows();
  int NumGlobalNonzerosB = B.NumGlobalNonzeros();

  Epetra_Vector X(A.RowMatrixRowMap());
  Epetra_Vector Y(A.RowMatrixRowMap());
  for (int i = 0 ; i < A.NumMyRows() ; ++i) 
    X[i] = 1.0* A.RowMatrixRowMap().GID(i);
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
    *(Ifpack_CreateOverlappingCrsMatrix(&A,OverlapLevel));
  if (Comm.MyPID() == 0)
    cout << "Time to create C = " << Time.ElapsedTime() << endl;

  // simple checks on global quantities
  int NumGlobalRowsC = C.NumGlobalRows();
  int NumGlobalNonzerosC = C.NumGlobalNonzeros();
  assert (NumGlobalRowsB == NumGlobalRowsC);
  assert (NumGlobalNonzerosB == NumGlobalNonzerosC);

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

  Ifpack_LocalFilter D(&B);

#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  if (Comm.MyPID() == 0)
    cout << "Test `TestOverlappingRowMatrix.exe' passed!" << endl;
  exit(EXIT_SUCCESS);
}

#else

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  if (Comm.MyPID() == 0) {
    puts("please configure IFPACK with --eanble-aztecoo --enable-teuchos");
    puts("--enable-amesos to run this test");
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return(EXIT_SUCCESS);
}

#endif
