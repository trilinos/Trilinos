//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "HYPRE_IJ_mv.h"
#include "EpetraExt_HypreIJMatrix.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "hypre_Helpers.hpp"
#include "Teuchos_Array.hpp"
#include <string>
#include <stdio.h>
#include <map>

namespace EpetraExt {

const int N = 10;
const int MatType = 3; //0 -> Unit diagonal, 1 -> Dense, val=col, 2 -> Random Dense, 3 -> Random Sparse
const double tol = 1E-6;

TEUCHOS_UNIT_TEST( EpetraExt_hypre, Construct ) {

  int ierr = 0;
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, 0);
  
  TEST_EQUALITY(Matrix.Filled(), true);
  
  for(int i = 0; i < Matrix.NumMyRows(); i++){
    int entries;
    ierr += Matrix.NumMyRowEntries(i, entries);
    TEST_EQUALITY(entries, 1);
    int numentries;
    Teuchos::Array<double> Values; Values.resize(entries);
    Teuchos::Array<int> Indices; Indices.resize(entries);
    ierr += Matrix.ExtractMyRowCopy(i, entries, numentries, &Values[0], &Indices[0]);
    TEST_EQUALITY(ierr, 0);
    TEST_EQUALITY(numentries,1);
    for(int j = 0; j < numentries; j++){
      TEST_FLOATING_EQUALITY(Values[j],1.0,tol);
      TEST_EQUALITY(Indices[j],i);
    }
  }
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, MatVec ) {
  int ierr = 0;
  int num_vectors = 5;
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, 0);
  
  Epetra_MultiVector X(Matrix.RowMatrixRowMap(), num_vectors, true);
  ierr += X.Random();
  Epetra_MultiVector Y(Matrix.RowMatrixRowMap(), num_vectors, true);
  
  ierr += Matrix.Multiply(false, X, Y);
  
  TEST_EQUALITY(EquivalentVectors(X,Y,tol),true);
  TEST_EQUALITY(ierr, 0);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, BetterMatVec ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  //TestMat.Print(std::cout);
  int num_vectors = 5;
  
  Epetra_MultiVector X(Matrix.RowMatrixRowMap(), num_vectors, true);
  ierr += X.Random();
  Epetra_MultiVector Y1(Matrix.RowMatrixRowMap(), num_vectors, true);
  Epetra_MultiVector Y2(Matrix.RowMatrixRowMap(), num_vectors, true);
  
  ierr += Matrix.Multiply(false, X, Y1);
  ierr += TestMat.Multiply(false, X, Y2);
  
  TEST_EQUALITY(EquivalentVectors(Y1,Y2,tol),true);
  TEST_EQUALITY_CONST(ierr, 0);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, TransposeMatVec ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  int num_vectors = 5;
  
  Epetra_MultiVector X(Matrix.RowMatrixRowMap(), num_vectors, true);
  ierr += X.Random();
  Epetra_MultiVector Y1(Matrix.RowMatrixRowMap(), num_vectors, true);
  Epetra_MultiVector Y2(Matrix.RowMatrixRowMap(), num_vectors, true);
  
  ierr += Matrix.Multiply(true, X, Y1);
  ierr += TestMat.Multiply(true, X, Y2);
  
  TEST_EQUALITY(EquivalentVectors(Y1,Y2,tol),true);
  TEST_EQUALITY(ierr, 0);
  
}

/*
TEUCHOS_UNIT_TEST( EpetraExt_hypre, ExtractEntry ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  
  int CurrEntry = 14;
  
  double *Value;
  int RowIndex;
  int ColIndex;
  ierr += Matrix.ExtractMyEntryView(CurrEntry, Value, RowIndex, ColIndex);
  printf("CurrIndex[%d] => [%d,%d] = %f", CurrEntry, RowIndex, ColIndex, Value[0]);
  TEST_EQUALITY(ierr, 0);
}
*/
TEUCHOS_UNIT_TEST( EpetraExt_hypre, LeftScale ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  //EpetraExt_HypreIJMatrix BackUp = EpetraExt_HypreIJMatrix(Matrix);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  Epetra_Vector X(Matrix.RowMatrixRowMap(), true);
  ierr += X.Random();
  Matrix.NumMyNonzeros();
  ierr += Matrix.LeftScale(X);
  ierr += TestMat.LeftScale(X);
  TEST_EQUALITY(EquivalentMatrices(Matrix,TestMat,tol), true);
  TEST_EQUALITY(ierr, 0);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, RightScale ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  Epetra_Vector X(Matrix.RowMatrixRowMap(), true);
  ierr += X.Random();
  Matrix.NumMyNonzeros();
  ierr += Matrix.RightScale(X);
  ierr += TestMat.RightScale(X);
  
  TEST_EQUALITY(EquivalentMatrices(Matrix,TestMat,tol), true);
  TEST_EQUALITY(ierr, 0);
  
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, ExtractDiagonalCopy ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  Epetra_Vector X(Matrix.RowMatrixRowMap(), true);
  Epetra_Vector Y(TestMat.RowMatrixRowMap(),true);
  
  ierr += Matrix.ExtractDiagonalCopy(X);
  ierr += TestMat.ExtractDiagonalCopy(Y);
  
  TEST_EQUALITY(EquivalentVectors(X,Y,tol), true);
  TEST_EQUALITY(ierr, 0);
  
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, InvRowSums ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  Epetra_Vector X(Matrix.RowMatrixRowMap(), true);
  Epetra_Vector Y(TestMat.RowMatrixRowMap(),true);
  
  ierr += Matrix.InvRowSums(X);
  ierr += TestMat.InvRowSums(Y);
  
  TEST_EQUALITY(EquivalentVectors(X,Y,tol), true);
  TEST_EQUALITY(ierr, 0);
  
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, InvColSums ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  Epetra_Vector X(Matrix.RowMatrixColMap(), true);
  Epetra_Vector Y(TestMat.RowMatrixColMap(),true);
  
  ierr += Matrix.InvColSums(X);
  //ierr += TestMat.InvColSums(Y);
  
  //TEST_EQUALITY(EquivalentVectors(X,Y,tol), true);
  TEST_EQUALITY(ierr, 0);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NormInf ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  double norm1 = Matrix.NormInf();
  double norm2 = TestMat.NormInf();
  
  TEST_FLOATING_EQUALITY(norm1, norm2, tol);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NormOne ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  double norm1 = Matrix.NormOne();
  double norm2 = TestMat.NormOne();
  
  TEST_FLOATING_EQUALITY(norm1, norm2, tol);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumGlobalNonzeros ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  int nnz1 = Matrix.NumGlobalNonzeros();
  int nnz2 = TestMat.NumGlobalNonzeros();
  
  TEST_EQUALITY(nnz1, nnz2);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumGlobalRows ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  int rows1 = Matrix.NumGlobalRows();
  int rows2 = TestMat.NumGlobalRows();
  
  TEST_EQUALITY(rows1, rows2);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumGlobalCols ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  int cols1 = Matrix.NumGlobalCols();
  int cols2 = TestMat.NumGlobalCols();
  
  TEST_EQUALITY(cols1, cols2);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumGlobalDiagonals ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  int hdiag1 = Matrix.NumGlobalDiagonals();
  int Ediag2 = TestMat.NumGlobalDiagonals();
  
  TEST_EQUALITY(hdiag1, Ediag2);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumMyNonzeros ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  int nnz1 = Matrix.NumMyNonzeros();
  int nnz2 = TestMat.NumMyNonzeros();
  
  TEST_EQUALITY(nnz1, nnz2);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumMyRows ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  int rows1 = Matrix.NumMyRows();
  int rows2 = TestMat.NumMyRows();
  
  TEST_EQUALITY(rows1, rows2);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumMyCols ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  int cols1 = Matrix.NumMyCols();
  int cols2 = TestMat.NumMyCols();
  
  TEST_EQUALITY(cols1, cols2);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumMyDiagonals ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  int diag1 = Matrix.NumMyDiagonals();
  int diag2 = TestMat.NumMyDiagonals();
  
  TEST_EQUALITY(diag1, diag2);
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, MaxNumEntries ) {
  EpetraExt_HypreIJMatrix Matrix = MatrixConstructor(N, MatType);
  Epetra_CrsMatrix TestMat = GetCrsMatrix(Matrix);
  
  int ent1 = Matrix.MaxNumEntries();
  int ent2 = TestMat.MaxNumEntries();
  
  TEST_EQUALITY(ent1, ent2);
}
} // namespace EpetraExt
