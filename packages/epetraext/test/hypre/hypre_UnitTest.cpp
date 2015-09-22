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

#include "Teuchos_UnitTestHarness.hpp"
#include "HYPRE_IJ_mv.h"
#include "EpetraExt_HypreIJMatrix.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_Array.hpp"
#include <string>
#include <stdio.h>
#include <map>

namespace EpetraExt {

const int N = 100;
const int MatType = 4; //0 -> Unit diagonal, 1 -> Random diagonal, 2 -> Dense, val=col, 3 -> Random Dense, 4 -> Random Sparse
const double tol = 1E-6;

TEUCHOS_UNIT_TEST( EpetraExt_hypre, Construct ) {

  int ierr = 0;
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, 0);
  
  TEST_EQUALITY(Matrix->Filled(), true);
  
  for(int i = 0; i < Matrix->NumMyRows(); i++){
    int entries;
    ierr += Matrix->NumMyRowEntries(i, entries);
    TEST_EQUALITY(entries, 1);
    int numentries;
    Teuchos::Array<double> Values; Values.resize(entries);
    Teuchos::Array<int> Indices; Indices.resize(entries);
    ierr += Matrix->ExtractMyRowCopy(i, entries, numentries, &Values[0], &Indices[0]);
    TEST_EQUALITY(ierr, 0);
    TEST_EQUALITY(numentries,1);
    for(int j = 0; j < numentries; j++){
      TEST_FLOATING_EQUALITY(Values[j],1.0,tol);
      TEST_EQUALITY(Indices[j],i);
    }
  }
  delete Matrix;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, MatVec ) {
  int ierr = 0;
  int num_vectors = 5;
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, 0);
  
  Epetra_MultiVector X(Matrix->RowMatrixRowMap(), num_vectors, true);
  ierr += X.Random();
  Epetra_MultiVector Y(Matrix->RowMatrixRowMap(), num_vectors, true);
  
  ierr += Matrix->Multiply(false, X, Y);
  
  TEST_EQUALITY(EquivalentVectors(X,Y,tol),true);
  TEST_EQUALITY(ierr, 0);
  delete Matrix;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, BetterMatVec ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  //TestMat->Print(std::cout);
  int num_vectors = 5;
  
  Epetra_MultiVector X(Matrix->RowMatrixRowMap(), num_vectors, true);
  ierr += X.Random();
  Epetra_MultiVector Y1(Matrix->RowMatrixRowMap(), num_vectors, true);
  Epetra_MultiVector Y2(Matrix->RowMatrixRowMap(), num_vectors, true);
  
  ierr += Matrix->Multiply(false, X, Y1);
  ierr += TestMat->Multiply(false, X, Y2);
  
  TEST_EQUALITY(EquivalentVectors(Y1,Y2,tol),true);

  ierr += Matrix->Multiply(false, Y1, X);
  ierr += TestMat->Multiply(false, Y1, Y2);

  TEST_EQUALITY(EquivalentVectors(X,Y2,tol),true);

  ierr += Matrix->Multiply(false, Y2, X);
  ierr += TestMat->Multiply(false, Y2, Y1);

  TEST_EQUALITY(EquivalentVectors(X,Y1,tol),true);
  TEST_EQUALITY_CONST(ierr, 0);
  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, TransposeMatVec ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int num_vectors = 5;
  
  Epetra_MultiVector X(Matrix->RowMatrixRowMap(), num_vectors, true);
  ierr += X.Random();
  Epetra_MultiVector Y1(Matrix->RowMatrixRowMap(), num_vectors, true);
  Epetra_MultiVector Y2(Matrix->RowMatrixRowMap(), num_vectors, true);
  
  ierr += Matrix->Multiply(true, X, Y1);
  ierr += TestMat->Multiply(true, X, Y2);
  
  TEST_EQUALITY(EquivalentVectors(Y1,Y2,tol),true);
  TEST_EQUALITY(ierr, 0);
  
  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, LeftScale ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  //EpetraExt_HypreIJMatrix* BackUp = EpetraExt_HypreIJMatrix(Matrix);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  Epetra_Vector X(Matrix->RowMatrixRowMap(), true);
  ierr += X.Random();
  Matrix->NumMyNonzeros();
  ierr += Matrix->LeftScale(X);
  ierr += TestMat->LeftScale(X);
  TEST_EQUALITY(EquivalentMatrices(*Matrix, *TestMat,tol), true);
  TEST_EQUALITY(ierr, 0);
  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, RightScale ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  Epetra_Vector X(Matrix->RowMatrixRowMap(), true);
  ierr += X.Random();
  Matrix->NumMyNonzeros();
  ierr += Matrix->RightScale(X);
  ierr += TestMat->RightScale(X);
  
  TEST_EQUALITY(EquivalentMatrices(*Matrix,*TestMat,tol), true);
  TEST_EQUALITY(ierr, 0);
  
  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, ExtractDiagonalCopy ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  Epetra_Vector X(Matrix->RowMatrixRowMap(), true);
  Epetra_Vector Y(TestMat->RowMatrixRowMap(),true);
  
  ierr += Matrix->ExtractDiagonalCopy(X);
  ierr += TestMat->ExtractDiagonalCopy(Y);
  
  TEST_EQUALITY(EquivalentVectors(X,Y,tol), true);
  TEST_EQUALITY(ierr, 0);
  
  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, InvRowSums ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  Epetra_Vector X(Matrix->RowMatrixRowMap(), true);
  Epetra_Vector Y(TestMat->RowMatrixRowMap(),true);
  
  ierr += Matrix->InvRowSums(X);
  ierr += TestMat->InvRowSums(Y);
  
  TEST_EQUALITY(EquivalentVectors(X,Y,tol), true);
  TEST_EQUALITY(ierr, 0);
  
  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, InvColSums ) {
  int ierr = 0;
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  Epetra_Vector X(Matrix->RowMatrixColMap(), true);
  Epetra_Vector Y(TestMat->RowMatrixColMap(),true);
  
  ierr += Matrix->InvColSums(X);
  ierr += TestMat->InvColSums(Y);
  
  TEST_EQUALITY(EquivalentVectors(X,Y,tol), true);
  TEST_EQUALITY(ierr, 0);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NormInf ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  double norm1 = Matrix->NormInf();
  double norm2 = TestMat->NormInf();
  
  TEST_FLOATING_EQUALITY(norm1, norm2, tol);
  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NormOne ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  double norm1 = Matrix->NormOne();
  double norm2 = TestMat->NormOne();
  
  TEST_FLOATING_EQUALITY(norm1, norm2, tol);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumGlobalNonzeros ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int nnz1 = Matrix->NumGlobalNonzeros();
  int nnz2 = TestMat->NumGlobalNonzeros();
  
  TEST_EQUALITY(nnz1, nnz2);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumGlobalRows ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int rows1 = Matrix->NumGlobalRows();
  int rows2 = TestMat->NumGlobalRows();
  
  TEST_EQUALITY(rows1, rows2);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumGlobalCols ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int cols1 = Matrix->NumGlobalCols();
  int cols2 = TestMat->NumGlobalCols();
  
  TEST_EQUALITY(cols1, cols2);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumGlobalDiagonals ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int hdiag1 = Matrix->NumGlobalDiagonals();
  int Ediag2 = TestMat->NumGlobalDiagonals();
  
  TEST_EQUALITY(hdiag1, Ediag2);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumMyNonzeros ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int nnz1 = Matrix->NumMyNonzeros();
  int nnz2 = TestMat->NumMyNonzeros();
  
  TEST_EQUALITY(nnz1, nnz2);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumMyRows ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int rows1 = Matrix->NumMyRows();
  int rows2 = TestMat->NumMyRows();
  
  TEST_EQUALITY(rows1, rows2);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumMyCols ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int cols1 = Matrix->NumMyCols();
  int cols2 = TestMat->NumMyCols();
  
  TEST_EQUALITY(cols1, cols2);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, NumMyDiagonals ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int diag1 = Matrix->NumMyDiagonals();
  int diag2 = TestMat->NumMyDiagonals();
  
  TEST_EQUALITY(diag1, diag2);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, MaxNumEntries ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int ent1 = Matrix->MaxNumEntries();
  int ent2 = TestMat->MaxNumEntries();
  
  TEST_EQUALITY(ent1, ent2);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, ApplyInverse ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, 1);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);

  int num_vectors = 1;
  Epetra_MultiVector X(Matrix->RowMatrixRowMap(), num_vectors, true);
  X.Random();
  Epetra_MultiVector Y1(Matrix->RowMatrixRowMap(), num_vectors, true);
  Epetra_MultiVector Y2(Matrix->RowMatrixRowMap(), num_vectors, true);
  
  TestMat->ApplyInverse(X,Y2);
  Matrix->ApplyInverse(X,Y1);
  TEST_EQUALITY(EquivalentVectors(Y1,Y2,tol), true);

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, SameMatVec ) {
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, MatType);
  Epetra_CrsMatrix* TestMat = newCrsMatrix(*Matrix);
  
  int num_vectors = 3;
  Epetra_MultiVector X1(Matrix->RowMatrixRowMap(), num_vectors, true);
  X1.Random();
  Epetra_MultiVector X2(X1);
  
  Matrix->Multiply(false,X1,X1);
  TestMat->Multiply(false,X2,X2);

  TEST_EQUALITY(EquivalentVectors(X1,X2,tol), true); 

  delete Matrix;
  delete TestMat;
}

TEUCHOS_UNIT_TEST( EpetraExt_hypre, Solve ) {
  int num_vectors = 2;
  EpetraExt_HypreIJMatrix* Matrix = newHypreMatrix(N, 1);
  Epetra_MultiVector TrueX(Matrix->RowMatrixRowMap(), num_vectors, true);
  TrueX.Random();
  Epetra_MultiVector RHS(Matrix->RowMatrixRowMap(), num_vectors, true);
  Matrix->Multiply(false,TrueX, RHS);
  Epetra_MultiVector X(Matrix->RowMatrixRowMap(), num_vectors, true);
  
  Matrix->SetParameter(Solver, PCG);
  Matrix->SetParameter(Preconditioner, BoomerAMG);

  /* Set some parameters (See Reference Manual for more parameters) */
  Matrix->SetParameter(Solver, &HYPRE_PCGSetMaxIter, 1000); /* max iterations */
  Matrix->SetParameter(Solver, &HYPRE_PCGSetTol, 1e-7); /* conv. tolerance */
  Matrix->SetParameter(Solver, &HYPRE_PCGSetTwoNorm, 1); /* use the two norm as the stopping criteria */
  Matrix->SetParameter(Solver, &HYPRE_PCGSetPrintLevel, 0); /* print solve info */
  Matrix->SetParameter(Solver, &HYPRE_PCGSetLogging, 0); /* needed to get run info later */

  /* Now set up the AMG preconditioner and specify any parameters */
  Matrix->SetParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, 0); /* print amg solution info */
  Matrix->SetParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, 6);
  Matrix->SetParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, 6); /* Sym G.S./Jacobi hybrid */ 
  Matrix->SetParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, 1);
  Matrix->SetParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, 0.0); /* conv. tolerance zero */
  Matrix->SetParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, 10); /* do only one iteration! */

  Matrix->SetParameter(true);
  /* Now setup and solve! */
  Matrix->Solve(false, false, false, RHS, X);
  TEST_EQUALITY(EquivalentVectors(X,TrueX,tol), true);
  
  delete Matrix;
}
} // namespace EpetraExt
