// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_core_Workspace.h"

#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_LAPACK.h"

int Galeri::core::Workspace::numDimensions_ = 3;

const int Galeri::core::Workspace::MIN                  = -1;
const int Galeri::core::Workspace::MAX                  = -2;
const int Galeri::core::Workspace::UNINITIALIZED        = 1550;
const int Galeri::core::Workspace::INITIALIZED          = 1551;
const int Galeri::core::Workspace::CONNECTIVITY_FREEZED = 1552;
const int Galeri::core::Workspace::COORDINATES_FREEZED  = 1553;

// ============================================================================ 
void Galeri::core::Workspace::solve_LAPACK(Epetra_RowMatrix& matrix,
                                                  Epetra_MultiVector& LHS,
                                                  Epetra_MultiVector& RHS)
{
  TEUCHOS_TEST_FOR_EXCEPTION(matrix.Comm().NumProc() != 1, std::logic_error,
                     "solve_LAPACK() works only in serial");

  TEUCHOS_TEST_FOR_EXCEPTION(LHS.NumVectors() != RHS.NumVectors(), std::logic_error,
                     "number of vectors in multivectors not consistent");

  int n = matrix.NumGlobalRows();
  int NumVectors = LHS.NumVectors();

  Epetra_SerialDenseMatrix DenseMatrix;
  DenseMatrix.Shape(n, n);

  for (int i = 0 ; i < n ; ++i)
    for (int j = 0 ; j < n ; ++j)
      DenseMatrix(i,j) = 0.0;

  // allocate storage to extract matrix rows.
  int Length = matrix.MaxNumEntries();
  vector<double> Values(Length);
  vector<int>    Indices(Length);

  for (int j = 0 ; j < matrix.NumMyRows() ; ++j) 
  {
    int NumEntries;
    int ierr = matrix.ExtractMyRowCopy(j, Length, NumEntries,
                                       &Values[0], &Indices[0]);

    for (int k = 0 ; k < NumEntries ; ++k) 
      DenseMatrix(j,Indices[k]) = Values[k];
  }

  Epetra_SerialDenseMatrix DenseX(n, NumVectors);
  Epetra_SerialDenseMatrix DenseB(n, NumVectors);

  for (int i = 0 ; i < n ; ++i)
    for (int j = 0 ; j < NumVectors ; ++j)
      DenseB(i,j) = RHS[j][i];

  Epetra_SerialDenseSolver DenseSolver;

  DenseSolver.SetMatrix(DenseMatrix);
  DenseSolver.SetVectors(DenseX,DenseB);

  DenseSolver.Factor();
  DenseSolver.Solve();

  for (int i = 0 ; i < n ; ++i)
    for (int j = 0 ; j < NumVectors ; ++j)
      LHS[j][i] = DenseX(i,j);

}

// ============================================================================ 
Epetra_MultiVector* 
Galeri::core::Workspace::createMultiVectorComponent(const Epetra_MultiVector& input)
{
  const Epetra_Comm& comm = input.Comm();
  const Epetra_BlockMap& inputMap = input.Map();
  const int* myGlobalElements = inputMap.MyGlobalElements();
  const int numMyElements = inputMap.NumMyElements();
  const int indexBase = inputMap.IndexBase();

  Epetra_Map extractMap(-1, numMyElements, myGlobalElements, indexBase, comm);
  Epetra_MultiVector* output = new Epetra_MultiVector(extractMap, input.NumVectors());

  return(output);
}

// ============================================================================ 
void 
Galeri::core::Workspace::extractMultiVectorComponent(const Epetra_MultiVector& input,
                                                     const int equation,
                                                     Epetra_MultiVector& output)
{
  const Epetra_BlockMap& inputMap = input.Map();
  const int numMyElements = inputMap.NumMyElements();

  for (int i = 0; i < numMyElements; ++i)
  {
    int j = inputMap.FirstPointInElement(i) + equation;
    for (int k = 0; k < input.NumVectors(); ++k)
      output[k][i] = input[k][j];
  }
}
