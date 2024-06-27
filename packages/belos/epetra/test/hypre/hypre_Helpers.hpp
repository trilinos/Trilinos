// @HEADER
// *****************************************************************************
//                 Didasko: Tutorial Package
//
// Copyright 2005 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EpetraExt_HYPRE_HELPERS_HPP
#define EpetraExt_HYPRE_HELPERS_HPP

#include "HYPRE_IJ_mv.h"
#include "EpetraExt_HypreIJMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

#include <string>

EpetraExt_HypreIJMatrix* newHypreMatrix(int N);

Epetra_CrsMatrix* newCrsMatrix(int N);

Epetra_CrsMatrix* GetCrsMatrix(EpetraExt_HypreIJMatrix &Matrix);

bool EquivalentVectors(const Epetra_MultiVector &X, const Epetra_MultiVector &Y, double tol);

bool EquivalentMatrices(const Epetra_RowMatrix &HypreMatrix, const Epetra_RowMatrix &CrsMatrix,double tol);

#endif // EpetraExt_HYPRE_HELPERS_HPP

