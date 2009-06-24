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

#ifndef EpetraExt_HYPRE_HELPERS_HPP
#define EpetraExt_HYPRE_HELPERS_HPP

#include "HYPRE_IJ_mv.h"
#include "EpetraExt_HypreIJMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

#include <string>

EpetraExt_HypreIJMatrix::EpetraExt_HypreIJMatrix MatrixConstructor(int N, int type);

Epetra_CrsMatrix::Epetra_CrsMatrix GetCrsMatrix(EpetraExt_HypreIJMatrix &Matrix);

bool EquivalentVectors(Epetra_MultiVector &X, Epetra_MultiVector &Y, double tol);

bool EquivalentMatrices(Epetra_RowMatrix &HypreMatrix, Epetra_RowMatrix &CrsMatrix,double tol);

#endif // EpetraExt_HYPRE_HELPERS_HPP

