/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include "Ifpack_Jacobi.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

//==============================================================================
Ifpack_Jacobi::Ifpack_Jacobi(const Ifpack_OverlapGraph * OverlapGraph, bool UseReciprocal, 
			       int NumSteps) 
  : Epetra_Object("Ifpack::Jacobi"),
    Epetra_CompObject(),
    Ifpack_OverlapFactorObject(OverlapGraph),
    Ifpack_OverlapSolveObject(Epetra_Object::Label(), OverlapGraph->OverlapGraph().Comm()),
    UseReciprocal_(UseReciprocal),
    NumSteps_(NumSteps),
    DiagValues_(0)
{
}
//==============================================================================
Ifpack_Jacobi::Ifpack_Jacobi(const Epetra_RowMatrix * UserMatrix, bool UseReciprocal, 
			       int NumSteps) 
  : Epetra_Object("Ifpack::Jacobi"),
    Epetra_CompObject(),
    Ifpack_OverlapFactorObject(UserMatrix),
    Ifpack_OverlapSolveObject(Epetra_Object::Label(), UserMatrix->Comm()),
    UseReciprocal_(UseReciprocal),
    NumSteps_(NumSteps),
    DiagValues_(0)
{
}
//==============================================================================
Ifpack_Jacobi::Ifpack_Jacobi(const Ifpack_Jacobi & Source) 
  : Epetra_Object(Source),
    Epetra_CompObject(Source),
    Ifpack_OverlapFactorObject(Source),
    Ifpack_OverlapSolveObject(Source),
    UseReciprocal_(Source.UseReciprocal_),
    NumSteps_(Source.NumSteps_),
    DiagValues_(Source.DiagValues_)
{
  if (DiagValues_!=0) DiagValues_ = new Epetra_Vector(*DiagValues_);
}
//==============================================================================
int Ifpack_Jacobi::ProcessOverlapMatrix(const Epetra_RowMatrix &A)
{

  if (DiagValues_=0) DiagValues_ = new Epetra_Vector(A.RowMap()); // Allocate if necessary
  EPETRA_CHK_ERR(A.ExtraMyDiagonalCopy(*DiagValues_)); // Get Diagonal Values of A

  // Compute inverse of diagonal if specified
  if (UseReciprocal()) {EPETRA_CHK_ERR(DiagValues().Reciprocal(DiagValues()));}

  if (NumSteps!=1) EPETRA_CHK_ERR(-1); // NumSteps != 1 not supported yet.
  return(0);
}
//==============================================================================
int Ifpack_Jacobi::DerivedFactor()
{

  return(0);
}
