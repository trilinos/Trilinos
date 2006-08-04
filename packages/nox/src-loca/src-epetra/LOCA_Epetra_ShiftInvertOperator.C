//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2005) Sandia Corporation
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

#include "LOCA_Epetra_ShiftInvertOperator.H"
#include "LOCA_Epetra_Group.H"
#include "NOX_Epetra_MultiVector.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "Epetra_MultiVector.h"

//=============================================================================
LOCA::Epetra::ShiftInvertOperator::ShiftInvertOperator(
                   const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		   const Teuchos::RefCountPtr<LOCA::Epetra::Group>& grp, 
		   const Teuchos::RefCountPtr<const Epetra_Operator>& jac, 
		   double shift) :
  globalData(global_data),
  locagrp(grp),
  jacOper(jac),
  shift_(shift),
  Label_(0)
{
  Label_ = "LOCA::Epetra::ShiftInvertOperator";

  // Compute shifted matrix J - shift * M
  grp->computeShiftedMatrix(1.0, -shift_);
}

LOCA::Epetra::ShiftInvertOperator::~ShiftInvertOperator()
{
}

int
LOCA::Epetra::ShiftInvertOperator::SetUseTranspose(bool UseTranspose)
{
  if (UseTranspose == false)
    return 0;
  else {
    globalData->locaErrorCheck->throwError(
	  "LOCA::Epetra::ShiftInvert::SetUseTranspose",
	  "Operator does not support transpose");
    return -1;
  }
}

int 
LOCA::Epetra::ShiftInvertOperator::Apply(const Epetra_MultiVector& X, 
					 Epetra_MultiVector& Y) const
{
  // Create NOX multivectors out of X and Y (views)
  NOX::Epetra::MultiVector nox_X(
		       Teuchos::rcp(&const_cast<Epetra_MultiVector&>(X),false),
		       NOX::DeepCopy,
		       NOX::Epetra::MultiVector::CreateView);
  NOX::Epetra::MultiVector nox_Y(Teuchos::rcp(&Y,false),
				 NOX::DeepCopy,
				 NOX::Epetra::MultiVector::CreateView);

  // Apply shifted matrix
  NOX::Abstract::Group::ReturnType result = 
    locagrp->applyShiftedMatrixMultiVector(nox_X, nox_Y);
  
  if (result == NOX::Abstract::Group::Ok)
    return 0;
  else 
    return -1;
} 

int 
LOCA::Epetra::ShiftInvertOperator::ApplyInverse(const Epetra_MultiVector& X, 
						Epetra_MultiVector&Y) const
{
  globalData->locaErrorCheck->throwError(
	  "LOCA::Epetra::ShiftInvertOperator::ApplyInverse",
	  "Operator does not support ApplyInverse");
  return -1;
}

double
LOCA::Epetra::ShiftInvertOperator::NormInf() const
{
  globalData->locaErrorCheck->throwError(
	  "LOCA::Epetra::ShiftInvertOperator::NormInf",
	  "Operator does not support NormInf");
  return -1;
}

const char*
LOCA::Epetra::ShiftInvertOperator::Label() const
{
  return(Label_);
}

bool
LOCA::Epetra::ShiftInvertOperator::UseTranspose() const
{
  return false;
}

bool
LOCA::Epetra::ShiftInvertOperator::HasNormInf() const
{
  return false;
}

const Epetra_Comm&
LOCA::Epetra::ShiftInvertOperator::Comm() const
{
  return jacOper->Comm();
}

const Epetra_Map&
LOCA::Epetra::ShiftInvertOperator::OperatorDomainMap() const
{
  return jacOper->OperatorDomainMap();
}

const Epetra_Map&
LOCA::Epetra::ShiftInvertOperator::OperatorRangeMap() const
{
  return jacOper->OperatorRangeMap();
}
