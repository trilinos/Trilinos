//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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

//==============================================================================
LOCA::Epetra::ShiftInvertOperator::ShiftInvertOperator(const LOCA::EpetraNew::Group& grp, const Epetra_Operator& jac, const double& shift, bool hasMassMatrix) :
  locagrp(grp),
  jacOper(jac),
  shift_(shift),
  massMatrix(hasMassMatrix),
  Label_(0)
{
  Label_ = "LOCA::Epetra::ShiftInvertOperator";
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
    LOCA::ErrorCheck::throwError(
	  "LOCA::Epetra::HouseholderJacOp::SetUseTranspose",
	  "Operator does not support transpose");
    return -1;
  }
}

int 
LOCA::Epetra::ShiftInvertOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector&Y) const
{
  bool applyMass = massMatrix;
  double shift = shift_;
  for (int j=0; j < X.NumVectors(); j++) {
    const Epetra_Vector xvec(Copy,X,j);
    Epetra_Vector yvec(View,Y,j);
    const NOX::Epetra::Vector Xvec(xvec); 
    NOX::Epetra::Vector Yvec(yvec);
    locagrp.applyShiftedMatrix(Xvec,Yvec,shift,applyMass);
    Epetra_Vector& result = Yvec.getEpetraVector();
    Y.Update(1.0,result,0.0);
  }
  return 0;
} 

int 
LOCA::Epetra::ShiftInvertOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector&Y) const
{
  LOCA::ErrorCheck::throwError(
	  "LOCA::Epetra::ShiftInvertOperator::ApplyInverse",
	  "Operator does not support ApplyInverse");
    return -1;
}

double
LOCA::Epetra::ShiftInvertOperator::NormInf() const
{
  LOCA::ErrorCheck::throwError(
	  "LOCA::Epetra::ShiftInvertOperator::NormInf",
	  "Operator does not support NormInf");
    return -1;
}

char*
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
  return jacOper.Comm();
}

const Epetra_Map&
LOCA::Epetra::ShiftInvertOperator::OperatorDomainMap() const
{
  return jacOper.OperatorDomainMap();
}

const Epetra_Map&
LOCA::Epetra::ShiftInvertOperator::OperatorRangeMap() const
{
  return jacOper.OperatorRangeMap();
}
