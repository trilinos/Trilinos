//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_Epetra_ShiftInvertOperator.H"
#include "LOCA_Epetra_Group.H"
#include "NOX_Epetra_MultiVector.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "Epetra_MultiVector.h"

//=============================================================================
LOCA::Epetra::ShiftInvertOperator::ShiftInvertOperator(
                   const Teuchos::RCP<LOCA::GlobalData>& global_data,
           const Teuchos::RCP<LOCA::Epetra::Group>& grp,
           const Teuchos::RCP<const Epetra_Operator>& jac,
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
