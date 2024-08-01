// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
LOCA::Epetra::ShiftInvertOperator::ApplyInverse(const Epetra_MultiVector& /* X */,
                        Epetra_MultiVector&/* Y */) const
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
