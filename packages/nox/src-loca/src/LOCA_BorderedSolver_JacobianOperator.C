// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_BorderedSolver_JacobianOperator.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"

LOCA::BorderedSolver::JacobianOperator::
JacobianOperator(const Teuchos::RCP<const NOX::Abstract::Group>& grp) :
  grpPtr(grp)
{
}

LOCA::BorderedSolver::JacobianOperator::
~JacobianOperator()
{
}

Teuchos::RCP<const NOX::Abstract::Group>
LOCA::BorderedSolver::JacobianOperator::
getGroup() const
{
  return grpPtr;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::JacobianOperator::
apply(const NOX::Abstract::MultiVector& X,
        NOX::Abstract::MultiVector& Y) const
{
  return grpPtr->applyJacobianMultiVector(X, Y);
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::JacobianOperator::
applyTranspose(const NOX::Abstract::MultiVector& X,
           NOX::Abstract::MultiVector& Y) const
{
  return grpPtr->applyJacobianTransposeMultiVector(X, Y);
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::JacobianOperator::
applyInverse(Teuchos::ParameterList& params,
         const NOX::Abstract::MultiVector& B,
         NOX::Abstract::MultiVector& X) const
{
  return grpPtr->applyJacobianInverseMultiVector(params, B, X);
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::JacobianOperator::
applyInverseTranspose(Teuchos::ParameterList& params,
              const NOX::Abstract::MultiVector& B,
              NOX::Abstract::MultiVector& X) const
{
  Teuchos::RCP<const LOCA::Abstract::TransposeSolveGroup> tsgrp =
    Teuchos::rcp_dynamic_cast<const LOCA::Abstract::TransposeSolveGroup>(grpPtr);
  if (tsgrp != Teuchos::null)
    return tsgrp->applyJacobianTransposeInverseMultiVector(params, B, X);
  else
    return NOX::Abstract::Group::NotDefined;
}
