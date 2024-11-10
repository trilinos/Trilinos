// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Abstract_Group.H"

#include "NOX_Abstract_MultiVector.H"
#include "Teuchos_ParameterList.hpp"

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::computeJacobian()
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::computeGradient()
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::computeNewton(Teuchos::ParameterList& /* params */)
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobian(const NOX::Abstract::Vector& /* input */,
                    NOX::Abstract::Vector& /* result */) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianTranspose(const NOX::Abstract::Vector& /* input */,
                         NOX::Abstract::Vector& /* result */) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianInverse(Teuchos::ParameterList& /* params */,
                       const NOX::Abstract::Vector& /* input */,
                       NOX::Abstract::Vector& /* result */) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyRightPreconditioning(bool /* useTranspose */,
                        Teuchos::ParameterList& /* params */,
                        const NOX::Abstract::Vector& /* input */,
                        NOX::Abstract::Vector& /* result */
                        ) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianMultiVector(
                    const NOX::Abstract::MultiVector& input,
                    NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobian(input[i], result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus = Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianTransposeMultiVector(
                    const NOX::Abstract::MultiVector& input,
                    NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobianTranspose(input[i], result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus = Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianInverseMultiVector(
                                    Teuchos::ParameterList& params,
                    const NOX::Abstract::MultiVector& input,
                    NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobianInverse(params, input[i], result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus = Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyRightPreconditioningMultiVector(
                   bool useTranspose,
                   Teuchos::ParameterList& params,
                   const NOX::Abstract::MultiVector& input,
                   NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<input.numVectors(); i++) {
    status = applyRightPreconditioning(useTranspose, params, input[i],
                       result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus = Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

bool NOX::Abstract::Group::isJacobian() const
{
  return false;
}

bool NOX::Abstract::Group::isGradient() const
{
  return false;
}

bool NOX::Abstract::Group::isNewton() const
{
  return false;
}

void NOX::Abstract::Group::logLastLinearSolveStats(NOX::SolverStats& ) const {}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::getNormLastLinearSolveResidual(double& /* residual */) const
{
  return NOX::Abstract::Group::NotDefined;
}

Teuchos::RCP<const NOX::Abstract::Group>
NOX::Abstract::Group::getNestedGroup() const
{
  return Teuchos::null;
}

Teuchos::RCP<NOX::Abstract::Group>
NOX::Abstract::Group::getNestedGroup()
{
  return Teuchos::null;
}
