// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Continuation_AbstractGroup.H"
#include "LOCA_SingularJacobianSolve_NicDay.H"
#include "LOCA_ErrorCheck.H"

LOCA::SingularJacobianSolve::NicDay::NicDay(Teuchos::ParameterList& params)
{
  reset(params);
}

LOCA::SingularJacobianSolve::NicDay::NicDay(
              const LOCA::SingularJacobianSolve::NicDay& source)
{
}

LOCA::SingularJacobianSolve::NicDay::~NicDay()
{
}

LOCA::SingularJacobianSolve::Generic*
LOCA::SingularJacobianSolve::NicDay::clone() const
{
  return new NicDay(*this);
}

LOCA::SingularJacobianSolve::Generic&
LOCA::SingularJacobianSolve::NicDay::operator=(
              const LOCA::SingularJacobianSolve::Generic& source)
{
  return operator=(dynamic_cast<const LOCA::SingularJacobianSolve::NicDay&>(source));
}

LOCA::SingularJacobianSolve::NicDay&
LOCA::SingularJacobianSolve::NicDay::operator=(
              const LOCA::SingularJacobianSolve::NicDay& source)
{
  return *this;
}

NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::NicDay::reset(Teuchos::ParameterList& params)
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::NicDay::compute(
                Teuchos::ParameterList& params,
                LOCA::Continuation::AbstractGroup& grp,
                const NOX::Abstract::Vector& input,
                    const NOX::Abstract::Vector& approxNullVec,
                const NOX::Abstract::Vector& jacApproxNullVec,
                NOX::Abstract::Vector& result)
{
  std::string callingFunction =
    "LOCA::SingularJacobianSolve::NicDay::compute()";
  NOX::Abstract::Group::ReturnType finalStatus;

  double alpha = jacApproxNullVec.innerProduct(input)
               / jacApproxNullVec.innerProduct(jacApproxNullVec);

  NOX::Abstract::Vector* tmpInput  = input.clone(NOX::DeepCopy);
  tmpInput->update(-alpha, jacApproxNullVec, 1.0);

  finalStatus = grp.applyJacobianInverse(params, *tmpInput, result);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  delete tmpInput;

  result.update(alpha, approxNullVec, 1.0);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::NicDay::computeMulti(
                Teuchos::ParameterList& params,
                LOCA::Continuation::AbstractGroup& grp,
                const NOX::Abstract::Vector*const* inputs,
                const NOX::Abstract::Vector& approxNullVec,
                const NOX::Abstract::Vector& jacApproxNullVec,
                NOX::Abstract::Vector** results,
                int nVecs)
{
  std::string callingFunction =
    "LOCA::SingularJacobianSolve::NicDay::computeMulti()";
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  double denom = jacApproxNullVec.innerProduct(jacApproxNullVec);

  double* alphas = new double[nVecs];
  NOX::Abstract::Vector** tmpInputs  = new NOX::Abstract::Vector*[nVecs];

  for (int i=0; i<nVecs; i++) {
    alphas[i] = jacApproxNullVec.innerProduct(*(inputs[i]))/denom;
    tmpInputs[i] = inputs[i]->clone(NOX::DeepCopy);
    tmpInputs[i]->update(-alphas[i], jacApproxNullVec, 1.0);
  }

  status = grp.applyJacobianInverseMulti(params, tmpInputs, results, nVecs);
  finalStatus =
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
                         callingFunction);

  for (int i=0; i<nVecs; i++) {
    results[i]->update(alphas[i], approxNullVec, 1.0);
    delete tmpInputs[i];
  }

  delete [] tmpInputs;
  delete [] alphas;

  return finalStatus;
}
