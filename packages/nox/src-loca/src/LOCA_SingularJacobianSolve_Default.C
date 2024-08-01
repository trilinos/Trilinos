// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Continuation_AbstractGroup.H"
#include "LOCA_SingularJacobianSolve_Default.H"

LOCA::SingularJacobianSolve::Default::Default(Teuchos::ParameterList& params)
{
  reset(params);
}

LOCA::SingularJacobianSolve::Default::Default(
              const LOCA::SingularJacobianSolve::Default& source)
{
}

LOCA::SingularJacobianSolve::Default::~Default()
{
}

LOCA::SingularJacobianSolve::Generic*
LOCA::SingularJacobianSolve::Default::clone() const
{
  return new Default(*this);
}

LOCA::SingularJacobianSolve::Generic&
LOCA::SingularJacobianSolve::Default::operator=(
              const LOCA::SingularJacobianSolve::Generic& source)
{
  return operator=(dynamic_cast<const LOCA::SingularJacobianSolve::Default&>(source));
}

LOCA::SingularJacobianSolve::Default&
LOCA::SingularJacobianSolve::Default::operator=(
              const LOCA::SingularJacobianSolve::Default& source)
{
  return *this;
}

NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Default::reset(Teuchos::ParameterList& params)
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Default::compute(
                Teuchos::ParameterList& params,
                LOCA::Continuation::AbstractGroup& grp,
                const NOX::Abstract::Vector& input,
                    const NOX::Abstract::Vector& approxNullVec,
                const NOX::Abstract::Vector& jacApproxNullVec,
                NOX::Abstract::Vector& result)
{
  return grp.applyJacobianInverse(params, input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Default::computeMulti(
                Teuchos::ParameterList& params,
                LOCA::Continuation::AbstractGroup& grp,
                const NOX::Abstract::Vector*const* inputs,
                const NOX::Abstract::Vector& approxNullVec,
                const NOX::Abstract::Vector& jacApproxNullVec,
                NOX::Abstract::Vector** results,
                int nVecs)
{
  return grp.applyJacobianInverseMulti(params, inputs, results, nVecs);
}
