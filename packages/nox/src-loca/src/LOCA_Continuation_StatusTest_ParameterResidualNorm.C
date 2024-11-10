// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

LOCA::Continuation::StatusTest::ParameterResidualNorm::ParameterResidualNorm(
                                  double rtol_,
                                  double atol_,
                                  double tol_) :
  rtol(rtol_),
  atol(atol_),
  tol(tol_),
  paramResidualNorm(0.0),
  status(NOX::StatusTest::Unconverged)
{
}


LOCA::Continuation::StatusTest::ParameterResidualNorm::~ParameterResidualNorm()
{
}

NOX::StatusTest::StatusType
LOCA::Continuation::StatusTest::ParameterResidualNorm::checkStatus(
                     const NOX::Solver::Generic& problem)
{
  // Get solution groups from solver
  const NOX::Abstract::Group& soln = problem.getSolutionGroup();

  // Cast soln group to continuation group (for parameter step)
  const LOCA::Continuation::ExtendedGroup* conGroupPtr =
    dynamic_cast<const LOCA::Continuation::ExtendedGroup*>(&soln);

  // Check that group is a continuation group, return converged if not
  if (conGroupPtr == NULL) {
    paramResidualNorm = 0.0;
    return NOX::StatusTest::Converged;
  }

  // Get residual vector
  const LOCA::Continuation::ExtendedVector& f =
    dynamic_cast<const LOCA::Continuation::ExtendedVector&>(soln.getF());

  paramResidualNorm =
    fabs(f.getParam()) / (rtol*fabs(conGroupPtr->getStepSize()) + atol);

  if (paramResidualNorm < tol)
    status = NOX::StatusTest::Converged;
  else
    status = NOX::StatusTest::Unconverged;

  return status;
}

NOX::StatusTest::StatusType
LOCA::Continuation::StatusTest::ParameterResidualNorm::getStatus() const
{
  return status;
}


ostream&
LOCA::Continuation::StatusTest::ParameterResidualNorm::print(ostream& stream,
                               int indent) const
{
  for (int j = 0; j < indent; j++)
    stream << ' ';
  stream << status;
  stream << "Continuation Scaled Parameter Residual = " << NOX::Utils::sciformat(paramResidualNorm, 3) << " < " << tol;
  stream << std::endl;

  return stream;
}

double
LOCA::Continuation::StatusTest::ParameterResidualNorm::getResidualNorm() const
{
  return paramResidualNorm;
}

double
LOCA::Continuation::StatusTest::ParameterResidualNorm::getRTOL() const
{
  return rtol;
}

double
LOCA::Continuation::StatusTest::ParameterResidualNorm::getATOL() const
{
  return atol;
}

double
LOCA::Continuation::StatusTest::ParameterResidualNorm::getTOL() const
{
  return tol;
}

