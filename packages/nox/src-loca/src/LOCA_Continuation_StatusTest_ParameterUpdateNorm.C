// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

LOCA::Continuation::StatusTest::ParameterUpdateNorm::ParameterUpdateNorm(
                                   double rtol_,
                                   double atol_,
                                   double tol_) :
  rtol(rtol_),
  atol(atol_),
  tol(tol_),
  paramUpdateNorm(0.0),
  status(NOX::StatusTest::Unconverged)
{
}


LOCA::Continuation::StatusTest::ParameterUpdateNorm::~ParameterUpdateNorm()
{
}

NOX::StatusTest::StatusType
LOCA::Continuation::StatusTest::ParameterUpdateNorm::checkStatus(
                     const NOX::Solver::Generic& problem)
{
  // Get solution groups from solver
  const NOX::Abstract::Group& soln = problem.getSolutionGroup();
  const NOX::Abstract::Group& oldsoln = problem.getPreviousSolutionGroup();

  // Cast soln group to continuation group
  const LOCA::Continuation::ExtendedGroup* conGroupPtr =
    dynamic_cast<const LOCA::Continuation::ExtendedGroup*>(&soln);

  // Check that group is a continuation group, return converged if not
  if (conGroupPtr == NULL) {
    paramUpdateNorm = 0.0;
    return NOX::StatusTest::Converged;
  }

  // Get solution vectors
  const LOCA::Continuation::ExtendedVector& x =
    dynamic_cast<const LOCA::Continuation::ExtendedVector&>(soln.getX());
  const LOCA::Continuation::ExtendedVector& xold =
    dynamic_cast<const LOCA::Continuation::ExtendedVector&>(oldsoln.getX());

  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0)
  {
    paramUpdateNorm = 1.0e+12;
    status = NOX::StatusTest::Unconverged;
    return status;
  }

  paramUpdateNorm =
    fabs(x.getParam() - xold.getParam()) / (rtol*fabs(xold.getParam()) + atol);

  if (paramUpdateNorm < tol)
    status = NOX::StatusTest::Converged;
  else
    status = NOX::StatusTest::Unconverged;

  return status;
}

NOX::StatusTest::StatusType
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getStatus() const
{
  return status;
}


ostream&
LOCA::Continuation::StatusTest::ParameterUpdateNorm::print(ostream& stream,
                               int indent) const
{
  for (int j = 0; j < indent; j++)
    stream << ' ';
  stream << status;
  stream << "Continuation Scaled Parameter Update = " << NOX::Utils::sciformat(paramUpdateNorm, 3) << " < " << tol;
  stream << std::endl;

  return stream;
}


double
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getUpdateNorm() const
{
  return paramUpdateNorm;
}

double
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getRTOL() const
{
  return rtol;
}

double
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getATOL() const
{
  return atol;
}

double
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getTOL() const
{
  return tol;
}

