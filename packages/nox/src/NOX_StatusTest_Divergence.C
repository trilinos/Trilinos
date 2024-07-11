// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_StatusTest_Divergence.H" // class definition
#include "NOX_Common.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Abstract_Group.H"

NOX::StatusTest::Divergence::Divergence(double threshold_, int maxSteps_) :
  maxSteps(maxSteps_),
  numSteps(0),
  lastIteration(-1),
  threshold(threshold_),
  status(NOX::StatusTest::Unevaluated)
{
}

NOX::StatusTest::Divergence::~Divergence()
{
}

NOX::StatusTest::StatusType NOX::StatusTest::Divergence::
checkStatus(const Solver::Generic& problem,
        NOX::StatusTest::CheckType /* checkType */)
{
  status = Unconverged;

  // This test should ignore the checkType!  This test must be run
  // each iteration because it triggers after a set number of
  // iterations.

  // First time through we don't do anything
  int niters = problem.getNumIterations();
  if (niters == 0) {
    lastIteration = 0;
    numSteps = 0;
    return Unconverged;
  }

  // Make sure we have not already counted the last nonlinear iteration.
  // This protects against multiple calls to checkStatus() in between
  // nonlinear iterations.
  bool isCounted = false;
  if (niters == lastIteration) {
    isCounted = true;
  }
  else
    lastIteration = niters;

  // Check the norm and see if it exceeds threshold
  if (!isCounted) {

    bool isOver = ( problem.getSolutionGroup().getNormF() > threshold );

    if ( isOver )
      numSteps ++;
    else
      numSteps = 0;

  }

  if (numSteps >= maxSteps)
    status = Failed;

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::Divergence::getStatus() const
{
  return status;
}

std::ostream& NOX::StatusTest::Divergence::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Divergence Count = " << numSteps << " < " << maxSteps << "\n";

  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << "             (max F-norm threshold = " << threshold << ")";
  stream << std::endl;
 return stream;
}


int NOX::StatusTest::Divergence::getMaxNumSteps() const
{
  return maxSteps;
}

int NOX::StatusTest::Divergence::getCurrentNumSteps() const
{
  return numSteps;
}

double NOX::StatusTest::Divergence::getThreshold() const
{
  return threshold;
}
