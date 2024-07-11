// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_StatusTest_Stagnation.H" // class definition
#include "NOX_Common.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Abstract_Group.H"

NOX::StatusTest::Stagnation::Stagnation(int maxSteps_, double tolerance_) :
  maxSteps(maxSteps_),
  numSteps(0),
  lastIteration(-1),
  tolerance(tolerance_),
  convRate(1.0),
  status(NOX::StatusTest::Unevaluated)
{

}

NOX::StatusTest::Stagnation::~Stagnation()
{
}

NOX::StatusTest::StatusType
NOX::StatusTest::Stagnation::
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

  // Compute the convergence rate and set counter appropriately
  if (!isCounted) {

    convRate = problem.getSolutionGroup().getNormF() /
               problem.getPreviousSolutionGroup().getNormF();

    if (convRate >= tolerance)
      numSteps ++;
    else
      numSteps = 0;

  }

  if (numSteps >= maxSteps)
    status = Failed;

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::Stagnation::getStatus() const
{
  return status;
}

std::ostream& NOX::StatusTest::Stagnation::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Stagnation Count = " << numSteps << " < " << maxSteps << "\n";

  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << "             (convergence rate = " << convRate << ")";
  stream << std::endl;
 return stream;
}


int NOX::StatusTest::Stagnation::getMaxNumSteps() const
{
  return maxSteps;
}

int NOX::StatusTest::Stagnation::getCurrentNumSteps() const
{
  return numSteps;
}

double NOX::StatusTest::Stagnation::getTolerance() const
{
  return tolerance;
}

double NOX::StatusTest::Stagnation::getConvRate() const
{
  return convRate;
}
