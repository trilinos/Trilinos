// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_StatusTest_MaxIters.H" // class definition
#include "NOX_Common.H"
#include "NOX_Solver_Generic.H"

NOX::StatusTest::MaxIters::
MaxIters(int maxIterations, const NOX::Utils* u) :
  maxiters(maxIterations),
  niters(0),
  status(Unevaluated)
{
  if (u != NULL)
    utils = *u;

  if (maxiters < 1)
  {
    utils.err() << "NOX::StatusTest::MaxIters - must choose a number greater than zero" << std::endl;
    throw std::runtime_error("NOX Error");
  }
}

NOX::StatusTest::MaxIters::~MaxIters()
{
}

NOX::StatusTest::StatusType NOX::StatusTest::MaxIters::
checkStatus(const Solver::Generic& problem,
        NOX::StatusTest::CheckType checkType)
{
  switch (checkType)
  {
  case NOX::StatusTest::Complete:
  case NOX::StatusTest::Minimal:
    niters = problem.getNumIterations();
    status = (niters >= maxiters) ? Failed : Unconverged;
    break;

  case NOX::StatusTest::None:
  default:
    niters = -1;
    status = Unevaluated;
    break;
  }

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::MaxIters::getStatus() const
{
  return status;
}

std::ostream& NOX::StatusTest::MaxIters::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Number of Iterations = " << niters << " < " << maxiters;
  stream << std::endl;
 return stream;
}

int NOX::StatusTest::MaxIters::getMaxIters() const
{
  return maxiters;
}

int NOX::StatusTest::MaxIters::getNumIters() const
{
  return niters;
}
