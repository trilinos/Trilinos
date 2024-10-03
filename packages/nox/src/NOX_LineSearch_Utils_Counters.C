// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_LineSearch_Utils_Counters.H"

#include "Teuchos_ParameterList.hpp"

NOX::LineSearchCounters::LineSearchCounters()
{
  reset();
}

NOX::LineSearchCounters::~LineSearchCounters()
{

}

void NOX::LineSearchCounters::reset()
{
  totalNumLineSearchCalls = 0;
  totalNumNonTrivialLineSearches = 0;
  totalNumFailedLineSearches = 0;
  totalNumIterations = 0;
}

bool NOX::LineSearchCounters::setValues(Teuchos::ParameterList& lineSearchParams)
{
  Teuchos::ParameterList& outputList = lineSearchParams.sublist("Output");
  outputList.set("Total Number of Line Search Calls", totalNumLineSearchCalls);
  outputList.set("Total Number of Non-trivial Line Searches", totalNumNonTrivialLineSearches);
  outputList.set("Total Number of Failed Line Searches", totalNumFailedLineSearches);
  outputList.set("Total Number of Line Search Inner Iterations", totalNumIterations);
  return true;
}

void NOX::LineSearchCounters::incrementNumLineSearches(int numLS)
{
  totalNumLineSearchCalls += numLS;
}

void NOX::LineSearchCounters::incrementNumNonTrivialLineSearches(int numNonTrivialLS)
{
  totalNumNonTrivialLineSearches += numNonTrivialLS;
}

void NOX::LineSearchCounters::incrementNumFailedLineSearches(int numFailedLS)
{
  totalNumFailedLineSearches += numFailedLS;
}

void NOX::LineSearchCounters::incrementNumIterations(int numInnerIters)
{
  totalNumIterations += numInnerIters;
}

int NOX::LineSearchCounters::getNumLineSearches() const
{
  return totalNumLineSearchCalls;
}

int NOX::LineSearchCounters::getNumNonTrivialLineSearches() const
{
  return totalNumNonTrivialLineSearches;
}

int NOX::LineSearchCounters::getNumFailedLineSearches() const
{
  return totalNumFailedLineSearches;
}

int NOX::LineSearchCounters::getNumIterations() const
{
  return totalNumIterations;
}
