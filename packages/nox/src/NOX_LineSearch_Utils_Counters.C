// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_LineSearch_Utils_Counters.H"

#include "NOX_Parameter_List.H"

NOX::LineSearch::Utils::Counters::Counters()
{
  reset();
}

NOX::LineSearch::Utils::Counters::~Counters()
{

}

void NOX::LineSearch::Utils::Counters::reset()
{
  totalNumLineSearchCalls = 0;
  totalNumNonTrivialLineSearches = 0;
  totalNumFailedLineSearches = 0;
  totalNumIterations = 0;
}

bool NOX::LineSearch::Utils::Counters::setValues(NOX::Parameter::List& lineSearchParams) 
{
  NOX::Parameter::List& outputList = lineSearchParams.sublist("Output");
  outputList.setParameter("Total Number of Line Search Calls", totalNumLineSearchCalls);
  outputList.setParameter("Total Number of Non-trivial Line Searches", totalNumNonTrivialLineSearches);
  outputList.setParameter("Total Number of Failed Line Searches", totalNumFailedLineSearches);
  outputList.setParameter("Total Number of Line Search Inner Iterations", totalNumIterations);
  return true;
}

void NOX::LineSearch::Utils::Counters::increaseNumLineSearches(int numLS)
{
  totalNumLineSearchCalls += numLS;
}

void NOX::LineSearch::Utils::Counters::increaseNumNonTrivialLineSearches(int numNonTrivialLS)
{
  totalNumNonTrivialLineSearches += numNonTrivialLS;
}
 
void NOX::LineSearch::Utils::Counters::increaseNumFailedLineSearches(int numFailedLS)
{
  totalNumFailedLineSearches += numFailedLS;
}

void NOX::LineSearch::Utils::Counters::increaseNumIterations(int numInnerIters)
{
  totalNumIterations += numInnerIters;
}

int NOX::LineSearch::Utils::Counters::getNumLineSearches() const
{
  return totalNumLineSearchCalls;
}
 
int NOX::LineSearch::Utils::Counters::getNumNonTrivialLineSearches() const
{
  return totalNumNonTrivialLineSearches;
}

int NOX::LineSearch::Utils::Counters::getNumFailedLineSearches() const
{
  return totalNumFailedLineSearches;
}

int NOX::LineSearch::Utils::Counters::getNumIterations() const
{
  return totalNumIterations;
}
