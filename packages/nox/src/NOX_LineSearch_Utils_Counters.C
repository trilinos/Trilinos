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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
