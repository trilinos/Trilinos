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
    throw "NOX Error";
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
