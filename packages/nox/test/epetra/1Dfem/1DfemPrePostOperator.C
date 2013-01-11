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

#include "1DfemPrePostOperator.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Solver_Generic.H"


UserPrePostOperator::UserPrePostOperator(const NOX::Utils& u) :
  numRunPreIterate(0),
  numRunPostIterate(0),
  numRunPreSolve(0),
  numRunPostSolve(0)
{ 
  utils = u;
}

UserPrePostOperator::~UserPrePostOperator()
{ 

}

void UserPrePostOperator::
runPreIterate(const NOX::Solver::Generic& solver)
{
  ++numRunPreIterate;
  utils.out(NOX::Utils::Details) << 
    "1Dfem's runPreIterate() routine called!" << std::endl;
}

void UserPrePostOperator::
runPostIterate(const NOX::Solver::Generic& solver)
{
  ++numRunPostIterate;
  utils.out(NOX::Utils::Details) 
    << "1Dfem's runPostIterate() routine called!" << std::endl;
}

void UserPrePostOperator::
runPreSolve(const NOX::Solver::Generic& solver)
{
  ++numRunPreSolve;
  utils.out(NOX::Utils::Details) 
    << "1Dfem's runPreSolve() routine called!" << std::endl;
}

void UserPrePostOperator::
runPostSolve(const NOX::Solver::Generic& solver)
{
  ++numRunPostSolve;
  utils.out(NOX::Utils::Details) 
    << "1Dfem's runPostSolve() routine called!" << std::endl;
}
