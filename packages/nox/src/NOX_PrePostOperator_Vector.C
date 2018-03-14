
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

#include "NOX_PrePostOperator_Vector.H"

void NOX::PrePostOperatorVector::runPreIterate(const NOX::Solver::Generic& solver)
{
  for (it i=ppop_vec_.begin(); i != ppop_vec_.end(); ++i)
    (*i)->runPreIterate(solver);
}

void NOX::PrePostOperatorVector::runPostIterate(const NOX::Solver::Generic& solver)
{
  for (it i=ppop_vec_.begin(); i != ppop_vec_.end(); ++i)
    (*i)->runPostIterate(solver);
}

void NOX::PrePostOperatorVector::runPreSolve(const NOX::Solver::Generic& solver)
{
  for (it i=ppop_vec_.begin(); i != ppop_vec_.end(); ++i)
    (*i)->runPreSolve(solver);
}

void NOX::PrePostOperatorVector::runPostSolve(const NOX::Solver::Generic& solver)
{
  for (it i=ppop_vec_.begin(); i != ppop_vec_.end(); ++i)
    (*i)->runPostSolve(solver);
}

void NOX::PrePostOperatorVector::runPreLineSearch(const NOX::Solver::Generic& solver)
{
  for (it i=ppop_vec_.begin(); i != ppop_vec_.end(); ++i)
    (*i)->runPreLineSearch(solver);
}

void NOX::PrePostOperatorVector::runPostLineSearch(const NOX::Solver::Generic& solver)
{
  for (it i=ppop_vec_.begin(); i != ppop_vec_.end(); ++i)
    (*i)->runPostLineSearch(solver);
}

void NOX::PrePostOperatorVector::pushBack(const Teuchos::RCP<NOX::Abstract::PrePostOperator>& ppop)
{
  ppop_vec_.push_back(ppop);
}

void NOX::PrePostOperatorVector::popBack()
{
  ppop_vec_.pop_back();
}

void NOX::PrePostOperatorVector::clear()
{
  ppop_vec_.clear();
}
