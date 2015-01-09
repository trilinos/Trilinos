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

#include "NOX_Solver_PrePostOperator.H"

#include "NOX_Utils.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Solver_Generic.H"


// Disallowed
NOX::Solver::PrePostOperator::PrePostOperator():
  havePrePostOperator(false)
{ }

// Disallowed
NOX::Solver::PrePostOperator::PrePostOperator(const PrePostOperator& p):
  havePrePostOperator(false)
{ }

// Disallowed
NOX::Solver::PrePostOperator& NOX::Solver::PrePostOperator::
operator=(const PrePostOperator& p)
{ return *this; }

NOX::Solver::PrePostOperator::
PrePostOperator(const Teuchos::RCP<NOX::Utils>& utils,
        Teuchos::ParameterList& p) :
  havePrePostOperator(false)
{ reset(utils, p); }

NOX::Solver::PrePostOperator::~PrePostOperator()
{ }

void NOX::Solver::PrePostOperator::
reset(const Teuchos::RCP<NOX::Utils>& utils, Teuchos::ParameterList& p)
{
  havePrePostOperator = false;

  if (p.INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RCP<NOX::Abstract::PrePostOperator> >
      ("User Defined Pre/Post Operator"))
  {
    prePostOperatorPtr = p.INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<NOX::Abstract::PrePostOperator> >
      ("User Defined Pre/Post Operator");
    havePrePostOperator = true;
  }
}
