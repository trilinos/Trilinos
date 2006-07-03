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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "NOX_Solver_PrePostOperator.H"

#include "NOX_Utils.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Solver_Generic.H"


// Disallowed
NOX::Solver::PrePostOperator::PrePostOperator()
{ }

// Disallowed
NOX::Solver::PrePostOperator::PrePostOperator(const PrePostOperator& p)
{ }

// Disallowed
NOX::Solver::PrePostOperator& NOX::Solver::PrePostOperator::
operator=(const PrePostOperator& p)
{ return *this; }

NOX::Solver::PrePostOperator::
PrePostOperator(const Teuchos::RefCountPtr<NOX::Utils>& utils,
		Teuchos::ParameterList& p) :
  havePrePostOperator(false)
{ reset(utils, p); }

NOX::Solver::PrePostOperator::~PrePostOperator()
{ }

void NOX::Solver::PrePostOperator::
reset(const Teuchos::RefCountPtr<NOX::Utils>& utils, Teuchos::ParameterList& p)
{
  havePrePostOperator = false;

  if (p.INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RefCountPtr<NOX::Abstract::PrePostOperator> >
      ("User Defined Pre/Post Operator")) 
  {
    prePostOperatorPtr = p.INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RefCountPtr<NOX::Abstract::PrePostOperator> >
      ("User Defined Pre/Post Operator");
    havePrePostOperator = true;
  }
}
