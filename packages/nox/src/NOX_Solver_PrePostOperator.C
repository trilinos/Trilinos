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
#include "NOX_Parameter_List.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Parameter_PrePostOperator.H"


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

NOX::Solver::PrePostOperator::PrePostOperator(NOX::Utils& utils, 
					      NOX::Parameter::List& p) :
  havePrePostOperator(false),
  prePostOperatorPtr(0)
{ 
  reset(utils, p);
}

NOX::Solver::PrePostOperator::~PrePostOperator()
{ 
  delete prePostOperatorPtr;
}

void NOX::Solver::PrePostOperator::
reset(NOX::Utils& utils, NOX::Parameter::List& p)
{
  //NOX::Parameter::List& p = params.sublist("Solver Options");

  if (prePostOperatorPtr != 0)
    delete prePostOperatorPtr;
  prePostOperatorPtr = 0;
  havePrePostOperator = false;

  if (p.isParameter("User Defined Pre/Post Operator")) {
    if (p.isParameterArbitrary("User Defined Pre/Post Operator")) {
      prePostOperatorPtr = dynamic_cast<NOX::Parameter::PrePostOperator*>
	(p.getArbitraryParameter("User Defined Pre/Post Operator").clone());
      if (prePostOperatorPtr != 0)
	havePrePostOperator = true;
      else
	if (utils.isPrintProcessAndType(NOX::Utils::Warning))
	  cout << "Warning: NOX::Solver::LineSearchBased::init() - " 
	       << "\"User Defined Pre/Post Operator\" not derived from " 
	       << "NOX::Parameter::PrePostOperator class!\n" 
	       << "Ignoring this flag!"<< endl;
    }
    else {
      cout << "ERROR: NOX::Solver::LineSearchBased::init() - the parameter "
	   << "\"User Defined Pre/Post Operator\" must be derived from an"
	   << "arbitrary parameter!" << endl;
      throw "NOX Error";
    }
  }
}

inline void NOX::Solver::PrePostOperator::
runPreIterate(const NOX::Solver::Generic& solver)
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPreIterate(solver);
}

inline void NOX::Solver::PrePostOperator::
runPostIterate(const NOX::Solver::Generic& solver)
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPostIterate(solver);
}

inline void NOX::Solver::PrePostOperator::
runPreSolve(const NOX::Solver::Generic& solver)
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPreSolve(solver);
}

inline void NOX::Solver::PrePostOperator::
runPostSolve(const NOX::Solver::Generic& solver)
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPostSolve(solver);
}
