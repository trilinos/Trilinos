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

#include "LOCA_Predictor_Manager.H"
#include "LOCA_Continuation_Group.H"
#include "LOCA_Predictor_Constant.H"
#include "LOCA_Predictor_Tangent.H"
#include "LOCA_Predictor_Secant.H"

LOCA::Predictor::Manager::Manager(NOX::Parameter::List& params) :
  method(),
  predictorPtr(NULL)
{
  reset(params);
}

LOCA::Predictor::Manager::~Manager()
{
  delete predictorPtr;
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Manager::reset(NOX::Parameter::List& params) 
{
  string newmethod = params.getParameter("Method", "Constant");

  if (method != newmethod) {
    delete predictorPtr;

    method = newmethod;

    if (method == "Constant")
      predictorPtr = new LOCA::Predictor::Constant(params);
    else if (method == "Tangent")
      predictorPtr = new LOCA::Predictor::Tangent(params);
    else if (method == "Secant")
      predictorPtr = new LOCA::Predictor::Secant(params);
    else {
      cerr << "LOCA::Predictor::Manager::reset() - invalid choice (" 
	   << method << ") for predictor method " << endl;
      return NOX::Abstract::Group::Failed;
    }
  }

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Manager::compute(LOCA::Continuation::Group& prevGroup,
				  LOCA::Continuation::Group& curGroup,
				  LOCA::Continuation::Vector& result) 
{
  if (predictorPtr == NULL) {
    cerr << "LOCA::Predictor::Manager::compute - Null pointer error" << endl;
    return NOX::Abstract::Group::Failed;
  }

  return predictorPtr->compute(prevGroup, curGroup, result);
}

const string&
LOCA::Predictor::Manager::getMethod() const 
{
  return method;
}
