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

#include "LOCA_Predictor_Manager.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Predictor_Constant.H"
#include "LOCA_Predictor_Tangent.H"
#include "LOCA_Predictor_Secant.H"
#include "LOCA_Predictor_Random.H"
#include "LOCA_Utils.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"

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
    else if (method == "Random")
      predictorPtr = new LOCA::Predictor::Random(params);
    else {
      if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
	cout << "LOCA::Predictor::Manager::reset() - invalid choice (" 
	     << method << ") for predictor method " << endl;
      }
      return NOX::Abstract::Group::Failed;
    }
  }

  return LOCA::Predictor::Generic::reset(params);
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Manager::compute(bool baseOnSecant, double stepSize,
				  LOCA::Continuation::ExtendedGroup& prevGroup,
				  LOCA::Continuation::ExtendedGroup& curGroup,
				  LOCA::Continuation::ExtendedVector& result) 
{
  if (predictorPtr == NULL) {
    if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
      cout << "LOCA::Predictor::Manager::compute - Null pointer error" << endl;
    }
    return NOX::Abstract::Group::Failed;
  }

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails))
    cout << "\n\tCalling Predictor with method: " << method << endl;

  return predictorPtr->compute(baseOnSecant, stepSize, prevGroup, curGroup, 
			       result);
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Manager::compute(bool baseOnSecant, 
	      const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      LOCA::MultiContinuation::ExtendedMultiVector& prevXMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& xMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result)
{
  if (predictorPtr == NULL) {
    if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
      cout << "LOCA::Predictor::Manager::compute - Null pointer error" << endl;
    }
    return NOX::Abstract::Group::Failed;
  }

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails))
    cout << "\n\tCalling Predictor with method: " << method << endl;

  return predictorPtr->compute(baseOnSecant, stepSize, grp, prevXMultiVec,
			       xMultiVec, result);
}

const string&
LOCA::Predictor::Manager::getMethod() const 
{
  return method;
}
