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

#include "LOCA_Predictor_Secant.H"
#include "LOCA_Predictor_Manager.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Utils.H"

LOCA::Predictor::Secant::Secant(NOX::Parameter::List& params) :
  firstStepPredictorPtr(NULL)
{
  reset(params);
}

LOCA::Predictor::Secant::~Secant()
{
  delete firstStepPredictorPtr;
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Secant::reset(NOX::Parameter::List& params) 
{
  NOX::Parameter::List& firstPredictorList 
    = LOCA::Utils::getSublist("First Step Predictor");

  delete firstStepPredictorPtr;
  firstStepPredictorPtr = new LOCA::Predictor::Manager(firstPredictorList);

  isFirstStep = true;

  return LOCA::Predictor::Generic::reset(params);
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Secant::compute(bool baseOnSecant, double stepSize,
				 LOCA::Continuation::ExtendedGroup& prevGroup,
				 LOCA::Continuation::ExtendedGroup& curGroup,
				 LOCA::Continuation::ExtendedVector& result) 
{
  NOX::Abstract::Group::ReturnType res;

  if (isFirstStep) {
    res = firstStepPredictorPtr->compute(baseOnSecant, stepSize, prevGroup, 
					 curGroup, result);
    isFirstStep = false;
  }
  else {

    // Compute x - xold
    result.update(1.0, curGroup.getX(), -1.0, prevGroup.getX(), 0.0);

    // Rescale so parameter component = 1
    result.scale(1.0/fabs(result.getParam()));

    // Rescale predictor
    curGroup.scalePredictor(result);

    // Set orientation based on parameter change
    setPredictorOrientation(baseOnSecant, stepSize, prevGroup, curGroup, 
			    result);
  
    // Set predictor in continuation group
    curGroup.setPredictorDirection(result);

    res = NOX::Abstract::Group::Ok;
  }

  return res;
}
