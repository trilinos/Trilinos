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

#include "LOCA_Predictor_Secant.H"
#include "LOCA_Predictor_Manager.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Utils.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"

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

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Secant::compute(
	      bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      LOCA::MultiContinuation::ExtendedMultiVector& prevXMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& xMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result)
{

  // Compute x - xold
  result[0].update(1.0, xMultiVec[0], -1.0, prevXMultiVec[0], 0.0);

  for (int i=0; i<result.numVectors(); i++) {

    result[i] = result[0];

    // Rescale so parameter component = 1
    result[i].scale(1.0/fabs(result.getScalar(i,i)));

    // Set off-diagonal elements to 0
    for (int j=0; j<result.numVectors(); j++)
      if (i != j)
	result.getScalar(i,j) = 0.0;

  }

  // Set orientation based on parameter change
  setPredictorOrientation(baseOnSecant, stepSize, grp, prevXMultiVec, 
			  xMultiVec, result);

  return NOX::Abstract::Group::Ok;
}
