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

#include "LOCA_Predictor_Restart.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_ErrorCheck.H"

LOCA::Predictor::Restart::Restart(NOX::Parameter::List& params) 
  : v(NULL)
{
  reset(params);
}

LOCA::Predictor::Restart::~Restart()
{
  if (v != NULL)
    delete v;
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Restart::reset(NOX::Parameter::List& params) 
{
  const char *func = "LOCA::Predictor::Restart::reset()";

  if (!params.isParameter("Solution Component"))
    LOCA::ErrorCheck::throwError(func, "\"Solution Component\" is not set!");
  const NOX::Abstract::Vector* v_x = 
    params.getAnyPtrParameter<NOX::Abstract::Vector>("Solution Component");

  if (!params.isParameter("Parameter Component"))
    LOCA::ErrorCheck::throwError(func, "\"Parameter Component\" is not set!");
  double v_p = params.getParameter("Parameter Component", 0.0);

  if (v != NULL)
    delete v;

  v = new LOCA::Continuation::ExtendedVector(*v_x, v_p);

  return LOCA::Predictor::Generic::reset(params);  
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Restart::compute(
				bool baseOnSecant, 
				double stepSize,
				LOCA::Continuation::ExtendedGroup& prevGroup,
				LOCA::Continuation::ExtendedGroup& curGroup,
				LOCA::Continuation::ExtendedVector& result) 
{
  // Set predictor based on stored predictor direction
  result = *v;

  // Set orientation based on parameter change
  setPredictorOrientation(baseOnSecant, stepSize, prevGroup, curGroup, result);

  curGroup.setPredictorDirection(result);
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Restart::compute(
	      bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      LOCA::MultiContinuation::ExtendedMultiVector& prevXMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& xMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result)
{
  LOCA::ErrorCheck::throwError("LOCA::Predictor::Restart::compute()", 
			       "\"Multivector version not implemented!");
  return NOX::Abstract::Group::Failed;
}
