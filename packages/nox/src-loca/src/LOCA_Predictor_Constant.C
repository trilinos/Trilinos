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

#include "LOCA_Predictor_Constant.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"

LOCA::Predictor::Constant::Constant(NOX::Parameter::List& params) 
{
}

LOCA::Predictor::Constant::~Constant()
{
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Constant::reset(NOX::Parameter::List& params) 
{
  return LOCA::Predictor::Generic::reset(params);
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Constant::compute(
				bool baseOnSecant, 
				double stepSize,
				LOCA::Continuation::ExtendedGroup& prevGroup,
				LOCA::Continuation::ExtendedGroup& curGroup,
				LOCA::Continuation::ExtendedVector& result) 
{
  result.init(0.0);
  result.getParam() = 1.0;

  // Set orientation based on parameter change
  setPredictorOrientation(baseOnSecant, stepSize, prevGroup, curGroup, result);

  curGroup.setPredictorDirection(result);
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::Predictor::Constant::compute(
	      bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      LOCA::MultiContinuation::ExtendedMultiVector& prevXMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& xMultiVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result)
{
  result.init(0.0);
  for (int i=0; i<result.numVectors(); i++)
    result.getScalar(i,i) = 1.0;

  // Set orientation based on parameter change
  setPredictorOrientation(baseOnSecant, stepSize, grp, prevXMultiVec, 
			  xMultiVec, result);

  return NOX::Abstract::Group::Ok;
}
