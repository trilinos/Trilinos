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

#include "LOCA_StepSize_Manager.H"
#include "LOCA_Continuation_Group.H"
#include "LOCA_Stepper.H"
#include "NOX_Solver_Generic.H"
#include "LOCA_StepSize_Constant.H"
#include "LOCA_StepSize_Adaptive.H"

LOCA::StepSize::Manager::Manager(NOX::Parameter::List& params) :
  method(),
  stepSizePtr(NULL)
{
  reset(params);
}

LOCA::StepSize::Manager::~Manager()
{
  delete stepSizePtr;
}

NOX::Abstract::Group::ReturnType 
LOCA::StepSize::Manager::reset(NOX::Parameter::List& params) 
{
  string newmethod = params.getParameter("Method", "Constant");

  if (method != newmethod) {
    delete stepSizePtr;

    method = newmethod;

    if (method == "Constant")
      stepSizePtr = new LOCA::StepSize::Constant(params);
    else if (method == "Adaptive")
      stepSizePtr = new LOCA::StepSize::Adaptive(params);
    else {
      cerr << "LOCA::StepSize::Manager::reset() - invalid choice (" 
	   << method << ") for step size method " << endl;
      return NOX::Abstract::Group::Failed;
    }
  }

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::StepSize::Manager::compute(LOCA::Continuation::Group& curGroup,
				 const LOCA::Continuation::Vector& predictor,
				 const NOX::Solver::Generic& solver,
				 const LOCA::Abstract::Iterator::StepStatus& stepStatus,
				 const LOCA::Stepper& stepper,
				 double& stepSize) 
{
  if (stepSizePtr == NULL) {
    cerr << "LOCA::StepSize::Manager::compute - Null pointer error" << endl;
    return NOX::Abstract::Group::Failed;
  }

  return stepSizePtr->compute(curGroup, predictor, solver, stepStatus, 
			      stepper, stepSize);
}

const string&
LOCA::StepSize::Manager::getMethod() const 
{
  return method;
}

double
LOCA::StepSize::Manager::getPrevStepSize() const {
  return stepSizePtr->getPrevStepSize();
}
