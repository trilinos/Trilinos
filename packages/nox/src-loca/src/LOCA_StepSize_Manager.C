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

#include "LOCA_StepSize_Manager.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Stepper.H"
#include "NOX_Solver_Generic.H"
#include "LOCA_StepSize_Constant.H"
#include "LOCA_StepSize_Adaptive.H"
#include "LOCA_Utils.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_NewStepper.H"

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
      if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
	cout << "LOCA::StepSize::Manager::reset() - invalid choice (" 
	     << method << ") for step size method " << endl;
      }
      return NOX::Abstract::Group::Failed;
    }
  }

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::StepSize::Manager::compute(
			LOCA::Continuation::ExtendedGroup& curGroup,
			const LOCA::Continuation::ExtendedVector& predictor,
			const NOX::Solver::Generic& solver,
			const LOCA::Abstract::Iterator::StepStatus& stepStatus,
			const LOCA::Stepper& stepper,
			double& stepSize) 
{
  if (stepSizePtr == NULL) {
    if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
      cout << "LOCA::StepSize::Manager::compute - Null pointer error" << endl;
    }
    return NOX::Abstract::Group::Failed;
  }

  return stepSizePtr->compute(curGroup, predictor, solver, stepStatus, 
			      stepper, stepSize);
}

NOX::Abstract::Group::ReturnType 
LOCA::StepSize::Manager::compute(
		      LOCA::MultiContinuation::ExtendedGroup& curGroup,
		      const LOCA::MultiContinuation::ExtendedVector& predictor,
		      const NOX::Solver::Generic& solver,
		      const LOCA::Abstract::Iterator::StepStatus& stepStatus,
		      const LOCA::NewStepper& stepper,
		      double& stepSize) 
{
  if (stepSizePtr == NULL) {
    if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
      cout << "LOCA::StepSize::Manager::compute - Null pointer error" << endl;
    }
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

double
LOCA::StepSize::Manager::getStartStepSize() const {
  return stepSizePtr->getStartStepSize();
}
