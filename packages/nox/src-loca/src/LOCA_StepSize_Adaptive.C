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

#include "LOCA_StepSize_Adaptive.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Stepper.H"
#include "NOX_Solver_Generic.H"
#include "LOCA_Utils.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_NewStepper.H"

LOCA::StepSize::Adaptive::Adaptive(NOX::Parameter::List& params) :
  LOCA::StepSize::Constant(params),
  agrValue(0.0)
{
  reset(params);
}

LOCA::StepSize::Adaptive::~Adaptive()
{
}

NOX::Abstract::Group::ReturnType 
LOCA::StepSize::Adaptive::reset(NOX::Parameter::List& params) 
{
  NOX::Abstract::Group::ReturnType res = 
    LOCA::StepSize::Constant::reset(params);
  
  agrValue = params.getParameter("Aggressiveness", 0.0);

  return res;
}

NOX::Abstract::Group::ReturnType 
LOCA::StepSize::Adaptive::compute(
		       LOCA::Continuation::ExtendedGroup& curGroup,
		       const LOCA::Continuation::ExtendedVector& predictor,
		       const NOX::Solver::Generic& solver,
		       const LOCA::Abstract::Iterator::StepStatus& stepStatus,
		       const LOCA::Stepper& stepper,
		       double& stepSize) 
{
  // If this is the first step, set step size to initial value
  if (isFirstStep) {
    double dpds = predictor.getParam();
    LOCA::StepSize::Constant::startStepSize /= dpds;
    LOCA::StepSize::Constant::maxStepSize /= dpds;
    LOCA::StepSize::Constant::minStepSize /= dpds;
    LOCA::StepSize::Constant::isFirstStep = false;
    stepSize = LOCA::StepSize::Constant::startStepSize;
    prevStepSize = 0.0;
  }
  else {
  
    // A failed nonlinear solve cuts the step size in half
    if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful) {
      stepSize *= LOCA::StepSize::Constant::failedFactor;    
    }
    else {

      double ds_ratio = curGroup.getStepSizeScaleFactor();
      LOCA::StepSize::Constant::startStepSize *= ds_ratio;
      LOCA::StepSize::Constant::maxStepSize *= ds_ratio;
      LOCA::StepSize::Constant::minStepSize *= ds_ratio;

      // Get maximum number of nonlinear iterations from stepper parameters
      const NOX::Parameter::List& p = LOCA::Utils::getSublist("Stepper");
      double maxNonlinearSteps 
	= static_cast<double>(p.getParameter("Max Nonlinear Iterations", 15));
      
      // Get number of nonlinear iterations in last step
      double numNonlinearSteps = 
	static_cast<double>(solver.getNumIterations());

      // Save successful stepsize as previous
      prevStepSize = stepSize;

      // adapive step size control
      double factor = (maxNonlinearSteps - numNonlinearSteps) 
               	      / (maxNonlinearSteps);

      stepSize *= (1.0 + agrValue * factor * factor);

      stepSize *= ds_ratio;
    } 
  }

  // Clip step size to be within prescribed bounds
  NOX::Abstract::Group::ReturnType res = clipStepSize(stepSize);

  return res;
}

NOX::Abstract::Group::ReturnType 
LOCA::StepSize::Adaptive::compute(
		     LOCA::MultiContinuation::ExtendedGroup& curGroup,
		     const LOCA::MultiContinuation::ExtendedVector& predictor,
		     const NOX::Solver::Generic& solver,
		     const LOCA::Abstract::Iterator::StepStatus& stepStatus,
		     const LOCA::NewStepper& stepper,
		     double& stepSize) 
{
  // If this is the first step, set step size to initial value
  if (isFirstStep) {
    double dpds = predictor.getScalar(0);
    LOCA::StepSize::Constant::startStepSize /= dpds;
    LOCA::StepSize::Constant::maxStepSize /= dpds;
    LOCA::StepSize::Constant::minStepSize /= dpds;
    LOCA::StepSize::Constant::isFirstStep = false;
    stepSize = LOCA::StepSize::Constant::startStepSize;
    prevStepSize = 0.0;
  }
  else {
  
    // A failed nonlinear solve cuts the step size in half
    if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful) {
      stepSize *= LOCA::StepSize::Constant::failedFactor;    
    }
    else {

      double ds_ratio = curGroup.getStepSizeScaleFactor();
      LOCA::StepSize::Constant::startStepSize *= ds_ratio;
      LOCA::StepSize::Constant::maxStepSize *= ds_ratio;
      LOCA::StepSize::Constant::minStepSize *= ds_ratio;

      // Get maximum number of nonlinear iterations from stepper parameters
      const NOX::Parameter::List& p = LOCA::Utils::getSublist("Stepper");
      double maxNonlinearSteps 
	= static_cast<double>(p.getParameter("Max Nonlinear Iterations", 15));
      
      // Get number of nonlinear iterations in last step
      double numNonlinearSteps = 
	static_cast<double>(solver.getNumIterations());

      // Save successful stepsize as previous
      prevStepSize = stepSize;

      // adapive step size control
      double factor = (maxNonlinearSteps - numNonlinearSteps) 
               	      / (maxNonlinearSteps);

      stepSize *= (1.0 + agrValue * factor * factor);

      stepSize *= ds_ratio;
    } 
  }

  // Clip step size to be within prescribed bounds
  NOX::Abstract::Group::ReturnType res = clipStepSize(stepSize);

  return res;
}

