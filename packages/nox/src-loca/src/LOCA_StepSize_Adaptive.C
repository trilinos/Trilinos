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

#include "LOCA_StepSize_Adaptive.H"
#include "LOCA_Continuation_Group.H"
#include "LOCA_Stepper.H"
#include "NOX_Solver_Generic.H"

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
LOCA::StepSize::Adaptive::compute(LOCA::Continuation::Group& curGroup,
				  const LOCA::Continuation::Vector& predictor,
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
  }
  else {
    double ds_ratio = curGroup.getStepSizeScaleFactor();
    LOCA::StepSize::Constant::startStepSize *= ds_ratio;
    LOCA::StepSize::Constant::maxStepSize *= ds_ratio;
    LOCA::StepSize::Constant::minStepSize *= ds_ratio;

    // Get maximum number of nonlinear iterations from stepper parameters
    const NOX::Parameter::List& locaParams = 
      stepper.getParameterList().sublist("LOCA");
    const NOX::Parameter::List& p = locaParams.sublist("Stepper");
    double maxNonlinearSteps 
      = static_cast<double>(p.getParameter("Max Nonlinear Iterations", 15));

    // Get number of nonlinear iterations in last step
    double numNonlinearSteps = static_cast<double>(solver.getNumIterations());
  
    // A failed nonlinear solve cuts the step size in half
    if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful) {
      stepSize = prevStepSize * 0.5;    
    }
    else {

      // adapive step size control
      double factor = (maxNonlinearSteps - numNonlinearSteps) 
	/ (maxNonlinearSteps - 1.0);
      if (agrValue != 0.0) {
	stepSize = prevStepSize * (1.0 + agrValue * factor * factor);
      }
      // if constant step size (agrValue = 0.0), the step size may still be 
      // reduced by a solver failure.  We should then slowly bring the step 
      // size back towards its constant value using agrValue = 0.5.
      else{
	if (prevStepSize != startStepSize) {  
	  stepSize = prevStepSize * (1.0 + 0.5 * factor * factor);

	  if (startStepSize > 0.0)
	    stepSize = min(stepSize, startStepSize);
	  else 
	    stepSize = max(stepSize, startStepSize);
	}
      }
    } 

    stepSize *= ds_ratio;
  }

  // Clip step size to be within prescribed bounds
  NOX::Abstract::Group::ReturnType res = clipStepSize(stepSize);

  // Save stepsize
  prevStepSize = stepSize;

  return res;
}

