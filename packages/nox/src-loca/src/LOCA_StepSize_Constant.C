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

#include "LOCA_StepSize_Constant.H"
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Stepper.H"
#include "NOX_Solver_Generic.H"
#include "LOCA_Utils.H"

LOCA::StepSize::Constant::Constant(NOX::Parameter::List& params) :
  maxStepSize(0.0),
  minStepSize(0.0)
{
  reset(params);
}

LOCA::StepSize::Constant::~Constant()
{
}

NOX::Abstract::Group::ReturnType 
LOCA::StepSize::Constant::reset(NOX::Parameter::List& params) 
{
  maxStepSize = params.getParameter("Max Step Size", 1.0e+12);
  minStepSize = params.getParameter("Min Step Size", 1.0e-12);
  startStepSize = params.getParameter("Initial Step Size", 1.0);
  failedFactor = params.getParameter("Failed Step Reduction Factor", 0.5);
  successFactor = params.getParameter("Successful Step Increase Factor", 1.26);
  prevStepSize = 0.0;
  isFirstStep = true;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::StepSize::Constant::compute(
			LOCA::Continuation::ExtendedGroup& curGroup,
			const LOCA::Continuation::ExtendedVector& predictor,
			const NOX::Solver::Generic& solver,
			const LOCA::Abstract::Iterator::StepStatus& stepStatus,
			const LOCA::Stepper& stepper,
			double& stepSize) 
{

  // If this is the first step, set step size to initial value adjusted 
  // to predicted change in parameter
  if (isFirstStep) {
    double dpds = predictor.getParam();
    startStepSize /= dpds;
    maxStepSize /= dpds;
    minStepSize /= dpds;
    stepSize = startStepSize;
    isFirstStep = false;
    prevStepSize = 0.0;
  }
  else {

    // Step size remains constant, unless...
    // A failed nonlinear solve cuts the step size by failedFactor
    if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful) {
      stepSize *= failedFactor;
    }
    else {

      double ds_ratio = curGroup.getStepSizeScaleFactor();
      startStepSize *= ds_ratio;
      maxStepSize *= ds_ratio;
      minStepSize *= ds_ratio;

      prevStepSize = stepSize;
      stepSize *= ds_ratio;

      // For constant step size, the step size may still have been
      // reduced by a solver failure.  We then increase the step size
      // by a factor of cube-root-2 until back to the original step size

      if (stepSize != startStepSize) {

        stepSize *= successFactor;

        if (startStepSize > 0.0)
          stepSize = min(stepSize, startStepSize);
        else
          stepSize = max(stepSize, startStepSize);
     }
    }
  }

  // Clip step size to be within prescribed bounds
  NOX::Abstract::Group::ReturnType res = clipStepSize(stepSize);

  return res;
}

NOX::Abstract::Group::ReturnType 
LOCA::StepSize::Constant::clipStepSize(double& stepSize)
{
  NOX::Abstract::Group::ReturnType res = NOX::Abstract::Group::Ok;

  // Compute sign of step size
  double signStep = 1.0;
  if (stepSize < 0.0)
    signStep = -1.0;

  // Clip the step size if above the bounds
  if (fabs(stepSize) > maxStepSize)
     stepSize = signStep*maxStepSize;

  // Clip step size at minimum, signal for failed run
  if (fabs(stepSize) < minStepSize) {
    res = NOX::Abstract::Group::Failed;
    stepSize =  signStep*minStepSize;
    if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
      cout << "\n\tStep size reached minimum step size bound" 
	   << endl;
    }
  }

  return res;
}

double
LOCA::StepSize::Constant::getPrevStepSize() const {
  return prevStepSize;
}

double
LOCA::StepSize::Constant::getStartStepSize() const {
  return startStepSize;
}
