// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_StepSize_Adaptive.H"
#include "NOX_Solver_Generic.H"
#include "LOCA_MultiContinuation_AbstractStrategy.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
//#include "LOCA_Stepper.H"
#include "LOCA_Abstract_Iterator.H"
#include "LOCA_Parameter_SublistParser.H"

LOCA::StepSize::Adaptive::Adaptive(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams) :
  LOCA::StepSize::Constant(global_data, topParams, stepsizeParams),
  agrValue(0.0),
  maxNonlinearSteps(0.0)
{
  agrValue = stepsizeParams->get("Aggressiveness", 0.5);

  // Get maximum number of nonlinear iterations from stepper parameters
  Teuchos::RCP<Teuchos::ParameterList> p =
    topParams->getSublist("Stepper");
  maxNonlinearSteps =
    static_cast<double>(p->get("Max Nonlinear Iterations", 15));
}

LOCA::StepSize::Adaptive::~Adaptive()
{
}

NOX::Abstract::Group::ReturnType
LOCA::StepSize::Adaptive::computeStepSize(
             LOCA::MultiContinuation::AbstractStrategy& curGroup,
             const LOCA::MultiContinuation::ExtendedVector& predictor,
             const NOX::Solver::Generic& solver,
             const LOCA::Abstract::Iterator::StepStatus& stepStatus,
//             const LOCA::Stepper& stepper,
             const LOCA::Abstract::Iterator& /* stepper */,
             double& stepSize)
{
  // If this is the first step, set step size to initial value
  if (isFirstStep) {
    double dpds = predictor.getScalar(0);
    if (dpds != 0.0) {
      LOCA::StepSize::Constant::startStepSize /= dpds;
      LOCA::StepSize::Constant::maxStepSize /= dpds;
      LOCA::StepSize::Constant::minStepSize /= dpds;
    }
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

