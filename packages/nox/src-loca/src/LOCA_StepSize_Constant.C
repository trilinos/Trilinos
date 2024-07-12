// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_StepSize_Constant.H"
#include "NOX_Solver_Generic.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"
#include "LOCA_MultiContinuation_AbstractStrategy.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
//#include "LOCA_Stepper.H"
#include "LOCA_Abstract_Iterator.H"

LOCA::StepSize::Constant::Constant(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams) :
  globalData(global_data),
  maxStepSize(1.0e+12),
  minStepSize(1.0e-12),
  startStepSize(1.0),
  failedFactor(0.5),
  successFactor(1.26),
  prevStepSize(0.0),
  isFirstStep(true)
{
  maxStepSize = stepsizeParams->get("Max Step Size", 1.0e+12);
  minStepSize = stepsizeParams->get("Min Step Size", 1.0e-12);
  startStepSize = stepsizeParams->get("Initial Step Size", 1.0);
  failedFactor =
    stepsizeParams->get("Failed Step Reduction Factor", 0.5);
  successFactor =
    stepsizeParams->get("Successful Step Increase Factor", 1.26);
}

LOCA::StepSize::Constant::~Constant()
{
}

NOX::Abstract::Group::ReturnType
LOCA::StepSize::Constant::computeStepSize(
             LOCA::MultiContinuation::AbstractStrategy& curGroup,
             const LOCA::MultiContinuation::ExtendedVector& predictor,
             const NOX::Solver::Generic& /* solver */,
             const LOCA::Abstract::Iterator::StepStatus& stepStatus,
//             const LOCA::Stepper& stepper,
             const LOCA::Abstract::Iterator& /* stepper */,
             double& stepSize)
{

  // If this is the first step, set step size to initial value adjusted
  // to predicted change in parameter
  if (isFirstStep) {
    double dpds = predictor.getScalar(0);
    if (dpds != 0.0) {
      startStepSize /= dpds;
      maxStepSize /= dpds;
      minStepSize /= dpds;
    }
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
          stepSize = NOX_MIN(stepSize, startStepSize);
        else
          stepSize = NOX_MAX(stepSize, startStepSize);
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
    if (globalData->locaUtils->isPrintType(NOX::Utils::Error)) {
      globalData->locaUtils->err() <<
    "\n\tStep size reached minimum step size bound" << std::endl;
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
