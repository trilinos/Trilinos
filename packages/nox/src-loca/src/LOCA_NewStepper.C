// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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
#include "LOCA_NewStepper.H"    // class definition
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "LOCA_StatusTest_Wrapper.H"

// LOCA Includes
#include "NOX_Utils.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_MultiPredictor_AbstractStrategy.H"
#include "LOCA_MultiContinuation_AbstractStrategy.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_StepSize_AbstractStrategy.H"
#include "LOCA_Eigensolver_AbstractStrategy.H"
#include "LOCA_SaveEigenData_AbstractStrategy.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "LOCA_MultiContinuation_ConstraintInterface.H"

LOCA::NewStepper::NewStepper(
                     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		     const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
		     const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t,
		     const Teuchos::RefCountPtr<Teuchos::ParameterList>& p) :
  LOCA::Abstract::Iterator(),
  globalData(),
  parsedParams(),
  predictor(),
  curGroupPtr(),
  prevGroupPtr(),
  eigensolver(),
  saveEigenData(),
  bifGroupPtr(),
  statusTestPtr(),
  paramListPtr(),
  stepperList(),
  solverPtr(),
  curPredictorPtr(),
  prevPredictorPtr(),
  stepSizeStrategyPtr(),
  conParamName(),
  conParamIDs(1),
  startValue(0.0),
  maxValue(0.0),
  minValue(0.0),
  stepSize(0.0),
  maxNonlinearSteps(15),
  targetValue(0.0),
  isTargetStep(false),
  doTangentFactorScaling(false),
  tangentFactor(1.0),
  minTangentFactor(0.1),
  tangentFactorExponent(1.0),
  calcEigenvalues(false)

{
  reset(global_data, initialGuess, t, p);
}

LOCA::NewStepper::~NewStepper()
{
}

bool
LOCA::NewStepper::reset(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
		    const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t,
		    const Teuchos::RefCountPtr<Teuchos::ParameterList>& p)
{
  globalData = global_data;
  paramListPtr = p;
  statusTestPtr = t;

  // Parse parameter list
  parsedParams = Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));
  parsedParams->parseSublists(paramListPtr);

  // Get stepper sublist
  stepperList = parsedParams->getSublist("Stepper");

  // Reset base class
  LOCA::Abstract::Iterator::resetIterator(*stepperList);

  // Create predictor strategy
  Teuchos::RefCountPtr<Teuchos::ParameterList> predictorParams = 
    parsedParams->getSublist("Predictor");
  predictor = globalData->locaFactory->createPredictorStrategy(
							      parsedParams,
							      predictorParams);

  // Create eigensolver
  Teuchos::RefCountPtr<Teuchos::ParameterList> eigenParams = 
    parsedParams->getSublist("Eigensolver");
  eigensolver = globalData->locaFactory->createEigensolverStrategy(
								parsedParams,
								eigenParams);

  // Create strategy to save eigenvectors/values
  saveEigenData = globalData->locaFactory->createSaveEigenDataStrategy(
								parsedParams,
								eigenParams);

  // Create step size strategy
  Teuchos::RefCountPtr<Teuchos::ParameterList> stepsizeParams = 
    parsedParams->getSublist("Step Size");
   stepSizeStrategyPtr = globalData->locaFactory->createStepSizeStrategy(
							     parsedParams,
							     stepsizeParams);

  // Get the continuation parameter starting value
  if (stepperList->isParameter("Initial Value"))
    startValue = stepperList->get("Initial Value", 0.0);
  else {
    globalData->locaErrorCheck->throwError(
		   "LOCA::Stepper::reset()",
		   "\"Initial Value\" of continuation parameter is not set!");
  }

  // Get the continuation parameter name
  if (stepperList->isParameter("Continuation Parameter")) {
    conParamName = stepperList->get("Continuation Parameter","None");
    initialGuess->setParam(conParamName,startValue);
  }
  else {
     globalData->locaErrorCheck->throwError(
			      "LOCA::Stepper::reset()",
			      "\"Continuation Parameter\" name is not set!");
  }

  // Get the continuation parameter index
  const LOCA::ParameterVector& pv = initialGuess->getParams();
  conParamIDs[0] = pv.getIndex(conParamName);

  // Get the max and min values of the continuation parameter
  if (stepperList->isParameter("Max Value"))
    maxValue = stepperList->get("Max Value", 0.0);
  else {
     globalData->locaErrorCheck->throwError(
		   "LOCA::Stepper::reset()",
		   "\"Maximum Value\" of continuation parameter is not set!");
  }
  if (stepperList->isParameter("Min Value"))
    minValue = stepperList->get("Min Value", 0.0);
  else {
    globalData->locaErrorCheck->throwError(
		   "LOCA::Stepper::reset()",
		   "\"Minimum Value\" of continuation parameter is not set!");
  }


  // Get the initial values or use their defaults
  stepSize = stepSizeStrategyPtr->getStartStepSize();
  maxNonlinearSteps = 
    stepperList->get("Max Nonlinear Iterations", 15);

  targetValue = 0.0;
  isTargetStep = false;
  tangentFactor = 1.0;
  doTangentFactorScaling =
    stepperList->get("Enable Tangent Factor Step Size Scaling", 
			      false);
  minTangentFactor = stepperList->get("Min Tangent Factor",0.1);
  tangentFactorExponent =
    stepperList->get("Tangent Factor Exponent",1.0);
  calcEigenvalues = stepperList->get("Compute Eigenvalues",false);

  // Make a copy of the parameter list, change continuation method to
  // natural
  Teuchos::RefCountPtr<Teuchos::ParameterList> firstStepperParams = 
    Teuchos::rcp(new Teuchos::ParameterList(*stepperList));
  firstStepperParams->set("Continuation Method", "Natural");

  // Create constrained group
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> constraintsGrp
    = buildConstrainedGroup(initialGuess);

  // Create bifurcation group
  Teuchos::RefCountPtr<Teuchos::ParameterList> bifurcationParams = 
    parsedParams->getSublist("Bifurcation");
  bifGroupPtr = globalData->locaFactory->createBifurcationStrategy(
						       parsedParams,
						       bifurcationParams,
						       constraintsGrp);

  // Create continuation strategy
  curGroupPtr = globalData->locaFactory->createContinuationStrategy(
							parsedParams,
							firstStepperParams,
							bifGroupPtr, predictor,
							conParamIDs);

  // Set step size
  curGroupPtr->setStepSize(0.0);

  // Set previous solution vector in current solution group
  curGroupPtr->setPrevX(curGroupPtr->getX());

  // Create solver using initial conditions
  solverPtr = Teuchos::rcp(new NOX::Solver::Manager(
					   curGroupPtr, 
					   statusTestPtr,
					   parsedParams->getSublist("NOX")));

  printInitializationInfo();

  if (globalData->locaUtils->isPrintType(NOX::Utils::Parameters))
    paramListPtr->print(globalData->locaUtils->out());

  return true;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::NewStepper::start() {
  NOX::StatusTest::StatusType solverStatus;
  string callingFunction = "LOCA::Stepper::start()";

  // Allow continuation group to preprocess the step
  curGroupPtr->preProcessContinuationStep(LOCA::Abstract::Iterator::Successful);

  printStartStep();

  // Perform solve of initial conditions
  solverStatus = solverPtr->solve();

  // Allow continuation group to postprocess the step
  if (solverStatus == NOX::StatusTest::Converged)
    curGroupPtr->postProcessContinuationStep(LOCA::Abstract::Iterator::Successful);
  else
    curGroupPtr->postProcessContinuationStep(LOCA::Abstract::Iterator::Unsuccessful);

  // Set up continuation groups
  const LOCA::MultiContinuation::ExtendedGroup& constSolnGrp =
    dynamic_cast<const LOCA::MultiContinuation::ExtendedGroup&>(
       solverPtr->getSolutionGroup());
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> underlyingGroup 
    = Teuchos::rcp_const_cast<LOCA::MultiContinuation::AbstractGroup>(constSolnGrp.getUnderlyingGroup());

  // Create continuation strategy
  curGroupPtr = globalData->locaFactory->createContinuationStrategy(
							parsedParams,
							stepperList,
							underlyingGroup, 
							predictor,
							conParamIDs);

  // Do printing (stepNumber==0 case) after continuation group set up
  if (solverStatus == NOX::StatusTest::Failed)
    printEndStep(LOCA::Abstract::Iterator::Unsuccessful);
  else
    printEndStep(LOCA::Abstract::Iterator::Successful);

  // Set the initial step size
  curGroupPtr->setStepSize(stepSize);

  prevGroupPtr = 
    Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::AbstractStrategy>(
	curGroupPtr->clone());

  // If nonlinear solve failed, return (this must be done after continuation
  // groups are created so Stepper::getSolutionGroup() functions correctly.
  if (solverStatus != NOX::StatusTest::Converged)
    return LOCA::Abstract::Iterator::Failed;

  // Save initial solution
  curGroupPtr->printSolution();

  // Compute eigenvalues/eigenvectors if requested
  if (calcEigenvalues) {
    Teuchos::RefCountPtr< std::vector<double> > evals_r;
    Teuchos::RefCountPtr< std::vector<double> > evals_i;
    Teuchos::RefCountPtr< NOX::Abstract::MultiVector > evecs_r;
    Teuchos::RefCountPtr< NOX::Abstract::MultiVector > evecs_i;
    eigensolver->computeEigenvalues(
				 *curGroupPtr->getBaseLevelUnderlyingGroup(),
				 evals_r, evals_i, evecs_r, evecs_i);

    saveEigenData->save(evals_r, evals_i, evecs_r, evecs_i);
  }

  // Compute predictor direction
  NOX::Abstract::Group::ReturnType predictorStatus =
    curGroupPtr->computePredictor();
  globalData->locaErrorCheck->checkReturnType(predictorStatus, 
					      callingFunction);
  curPredictorPtr =
    Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(curGroupPtr->getPredictorTangent()[0].clone(NOX::DeepCopy));
  prevPredictorPtr =
    Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(curGroupPtr->getPredictorTangent()[0].clone(NOX::ShapeCopy));

  // Create new solver using new continuation groups and combo status test
  solverPtr = Teuchos::rcp(new NOX::Solver::Manager(
					   curGroupPtr, statusTestPtr,
					   parsedParams->getSublist("NOX")));

  return LOCA::Abstract::Iterator::NotFinished;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::NewStepper::finish(LOCA::Abstract::Iterator::IteratorStatus itStatus)
{
  string callingFunction = "LOCA::Stepper::finish()";

  //
  // We don't need to check if the last step was successful since finish
  // is never called if it wasn't.  We might want to change that if there is
  // some post processing we want to do even if the run failed.
  //

  // Copy last solution
  curGroupPtr->copy(solverPtr->getSolutionGroup());

  // Return if iteration failed (reached max number of steps)
  if (itStatus == LOCA::Abstract::Iterator::Failed)
    return itStatus;

  // Do one additional step using natural continuation to hit target value
  double value = curGroupPtr->getContinuationParameter();

  if (fabs(value-targetValue) > 1.0e-15*(1.0 + fabs(targetValue))) {

    isTargetStep = true;

    // Save previous successful step information
    prevGroupPtr->copy(*curGroupPtr);

    // Get bifurcation group if there is one, or solution group if not
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> underlyingGrp
      = curGroupPtr->getUnderlyingGroup();

    // Create predictor strategy
    Teuchos::RefCountPtr<Teuchos::ParameterList> lastStepPredictorParams = 
      parsedParams->getSublist("Last Step Predictor");
    predictor = globalData->locaFactory->createPredictorStrategy(
						     parsedParams,
						     lastStepPredictorParams);

    // Make a copy of the parameter list, change continuation method to
    // natural
    Teuchos::RefCountPtr<Teuchos::ParameterList> lastStepperParams = 
      Teuchos::rcp(new Teuchos::ParameterList(*stepperList));
    lastStepperParams->set("Continuation Method", "Natural");

    // Create continuation strategy
    curGroupPtr = globalData->locaFactory->createContinuationStrategy(
							  parsedParams,
							  lastStepperParams,
							  underlyingGrp, 
							  predictor,
							  conParamIDs);

    // Set step size
    stepSize = targetValue - value;
    curGroupPtr->setStepSize(stepSize);

    // Get predictor direction
    NOX::Abstract::Group::ReturnType predictorStatus =
      curGroupPtr->computePredictor();
    globalData->locaErrorCheck->checkReturnType(predictorStatus, 
						callingFunction);
    *curPredictorPtr = curGroupPtr->getPredictorTangent()[0];

    // Set previous solution vector in current solution group
    curGroupPtr->setPrevX(curGroupPtr->getX());

    // Take step in predictor direction
    curGroupPtr->computeX(*curGroupPtr, *curPredictorPtr, stepSize);

    // Allow continuation group to preprocess the step
    curGroupPtr->preProcessContinuationStep(LOCA::Abstract::Iterator::Successful);

    printStartStep();

    // Create new solver
    solverPtr = Teuchos::rcp(new NOX::Solver::Manager(
					   curGroupPtr, statusTestPtr,
					   parsedParams->getSublist("NOX")));

    // Solve step
    NOX::StatusTest::StatusType solverStatus = solverPtr->solve();

    // Allow continuation group to postprocess the step
    if (solverStatus == NOX::StatusTest::Converged)
      curGroupPtr->postProcessContinuationStep(LOCA::Abstract::Iterator::Successful);
    else
      curGroupPtr->postProcessContinuationStep(LOCA::Abstract::Iterator::Unsuccessful);

    // Get solution
    curGroupPtr->copy(solverPtr->getSolutionGroup());

    if (solverStatus != NOX::StatusTest::Converged) {
      printEndStep(LOCA::Abstract::Iterator::Unsuccessful);
      return LOCA::Abstract::Iterator::Failed;
    }

    printEndStep(LOCA::Abstract::Iterator::Successful);

    curGroupPtr->printSolution();

  }

  return LOCA::Abstract::Iterator::Finished;
}

LOCA::Abstract::Iterator::StepStatus
LOCA::NewStepper::preprocess(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful) {

    // Restore previous step information
    curGroupPtr->copy(*prevGroupPtr);
  }
  else {

    // Save previous successful step information
    prevGroupPtr->copy(*curGroupPtr);
  }

  // Compute step size
  stepStatus = computeStepSize(stepStatus, stepSize);

  // Set step size in current solution group
  curGroupPtr->setStepSize(stepSize);

  // Set previous solution vector in current solution group
  curGroupPtr->setPrevX(prevGroupPtr->getX());

  // Take step in predictor direction
  curGroupPtr->computeX(*prevGroupPtr, *curPredictorPtr, stepSize);

  // Allow continuation group to preprocess the step
  curGroupPtr->preProcessContinuationStep(stepStatus);

  // Reset solver to compute new solution
//   solverPtr->reset(*curGroupPtr, *statusTestPtr,
// 		   parsedParams->getSublist("NOX"));

  solverPtr = Teuchos::rcp(new NOX::Solver::Manager(
					    curGroupPtr, statusTestPtr,
					    parsedParams->getSublist("NOX")));

  return stepStatus;
}

LOCA::Abstract::Iterator::StepStatus
LOCA::NewStepper::compute(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  NOX::StatusTest::StatusType solverStatus;

  // Print info for beginning of step
  printStartStep();

  // Compute next point on continuation curve
  solverStatus = solverPtr->solve();

  // Check solver status
  if (solverStatus == NOX::StatusTest::Failed) {
    printEndStep(LOCA::Abstract::Iterator::Unsuccessful);
    return LOCA::Abstract::Iterator::Unsuccessful;
  }

  // Copy solution out of solver
  curGroupPtr->copy(solverPtr->getSolutionGroup());

  // Print successful info for end of step
  printEndStep(LOCA::Abstract::Iterator::Successful);

  return LOCA::Abstract::Iterator::Successful;
}

LOCA::Abstract::Iterator::StepStatus
LOCA::NewStepper::postprocess(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  string callingFunction = "LOCA::Stepper::postprocess()";

  // Allow continuation group to postprocess the step
  curGroupPtr->postProcessContinuationStep(stepStatus);

  if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful)
    return stepStatus;

  *prevPredictorPtr = *curPredictorPtr;

  NOX::Abstract::Group::ReturnType predictorStatus =
    curGroupPtr->computePredictor();
  globalData->locaErrorCheck->checkReturnType(predictorStatus, 
					      callingFunction);
  *curPredictorPtr = curGroupPtr->getPredictorTangent()[0];

  if (doTangentFactorScaling && (getStepNumber() > 1)) {
    tangentFactor = curGroupPtr->computeScaledDotProduct(*curPredictorPtr,
							 *prevPredictorPtr) /
      sqrt(curGroupPtr->computeScaledDotProduct(*curPredictorPtr,
						*curPredictorPtr) *
	   curGroupPtr->computeScaledDotProduct(*prevPredictorPtr,
						 *prevPredictorPtr));

    if (tangentFactor < minTangentFactor) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
	globalData->locaUtils->out() 
	  << "\n\tTangent factor scaling:  Failing step!  Tangent factor " 
	  << "less than" << std::endl << "\t\tspecified bound: " 
	  << globalData->locaUtils->sciformat(tangentFactor) << " < " 
	  << globalData->locaUtils->sciformat(minTangentFactor) << std::endl;
      }
      return LOCA::Abstract::Iterator::Unsuccessful;
    }
  }

  // Print (save) solution
  curGroupPtr->printSolution();

  // Compute eigenvalues/eigenvectors
  if (calcEigenvalues) {
    Teuchos::RefCountPtr< std::vector<double> > evals_r;
    Teuchos::RefCountPtr< std::vector<double> > evals_i;
    Teuchos::RefCountPtr< NOX::Abstract::MultiVector > evecs_r;
    Teuchos::RefCountPtr< NOX::Abstract::MultiVector > evecs_i;
    eigensolver->computeEigenvalues(
				 *curGroupPtr->getBaseLevelUnderlyingGroup(),
				 evals_r, evals_i, evecs_r, evecs_i);

    saveEigenData->save(evals_r, evals_i, evecs_r, evecs_i);
  }

  return stepStatus;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::NewStepper::stop(LOCA::Abstract::Iterator::StepStatus stepStatus)
{

  // Check to see if max number of steps has been reached
  if (LOCA::Abstract::Iterator::numTotalSteps
        >= LOCA::Abstract::Iterator::maxSteps) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
      globalData->locaUtils->out() 
	<< "\n\tContinuation run stopping: reached maximum number of steps " 
	<< LOCA::Abstract::Iterator::maxSteps << std::endl;
    }
    return LOCA::Abstract::Iterator::Failed;
  }

  if (stepStatus == LOCA::Abstract::Iterator::Successful) {

    double value = curGroupPtr->getContinuationParameter();
    double paramStep = value - prevGroupPtr->getContinuationParameter();

    // See if we went past bounds for parameter
    if ( value >= maxValue*(1.0 - 1.0e-15) && paramStep > 0 ) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
	 globalData->locaUtils->out() 
	   << "\n\tContinuation run stopping: parameter reached bound of " 
	   << globalData->locaUtils->sciformat(maxValue) << std::endl;
      }
      targetValue = maxValue;
      return LOCA::Abstract::Iterator::Finished;
    }
    if ( value <= minValue*(1.0 + 1.0e-15) && paramStep < 0 ) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
	globalData->locaUtils->out() 
	  << "\n\tContinuation run stopping: parameter reached bound of " 
	  << globalData->locaUtils->sciformat(minValue) << std::endl;
      }
      targetValue = minValue;
      return LOCA::Abstract::Iterator::Finished;
    }

    // Check to see if arclength step was aimed to reach bound
    if (isLastIteration()) {

      // Check to see if continuation parameter is within threshold of bound
      if (withinThreshold()) {
	if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
	  globalData->locaUtils->out() 
	    << "\n\tContinuation run stopping: parameter stepped to bound" 
	    << std::endl;
	}
	return LOCA::Abstract::Iterator::Finished;
      }
      else
	return LOCA::Abstract::Iterator::NotFinished;
    }
  }
  else if (isLastIteration())  // Failed step did not reach bounds as predicted
    return LOCA::Abstract::Iterator::NotFinished;

  return LOCA::Abstract::Iterator::NotFinished;
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
LOCA::NewStepper::buildConstrainedGroup(
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp)
{
  // Get constraints sublist
  Teuchos::RefCountPtr<Teuchos::ParameterList> constraintsList =
    parsedParams->getSublist("Constraints");

  // If we don't have a constraint object, return original group
  if (!constraintsList->isParameter("Constraint Object"))
    return grp;

  string methodName = "LOCA::NewStepper::buildConstrainedGroup()";

  Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface> constraints;
  Teuchos::RefCountPtr< vector<string> > constraintParamNames;

  // Get constraint object
  if ((*constraintsList).INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface> >("Constraint Object"))
    constraints = (*constraintsList).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface> >("Constraint Object");
  else
    globalData->locaErrorCheck->throwError(methodName,
	  "\"Constraint Object\" parameter is not of type Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>!");

  // Get parameter names for constraints
  if ((*constraintsList).INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RefCountPtr< vector<string> > > ("Constraint Parameter Names"))
    constraintParamNames = (*constraintsList).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RefCountPtr< vector<string> > > ("Constraint Parameter Names");
  else
    globalData->locaErrorCheck->throwError(methodName,
	  "\"Constraint Parameter Names\" parameter is not of type Teuchos::RefCountPtr< vector<string> >!");

  // Convert names to integer IDs
  vector<int> constraintParamIDs(constraintParamNames->size());
  const LOCA::ParameterVector& pvec = grp->getParams();
  for (unsigned int i=0; i<constraintParamIDs.size(); i++)
    constraintParamIDs[i] = pvec.getIndex((*constraintParamNames)[i]);

  // Create constrained group
  return 
    Teuchos::rcp(new LOCA::MultiContinuation::ConstrainedGroup(
							globalData,
							parsedParams,
							constraintsList,
							grp,
							constraints,
							constraintParamIDs));
}

LOCA::Abstract::Iterator::StepStatus
LOCA::NewStepper::computeStepSize(LOCA::Abstract::Iterator::StepStatus stepStatus,
			       double& stepSz)
{
  NOX::Abstract::Group::ReturnType res =
    stepSizeStrategyPtr->computeStepSize(*curGroupPtr, *curPredictorPtr, 
					 *solverPtr,
					 stepStatus, *this, stepSz);

  if (res == NOX::Abstract::Group::Failed)
    return LOCA::Abstract::Iterator::Provisional;

  if (doTangentFactorScaling) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
      globalData->locaUtils->out() 
	<< "\n\tTangent factor scaling:  Rescaling step size by " 
	<< globalData->locaUtils->sciformat(pow(fabs(tangentFactor), 
						tangentFactorExponent)) 
	<< std::endl;
    }

    stepSz *= pow(fabs(tangentFactor), tangentFactorExponent);
  }

  // Cap the con parameter so we don't go past bounds
  double prevValue = curGroupPtr->getContinuationParameter();
  double dpds = curPredictorPtr->getScalar(0);
  if ( (prevValue+stepSz*dpds > maxValue*(1.0 - 1.0e-15)) ) {
    stepSz = (maxValue - prevValue)/dpds;
    targetValue = maxValue;
    setLastIteration();
  }
  if ( (prevValue+stepSz*dpds < minValue*(1.0 + 1.0e-15)) ) {
    stepSz = (minValue - prevValue)/dpds;
    targetValue = minValue;
    setLastIteration();
  }

  return LOCA::Abstract::Iterator::Successful;
}

Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup>
LOCA::NewStepper::getSolutionGroup() const
{
  return curGroupPtr->getBaseLevelUnderlyingGroup();
}

Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup>
LOCA::NewStepper::getBifurcationGroup() const
{
  return curGroupPtr->getUnderlyingGroup();
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
LOCA::NewStepper::getList() const
{
  return paramListPtr;
}

Teuchos::RefCountPtr<const NOX::Solver::Generic>
LOCA::NewStepper::getSolver() const
{
  if (solverPtr.get() == NULL) {
    globalData->locaErrorCheck->throwError(
				    "LOCA::Stepper::getSolver()",
				    "Solver has not been constructed yet!");
  }

  return solverPtr;
}

void
LOCA::NewStepper::printInitializationInfo()
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
    globalData->locaUtils->out() 
      << std::endl 
      << globalData->locaUtils->fill(72, '~') << std::endl;

     globalData->locaUtils->out() 
       << "Beginning Continuation Run \n"
       << "Stepper Method:             " 
       << stepperList->get("Continuation Method", "None")
       << "\n"
       << "Initial Parameter Value = " 
       << globalData->locaUtils->sciformat(startValue)
       << "\n"
       << "Maximum Parameter Value = " 
       << globalData->locaUtils->sciformat(maxValue) << "\n"
       << "Minimum Parameter Value = " 
       << globalData->locaUtils->sciformat(minValue) << "\n"
       << "Maximum Number of Continuation Steps = "
       << LOCA::Abstract::Iterator::maxSteps
       << std::endl;

     globalData->locaUtils->out() 
       << globalData->locaUtils->fill(72, '~') << std::endl << std::endl;
  }
}

void
LOCA::NewStepper::printStartStep()
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
    globalData->locaUtils->out() 
      << std::endl << globalData->locaUtils->fill(72, '~') << std::endl;

    globalData->locaUtils->out() 
      << "Start of Continuation Step " << stepNumber << " : ";
    if (stepNumber==0) {
      globalData->locaUtils->out() 
	<< "Attempting to converge initial guess at initial parameter "
	<< "values." << std::endl;
    }
    else if (isTargetStep) {
      globalData->locaUtils->out() 
	<< "Attempting to hit final target value "
	<< globalData->locaUtils->sciformat(targetValue) << std::endl;
    }
    else {
      globalData->locaUtils->out() 
	<< "Parameter: " << conParamName
	<< " = "
	<< globalData->locaUtils->sciformat(curGroupPtr->getContinuationParameter())
	<< " from "
	<< globalData->locaUtils->sciformat(prevGroupPtr->getContinuationParameter())
	<< std::endl;
      globalData->locaUtils->out() 
	<< "Continuation Method: " 
	<< stepperList->get("Continuation Method", "None")
	<< std::endl;
      globalData->locaUtils->out() 
	<< "Current step size  = " 
	<< globalData->locaUtils->sciformat(stepSize) << "   "
	<< "Previous step size = "
	<< globalData->locaUtils->sciformat(stepSizeStrategyPtr->getPrevStepSize()) 
	<< std::endl;
    }
    globalData->locaUtils->out() 
      << globalData->locaUtils->fill(72, '~') << std::endl << std::endl;
  }
}

void
LOCA::NewStepper::printEndStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Successful) {
    // Print results of successful continuation step
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
      globalData->locaUtils->out() 
	<< std::endl << globalData->locaUtils->fill(72, '~') << std::endl;
      globalData->locaUtils->out() 
	<< "End of Continuation Step " << stepNumber << " : "
	<< "Parameter: " << conParamName << " = "
	<< globalData->locaUtils->sciformat(curGroupPtr->getContinuationParameter());
      if (stepNumber != 0)
        globalData->locaUtils->out() 
	  << " from "
	  << globalData->locaUtils->sciformat(prevGroupPtr->getContinuationParameter());
      globalData->locaUtils->out()
	<< std::endl << "--> Step Converged in "
	<< solverPtr->getNumIterations()
	<<" Nonlinear Solver Iterations!\n";
      globalData->locaUtils->out() 
	<< globalData->locaUtils->fill(72, '~') << std::endl << std::endl;
    }
  }
  else {
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
      // RPP: We may not need this, the failure info should be
      // at the method level!
      globalData->locaUtils->out() 
	<< std::endl << globalData->locaUtils->fill(72, '~') << std::endl;
      globalData->locaUtils->out() 
	<< "Continuation Step Number " << stepNumber
	<< " experienced a convergence failure in\n"
	<< "the nonlinear solver after "<< solverPtr->getNumIterations()
	<<" Iterations\n";
      globalData->locaUtils->out() 
	<< "Value of continuation parameter at failed step = "
	<< globalData->locaUtils->sciformat(curGroupPtr->getContinuationParameter());
      if (stepNumber != 0)
        globalData->locaUtils->out() 
	  << " from "
	  << globalData->locaUtils->sciformat(prevGroupPtr->getContinuationParameter());
      globalData->locaUtils->out() 
	<< std::endl << globalData->locaUtils->fill(72, '~') << endl;
    }
  }
}

void
LOCA::NewStepper::printEndInfo()
{

}

bool
LOCA::NewStepper::withinThreshold()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> stepSizeList = 
    parsedParams->getSublist("Step Size");
  double relt = stepperList->get("Relative Stopping Threshold", 0.9);
  double initialStep = stepSizeList->get("Initial Step Size", 1.0);
  double conParam = curGroupPtr->getContinuationParameter();

  return (fabs(conParam-targetValue) < relt*fabs(initialStep));
}
