// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Epetra_AdaptiveStepper.H"    // class definition
#include "NOX_StatusTest_Generic.H"

// LOCA Includes
#include "NOX_Utils.H"
#include "NOX_Solver_Factory.H"
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

LOCA::Epetra::AdaptiveStepper::AdaptiveStepper(
          const Teuchos::RCP<Teuchos::ParameterList>& pList,
          const Teuchos::RCP<LOCA::Epetra::AdaptiveSolutionManager>& solnManager_,
          const Teuchos::RCP<LOCA::GlobalData>& global_data,
          const Teuchos::RCP<NOX::StatusTest::Generic>& nt) :

  LOCA::Abstract::Iterator(),
  mgr(solnManager_),
  globalData(global_data),
  parsedParams(),
  predictor(),
  curGroupPtr(),
  prevGroupPtr(),
  eigensolver(),
  saveEigenData(),
  bifGroupPtr(),
  noxStatusTestPtr(nt),
  paramListPtr(pList),
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
  calcEigenvalues(false),
  return_failed_on_max_steps(true),
  max_steps_exceeded(false)
{

  // Parse parameter list
  parsedParams = Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));

  parsedParams->parseSublists(paramListPtr);

  // Get stepper sublist
  stepperList = parsedParams->getSublist("Stepper");

  // Reset base class
  LOCA::Abstract::Iterator::resetIterator(*stepperList);

  // Build the LOCA Factory needed to create the various strategies for stepping
  buildLOCAFactory();

  // Get the continuation parameter starting value
  if (stepperList->isParameter("Initial Value"))
    startValue = stepperList->get("Initial Value", 0.0);
  else {
    globalData->locaErrorCheck->throwError(
           "LOCA::Epetra::AdaptiveStepper::reset()",
           "\"Initial Value\" of continuation parameter is not set!");
  }

  // Get the continuation parameter name
  if (stepperList->isParameter("Continuation Parameter")) {
    conParamName = stepperList->get("Continuation Parameter","None");
  }
  else {
     globalData->locaErrorCheck->throwError(
                  "LOCA::Epetra::AdaptiveStepper::reset()",
                  "\"Continuation Parameter\" name is not set!");
  }

  // Get the max and min values of the continuation parameter
  if (stepperList->isParameter("Max Value"))
    maxValue = stepperList->get("Max Value", 0.0);
  else {
     globalData->locaErrorCheck->throwError(
           "LOCA::Epetra::AdaptiveStepper::reset()",
           "\"Maximum Value\" of continuation parameter is not set!");
  }
  if (stepperList->isParameter("Min Value"))
    minValue = stepperList->get("Min Value", 0.0);
  else {
    globalData->locaErrorCheck->throwError(
           "LOCA::Epetra::AdaptiveStepper::reset()",
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

  return_failed_on_max_steps =
    stepperList->get("Return Failed on Reaching Max Steps", true);

  // Make a copy of the parameter list, change continuation method to
  // natural
  firstStepperParams =
    Teuchos::rcp(new Teuchos::ParameterList(*stepperList));
  firstStepperParams->set("Continuation Method", "Natural");

  bifurcationParams =
    parsedParams->getSublist("Bifurcation");

  setSolutionGroup(mgr->buildSolutionGroup(), startValue);

  printInitializationInfo();

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters))
    paramListPtr->print(globalData->locaUtils->out());

}


LOCA::Epetra::AdaptiveStepper::~AdaptiveStepper()
{
}


void
LOCA::Epetra::AdaptiveStepper::buildLOCAFactory(){

  // Create predictor strategy
  Teuchos::RCP<Teuchos::ParameterList> predictorParams =
    parsedParams->getSublist("Predictor");
  predictor = globalData->locaFactory->createPredictorStrategy(
                                  parsedParams,
                                  predictorParams);

  // Create eigensolver
  Teuchos::RCP<Teuchos::ParameterList> eigenParams =
    parsedParams->getSublist("Eigensolver");
  eigensolver = globalData->locaFactory->createEigensolverStrategy(
                                parsedParams,
                                eigenParams);

  // Create strategy to save eigenvectors/values
  saveEigenData = globalData->locaFactory->createSaveEigenDataStrategy(
                                parsedParams,
                                eigenParams);

  // Create step size strategy
  Teuchos::RCP<Teuchos::ParameterList> stepsizeParams =
    parsedParams->getSublist("Step Size");

   if(Teuchos::is_null(stepSizeStrategyPtr))
     stepSizeStrategyPtr = globalData->locaFactory->createStepSizeStrategy(
                                 parsedParams,
                                 stepsizeParams);

}

void
LOCA::Epetra::AdaptiveStepper::setSolutionGroup(
  const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initialGuess, const double curTimeValue){

  initialGuess->setParam(conParamName, curTimeValue);

  // Get the continuation parameter index
  const LOCA::ParameterVector& pv = initialGuess->getParams();
  conParamIDs[0] = pv.getIndex(conParamName);

  // Create constrained group
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> constraintsGrp
    = buildConstrainedGroup(initialGuess);

  // Create bifurcation group
  bifGroupPtr = globalData->locaFactory->createBifurcationStrategy(
                               parsedParams,
                               bifurcationParams,
                               constraintsGrp);

  // Create continuation strategy
  curGroupPtr = globalData->locaFactory->createContinuationStrategy(
                            parsedParams,
                            firstStepperParams,
                            bifGroupPtr,
                            predictor,
                            conParamIDs);

  // Set step size
  curGroupPtr->setStepSize(0.0);

  // Set previous solution vector in current solution group
  curGroupPtr->setPrevX(curGroupPtr->getX());

  // Create solver using initial conditions
  solverPtr = NOX::Solver::buildSolver(curGroupPtr, noxStatusTestPtr,
                       parsedParams->getSublist("NOX"));

  return;

}

bool
LOCA::Epetra::AdaptiveStepper::eigensolverReset( Teuchos::RCP<Teuchos::ParameterList> & newEigensolverList ) {

   // overwrite the eigensolver parameter list
   const Teuchos::RCP<Teuchos::ParameterList> & eigenParams = parsedParams->getSublist("Eigensolver");
   *eigenParams = *newEigensolverList;
   // recreate the eigensolver
   eigensolver = globalData->locaFactory->createEigensolverStrategy( parsedParams,
                                                                     eigenParams  );
   return true;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Epetra::AdaptiveStepper::start() {

  // This is the relaxation (equilibration) step

  NOX::StatusTest::StatusType solverStatus;
  std::string callingFunction = "LOCA_AdaptiveStepper::start()";

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
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> underlyingGroup
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

  // If the equilibration nonlinear solve failed, return failure, as it makes little sense to proceed.
  //  (this must be done after continuation
  // groups are created so AdaptiveStepper::getSolutionGroup() functions correctly.)

  if (solverStatus != NOX::StatusTest::Converged)

    return LOCA::Abstract::Iterator::Failed;

  // Save initial solution
  curGroupPtr->printSolution();

  // Compute eigenvalues/eigenvectors if requested
  if (calcEigenvalues) {
    Teuchos::RCP< std::vector<double> > evals_r;
    Teuchos::RCP< std::vector<double> > evals_i;
    Teuchos::RCP< NOX::Abstract::MultiVector > evecs_r;
    Teuchos::RCP< NOX::Abstract::MultiVector > evecs_i;
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
    Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(
      curGroupPtr->getPredictorTangent()[0].clone(NOX::DeepCopy));

  prevPredictorPtr =
    Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(
      curGroupPtr->getPredictorTangent()[0].clone(NOX::ShapeCopy));

  // Create new solver using new continuation groups and combo status test
  solverPtr = NOX::Solver::buildSolver(curGroupPtr, noxStatusTestPtr,
                       parsedParams->getSublist("NOX"));

  // We're not done yet, just finished the equilibration step
  return LOCA::Abstract::Iterator::NotFinished;

}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Epetra::AdaptiveStepper::finish(LOCA::Abstract::Iterator::IteratorStatus /* itStatus */)
{
  std::string callingFunction = "LOCA_AdaptiveStepper::finish()";

  //
  // We don't need to check if the last step was successful since finish
  // is never called if it wasn't.  We might want to change that if there is
  // some post processing we want to do even if the run failed.
  //

  // Copy last solution
  curGroupPtr->copy(solverPtr->getSolutionGroup());

  bool do_target = stepperList->get("Hit Continuation Bound", true);

  // If we do not need to hit the final value exactly, return Finished
  if (!do_target)

    return LOCA::Abstract::Iterator::Finished;

  // Do one additional step using natural continuation to hit target value
  double value = curGroupPtr->getContinuationParameter();

  if (fabs(value-targetValue) > 1.0e-15*(1.0 + fabs(targetValue))) {

    isTargetStep = true;

    // Save previous successful step information
    prevGroupPtr->copy(*curGroupPtr);

    // Get bifurcation group if there is one, or solution group if not
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> underlyingGrp
      = curGroupPtr->getUnderlyingGroup();

    // Create predictor strategy
    Teuchos::RCP<Teuchos::ParameterList> lastStepPredictorParams =
      parsedParams->getSublist("Last Step Predictor");
    // change default method to constant to avoid infinite stack recursion
    lastStepPredictorParams->get("Method", "Constant");
    predictor = globalData->locaFactory->createPredictorStrategy(
                             parsedParams,
                             lastStepPredictorParams);

    // Make a copy of the parameter list, change continuation method to
    // natural
    Teuchos::RCP<Teuchos::ParameterList> lastStepperParams =
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
    solverPtr = NOX::Solver::buildSolver(curGroupPtr, noxStatusTestPtr,
                     parsedParams->getSublist("NOX"));

    // Solve step
    NOX::StatusTest::StatusType solverStatus = solverPtr->solve();

    // Allow continuation group to postprocess the step
    if (solverStatus == NOX::StatusTest::Converged)
      curGroupPtr->postProcessContinuationStep(LOCA::Abstract::Iterator::Successful);
    else
      curGroupPtr->postProcessContinuationStep(LOCA::Abstract::Iterator::Unsuccessful);

    // Get solution
    curGroupPtr->copy(solverPtr->getSolutionGroup());

    // Failed to converge on the last step
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
LOCA::Epetra::AdaptiveStepper::adapt(LOCA::Abstract::Iterator::StepStatus stepStatus)
{

  NOX::StatusTest::StatusType solverStatus;

  // get current value of continuation parameter
  double value = curGroupPtr->getContinuationParameter();

// Project the current solution into the solution vector corresponding to the new mesh
// This must be called prior to the call of buildSolutionGroup() below, or the current
// solution will be lost.

  mgr->projectCurrentSolution();

// The new solution group is created by the solution manager, using the current discretization

  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> newSolnGroup = mgr->buildSolutionGroup();

// Re-build the LOCA factory that creates the stepping criteria for the adapted problem

  buildLOCAFactory();

// Re-build the solution group at the larger size needed for the adapted mesh

  setSolutionGroup(newSolnGroup, value);

  // Allow continuation group to preprocess the step
  curGroupPtr->preProcessContinuationStep(LOCA::Abstract::Iterator::Successful);

  printRelaxationStep();

  // Perform solve of newSolnGroup conditions - this is the "relaxation" (equilibration) solve
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
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> underlyingGroup
    = Teuchos::rcp_const_cast<LOCA::MultiContinuation::AbstractGroup>(constSolnGrp.getUnderlyingGroup());

  // Create continuation strategy
  curGroupPtr = globalData->locaFactory->createContinuationStrategy(
                            parsedParams,
                            stepperList,
                            underlyingGroup,
                            predictor,
                            conParamIDs);

  // Do printing (relaxation case) after continuation group set up
  if (solverStatus == NOX::StatusTest::Failed)
    printRelaxationEndStep(LOCA::Abstract::Iterator::Unsuccessful);
  else
    printRelaxationEndStep(LOCA::Abstract::Iterator::Successful);

  prevGroupPtr =
    Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::AbstractStrategy>(
    curGroupPtr->clone());

  // If relaxation solve failed, return (this must be done after continuation
  // groups are created so AdaptiveStepper::getSolutionGroup() functions correctly.

  if (solverStatus != NOX::StatusTest::Converged)
    return LOCA::Abstract::Iterator::Unsuccessful;

  if(mgr->getAdaptParamsNonConst()->get<bool>("Print Relaxation Solution", false)){

    // Save relaxation solution to the output file
    curGroupPtr->printSolution();

    // Compute eigenvalues/eigenvectors if requested
    if (calcEigenvalues) {
      Teuchos::RCP< std::vector<double> > evals_r;
      Teuchos::RCP< std::vector<double> > evals_i;
      Teuchos::RCP< NOX::Abstract::MultiVector > evecs_r;
      Teuchos::RCP< NOX::Abstract::MultiVector > evecs_i;
      eigensolver->computeEigenvalues(
                   *curGroupPtr->getBaseLevelUnderlyingGroup(),
                   evals_r, evals_i, evecs_r, evecs_i);

      saveEigenData->save(evals_r, evals_i, evecs_r, evecs_i);
    }
  }

  // We have successfully relaxed the solution from the remesh, now prepare to resume stepping

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
  solverPtr = NOX::Solver::buildSolver(curGroupPtr, noxStatusTestPtr,
                       parsedParams->getSublist("NOX"));

  return stepStatus;

}

LOCA::Abstract::Iterator::StepStatus
LOCA::Epetra::AdaptiveStepper::preprocess(LOCA::Abstract::Iterator::StepStatus stepStatus)
{

  if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful) { // Previous step unsuccessful

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
  solverPtr = NOX::Solver::buildSolver(curGroupPtr, noxStatusTestPtr,
                       parsedParams->getSublist("NOX"));

  return stepStatus;

}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Epetra::AdaptiveStepper::run()
{
  // We return one of two successful states to Piro:
  //   LOCA::Abstract::Iterator::Finished || LOCA::Abstract::Iterator::NotFinished
  // Anything else indicates a fatal error and we are giving up

  iteratorStatus = start();

  // equilibration step failed, no sense in going on
  if (iteratorStatus == LOCA::Abstract::Iterator::Failed)

    return LOCA::Abstract::Iterator::Failed;

  stepNumber++;
  mgr->getAdaptManager()->setIteration(stepNumber);
  mgr->getAdaptManager()->setTime(getContinuationParameter());

  iteratorStatus = iterate();

  if (iteratorStatus == LOCA::Abstract::Iterator::Failed){
    if(max_steps_exceeded){
      // not a true failure, the last iteration converged but we exceeded the maximum iterations specified by the user
      max_steps_exceeded = false;
      return LOCA::Abstract::Iterator::NotFinished;
    }
    else
      // we failed in the LOCA iteration sequence and cannot recover, bail out
      return LOCA::Abstract::Iterator::Failed;
  }
  else if(iteratorStatus == LOCA::Abstract::Iterator::Finished){
      max_steps_exceeded = false;
      return LOCA::Abstract::Iterator::Finished;
  }

// This should only be called if the last Newton solve was successful but the status is NotFinished
  iteratorStatus = finish(iteratorStatus);

  return iteratorStatus;

}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Epetra::AdaptiveStepper::iterate()
{

  LOCA::Abstract::Iterator::StepStatus stepStatus =
    LOCA::Abstract::Iterator::Successful;
  LOCA::Abstract::Iterator::StepStatus preStatus;
  LOCA::Abstract::Iterator::StepStatus compStatus;
  LOCA::Abstract::Iterator::StepStatus postStatus;

  iteratorStatus = stop(stepStatus);

  // Loop until we finish the LOCA stepping trajectory, or get stuck and cannot recover.

  bool lastStepFailed = false;

  while (iteratorStatus == LOCA::Abstract::Iterator::NotFinished) {

    Teuchos::RCP<NOX::Epetra::AdaptManager> adaptManager = mgr->getAdaptManager();

    // Adapt the mesh and move forward
    if(!adaptManager.is_null()
      && stepStatus != LOCA::Abstract::Iterator::Unsuccessful
      && adaptManager->queryAdaptationCriteria() && !lastStepFailed) {

      // Adapt the mesh and problem size
      if(!mgr->adaptProblem())

        return LOCA::Abstract::Iterator::Failed; // Abort if mesh manipulation fails - this is Fatal

      // Project (remap the physics) and relax the solution on the new mesh

      preStatus = adapt(stepStatus);

      // Abort if this fails

      if(preStatus == LOCA::Abstract::Iterator::Unsuccessful)

        // Bail out if the projection cannot be done
        return LOCA::Abstract::Iterator::Failed;

    }
    else { // do not adapt the mesh

      // We skip adaptation and go directly here if the nonlinear solver just failed
      preStatus = preprocess(stepStatus);

    }

    compStatus = compute(preStatus);

    postStatus = postprocess(compStatus);

    stepStatus = computeStepStatus(preStatus, compStatus, postStatus);

    ++numTotalSteps;

    if (stepStatus ==  LOCA::Abstract::Iterator::Successful){
      // we sucessfully stepped, increment the step number
      ++stepNumber;
      lastStepFailed = false;
    }
    else {
      // we didn't, increment the failed step number
      ++numFailedSteps;
      lastStepFailed = true;
      continue; // repeat the failed step with a smaller LOCA step increment

    }

    mgr->getAdaptManager()->setIteration(stepNumber);
    mgr->getAdaptManager()->setTime(getContinuationParameter());

    if (iteratorStatus != LOCA::Abstract::Iterator::Failed)
      iteratorStatus = stop(stepStatus);

  }

  return iteratorStatus;
}

LOCA::Abstract::Iterator::StepStatus
LOCA::Epetra::AdaptiveStepper::compute(LOCA::Abstract::Iterator::StepStatus /* stepStatus */)
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
LOCA::Epetra::AdaptiveStepper::postprocess(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  std::string callingFunction = "LOCA_AdaptiveStepper::postprocess()";

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
    Teuchos::RCP< std::vector<double> > evals_r;
    Teuchos::RCP< std::vector<double> > evals_i;
    Teuchos::RCP< NOX::Abstract::MultiVector > evecs_r;
    Teuchos::RCP< NOX::Abstract::MultiVector > evecs_i;
    eigensolver->computeEigenvalues(
                 *curGroupPtr->getBaseLevelUnderlyingGroup(),
                 evals_r, evals_i, evecs_r, evecs_i);

    saveEigenData->save(evals_r, evals_i, evecs_r, evecs_i);
  }

  return stepStatus;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Epetra::AdaptiveStepper::stop(LOCA::Abstract::Iterator::StepStatus stepStatus)
{

  // Check to see if max number of steps has been reached

  if (LOCA::Abstract::Iterator::numTotalSteps
        >= LOCA::Abstract::Iterator::maxSteps) {

    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
      globalData->locaUtils->out()
        << "\n\tContinuation run stopping: reached maximum number of steps "
        << LOCA::Abstract::Iterator::maxSteps << std::endl;
    }

    if (return_failed_on_max_steps){
      max_steps_exceeded = true;
      return LOCA::Abstract::Iterator::Failed;
    }
    else {
      max_steps_exceeded = true;
      return LOCA::Abstract::Iterator::Finished;
    }

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

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::Epetra::AdaptiveStepper::buildConstrainedGroup(
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp)
{
  // Get constraints sublist
  Teuchos::RCP<Teuchos::ParameterList> constraintsList =
    parsedParams->getSublist("Constraints");

  // If we don't have a constraint object, return original group
  if (!constraintsList->isParameter("Constraint Object"))
    return grp;

  std::string methodName = "LOCA::Epetra::AdaptiveStepper::buildConstrainedGroup()";

  Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> constraints;
  Teuchos::RCP< std::vector<std::string> > constraintParamNames;

  // Get constraint object
  if ((*constraintsList).INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> >("Constraint Object"))
    constraints = (*constraintsList).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> >("Constraint Object");
  else
    globalData->locaErrorCheck->throwError(methodName,
      "\"Constraint Object\" parameter is not of type Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>!");

  // Get parameter names for constraints
  if ((*constraintsList).INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RCP< std::vector<std::string> > > ("Constraint Parameter Names"))
    constraintParamNames = (*constraintsList).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP< std::vector<std::string> > > ("Constraint Parameter Names");
  else
    globalData->locaErrorCheck->throwError(methodName,
      "\"Constraint Parameter Names\" parameter is not of type Teuchos::RCP< std::vector<std::string> >!");

  // Convert names to integer IDs
  std::vector<int> constraintParamIDs(constraintParamNames->size());
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
LOCA::Epetra::AdaptiveStepper::computeStepSize(LOCA::Abstract::Iterator::StepStatus stepStatus,
                   double& stepSz)
{
  NOX::Abstract::Group::ReturnType res =
    stepSizeStrategyPtr->computeStepSize(*curGroupPtr, *curPredictorPtr,
                     *solverPtr, stepStatus, *this, stepSz);

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

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::Epetra::AdaptiveStepper::getSolutionGroup() const
{
  return curGroupPtr->getBaseLevelUnderlyingGroup();
}

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::Epetra::AdaptiveStepper::getBifurcationGroup() const
{
  return curGroupPtr->getUnderlyingGroup();
}

Teuchos::RCP<const Teuchos::ParameterList>
LOCA::Epetra::AdaptiveStepper::getList() const
{
  return paramListPtr;
}

Teuchos::RCP<const NOX::Solver::Generic>
LOCA::Epetra::AdaptiveStepper::getSolver() const
{
  if (solverPtr.get() == NULL) {
    globalData->locaErrorCheck->throwError(
                    "LOCA_AdaptiveStepper::getSolver()",
                    "Solver has not been constructed yet!");
  }

  return solverPtr;
}

double
LOCA::Epetra::AdaptiveStepper::getContinuationParameter() const
{
  return curGroupPtr->getContinuationParameter();
}

void
LOCA::Epetra::AdaptiveStepper::printInitializationInfo()
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
    globalData->locaUtils->out()
      << std::endl
      << globalData->locaUtils->fill(72, '~') << std::endl;

     globalData->locaUtils->out()
       << "Beginning Continuation Run \n"
       << "AdaptiveStepper Method:             "
       << stepperList->get("Continuation Method", "Arc Length")
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
LOCA::Epetra::AdaptiveStepper::printStartStep()
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
    << stepperList->get("Continuation Method", "Arc Length")
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
LOCA::Epetra::AdaptiveStepper::printRelaxationStep()
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {

    globalData->locaUtils->out()
      << std::endl << globalData->locaUtils->fill(72, '~') << std::endl;

    globalData->locaUtils->out()
      << "Start of Continuation Step " << stepNumber << " : ";

    globalData->locaUtils->out()
    << "Attempting to converge the remeshed solution at current parameter "
    << "values." << std::endl;

    globalData->locaUtils->out()
      << globalData->locaUtils->fill(72, '~') << std::endl << std::endl;
  }
}

void
LOCA::Epetra::AdaptiveStepper::printEndStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
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
    << std::endl << globalData->locaUtils->fill(72, '~') << std::endl;
    }
  }
}

void
LOCA::Epetra::AdaptiveStepper::printRelaxationEndStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Successful) {
    // Print results of successful continuation step
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
      globalData->locaUtils->out()
    << std::endl << globalData->locaUtils->fill(72, '~') << std::endl;
      globalData->locaUtils->out()
    << "End of Relaxation Step " << stepNumber << " : "
    << "Parameter: " << conParamName << " = "
    << globalData->locaUtils->sciformat(curGroupPtr->getContinuationParameter());

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
    << "Relaxation Step Number " << stepNumber
    << " experienced a convergence failure in\n"
    << "the nonlinear solver after "<< solverPtr->getNumIterations()
    <<" Iterations\n";
      globalData->locaUtils->out()
    << "Value of continuation parameter at failed step = "
    << globalData->locaUtils->sciformat(curGroupPtr->getContinuationParameter());
     globalData->locaUtils->out()
    << std::endl << globalData->locaUtils->fill(72, '~') << std::endl;
    }
  }
}

void
LOCA::Epetra::AdaptiveStepper::printEndInfo()
{

}

bool
LOCA::Epetra::AdaptiveStepper::withinThreshold()
{
  Teuchos::RCP<Teuchos::ParameterList> stepSizeList =
    parsedParams->getSublist("Step Size");
  double relt = stepperList->get("Relative Stopping Threshold", 0.9);
  double initialStep = stepSizeList->get("Initial Step Size", 1.0);
  double conParam = curGroupPtr->getContinuationParameter();

  return (fabs(conParam-targetValue) < relt*fabs(initialStep));
}

