// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2012) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact 
// Glen Hansen (gahanse@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
          Teuchos::RCP<LOCA::Epetra::AdaptiveSolutionManager>& solnManager_,
          const Teuchos::RCP<NOX::StatusTest::Generic>& nt) :
  LOCA::Abstract::Iterator(),
  mgr(solnManager_),
  globalData(),
  parsedParams(),
  predictor(),
  curGroupPtr(),
  prevGroupPtr(),
  eigensolver(),
  saveEigenData(),
  bifGroupPtr(),
  noxStatusTestPtr(),
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
  calcEigenvalues(false),
  return_failed_on_max_steps(true)
{

  reset(solnManager_->getGlobalDataNonConst(), 
              nt, solnManager_->getPiroParamsNonConst());
  reSize();

}


LOCA::Epetra::AdaptiveStepper::~AdaptiveStepper()
{
}

bool
LOCA::Epetra::AdaptiveStepper::reset(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const Teuchos::RCP<NOX::StatusTest::Generic>& nt,
		    const Teuchos::RCP<Teuchos::ParameterList>& p )
{
  globalData = global_data;
  paramListPtr = p;
  noxStatusTestPtr = nt;

  // Parse parameter list
  parsedParams = Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));
  parsedParams->parseSublists(paramListPtr);

  // Get stepper sublist
  stepperList = parsedParams->getSublist("Stepper");

  // Reset base class
  LOCA::Abstract::Iterator::resetIterator(*stepperList);

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
   stepSizeStrategyPtr = globalData->locaFactory->createStepSizeStrategy(
							     parsedParams,
							     stepsizeParams);

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

  // TODO Deprecated as moved to LOCA::StatusTest::MaxIters
  return_failed_on_max_steps =
    stepperList->get("Return Failed on Reaching Max Steps", true);

  // Make a copy of the parameter list, change continuation method to
  // natural
  firstStepperParams = 
    Teuchos::rcp(new Teuchos::ParameterList(*stepperList));
  firstStepperParams->set("Continuation Method", "Natural");

  bifurcationParams = 
    parsedParams->getSublist("Bifurcation");

  return true;
}

void
LOCA::Epetra::AdaptiveStepper::reSize(){

  setSolutionGroup(mgr->buildSolutionGroup());

}

void
LOCA::Epetra::AdaptiveStepper::setSolutionGroup(
  const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initialGuess){

  initialGuess->setParam(conParamName, startValue);

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

  printInitializationInfo();

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters))
    paramListPtr->print(globalData->locaUtils->out());

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
  NOX::StatusTest::StatusType solverStatus;
  string callingFunction = "LOCA_AdaptiveStepper::start()";

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

  // If nonlinear solve failed, return (this must be done after continuation
  // groups are created so AdaptiveStepper::getSolutionGroup() functions correctly.
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

  return LOCA::Abstract::Iterator::NotFinished;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Epetra::AdaptiveStepper::finish(LOCA::Abstract::Iterator::IteratorStatus itStatus)
{
  string callingFunction = "LOCA_AdaptiveStepper::finish()";

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

  bool do_target = stepperList->get("Hit Continuation Bound", true);
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

//  mgr->projectCurrentSolution(curGroupPtr->getX());
//  mgr->projectCurrentSolution(curGroupPtr);
  mgr->projectCurrentSolution();

// The new solution group is created by the solution manager, using the current discretization

  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> newSolnGroup = mgr->buildSolutionGroup();

  newSolnGroup->setParam(conParamName, value); 

  // Get the continuation parameter index
  const LOCA::ParameterVector& pv = newSolnGroup->getParams();
  conParamIDs[0] = pv.getIndex(conParamName);

  // Rebuild constrained group
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> constraintsGrp
    = buildConstrainedGroup(newSolnGroup);

  // Rebuild bifurcation group
  bifGroupPtr = globalData->locaFactory->createBifurcationStrategy(
						       parsedParams,
						       bifurcationParams,
					       constraintsGrp);

  // Rebuild continuation strategy
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

  // Create solver using newSolnGroup conditions
  solverPtr = NOX::Solver::buildSolver(curGroupPtr, noxStatusTestPtr,
				       parsedParams->getSublist("NOX"));

  // Allow continuation group to preprocess the step
  curGroupPtr->preProcessContinuationStep(LOCA::Abstract::Iterator::Successful);

  printRelaxationStep();

  // Perform solve of newSolnGroup conditions
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

  if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful) { // Previous step unsuccessful, stop

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
  iteratorStatus = start();
  if (iteratorStatus == LOCA::Abstract::Iterator::Failed)
    return iteratorStatus;

  stepNumber++;
  mgr->getAdaptManager()->setIteration(stepNumber);
  mgr->getAdaptManager()->setTime(getContinuationParameter());

  iteratorStatus = iterate();

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

  while (iteratorStatus == LOCA::Abstract::Iterator::NotFinished) {

    Teuchos::RCP<NOX::Epetra::AdaptManager> adaptManager = mgr->getAdaptManager();

    if(!adaptManager.is_null() 
      && stepStatus != LOCA::Abstract::Iterator::Unsuccessful
      && adaptManager->queryAdaptationCriteria()) {

      // Start adapting the mesh
      if(!adaptManager->startMeshAdaptation())

        return LOCA::Abstract::Iterator::Failed; // Abort if mesh manipulation fails

      // Project (remap the physics) and relax the solution on the new mesh
    
      preStatus = adapt(stepStatus);

      // Abort if this fails

      if(preStatus == LOCA::Abstract::Iterator::Unsuccessful)
        return LOCA::Abstract::Iterator::Failed;

      // Commit any remaining changes and move forward

      adaptManager->finalizeMeshAdaptation();

    }
    else

      preStatus = preprocess(stepStatus);

    compStatus = compute(preStatus);

    postStatus = postprocess(compStatus);

    stepStatus = computeStepStatus(preStatus, compStatus, postStatus);

    ++numTotalSteps;
    if (stepStatus ==  LOCA::Abstract::Iterator::Successful)
      ++stepNumber;
    else
      ++numFailedSteps;

    mgr->getAdaptManager()->setIteration(stepNumber);
    mgr->getAdaptManager()->setTime(getContinuationParameter());

    if (iteratorStatus != LOCA::Abstract::Iterator::Failed)
      iteratorStatus = stop(stepStatus);
  }

  return iteratorStatus;
}

LOCA::Abstract::Iterator::StepStatus
LOCA::Epetra::AdaptiveStepper::compute(LOCA::Abstract::Iterator::StepStatus stepStatus)
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
  string callingFunction = "LOCA_AdaptiveStepper::postprocess()";

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
    if (return_failed_on_max_steps)
      return LOCA::Abstract::Iterator::Failed;
    else
      return LOCA::Abstract::Iterator::Finished;
      
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

  string methodName = "LOCA::Epetra::AdaptiveStepper::buildConstrainedGroup()";

  Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> constraints;
  Teuchos::RCP< vector<string> > constraintParamNames;

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
      isType< Teuchos::RCP< vector<string> > > ("Constraint Parameter Names"))
    constraintParamNames = (*constraintsList).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP< vector<string> > > ("Constraint Parameter Names");
  else
    globalData->locaErrorCheck->throwError(methodName,
	  "\"Constraint Parameter Names\" parameter is not of type Teuchos::RCP< vector<string> >!");

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
	<< std::endl << globalData->locaUtils->fill(72, '~') << endl;
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
	<< std::endl << globalData->locaUtils->fill(72, '~') << endl;
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

