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
#include "LOCA_Stepper.H"    // class definition
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "LOCA_StatusTest_Wrapper.H"

// LOCA Includes
#include "LOCA_Utils.H"		                // for static function doPrint
#include "LOCA_Continuation_AbstractGroup.H"   // class data element
#include "LOCA_Continuation_ExtendedGroup.H"
#include "LOCA_Continuation_NaturalGroup.H"
#include "LOCA_Continuation_ArcLengthGroup.H"   //

/* Some compilers (in particular the SGI and ASCI Red - TFLOP) 
 * fail to find the max and min function.  Therfore we redefine them 
 * here. 
 */ 
#ifdef max
#undef max
#endif
#define max(a,b) ((a)>(b)) ? (a) : (b);

#ifdef min
#undef min
#endif
#define min(a,b) ((a)<(b)) ? (a) : (b);


LOCA::Stepper::Stepper(LOCA::Continuation::AbstractGroup& initialGuess, 
		       NOX::StatusTest::Generic& t,
		       NOX::Parameter::List& p) :
  LOCA::Abstract::Iterator(),
  conGroupManagerPtr(NULL),
  curGroupPtr(NULL),
  prevGroupPtr(NULL),
  statusTestPtr(NULL),
  paramListPtr(NULL),
  solverPtr(NULL),
  predictorManagerPtr(NULL),
  curPredictorPtr(NULL),
  prevPredictorPtr(NULL),
  stepSizeManagerPtr(NULL)
  
{
  reset(initialGuess, t, p);
}

LOCA::Stepper::Stepper(const LOCA::Stepper& s) :
  LOCA::Abstract::Iterator(s),
  conGroupManagerPtr(NULL),
  curGroupPtr(NULL),
  prevGroupPtr(NULL),
  statusTestPtr(s.statusTestPtr),
  paramListPtr(s.paramListPtr),
  solverPtr(NULL),
  predictorManagerPtr(NULL),
  curPredictorPtr(NULL),
  prevPredictorPtr(NULL),
  stepSizeManagerPtr(NULL),
  startValue(s.startValue),
  maxValue(s.maxValue),
  minValue(s.minValue),
  stepSize(s.stepSize),
  maxNonlinearSteps(s.maxNonlinearSteps),
  targetValue(s.targetValue),
  isTargetStep(s.isTargetStep),
  tangentFactor(s.tangentFactor),
  minTangentFactor(s.minTangentFactor),
  tangentFactorExponent(s.tangentFactorExponent),
  calcEigenvalues(s.calcEigenvalues)
{ 
  conGroupManagerPtr = 
    new LOCA::Continuation::Manager(*s.conGroupManagerPtr);
  curGroupPtr = 
    dynamic_cast<LOCA::Continuation::ExtendedGroup*>(s.curGroupPtr->clone());
  prevGroupPtr = 
    dynamic_cast<LOCA::Continuation::ExtendedGroup*>(s.prevGroupPtr->clone());
  predictorManagerPtr = 
    new LOCA::Predictor::Manager(*s.predictorManagerPtr);
  curPredictorPtr = 
    dynamic_cast<LOCA::Continuation::ExtendedVector*>(s.curPredictorPtr->clone());
  prevPredictorPtr = 
    dynamic_cast<LOCA::Continuation::ExtendedVector*>(s.prevPredictorPtr->clone());
  stepSizeManagerPtr = 
    new LOCA::StepSize::Manager(*s.stepSizeManagerPtr);

  // Right now this doesn't work because we can't copy the solver
}

LOCA::Stepper::~Stepper() 
{ 
  delete conGroupManagerPtr;
  delete curGroupPtr;
  delete prevGroupPtr;
  delete predictorManagerPtr;
  delete curPredictorPtr;
  delete prevPredictorPtr;
  delete stepSizeManagerPtr;
  delete solverPtr;
}

bool 
LOCA::Stepper::reset(LOCA::Continuation::AbstractGroup& initialGuess,
		     NOX::StatusTest::Generic& t,
		     NOX::Parameter::List& p) 
{
  delete curGroupPtr;
  delete prevGroupPtr;
  delete curPredictorPtr;
  delete prevPredictorPtr;
  delete conGroupManagerPtr;
  delete predictorManagerPtr;
  delete stepSizeManagerPtr;
  delete solverPtr;

  paramListPtr = &p;
  statusTestPtr = &t;

  // Initialize the utilities
  LOCA::Utils::setUtils(paramListPtr->sublist("LOCA").sublist("Utilities"));

  // Get LOCA sublist
  NOX::Parameter::List& locaList = paramListPtr->sublist("LOCA");

  // Get stepper sublist
  NOX::Parameter::List& stepperList = locaList.sublist("Stepper");

  // Reset base class
  LOCA::Abstract::Iterator::reset(stepperList);

  // Reset group, predictor, step-size managers
  conGroupManagerPtr = 
    new LOCA::Continuation::Manager(stepperList);
  predictorManagerPtr = 
    new LOCA::Predictor::Manager(locaList.sublist("Predictor"));
  stepSizeManagerPtr = 
    new LOCA::StepSize::Manager(locaList.sublist("Step Size"));

  // Get the continuation parameter starting value
  if (stepperList.isParameter("Initial Value"))
    startValue = stepperList.getParameter("Initial Value", 0.0);
  else {
    errorCheck.throwError("LOCA::Stepper::reset()",
				 "\"Initial Value\" of continuation parameter is not set!");
  }

  // Get the continuation parameter name
  if (stepperList.isParameter("Continuation Parameter"))
    initialGuess.setParam(stepperList.getParameter("Continuation Parameter", 
						   "None"), 
			  startValue);
  else {
     errorCheck.throwError("LOCA::Stepper::reset()",
				  "\"Continuation Parameter\" name is not set!");
  }
  
  // Get the max and min values of the continuation parameter
  if (stepperList.isParameter("Max Value"))
    maxValue = stepperList.getParameter("Max Value", 0.0);
  else {
     errorCheck.throwError("LOCA::Stepper::reset()",
				  "\"Maximum Value\" of continuation parameter is not set!");
  }
  if (stepperList.isParameter("Min Value"))
    minValue = stepperList.getParameter("Min Value", 0.0);
  else {
    errorCheck.throwError("LOCA::Stepper::reset()",
				 "\"Minimum Value\" of continuation parameter is not set!");
  }
  

  // Get the initial values or use their defaults
  stepSize = stepSizeManagerPtr->getStartStepSize();
  maxNonlinearSteps = stepperList.getParameter("Max Nonlinear Iterations", 15);

  targetValue = 0.0;
  isTargetStep = false;
  tangentFactor = 1.0;
  minTangentFactor = stepperList.getParameter("Min Tangent Factor",0.1);
  tangentFactorExponent = 
    stepperList.getParameter("Tangent Factor Exponent",1.0);
  calcEigenvalues = stepperList.getParameter("Compute Eigenvalues",false);

  // Make a copy of the parameter list, change continuation method to 
  // natural
  NOX::Parameter::List firstStepParams(*paramListPtr);
  NOX::Parameter::List& firstStepperParams 
      = firstStepParams.sublist("LOCA").sublist("Stepper");
  firstStepperParams.setParameter("Continuation Method", "Natural");

  // Reset continuation manager
  conGroupManagerPtr->reset(firstStepperParams);

  // Create continuation group
  curGroupPtr = conGroupManagerPtr->createContinuationGroup(initialGuess, firstStepParams.sublist("NOX").sublist("Direction").sublist("Newton").sublist("Linear Solver"));
      
  // Set step size
  curGroupPtr->setStepSize(0.0);
  
  // Set previous solution vector in current solution group
  curGroupPtr->setPrevX(curGroupPtr->getX());

  // Create solver using initial conditions
  solverPtr = new NOX::Solver::Manager(*curGroupPtr, *statusTestPtr, 
				       paramListPtr->sublist("NOX"));

  printInitializationInfo();

  if (LOCA::Utils::doPrint(LOCA::Utils::Parameters))
    paramListPtr->print(cout);
  
  return true;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Stepper::start() {
  NOX::StatusTest::StatusType solverStatus;

  printStartStep();

  // Perform solve of initial conditions
  solverStatus = solverPtr->solve();

  // Reset continuation manager
  conGroupManagerPtr->reset(paramListPtr->sublist("LOCA").sublist("Stepper"));

  // Set up continuation groups
  const LOCA::Continuation::ExtendedGroup& constSolnGrp = 
    dynamic_cast<const LOCA::Continuation::ExtendedGroup&>(solverPtr->getSolutionGroup());
  LOCA::Continuation::AbstractGroup& solnGrp = 
    const_cast<LOCA::Continuation::AbstractGroup&>(constSolnGrp.getUnderlyingGroup());
  delete curGroupPtr;
  curGroupPtr = 
    conGroupManagerPtr->createContinuationGroup(solnGrp, paramListPtr->sublist("NOX").sublist("Direction").sublist("Newton").sublist("Linear Solver"));
  
  // Do printing (stepNumber==0 case) after continuation group set up
  if (solverStatus == NOX::StatusTest::Failed) 
    printEndStep(LOCA::Abstract::Iterator::Unsuccessful);
  else
    printEndStep(LOCA::Abstract::Iterator::Successful);

  // Set the initial step size
  curGroupPtr->setStepSize(stepSize);

  prevGroupPtr = 
    dynamic_cast<LOCA::Continuation::ExtendedGroup*>(curGroupPtr->clone());

  // If nonlinear solve failed, return (this must be done after continuation 
  // groups are created so Stepper::getSolutionGroup() functions correctly.
  if (solverStatus != NOX::StatusTest::Converged)
    return LOCA::Abstract::Iterator::Failed;

  curGroupPtr->printSolution();

  // Initialize predictor direction
  curPredictorPtr = 
    dynamic_cast<LOCA::Continuation::ExtendedVector*>(curGroupPtr->getX().clone(NOX::ShapeCopy));

  // Compute predictor direction
  predictorManagerPtr->compute(*prevGroupPtr, *curGroupPtr, *curPredictorPtr);

  prevPredictorPtr = 
    dynamic_cast<LOCA::Continuation::ExtendedVector*>(curPredictorPtr->clone());

  // Create new solver using new continuation groups and combo status test
  delete solverPtr;
  solverPtr = new NOX::Solver::Manager(*curGroupPtr, *statusTestPtr, 
				       paramListPtr->sublist("NOX"));

  return LOCA::Abstract::Iterator::NotFinished;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Stepper::finish(LOCA::Abstract::Iterator::IteratorStatus iteratorStatus)
{

  //
  // We don't need to check if the last step was successful since finish
  // is never called if it wasn't.  We might want to change that if there is
  // some post processing we want to do even if the run failed.
  //

  // Copy last solution
  *curGroupPtr = solverPtr->getSolutionGroup();

  // Return if iteration failed (reached max number of steps)
  if (iteratorStatus == LOCA::Abstract::Iterator::Failed)
    return iteratorStatus;

  // Do one additional step using natural continuation to hit target value
  double value = curGroupPtr->getContinuationParameter();

  if (fabs(value-targetValue) > 1.0e-15*(1.0 + fabs(targetValue))) {
      
    isTargetStep = true;

    // Save previous successful step information
    *prevGroupPtr = *curGroupPtr;

    // Get underyling solution group
    //LOCA::Continuation::AbstractGroup& underlyingGroup 
    //  = dynamic_cast<LOCA::Continuation::AbstractGroup&>(prevGroupPtr->getUnderlyingGroup());
    LOCA::Continuation::AbstractGroup& underlyingGroup 
      = dynamic_cast<LOCA::Continuation::AbstractGroup&>(getSolutionGroup());

    // Make a copy of the parameter list, change continuation method to 
    // natural, predictor method to constant
    NOX::Parameter::List lastStepParams(*paramListPtr);
    NOX::Parameter::List& lastStepperParams 
      = lastStepParams.sublist("LOCA").sublist("Stepper");
    lastStepperParams.setParameter("Continuation Method", "Natural");
    lastStepperParams.sublist("Predictor").setParameter("Method", "Random");

    // Reset continuation manager
    //conGroupManagerPtr->reset(lastStepperParams);
    delete conGroupManagerPtr;
    conGroupManagerPtr = new LOCA::Continuation::Manager(lastStepperParams);
      
    // Reset predictor manager
    predictorManagerPtr->reset(lastStepperParams.sublist("Predictor"));

    // Get new continuation group
    delete curGroupPtr;

    curGroupPtr = conGroupManagerPtr->createContinuationGroup(underlyingGroup, lastStepParams.sublist("NOX").sublist("Direction").sublist("Newton").sublist("Linear Solver"));
      
    // Set step size
    stepSize = targetValue - value;
    curGroupPtr->setStepSize(stepSize);

    // Get predictor direction
    predictorManagerPtr->compute(*curGroupPtr, *curGroupPtr, *curPredictorPtr);
      
    // Set previous solution vector in current solution group
    curGroupPtr->setPrevX(curGroupPtr->getX());

    // Take step in predictor direction
    curGroupPtr->computeX(*curGroupPtr, *curPredictorPtr, stepSize);

    printStartStep();
      
    // Create new solver
    delete solverPtr;
    solverPtr = new NOX::Solver::Manager(*curGroupPtr, *statusTestPtr, 
					 lastStepParams.sublist("NOX"));

    // Solve step
    NOX::StatusTest::StatusType solverStatus = solverPtr->solve();

    // Get solution
    *curGroupPtr 
      = dynamic_cast<const LOCA::Continuation::ExtendedGroup&>(solverPtr->getSolutionGroup());

    if (solverStatus == NOX::StatusTest::Failed) {
      printEndStep(LOCA::Abstract::Iterator::Unsuccessful);
      return LOCA::Abstract::Iterator::Failed;
    }

    printEndStep(LOCA::Abstract::Iterator::Successful);

    curGroupPtr->printSolution();

  }

  return LOCA::Abstract::Iterator::Finished;
}

LOCA::Abstract::Iterator::StepStatus
LOCA::Stepper::preprocess(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful) {

    // Restore previous step information
    *curGroupPtr = *prevGroupPtr;
  }
  else {

    // Save previous successful step information
    *prevGroupPtr = *curGroupPtr;
  }
  
  // Compute step size
  stepStatus = computeStepSize(stepStatus, stepSize);

  // Set step size in current solution group
  curGroupPtr->setStepSize(stepSize);

  // Set previous solution vector in current solution group
  curGroupPtr->setPrevX(prevGroupPtr->getX());

  // Take step in predictor direction
  curGroupPtr->computeX(*prevGroupPtr, *curPredictorPtr, stepSize);

  // Reset solver to compute new solution
  solverPtr->reset(*curGroupPtr, *statusTestPtr, 
		   paramListPtr->sublist("NOX"));

  return stepStatus;
}
  
LOCA::Abstract::Iterator::StepStatus
LOCA::Stepper::compute(LOCA::Abstract::Iterator::StepStatus stepStatus) 
{
  NOX::StatusTest::StatusType solverStatus;

  printStartStep();

  solverStatus = solverPtr->solve();

  if (solverStatus == NOX::StatusTest::Failed) {
    printEndStep(LOCA::Abstract::Iterator::Unsuccessful);
    return LOCA::Abstract::Iterator::Unsuccessful;
  }

  // Copy solution out of solver
  *curGroupPtr = solverPtr->getSolutionGroup();

  printEndStep(LOCA::Abstract::Iterator::Successful);

  curGroupPtr->printSolution();

  return LOCA::Abstract::Iterator::Successful;
}

LOCA::Abstract::Iterator::StepStatus
LOCA::Stepper::postprocess(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful)
    return stepStatus;

  // Compute eigenvalues/eigenvectors
  if (calcEigenvalues) {
    curGroupPtr->getBaseLevelUnderlyingGroup().computeEigenvalues(*paramListPtr);
  }

  *prevPredictorPtr = *curPredictorPtr;

  predictorManagerPtr->compute(*prevGroupPtr, *curGroupPtr, *curPredictorPtr);

  if (getStepNumber() > 1) {
    tangentFactor = curGroupPtr->computeScaledDotProduct(*curPredictorPtr, 
							 *prevPredictorPtr) / 
      sqrt(curGroupPtr->computeScaledDotProduct(*curPredictorPtr, 
						*curPredictorPtr) * 
	   curGroupPtr->computeScaledDotProduct(*prevPredictorPtr, 
						 *prevPredictorPtr));

    if (tangentFactor < minTangentFactor)
      return LOCA::Abstract::Iterator::Unsuccessful;

  }

  return stepStatus;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Stepper::stop(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Successful) {
    
    double value = curGroupPtr->getContinuationParameter();
    double paramStep = value - prevGroupPtr->getContinuationParameter();

    // See if we went past bounds for parameter
    if ( value >= maxValue*(1.0 - 1.0e-15) && paramStep > 0 ) {
      if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
	cout << "\n\tContinuation run stopping: parameter reached bound of "
	     << LOCA::Utils::sci(maxValue) << endl;
      }
      targetValue = maxValue;
      return LOCA::Abstract::Iterator::Finished;
    }
    if ( value <= minValue*(1.0 + 1.0e-15) && paramStep < 0 ) {
      if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
	cout << "\n\tContinuation run stopping: parameter reached bound of "
	     << LOCA::Utils::sci(minValue) << endl;
      }
      targetValue = minValue;
      return LOCA::Abstract::Iterator::Finished;
    }

    // Check to see if arclength step was aimed to reach bound (should be near)
    if (isLastIteration()) {
      if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
	cout << "\n\tContinuation run stopping: parameter stepped to bound" 
	     << endl;
      }
      return LOCA::Abstract::Iterator::Finished;
    }
  }
  else if (isLastIteration())  // Failed step did not reach bounds as predicted
    return LOCA::Abstract::Iterator::NotFinished;

  // Check to see if max number of steps has been reached
  if (LOCA::Abstract::Iterator::numTotalSteps
        >= LOCA::Abstract::Iterator::maxSteps) {
    if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
      cout << "\n\tContinuation run stopping: reached maximum number of steps "
	   << LOCA::Abstract::Iterator::maxSteps << endl;
    }
    return LOCA::Abstract::Iterator::Failed;
  }

  return LOCA::Abstract::Iterator::NotFinished;
}

LOCA::Abstract::Iterator::StepStatus
LOCA::Stepper::computeStepSize(LOCA::Abstract::Iterator::StepStatus stepStatus,
			       double& stepSize)
{
  NOX::Abstract::Group::ReturnType res = 
    stepSizeManagerPtr->compute(*curGroupPtr, *curPredictorPtr, *solverPtr, 
				stepStatus, *this, stepSize);

  if (res == NOX::Abstract::Group::Failed)
    return LOCA::Abstract::Iterator::Provisional;

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails)) {
    cout << "\n\tTangent factor scaling:  Rescaling step size by "
	 << LOCA::Utils::sci(pow(fabs(tangentFactor), tangentFactorExponent))
	 << endl;
  }

  stepSize *= pow(fabs(tangentFactor), tangentFactorExponent);
  
  // Cap the con parameter so we don't go past bounds
  double prevValue = curGroupPtr->getContinuationParameter();
  double dpds = curPredictorPtr->getParam();
  if ( (prevValue+stepSize*dpds > maxValue*(1.0 - 1.0e-15)) ) {
    stepSize = (maxValue - prevValue)/dpds;
    targetValue = maxValue;
    setLastIteration();
  }
  if ( (prevValue+stepSize*dpds < minValue*(1.0 + 1.0e-15)) ) {
    stepSize = (minValue - prevValue)/dpds;
    targetValue = minValue;
    setLastIteration();
  }

  return LOCA::Abstract::Iterator::Successful;
}

LOCA::Continuation::AbstractGroup& 
LOCA::Stepper::getSolutionGroup()
{
  return curGroupPtr->getUnderlyingGroup();
}

const NOX::Parameter::List& 
LOCA::Stepper::getParameterList() const
{
  return *paramListPtr;
}

void 
LOCA::Stepper::printInitializationInfo()
{  
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
    cout << endl << LOCA::Utils::fill(72, '~') << endl;
    cout << "Beginning Continuation Run \n" 
	 << "Stepper Method:             " << conGroupManagerPtr->getMethod() 
	 << "\n"
	 << "Initial Parameter Value = " << LOCA::Utils::sci(startValue) 
	 << "\n"
	 << "Maximum Parameter Value = " << LOCA::Utils::sci(maxValue) << "\n"
	 << "Minimum Parameter Value = " << LOCA::Utils::sci(minValue) << "\n"
	 << "Maximum Number of Continuation Steps = " 
	 << LOCA::Abstract::Iterator::maxSteps 
	 << endl;
    cout << LOCA::Utils::fill(72, '~') << endl << endl;
  }
}
 
void 
LOCA::Stepper::printStartStep()
{  
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
    cout << "\n" << LOCA::Utils::fill(72, '~') << "\n";
    cout << "Start of Continuation Step " << stepNumber <<" : ";
    if (stepNumber==0) {
      cout << "Attempting to converge initial guess at initial parameter"
	   << "values." << endl;
    }
    else if (isTargetStep) {
      cout << "Attempting to hit final target value " 
	   << LOCA::Utils::sci(targetValue) << endl;
    }
    else {
      cout << "Parameter: " << conGroupManagerPtr->getConParamID()
  	   << " = " 
	   << LOCA::Utils::sci(curGroupPtr->getContinuationParameter())
           << " from " 
	   << LOCA::Utils::sci(prevGroupPtr->getContinuationParameter()) 
	   << endl;
      cout << "Continuation Method: " << conGroupManagerPtr->getMethod() 
	   << endl;
      cout << "Current step size  = " << LOCA::Utils::sci(stepSize) << "   "
	   << "Previous step size = " 
	   << LOCA::Utils::sci(stepSizeManagerPtr->getPrevStepSize()) << endl;
    }
    cout << LOCA::Utils::fill(72, '~') << "\n" << endl;
  }
}

void 
LOCA::Stepper::printEndStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Successful) {
    // Print results of successful continuation step
    if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
      cout << "\n" << LOCA::Utils::fill(72, '~') << "\n";
      cout << "End of Continuation Step " << stepNumber << " : ";
      cout << "Parameter: " << conGroupManagerPtr->getConParamID()
	   << " = " 
	   << LOCA::Utils::sci(curGroupPtr->getContinuationParameter());
      if (stepNumber != 0) 
        cout << " from " 
	     << LOCA::Utils::sci(prevGroupPtr->getContinuationParameter());
      cout << endl << "--> Step Converged in "
           << solverPtr->getNumIterations() 
	   <<" Nonlinear Solver Iterations!\n";
      cout << LOCA::Utils::fill(72, '~') << "\n" << endl;
    }
  }
  else {
    if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
      // RPP: We may not need this, the failure info should be 
      // at the method level!
      cout << endl << LOCA::Utils::fill(72, '~') << endl;
      cout << "Continuation Step Number " << stepNumber 
           << " experienced a convergence failure in\n"
           << "the nonlinear solver after "<< solverPtr->getNumIterations() 
	   <<" Iterations\n";
      cout << "Value of continuation parameter at failed step = "
           << LOCA::Utils::sci(curGroupPtr->getContinuationParameter());
      if (stepNumber != 0) 
        cout << " from " 
	     << LOCA::Utils::sci(prevGroupPtr->getContinuationParameter());
      cout << endl << LOCA::Utils::fill(72, '~') << endl;
    }
  }
}

void 
LOCA::Stepper::printEndInfo()
{

}

