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

// LOCA Includes
#include "LOCA_Utils.H"		      // for static function doPrint
#include "LOCA_Abstract_Group.H"      // class data element
#include "LOCA_Continuation_Group.H"
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


LOCA::Stepper::Stepper(LOCA::Abstract::Group& initialGuess, 
		       NOX::StatusTest::Generic& t,
		       NOX::Parameter::List& p) :
  LOCA::Abstract::Iterator(),
  conGroupManagerPtr(NULL),
  curGroupPtr(NULL),
  prevGroupPtr(NULL),
  statusTestPtr(NULL),
  solverPtr(NULL),
  predictorManagerPtr(NULL),
  curPredictorPtr(NULL),
  prevPredictorPtr(NULL),
  stepSizeManagerPtr(NULL),
  paramList()
{
  // Initialize the utilities
  Utils::setUtils(p.sublist("LOCA").sublist("Utilities"));
  reset(initialGuess, t, p);
}

LOCA::Stepper::Stepper(const LOCA::Stepper& s) :
  LOCA::Abstract::Iterator(s),
  conGroupManagerPtr(NULL),
  curGroupPtr(NULL),
  prevGroupPtr(NULL),
  statusTestPtr(s.statusTestPtr),
  paramList(s.paramList),
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
  isLastStep(s.isLastStep),
  tangentFactor(s.tangentFactor),
  minTangentFactor(s.minTangentFactor),
  tangentFactorExponent(s.tangentFactorExponent)
{ 
  conGroupManagerPtr = 
    new LOCA::Continuation::Manager(*s.conGroupManagerPtr);
  curGroupPtr = 
    dynamic_cast<LOCA::Continuation::Group*>(s.curGroupPtr->clone());
  prevGroupPtr = 
    dynamic_cast<LOCA::Continuation::Group*>(s.prevGroupPtr->clone());
  predictorManagerPtr = 
    new LOCA::Predictor::Manager(*s.predictorManagerPtr);
  curPredictorPtr = 
    dynamic_cast<LOCA::Continuation::Vector*>(s.curPredictorPtr->clone());
  prevPredictorPtr = 
    dynamic_cast<LOCA::Continuation::Vector*>(s.prevPredictorPtr->clone());
  stepSizeManagerPtr = 
    new LOCA::StepSize::Manager(*s.stepSizeManagerPtr);
  // Right now this doesn't work because we can't copy the solver
}

LOCA::Stepper::~Stepper() 
{ 
  paramList.print(cout);

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
LOCA::Stepper::reset(LOCA::Abstract::Group& initialGuess,
		     NOX::StatusTest::Generic& t,
		     NOX::Parameter::List& p) 
{
  delete curGroupPtr;
  delete prevGroupPtr;
  delete curPredictorPtr;
  delete prevPredictorPtr;

  statusTestPtr = &t;
  paramList = p;

  // Get LOCA sublist
  NOX::Parameter::List& locaList = paramList.sublist("LOCA");

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
    cout << "ERROR: LOCA::Stepper::Stepper::resetStepperMembers() - "
	 << "\"Initial Value\" of continuation param is not set!" << endl;
    throw "LOCA Error";
  }
  
  // Get the max and min values of the continuation parameter
  if (stepperList.isParameter("Max Value"))
    maxValue = stepperList.getParameter("Max Value", 0.0);
  else {
    cout << "ERROR: LOCA::Stepper::Stepper::resetStepperMembers() - "
	 << "\"Maximum Value\" of continuation param is not set!" << endl;
    throw "LOCA Error";
  }
  if (stepperList.isParameter("Min Value"))
    minValue = stepperList.getParameter("Min Value", 0.0);
  else {
    cout << "ERROR: LOCA::Stepper::Stepper::resetStepperMembers() - "
	 << "\"Minimum Value\" of continuation param is not set!" << endl;
    throw "LOCA Error";
  }
  

  // Get the initial values or use their defaults
  stepSize = stepSizeManagerPtr->getStartStepSize();
  maxNonlinearSteps = stepperList.getParameter("Max Nonlinear Iterations", 15);

  isLastStep = false;
  tangentFactor = 1.0;
  minTangentFactor = stepperList.getParameter("Min Tangent Factor",0.1);
  tangentFactorExponent = 
    stepperList.getParameter("Tangent Factor Exponent",1.0);

  // Create solver using initial conditions
  solverPtr = new NOX::Solver::Manager(initialGuess, *statusTestPtr, 
				       paramList.sublist("NOX"));

  printInitializationInfo();

  //  if (Utils::doPrint(Utils::Parameters))
    paramList.print(cout);

  return true;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Stepper::start() {
  NOX::StatusTest::StatusType solverStatus;

  // Perform solve of initial conditions
  solverStatus = solverPtr->solve();

  // Set up continuation groups
  const LOCA::Abstract::Group& solnGrp = 
    dynamic_cast<const LOCA::Abstract::Group&>(solverPtr->getSolutionGroup());
  curGroupPtr = 
    conGroupManagerPtr->createContinuationGroup(solnGrp, paramList.sublist("NOX").sublist("Direction").sublist("Linear Solver"));

  // Set the initial step size
  curGroupPtr->setStepSize(stepSize);

  prevGroupPtr = 
    dynamic_cast<LOCA::Continuation::Group*>(curGroupPtr->clone());

  // If nonlinear solve failed, return (this must be done after continuation 
  // groups are created so Stepper::getSolutionGroup() functions correctly.
  if (solverStatus != NOX::StatusTest::Converged)
    return LOCA::Abstract::Iterator::Failed;

  curGroupPtr->printSolution();

  // Initialize predictor direction
  curPredictorPtr = 
    dynamic_cast<LOCA::Continuation::Vector*>(curGroupPtr->getX().clone(NOX::ShapeCopy));

  // Compute predictor direction
  predictorManagerPtr->compute(*prevGroupPtr, *curGroupPtr, *curPredictorPtr);

  prevPredictorPtr = 
    dynamic_cast<LOCA::Continuation::Vector*>(curPredictorPtr->clone());

  // Create new solver using new continuation groups
  delete solverPtr;
  solverPtr = new NOX::Solver::Manager(*curGroupPtr, *statusTestPtr, 
				       paramList.sublist("NOX"));

  return LOCA::Abstract::Iterator::NotFinished;
}

LOCA::Abstract::Iterator::IteratorStatus
LOCA::Stepper::finish()
{
  // Get last solution
  *curGroupPtr 
    = dynamic_cast<const LOCA::Continuation::Group&>(solverPtr->getSolutionGroup());

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
  solverPtr->reset(*curGroupPtr, *statusTestPtr, paramList.sublist("NOX"));

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

  // Get solution
  *curGroupPtr 
    = dynamic_cast<const LOCA::Continuation::Group&>(solverPtr->getSolutionGroup());

  printEndStep(LOCA::Abstract::Iterator::Successful);

  curGroupPtr->printSolution();

  return LOCA::Abstract::Iterator::Successful;
}

LOCA::Abstract::Iterator::StepStatus
LOCA::Stepper::postprocess(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful)
    return stepStatus;

  *prevPredictorPtr = *curPredictorPtr;

  predictorManagerPtr->compute(*prevGroupPtr, *curGroupPtr, *curPredictorPtr);

  if (getStepNumber() > 0) {
    //tangentFactor = curPredictorPtr->dot(*prevPredictorPtr) 
    //  / (curPredictorPtr->norm() * prevPredictorPtr->norm());

    //tangentFactor = curGroupPtr->scaledDotProduct(*curPredictorPtr, 
    //						*prevPredictorPtr) 
    // / (curPredictorPtr->norm() * prevPredictorPtr->norm());

    tangentFactor = curGroupPtr->scaledDotProduct(*curPredictorPtr, 
						  *prevPredictorPtr) / 
      sqrt(curGroupPtr->scaledDotProduct(*curPredictorPtr, *curPredictorPtr) * 
	   curGroupPtr->scaledDotProduct(*prevPredictorPtr, *prevPredictorPtr));

    cout << "LOCA::Stepper::postprocess():  tangentFactor = "
	 << tangentFactor << endl;

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

    // See if we went past bounds for parameter
    if ( value >= maxValue*(1.0 - 1.0e-15) || 
	 value <= minValue*(1.0 + 1.0e-15) )
    return LOCA::Abstract::Iterator::Finished;

    // Check to see if we hit final parameter value
    if (isLastStep)
      return LOCA::Abstract::Iterator::Finished;

  }

  // Check to see if max number of steps has been reached
  if (LOCA::Abstract::Iterator::numTotalSteps >= LOCA::Abstract::Iterator::maxSteps) 
    return LOCA::Abstract::Iterator::Failed;

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
    return LOCA::Abstract::Iterator::Unsuccessful;

  stepSize *= pow(tangentFactor, tangentFactorExponent);
  
  // Cap the con parameter so we don't go past bounds
  double prevValue = curGroupPtr->getContinuationParameter();
  double dpds = curPredictorPtr->getParam();
  if ( (prevValue+stepSize*dpds > maxValue*(1.0 - 1.0e-15)) ) {
    stepSize = (maxValue - prevValue)/dpds;
    isLastStep = true;
  }
  if ( (prevValue+stepSize*dpds < minValue*(1.0 + 1.0e-15)) ) {
    stepSize = (minValue - prevValue)/dpds;
    isLastStep = true;
  }

  return LOCA::Abstract::Iterator::Successful;
}

const LOCA::Abstract::Group& 
LOCA::Stepper::getSolutionGroup() const
{
  return curGroupPtr->getGroup();
}

const NOX::Parameter::List& 
LOCA::Stepper::getParameterList() const
{
  return paramList;
}

void 
LOCA::Stepper::printInitializationInfo()
{  
  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << endl << Utils::fill(72, '~') << endl;
    cout << "Beginning Continuation Run \n" 
	 << "Stepper Method:             " << conGroupManagerPtr->getMethod() << "\n"
	 << "Initial Parameter Value = " << startValue << "\n"
	 << "Maximum Parameter Value = " << maxValue << "\n"
	 << "Minimum Parameter Value = " << minValue << "\n"
	 << "Maximum Number of Continuation Steps = " << LOCA::Abstract::Iterator::maxSteps 
	 << endl;
    cout << Utils::fill(72, '~') << endl << endl;
  }
}
 
void 
LOCA::Stepper::printStartStep()
{  
  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << "\n" << Utils::fill(72, '~') << "\n";
    cout << "Start of Continuation Step " << stepNumber << endl;
    cout << "Continuation Method: " << conGroupManagerPtr->getMethod() << endl;
    cout << "Continuation Parameter: " << conGroupManagerPtr->getConParamID()
	 << " = " << curGroupPtr->getContinuationParameter() << " from " << prevGroupPtr->getContinuationParameter() << endl;
    cout << "Current step size  = " << stepSize << "   "
	 << "Previous step size = " << stepSizeManagerPtr->getPrevStepSize() << endl;
    cout << Utils::fill(72, '~') << "\n" << endl;
  }
}

void 
LOCA::Stepper::printEndStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Successful) {
    // Print results of successful continuation step
    if (Utils::doPrint(Utils::StepperIteration)) {
      cout << "\n" << Utils::fill(72, '~') << "\n";
      cout << "End of  Continuation Step " << stepNumber << endl;
      cout << "Continuation Parameter: " << conGroupManagerPtr->getConParamID()
	   << " = " << curGroupPtr->getContinuationParameter() << " from " << prevGroupPtr->getContinuationParameter() << endl;
      cout << "--> Step Converged in "
           << solverPtr->getNumIterations() <<" Nonlinear Solver Iterations!\n";
      cout << Utils::fill(72, '~') << "\n" << endl;
    }
  }
  else {
    if (Utils::doPrint(Utils::StepperIteration)) {
      // RPP: We may not need this, the failure info should be 
      // at the method level!
      cout << endl << Utils::fill(72, '~') << endl;
      cout << "Continuation Step Number " << stepNumber 
           << " experienced a convergence failure in\n"
           << "the nonlinear solver after "<< solverPtr->getNumIterations() <<" Iterations\n";
      cout << "Value of continuation parameter at failed step = " << curGroupPtr->getContinuationParameter()
	   << " from " << prevGroupPtr->getContinuationParameter() << endl;
      cout << Utils::fill(72, '~') << endl;
    }
  }
}

void 
LOCA::Stepper::printEndInfo()
{

}

