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

// LOCA Includes
#include "LOCA_Utils.H"		      // for static function doPrint
#include "LOCA_Abstract_Group.H"      // class data element
#include "LOCA_Continuation_Group.H"
#include "LOCA_Continuation_NaturalGroup.H"
#include "LOCA_Continuation_ArcLengthGroup.H"   //

using namespace LOCA;
using namespace NOX::StatusTest;

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


Stepper::Stepper(LOCA::Abstract::Group& initialGuess, 
		 NOX::StatusTest::Generic& t,
		 NOX::Parameter::List& p,
		 LOCA::Abstract::DataOutput& dataOut) :
  conGroupManager(p.sublist("Stepper")),
  curGroupPtr(NULL),
  prevGroupPtr(NULL),
  statusTestPtr(&t),
  solverPtr(NULL),
  predictor(p.sublist("Predictor")),
  predictorDirection(NULL),
  stepSizeManager(p.sublist("Step Size")),
  dataOutput(dataOut),
  paramList()
{
  // Initialize the utilities
  Utils::setUtils(paramList.sublist("Utilities"));
  paramList = p;
  init(initialGuess);
}

Stepper::Stepper(const Stepper& s) :
  dataOutput(s.dataOutput),
  solverPtr(s.solverPtr),
  conGroupManager(s.conGroupManager),
  predictor(s.predictor),
  stepSizeManager(s.stepSizeManager)
{ //Copy constructor declared private since it is not expected to be neeeded.
}

Stepper::~Stepper() 
{ 
  delete curGroupPtr;
  delete prevGroupPtr;
  delete predictorDirection;
  delete solverPtr;
}

bool Stepper::reset(LOCA::Abstract::Group& initialGuess,
		    NOX::StatusTest::Generic& t,
		    const NOX::Parameter::List& p) 
{
  delete curGroupPtr;
  delete prevGroupPtr;
  delete predictorDirection;
  statusTestPtr = &t;
  paramList = p;

  return init(initialGuess);
}

StatusType Stepper::getStatus()
{
  return stepperStatus;
}

StatusType Stepper::solve()
{
  // Perform solve of initial conditions
  solverStatus = Unconverged;
  solverStatus = solverPtr->solve();

  stepperStatus = checkStepperStatus();

  // Set up continuation groups
  const LOCA::Abstract::Group& solnGrp = 
    dynamic_cast<const LOCA::Abstract::Group&>(solverPtr->getSolutionGroup());
  curGroupPtr = 
    conGroupManager.createContinuationGroup(solnGrp, paramList.sublist("Solver").sublist("Direction").sublist("Linear Solver"));
  prevGroupPtr = 
    dynamic_cast<LOCA::Continuation::Group*>(curGroupPtr->clone());

  // If nonlinear solve failed, return (this must be done after continuation 
  // groups are created so Stepper::getSolutionGroup() functions correctly.
  if (stepperStatus != Unconverged) 
    return stepperStatus;

  stepNumber = 1;
  numTotalSteps = 1;
  dataOutput.saveGroupData(getSolutionGroup());

  // Initialize predictor direction
  predictorDirection = 
    dynamic_cast<LOCA::Continuation::Vector*>(curGroupPtr->getX().clone(NOX::ShapeCopy));

  // Create new solver using new continauation groups
  delete solverPtr;
  solverPtr = new NOX::Solver::Manager(*curGroupPtr, *statusTestPtr, 
				       paramList.sublist("Solver"));

  while (stepperStatus == Unconverged) {
    stepperStatus = step();
  }

  // Get last solution
  *curGroupPtr 
    = dynamic_cast<const LOCA::Continuation::Group&>(solverPtr->getSolutionGroup());

  return stepperStatus;
}

StatusType Stepper::step()
{
  stepperStatus = Unconverged;
 
  if (solverStatus == Failed) {

    *curGroupPtr = *prevGroupPtr;

    stepSize = computeStepSize(Failed);
    curGroupPtr->setStepSize(stepSize);

    // Set previous solution vector in current solution group
    curGroupPtr->setPrevX(prevGroupPtr->getX());

    // Take (reduced) step in predictor direction, which must be still valid
    curGroupPtr->computeX(*prevGroupPtr, *predictorDirection, stepSize);

  }
  else {

    // Set previous solution vector in current solution group
    if (stepNumber > 1)
      curGroupPtr->setPrevX(prevGroupPtr->getX());

    // Compute predictor direction
    predictor.compute(*prevGroupPtr, *curGroupPtr, *predictorDirection);

    // Compute step size
    stepSize = computeStepSize(solverStatus);

    //if (stepSize > 1.5e6)
    //  stepSize /= 2.0;

    curGroupPtr->setStepSize(stepSize);

    // Save previous successful step information
    *prevGroupPtr = *curGroupPtr;

    // Set previous solution vector in current solution group
    curGroupPtr->setPrevX(prevGroupPtr->getX());

    // Take step in predictor direction
    curGroupPtr->computeX(*prevGroupPtr, *predictorDirection, stepSize);
  }

  solverPtr->reset(*curGroupPtr, *statusTestPtr, paramList.sublist("Solver"));

  return nonlinearSolve();
}

StatusType Stepper::nonlinearSolve()
{
  printStartStep();
  
  solverStatus = Unconverged;
  solverStatus = solverPtr->solve();
  
  if (solverStatus == Failed) {
    numFailedSteps += 1;
    printEndStep(solverStatus);
  }
  else  {

    dataOutput.saveGroupData(getSolutionGroup());

    // Get solution
    *curGroupPtr 
       = dynamic_cast<const LOCA::Continuation::Group&>(solverPtr->getSolutionGroup());
    printEndStep(solverStatus);

    stepNumber += 1;

    // Recalculate scale factor
    //curGroupPtr->recalculateScaleFactor();
  }
  
  // See if we went past bounds for parameter
  if ( curGroupPtr->getContinuationParameter() >= maxValue*(1.0 - 1.0e-15) || 
       curGroupPtr->getContinuationParameter() <= minValue*(1.0 + 1.0e-15) )
    isLastStep = true;
  
  numTotalSteps += 1;

  stepperStatus = checkStepperStatus();

  return stepperStatus;
}

double Stepper::computeStepSize(StatusType solverStatus)
{
  NOX::Abstract::Group::ReturnType res = 
    stepSizeManager.compute(*curGroupPtr, *predictorDirection, *solverPtr, 
			    solverStatus, *this, stepSize);

  if (res == NOX::Abstract::Group::Failed)
    stepperStatus = Failed;
  
  // Cap the con parameter so we don't go past bounds
  double prevValue = curGroupPtr->getContinuationParameter();
  double dpds = predictorDirection->getParam();
  if ( (prevValue+stepSize*dpds > maxValue*(1.0 - 1.0e-15)) ) {
    stepSize = (maxValue - prevValue)/dpds;
    isLastStep = true;
  }
  if ( (prevValue+stepSize*dpds < minValue*(1.0 + 1.0e-15)) ) {
    stepSize = (minValue - prevValue)/dpds;
    isLastStep = true;
  }

  return stepSize;
}

const LOCA::Abstract::Group& Stepper::getSolutionGroup() const
{
  return curGroupPtr->getGroup();
}

const NOX::Parameter::List& Stepper::getParameterList() const
{
  return paramList;
}

bool Stepper::init(LOCA::Abstract::Group& initialGuess)
{ 
  // Get a reference to the parameter stepper sublist from the Solver.
  const NOX::Parameter::List& p = paramList.sublist("Stepper");

  // Get the continuation parameter starting value
  if (p.isParameter("Initial Value"))
    startValue = p.getParameter("Initial Value", 0.0);
  else {
    cout << "ERROR: LOCA::Stepper::Stepper::resetStepperMembers() - "
	 << "\"Initial Value\" of continuation param is not set!" << endl;
    throw "LOCA Error";
  }
  
  // Get the max and min values of the continuation parameter
  if (p.isParameter("Max Value"))
    maxValue = p.getParameter("Max Value", 0.0);
  else {
    cout << "ERROR: LOCA::Stepper::Stepper::resetStepperMembers() - "
	 << "\"Maximum Value\" of continuation param is not set!" << endl;
    throw "LOCA Error";
  }
  if (p.isParameter("Min Value"))
    minValue = p.getParameter("Min Value", 0.0);
  else {
    cout << "ERROR: LOCA::Stepper::Stepper::resetStepperMembers() - "
	 << "\"Minimum Value\" of continuation param is not set!" << endl;
    throw "LOCA Error";
  }
  

  // Get the initial values or use their defaults
  stepSize = 0.0;
  stepNumber = 0;                    
  numFailedSteps = 0;
  numTotalSteps = 0;
  maxNonlinearSteps = p.getParameter("Max Nonlinear Iterations", 15);
  maxConSteps = p.getParameter("Max Continuation Steps", 100);

  stepperStatus = Unconverged;
  solverStatus = Unconverged;
  isLastStep = false;

  // Create solver using initial conditions
  solverPtr = new NOX::Solver::Manager(initialGuess, *statusTestPtr, 
				       paramList.sublist("Solver"));

  printInitializationInfo();

  //  if (Utils::doPrint(Utils::Parameters))
    paramList.print(cout);

  return true;
}

int Stepper::getNumContinuationSteps() const
{
  return stepNumber;
}

int Stepper::getNumFailedSteps() const
{
  return numFailedSteps;
}

int Stepper::getNumTotalSteps() const
{
  return numTotalSteps;
}

void Stepper::printInitializationInfo()
{  
  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << endl << Utils::fill(72, '~') << endl;
    cout << "Beginning Continuation Run \n" 
	 << "Stepper Method:             " << conGroupManager.getMethod() << "\n"
	 << "Initial Parameter Value = " << startValue << "\n"
	 << "Maximum Parameter Value = " << maxValue << "\n"
	 << "Minimum Parameter Value = " << minValue << "\n"
	 << "Maximum Number of Continuation Steps = " << maxConSteps 
	 << endl;
    cout << Utils::fill(72, '~') << endl << endl;
  }
}
 
void Stepper::printStartStep()
{  
  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << "\n" << Utils::fill(72, '~') << "\n";
    cout << "Start of Continuation Step " << stepNumber << endl;
    cout << "Continuation Method: " << conGroupManager.getMethod() << endl;
    cout << "Continuation Parameter: " << conGroupManager.getConParamID()
	 << " = " << curGroupPtr->getContinuationParameter() << " from " << prevGroupPtr->getContinuationParameter() << endl;
    cout << "Current step size  = " << stepSize << "   "
	 << "Previous step size = " << stepSizeManager.getPrevStepSize() << endl;
    cout << Utils::fill(72, '~') << "\n" << endl;
  }
}

void Stepper::printEndStep(StatusType& solverStatus)
{
  if (solverStatus != Failed) {
    // Print results of successful continuation step
    if (Utils::doPrint(Utils::StepperIteration)) {
      cout << "\n" << Utils::fill(72, '~') << "\n";
      cout << "End of  Continuation Step " << stepNumber << endl;
      cout << "Continuation Parameter: " << conGroupManager.getConParamID()
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

void Stepper::printEndInfo()
{

}

StatusType Stepper::checkStepperStatus()
{
  if ((solverStatus != Converged) && (stepperStatus == Failed))
    return Failed;

  // Check to see if we hit final parameter value
  if ((solverStatus == Converged) && isLastStep)
    return Converged;

  // Check to see if max number of steps has been reached
  if (numTotalSteps == maxConSteps) 
    return Failed;

  // Check to see if the initial solve (step=0) failed
  // if (stepNumber == 0) 
  //  return Failed;

  return Unconverged;
}
