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
#include "LOCA_Abstract_Vector.H"     // class data element

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
		 const NOX::Parameter::List& p,
		 LOCA::Abstract::DataOutput& dataOut) :
  curGroupPtr(&initialGuess),
  prevGroupPtr(dynamic_cast<LOCA::Abstract::Group*>(initialGuess.clone())),
  statusTestPtr(&t),
  paramList(p),
  conParams(initialGuess.getParams()),
  solver(initialGuess, t, paramList.sublist("Solver")),
  dataOutput(dataOut)
{
  // Initialize the utilities
  Utils::setUtils(paramList.sublist("Utilities"));

  init();
}

Stepper::~Stepper() 
{ 
  delete prevGroupPtr;
}

bool Stepper::reset(LOCA::Abstract::Group& initialGuess,
		    NOX::StatusTest::Generic& t,
		    const NOX::Parameter::List& p) 
{
  curGroupPtr = &initialGuess;
  delete prevGroupPtr;
  prevGroupPtr = 0;
  prevGroupPtr = dynamic_cast<LOCA::Abstract::Group*>(initialGuess.clone());
  statusTestPtr = &t;
  paramList = p;
  conParams = initialGuess.getParams();

  return init();
}

StatusType Stepper::getStatus()
{
  return stepperStatus;
}

StatusType Stepper::solve()
{
  while (stepperStatus == Unconverged) {
    stepperStatus = step();
  }
  return stepperStatus;
}

StatusType Stepper::step()
{
  stepperStatus = Unconverged;

  if (stepNumber != 0) {
 
    if (solverStatus == Failed) {

      *curGroupPtr = *prevGroupPtr;

      curStepSize = computeStepSize(Failed);
      curValue = prevValue + curStepSize;

    }
    else {
      // Set initial guess for next step equal to last step's final solution 
      *curGroupPtr = dynamic_cast<const LOCA::Abstract::Group&>(solver.getSolutionGroup());

      prevStepSize = curStepSize;
      prevValue = curValue;
      curStepSize = computeStepSize(solverStatus);

      // Cap the con parameter so we don't go past the final value
    
      if ( (prevValue+curStepSize-finalValue)*(prevValue-finalValue) < 0.0)
        curStepSize = finalValue - prevValue;

      curValue += curStepSize;

      // PREDICTOR
    }


    conParams.setValue(conParamID, curValue);
    curGroupPtr->setParams(conParams);

    // is this computeF needed?
    curGroupPtr->computeF();
    solver.reset(*curGroupPtr, *statusTestPtr, paramList.sublist("Solver"));
    
  }      
  
  printStartStep();
  
  solverStatus = Unconverged;
  solverStatus = solver.solve();
  
  printEndStep(solverStatus);
  
  if (solverStatus == Failed) {
    numFailedSteps += 1;
  }
  else  {
    stepNumber += 1;
  }

  if (solverStatus != Failed) {
    dataOutput.saveGroupData(dynamic_cast<const LOCA::Abstract::Group&>(solver.getSolutionGroup()));
  }
  
  numTotalSteps += 1;

  stepperStatus = checkStepperStatus();

  return stepperStatus;
}

const Abstract::Group& Stepper::getSolutionGroup() const
{
  return dynamic_cast<const LOCA::Abstract::Group&>(solver.getSolutionGroup());
}

const Abstract::Group& Stepper::getPreviousSolutionGroup() const
{
  return *prevGroupPtr;
}

const NOX::Parameter::List& Stepper::getParameterList() const
{
  return paramList;
}

bool Stepper::init()
{ 
  // Get a reference to the parameter stepper sublist from the Solver.
  const NOX::Parameter::List& p = paramList.sublist("Stepper");
  
  // Get the stepper method
  if (p.isParameter("Stepper Method"))
    stepperMethod = p.getParameter("Stepper Method", "???");
  else {
    cout << "ERROR: LOCA::Stepper::Stepper::resetStepperMembers() - "
	 << "\"Stepper Method\" is not set in the input parameter list!" 
	 << endl;
    throw "LOCA Error";
  }
 
  // Get the solver method
  
  solverMethod = "<NOT IMPLEMENTED YET>";

  // Get the continuation parameter
  conParamID = p.getParameter("Continuation Parameter", "???");
  if (conParamID == "???") {
    cout << "ERROR: LOCA::Stepper::Stepper::resetStepperMembers() - "
	 << "\"Continuation Parameter\" is not set" << endl;
    throw "LOCA Error";
  }

  // Get the continuation parameter starting value
  if (p.isParameter("Initial Value"))
    startValue = p.getParameter("Initial Value", 0.0);
  else {
    cout << "ERROR: LOCA::Stepper::Stepper::resetStepperMembers() - "
	 << "\"Initial Value\" of continuation param is not set!" << endl;
    throw "LOCA Error";
  }
  
  // Get the final value of the continuation parameter
  if (p.isParameter("Final Value"))
    finalValue = p.getParameter("Final Value", 0.0);
  else {
    cout << "ERROR: LOCA::Stepper::Stepper::resetStepperMembers() - "
	 << "\"Final Value\" of continuation param is not set!" << endl;
    throw "LOCA Error";
  }

  // First step of a new continuation run - set con param to starting value
  curValue = startValue;
  prevValue = startValue;

  // Get the initial values or use their defaults
  startStepSize = p.getParameter("Initial Step Size", 1.0);
  minStepSize = p.getParameter("Min Step Size", 1.0e-12);
  maxStepSize = p.getParameter("Max Step Size", 1.0e+12);
  curStepSize = startStepSize;
  prevStepSize = 0.0;
  stepNumber = 0;                    
  numFailedSteps = 0;
  numTotalSteps = 0;
  agrValue = p.getParameter("Step Size Aggressiveness", 1.0e+6);
  maxNonlinearSteps = p.getParameter("Max Nonlinear Iterations", 15);
  maxConSteps = p.getParameter("Max Continuation Steps", 100);
  stepperStatus = Unconverged;
  solverStatus = Unconverged;

  printInitializationInfo();

  if (Utils::doPrint(Utils::Parameters))
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
	 << "Stepper Method:             " << stepperMethod << "\n"
	 << "Solver Method:              " << solverMethod << "\n\n"
	 << "Initial Parameter Value = " << startValue << "\n"
	 << "Final Parameter Value = " << finalValue << "\n"
	 << "Initial Step Size = " << startStepSize << "\n" 
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
    cout << "Continuation Method: " << stepperMethod << endl;
    cout << "Continuation Parameter: " << conParamID
	 << " = " << curValue << " from " << prevValue << endl;
    cout << "Current step size  = " << curStepSize << "   "
	 << "Previous step size = " << prevStepSize << endl;
    cout << "\n" << Utils::fill(72, '~') << "\n" << endl;
  }
}

void Stepper::printEndStep(StatusType& solverStatus)
{
  if (solverStatus != Failed) {
    // Print results of successful continuation step
    if (Utils::doPrint(Utils::StepperIteration)) {
      cout << "\n" << Utils::fill(72, '~') << "\n";
      cout << "End of  Continuation Step " << stepNumber << endl;
      cout << "Continuation Parameter: " << conParamID
	   << " = " << curValue << " from " << prevValue << endl;
      cout << "--> Step Converged!";
      cout << "\n" << Utils::fill(72, '~') << "\n" << endl;
    }
  }
  else {
    if (Utils::doPrint(Utils::StepperIteration)) {
      // RPP: We may not need this, the failure info should be 
      // at the method level!
      cout << endl << Utils::fill(72, '~') << endl;
      cout << "Continuation Step Number " << stepNumber << " experienced a "
	   << "convergence failure in the nonlinear solver!" << endl;
      cout << "Value of continuation parameter at failed step = " << curValue
	   << " from " << prevValue << endl;
      cout << endl << Utils::fill(72, '~') << endl;
    }
  }
}

void Stepper::printEndInfo()
{

}


double Stepper::computeStepSize(StatusType solverStatus)
{
  double tmpStepSize = curStepSize;
  double predStepSize = curStepSize;

  if ((solverStatus == Failed) || 
      (solverStatus == Unconverged)) {

    // A failed nonlinear solve cuts the current step size in half
    tmpStepSize = curStepSize * 0.5;    

  }
  else if (stepNumber > 1) {

    // adapive step size control
    if (agrValue != 0.0) {

      // Number of nonlinear iterations to reach convergence for last 
      // nonlinear solve -- used to pick next step size
      double numNonlinearSteps = ((double) solver.getNumIterations());

      double factor = ((double) maxNonlinearSteps - numNonlinearSteps)
                      / ((double) maxNonlinearSteps - 1.0);

      tmpStepSize = curStepSize * (1.0 + agrValue * factor * factor);
  
    }
    // if constant step size (agrValue = 0.0), the step size may still be 
    // reduced by a solver failure.  We should then slowly bring the step 
    // size back towards its constant value using agrValue = 0.5.
    else{
      if (curStepSize != startStepSize) {  
	double numNonlinearSteps = ((double) solver.getNumIterations());

        double factor = ((double) maxNonlinearSteps - numNonlinearSteps)
                        / ((double) maxNonlinearSteps - 1.0);

        tmpStepSize = curStepSize * (1.0 + 0.5 * factor * factor);

        if (startStepSize > 0.0) {
          tmpStepSize = min(tmpStepSize, startStepSize);
	}
        else {
	  tmpStepSize = max(tmpStepSize, startStepSize);
	}
      }
    }
  }  
  predStepSize = tmpStepSize;

  // Clip the step size if above the bounds
  if (fabs(tmpStepSize) > maxStepSize) {
     predStepSize = maxStepSize;
     if (tmpStepSize < 0.0) predStepSize *= -1.0;
  }

  // Abort run if step size below bounds
  if (fabs(tmpStepSize) < minStepSize) {
    stepperStatus = Failed;
    predStepSize =  minStepSize;
    if (tmpStepSize < 0.0) predStepSize *= -1.0;
  }
  
  return predStepSize;
}

StatusType Stepper::checkStepperStatus()
{
  if ((solverStatus != Converged) && (stepperStatus == Failed))
    return Failed;

  // Check to see if we hit final parameter value
  if ((solverStatus == Converged) && 
      (curValue == finalValue))
    return Converged;

  // Check to see if max number of steps has been reached
  if (numTotalSteps == maxConSteps) 
    return Failed;

  // Check to see if the initial solve (step=0) failed
  // if (stepNumber == 0) 
  //  return Failed;

  return Unconverged;
}
