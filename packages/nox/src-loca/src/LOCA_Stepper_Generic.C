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

#include "LOCA_Stepper_Generic.H"    // class definition

// LOCA Includes
#include "LOCA_Utils.H"		     // for static function doPrint
#include "LOCA_Solver_Generic.H"
#include "LOCA_Abstract_Group.H"
#include "LOCA_Abstract_Vector.H"

using namespace LOCA;
using namespace LOCA::Stepper;

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


Generic::Generic(Solver::Generic& s) :
  paramsPtr(0),
  conParams(s.getSolutionGroup().getParams()),
  solverPtr(0)
{
  resetGenericMembers(s);
}

Generic::~Generic() 
{ 

}


void Generic::printInitializationInfo()
{  
  if (Utils::doPrint(Utils::StepperIteration)) {
    cout << endl << Utils::fill(72, '~') << endl;
    cout << "Beginning Continuation Run \n" 
	 << "Stepper Method:     " << stepperMethod << "\n"
	 << "Solver Method: " << solverMethod << "\n\n"
	 << "Initial Parameter Value = " << startValue << "\n"
	 << "Final Parameter Value = " << finalValue << "\n"
	 << "Initial Step Size = " << startStepSize << "\n" 
	 << "Maximum Number of Continuation Steps = " << maxConSteps 
	 << endl;
    cout << Utils::fill(72, '~') << endl << endl;
  }
}
 
void Generic::printStartStep()
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

void Generic::printEndStep(NOX::StatusTest::StatusType& solverStatus)
{
  if (solverStatus != NOX::StatusTest::Failed) {
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

void Generic::printEndInfo()
{

}

void Generic::resetGenericMembers(Solver::Generic& s)
{
  // Set the solver pointer
  if (&s == 0) {
    cout << "ERROR: LOCA::Stepper::Generic::resetGenericMembers() - Solver is NULL!"
	 << endl;
    throw "LOCA Error";
  }
  else 
    solverPtr = &s;
  
  // Get a reference to the parameter stepper sublist from the Solver.
  const NOX::Parameter::List& p = s.getParameterList().sublist("Stepper");
  
  // Get the stepper method
  if (p.isParameter("Stepper Method"))
    stepperMethod = p.getParameter("Stepper Method", "???");
  else {
    if(Utils::doPrint(Utils::StepperDetails));
    cout << "ERROR: LOCA::Stepper::Generic::resetGenericMembers() - "
	 << "\"Stepper Method\" is not set in the input parameter list!" 
	 << endl;
    throw "LOCA Error";
  }
 
  // Get the solver method
  solverMethod = solverPtr->getLabel();

  // Get the continuation parameter
  conParamID = p.getParameter("Continuation Parameter", "???");
  if (conParamID == "???") {
    cout << "ERROR: LOCA::Stepper::Generic::resetGenericMembers() - "
	 << "\"Continuation Parameter\" is not set" << endl;
    throw "LOCA Error";
  }

  // Get the continuation parameter starting value
  if (p.isParameter("Initial Value"))
    startValue = p.getParameter("Initial Value", 0.0);
  else {
    cout << "ERROR: LOCA::Stepper::Generic::resetGenericMembers() - "
	 << "\"Initial Value\" of continuation param is not set!" << endl;
    throw "LOCA Error";
  }
  
  // Get the final value of the continuation parameter
  if (p.isParameter("Final Value"))
    finalValue = p.getParameter("Final Value", 0.0);
  else {
    cout << "ERROR: LOCA::Stepper::Generic::resetGenericMembers() - "
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
  status = NOX::StatusTest::Unconverged;
  solverStatus = NOX::StatusTest::Unconverged;
}

double Generic::computeStepSize(NOX::StatusTest::StatusType solverStatus)
{
  double tmpStepSize = curStepSize;
  double predStepSize = curStepSize;

  if ((solverStatus == NOX::StatusTest::Failed) || 
      (solverStatus == NOX::StatusTest::Unconverged)) {

    // A failed nonlinear solve cuts the current step size in half
    tmpStepSize = curStepSize * 0.5;    

  }
  else  {

    // adapive step size control
    if (agrValue != 0.0) {

      // Number of nonlinear iterations to reach convergence for last 
      // nonlinear solve -- used to pick next step size
      double numNonlinearSteps = ((double) solverPtr->getNumIterations());
      
      tmpStepSize = curStepSize * (1.0 + agrValue * (numNonlinearSteps/
					   ((double) maxNonlinearSteps)));
  
    }
    // if constant step size (agrValue = 0.0), the step size may still be 
    // reduced by a solver failure.  We should then slowly bring the step 
    // size back towards its constant value using agrValue = 0.5.
    else{
      if (curStepSize != startStepSize) {  
	double numNonlinearSteps = ((double) solverPtr->getNumIterations());

	tmpStepSize = prevStepSize * (1.0 + 0.5 * (numNonlinearSteps/
					      ((double) maxNonlinearSteps)));
	tmpStepSize = min(tmpStepSize, startStepSize);
      }
    }
  }  
  
  // Clip the step size if outside the bounds
  predStepSize = min(tmpStepSize, maxStepSize);
  predStepSize = min(tmpStepSize, maxStepSize);
  
  return predStepSize;
}

NOX::StatusTest::StatusType Generic::checkStepperStatus()
{
  // Check to see if we hit final parameter value
  if ((solverStatus == NOX::StatusTest::Converged) && 
      (curValue == finalValue))
    return NOX::StatusTest::Converged;

  // Check to see if max number of steps has been reached
  if (numTotalSteps == maxConSteps) 
    return NOX::StatusTest::Failed;

  return NOX::StatusTest::Unconverged;
}
