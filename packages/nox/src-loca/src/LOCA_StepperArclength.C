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

#include "LOCA_StepperArclength.H"    // class definition

// LOCA Includes
#include "LOCA_Utils.H"		      // for static function doPrint
#include "LOCA_Abstract_Group.H"      // class data element
#include "LOCA_Abstract_Vector.H"     // class data element
#include "LOCA_Bifurcation_ArcLengthGroup.H"   //

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


StepperArclength::StepperArclength(const Stepper* sPtr,
                      LOCA::Bifurcation::ArcLengthGroup& alGroup) :
   LOCA::Stepper(*sPtr, alGroup),
   isFirstArclengthStep(true),
   originalParamFinalValue(finalValue)
{
  curGroupPtr = dynamic_cast<LOCA::Bifurcation::ArcLengthGroup*>(Stepper::curGroupPtr);
  prevGroupPtr = dynamic_cast<LOCA::Bifurcation::ArcLengthGroup*>(Stepper::prevGroupPtr);

  //Change member data here that is different for arclength than the original stepper
  if (!doFirstOrderPredictor) {
    cout << "First Order Predictor currently required for Arclength"
         << "continuation:\n    setting doFirstOrderPredictor=true" << endl;
    doFirstOrderPredictor = true;
  }

  stepperMethod = "Arclength Continuation";
  conParamID = "arclengthXXX"; // Arclength is parameter: no ID should be needed
  startValue = 0.0; // Arclength variable can be defined from anywhere
  curValue = startValue;
  prevValue = startValue;
  prevStepSize = 0.0;
}

StepperArclength::StepperArclength(const StepperArclength& s) :
  Stepper(s)
{ //Should not be called: declared private for safety
}

StepperArclength::~StepperArclength() 
{ 
  //delete prevGroupPtr;
  cout <<" REMINDER: StepperArclength::~StepperArclength():\n"
       <<" COPY stepper data (e.g. numFailedSteps) from arclengthtepper back to original stepper"<< endl;
}

bool StepperArclength::reset(LOCA::Abstract::Group& initialGuess,
		    NOX::StatusTest::Generic& t,
		    const NOX::Parameter::List& p) 
{
  return false;
}

StatusType StepperArclength::solve()
{
  while (stepperStatus == Unconverged) {
    stepperStatus = step();
  }
  return stepperStatus;
}

StatusType StepperArclength::step()
{
  stepperStatus = Unconverged;
 
  if (solverStatus == Failed) {

    *curGroupPtr = *prevGroupPtr;

    curStepSize = computeStepSize(Failed);
    curValue = prevValue + curStepSize;

  }
  else {

    // Save previous successful step information
    *prevGroupPtr = *curGroupPtr;
     
    // Set initial guess for next step equal to last step's final solution 
    *curGroupPtr = getSolutionGroup();
    curGroupPtr->setPrevX(prevGroupPtr->getX()); //New for arclength
 
    // Compute tangent values for first order predictor
    
    if (doFirstOrderPredictor) {
      curGroupPtr->computeTangent(paramList.sublist("Solver").sublist("Direction").sublist("Linear Solver"), -1); 

      if (isFirstArclengthStep) { //New for arclength
         double dp_ds = curGroupPtr->getTangent().getArcParam();
         cout << "Origin dp_ds = " << dp_ds << endl;
         //AGS: Add absolute value on whole rhs or last term?
         startStepSize /=  dp_ds;
         minStepSize /= fabs(dp_ds);
         maxStepSize /= fabs(dp_ds);
         isFirstArclengthStep = false;
      }
      // finalValue = originalParamFinalValue * ????WRITE THIS;
    }

    prevStepSize = curStepSize;
    prevValue = curValue;
    curStepSize = computeStepSize(solverStatus);

    curValue += curStepSize;

    if (doFirstOrderPredictor)
      curGroupPtr->computeX(*curGroupPtr, curGroupPtr->getTangent(), -curStepSize);
  }

  curGroupPtr->setArclengthStep(curStepSize); //New for arclength
  //conParams.setValue(conParamID, curValue);
  //curGroupPtr->setParams(conParams);

  solver.reset(*curGroupPtr, *statusTestPtr, paramList.sublist("Solver"));

  return nonlinearSolve();
}

double StepperArclength::computeStepSize(StatusType solverStatus)
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

  // Clip step size at minimum, signal for failed run
  if (fabs(tmpStepSize) < minStepSize) {
    stepperStatus = Failed;
    predStepSize =  minStepSize;
    if (tmpStepSize < 0.0) predStepSize *= -1.0;
  }
  
  // Cap the con parameter so we don't go past the final value
    
  if ( (prevValue+predStepSize-finalValue)*(prevValue-finalValue) < 0.0)
    predStepSize = finalValue - prevValue;

  return predStepSize;
}

bool StepperArclength::init()
// 
{ 
  // This function is here to prevent Stepper::init() from being called again, and thereby rereading parameter values set in constructor
  throw "LOCA Error";
}
