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

#include "LOCA_Stepper.H" // class definition

// NOX Support
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Manager.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

// LOCA Methods
#include "LOCA_Solver_ZeroOrder.H"

using namespace NOX;
using namespace NOX::Solver;
using namespace LOCA;

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

Stepper::Stepper(const Parameter::List& p,
		 Status::Test& test,
		 Abstract::Group& initialGuess, 
		 Abstract::Vector& nullVectorX,
		 Abstract::Vector& nullVectorY,
		 Abstract::Group& massMatrix) :
  solverPtr(0),
  iparams(p),                        // copy p
  testPtr(&test),                  // pointer to nonlinear status test 
  solnPtr(&initialGuess),            // pointer to initial guess and Jacobian
  oldSolnPtr(initialGuess.clone(DeepCopy)),  // create via clone
  oldSoln(*oldSolnPtr),		     // reference to just-created pointer
  nullVectorXPtr(&nullVectorX),      // pointer to null vector X
  nullVectorYPtr(&nullVectorY),      // pointer to null vector Y
  massMatrixPtr(&massMatrix),        // pointer to mass matrix
  method(""),
  conParamLabel("")
{
  // Set up utilities (i.e., set print processor, etc)
  Utils::setUtils(iparams.sublist("Nonlinear Solver"));
  
  // Print LOCA Copyright
  if (Utils::doPrint(Utils::Loca))
    cout << endl << "LOCA v2.0, Copyright 2002 Sandia Corporation" << endl;
  
  init();
}

// Protected
void Stepper::init()
{
  status = Status::Unconverged;
  
  // Set the continuation parameters  
  startValue = iparams.getParameter("Initial Value", 0.0);
  oldValue = 0.0;
  curValue = startValue;
  finalValue = iparams.getParameter("Final Value", 1.0);
  initStepSize = iparams.getParameter("Initial Step Size", 1.0);
  minStepSize = iparams.getParameter("Min Step Size", 1.0e-12);
  maxStepSize = iparams.getParameter("Max Step Size", 1.0e+12);
  curStepSize = iparams.getParameter("Initial Step Size", 1.0);
  oldStepSize = 0.0;
  stepNumber = 0;                    
  agrValue = iparams.getParameter("Step Size Aggressiveness", 0.0);
  // RPP add accessor method to maxiters test to get the following value
  maxNonlinearIters = iparams.getParameter("Max Nonlinear Iterations", 20);
  maxConSteps = iparams.getParameter("Max Continuation Steps", 1);
  order = iparams.getParameter("Order of Continuation", ZERO);

  // Print out parameter list information
  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "\n-- User Defined Parameters Passed to LOCA --\n\n";
    iparams.print(cout,5);
  }

  // Print the parameters if requested
  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n-- Status Test Passed to LOCA for nonlinear solves --\n\n";
    testPtr->print(cout, 5);
    cout <<"\n" << Utils::fill(72) << "\n";
  }

  // Get the continuation parameter label
  if (iparams.isParameter("Continuation Parameter") == true) {
    conParamLabel = iparams.getParameter("Continuation Parameter", "");
    
    // Make sure the label is valid
    bool isValid = solnPtr->setParameter(conParamLabel, curValue);
    if (isValid == false) {
      cout << "ERROR: LOCA::Stepper::init() - Continuation Parameter named \""
	   << conParamLabel << "\" has failed in call to setParameter()" 
	   << endl;
    }
  }
  else {
    cout << "ERROR: LOCA::Stepper::init() - \"Continuation Parameter\" "
	 << "is not set in the parameter list!" << endl;
    throw "NOX Error";
  }

  //Create the solver type, throw error if invalid method is chosen
  setSolver();

  // Print out initial information
  printStart();
}

bool Stepper::reset(Abstract::Group& xgrp, Status::Test& t, 
		 const Parameter::List& p) 
{
  solnPtr = &xgrp;
  testPtr = &t;
  iparams = p;	
  init();
  return true;
}

Stepper::~Stepper() 
{
  delete oldSolnPtr;
}


Status::StatusType Stepper::getStatus()
{
  return status;
}

Status::StatusType Stepper::iterate()
{
  // Copy pointers into temporary references
  Abstract::Group& soln = *solnPtr;
  Solver::Generic& solver = *solverPtr;

  // Print out initial step information
  printStartStep();

  // Copy current soln to the old soln.
  oldSoln = soln;

  // Set the continuation parameter value in the application code
  // Must be called before we reset the solver!!!
  soln.setParameter(conParamLabel, curValue);

  // Reset nonlinear status test and solver 
  // RPP: In solver.reset() the status of the solver can be changed to
  // converged.  So we must set the parameter before we reset the solver
  // otherwise the solve won't even try to solve things 
  Status::StatusType nlSolverStatus = Status::Unconverged;
  solver.reset(soln, *testPtr, iparams);

  // Do a complete nonlinear solve
  nlSolverStatus = solver.solve();
  
  // Return the results of this step
  printEndStep(nlSolverStatus);

  // If the final step length of this step is on or past the finalValue
  // then flag it as converged
  if (curValue >= finalValue) {
    
    return Status::Converged;
  }

  // If the nonlinear solver failed, cut the step length in half and try again.
  // RPP: We will move this into the solver (0,1,arc) level later.
  if (nlSolverStatus != Status::Converged) {
    soln = oldSoln;
    curStepSize *= 0.5;
    curValue = oldValue + curStepSize;
    return status;
  }

  // Get the step length for the next continuation step
  oldStepSize = curStepSize;
  if (stepNumber != 0)
    solver.getStepLength(curStepSize);
  
  // Step the continuation parameter
  oldValue = curValue;
  curValue += curStepSize;

  // Update iteration count.
  stepNumber ++;
  
  return status;
}

Status::StatusType Stepper::solve()
{
  printStartStep();

  // For continuation step 0 (the first time through), check to make 
  // sure the initial guess is on a steady state solution - don't
  // do an actual continuation step, just a pure steady state solve.
  if (Utils::doPrint(Utils::Loca))
    cout << "    Forcing the initial guess to a steady state." << endl;

  NOX::Solver::Manager* noxSolver = new NOX::Solver::Manager(*solnPtr, *testPtr, iparams.sublist("Nonlinear Solver"));
  status = noxSolver->solve();

  // If it failed, exit program.
  if (status != Status::Converged) {
    cout << "ERROR: LOCA::Stepper::Solve() - the initial guess failed "
	 << "to converge to a steady state!  Aborting continuation run!"
	 << endl;
    throw "NOX Error";
  }

  if (Utils::doPrint(Utils::Loca))
    cout << endl << "    Initial guess converged to a steady state!" << endl;
  oldSoln = *solnPtr;
  delete noxSolver;
  printEndStep(status);
  solverPtr->getStepLength(curStepSize);
  oldValue = curValue;
  curValue += curStepSize;
  stepNumber ++;
  
  status = Status::Unconverged;

  // Iterate until converged or failed
  // LOCA does not use the normal convergence tests!  For continuation
  // we need to either hit the max iterations or the final parameter value.
  while ((status == Status::Unconverged) && 
	 (stepNumber <= maxConSteps)) {
    status = iterate();
  }

  return status;
}

const Abstract::Group& Stepper::getSolutionGroup() const
{
  return *solnPtr;
}

const Abstract::Group& Stepper::getPreviousSolutionGroup() const
{
  return oldSoln;
}

int Stepper::getNumIterations() const
{
  return stepNumber;
}

const Parameter::List& Stepper::getOutputParameters() const
{
  oparams.setParameter("Continuation Steps", stepNumber);
  oparams.setParameter("2-Norm of Residual", solnPtr->getNormRHS());
  return oparams;
}

// protected
void Stepper::printStart() 
{
  if (Utils::doPrint(Utils::Loca)) {
    cout << endl << Utils::fill(72, '~') << endl;
    cout << "Beginning LOCA Continuation Run: " 
	 << iparams.getParameter("LOCA Method","???") << endl;
    cout << "Initial Step Size = " << initStepSize << endl;
    cout << "Maximum Number of Continuation Steps = " << maxConSteps 
	 << endl;
    cout << "Final Parameter Value = " << finalValue << endl;
    cout << Utils::fill(72, '~') << endl << endl;
  }
}

// protected
void Stepper::printStartStep() 
{
  if (Utils::doPrint(Utils::Loca)) {
    cout << "\n" << Utils::fill(72, '~') << "\n";
    cout << "Start of Continuation Step " << stepNumber << endl;
    cout << "Continuation Method: " << method << endl;
    cout << "Continuation Parameter: " << conParamLabel
	 << " = " << curValue << " from " << oldValue << endl;
    cout << "Current step size  = " << curStepSize << "   "
	 << "Previous step size = " << oldStepSize << endl;
    cout << "\n" << Utils::fill(72, '~') << "\n" << endl;
  }
}

// protected
void Stepper::printEndStep(NOX::Status::StatusType& solverStatus) 
{
  if (solverStatus != Status::Failed) {
    // Print results of successful continuation step
    if (Utils::doPrint(Utils::Loca)) {
      cout << "\n" << Utils::fill(72, '~') << "\n";
      cout << "End of  Continuation Step " << stepNumber << endl;
      cout << "Continuation Parameter: " << conParamLabel
	   << " = " << curValue << " from " << oldValue << endl;
      cout << "--> Step Converged!";
      cout << "\n" << Utils::fill(72, '~') << "\n" << endl;
    }
  }
  else {
    if (Utils::doPrint(Utils::Loca)) {
      // RPP: We may not need this, the failure info should be 
      // at the method level!
      cout << endl << Utils::fill(72, '~') << endl;
      cout << "Continuation Step Number " << stepNumber << " experienced a "
	   << "convergence failure in the nonlinear solver!" << endl;
      cout << "Value of continuation parameter at failed step = " << curValue
	   << " from " << oldValue << endl;
      cout << endl << Utils::fill(72, '~') << endl;
    }
  }
}

bool Stepper::setSolver()
{
  // Check to make sure "LOCA Method" is set.
  if (iparams.isParameter("LOCA Method") == false) {
    cout << "ERROR: LOCA::Stepper::setSolver() - The parameter \"LOCA Method\""
	 << " is not set in the parameter list!" << endl;
    throw "NOX Error";
  }

  string newMethod = iparams.getParameter("LOCA Method", "Newton");

  if (method != newMethod) {
    
    method = newMethod;

    delete solverPtr;
    solverPtr = 0;
    
    if (method == "Zero Order Continuation") {
      solverPtr = new LOCA::Solver::ZeroOrder(*solnPtr, *testPtr, iparams);
    } 
    else if (method == "First Order Continuation") {
      //ptr = new LineSearch(grp, tests, params);
    } 
    else if (method == "Arc Length Continuation") {
      //ptr = new NonlinearCG(grp, tests, params);
    } 
    else if (method == "Turning Point") {
      //ptr = new TrustRegion(grp, tests, params);
    } 
    else if (method == "Pitchfork") {
      //ptr = new TrustRegion(grp, tests, params);
    } 
    else if (method == "Hopf") {
      //ptr = new TrustRegion(grp, tests, params);
    } 
    else if (method == "Homotopy") {
      //ptr = new TrustRegion(grp, tests, params);
    } 
    else if (method == "Phase Transition") {
      //ptr = new TrustRegion(grp, tests, params);
    } 
    else {
      cout << "ERROR: LOCA::Stepper - Invalid solver choice for \"" 
	   << "LOCA Method\" parameter!" << endl;
      throw "LOCA Error";
    }

    if (solverPtr == NULL) {
      cerr << "LOCA::Stepper::setSolver() - failed to allocate memory "
	   << "for solver!" << endl;
      return false;
    }

    return true;
  }
  else {

    if (solverPtr == NULL) {
      cerr << "LOCA::Stepper::setSolver() - Null pointer error for solver" 
	   << endl;
      return false;
    }

    return solverPtr->reset(*solnPtr, *testPtr, iparams);
  }
}

bool Stepper::perturb()
{
  /*
  Epetra_Vector* tmpVec1 = new Epetra_Vector(this->getEpetraVector());
  Epetra_Vector* tmpVec2 = new Epetra_Vector(this->getEpetraVector());
  tmpVec1->PutScalar(1.0);
  tmpVec2->Random();
  tmpVec2->Abs(this->getEpetraVector());
  tmpVec2->Update(-0.5, *tmpVec1, 1.0);
  tmpVec2->Scale(1.0e-4);
  tmpVec2->Update(1.0, *tmpVec1, 1.0);
  (this->getEpetraVector()).Multiply(1.0, this->getEpetraVector(), 
				      *tmpVec2, 0.0);
  delete tmpVec1;
  delete tmpVec2;
  */
  return true;
}
