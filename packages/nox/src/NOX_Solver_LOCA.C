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

#include "NOX_Solver_LOCA.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::Solver;

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

LOCA::LOCA(Abstract::Group& xgrp, Status::Test& continuationTest, 
		       Status::Test& nonlinearSolverTest, 
		       const Parameter::List& p) :
  solnPtr(&xgrp),		     // pointer to xgrp
  oldSolnPtr(xgrp.clone(DeepCopy)),  // create via clone
  oldSoln(*oldSolnPtr),		     // reference to just-created pointer
  testConPtr(&continuationTest),     // pointer to status test 
  iparams(p),			     // copy p
  method("")
{
  init();
}

// Protected
void LOCA::init()
{
  // Print LOCA Copyright
  if (Utils::doPrint(Utils::LocaOuterIteration))
    cout << "LOCA v2.0, Copyright 2002 Sandia Corporation";

  // Get the continuation sublist
  Parameter::List& conParam = iparams.sublist("Continuation");

  status = Status::Unconverged;
  
  // Set the continuation parameters  
  startValue = conParam.getParameter("Initial Value", 0.0);
  curValue = startValue;
  finalValue = conParam.getParameter("Final Value", 1.0);
  initStepSize = conParam.getParameter("Initial Step Size", 1.0);
  minStepSize = conParam.getParameter("Min Step Size", 1.0e-6);
  maxStepSize = conParam.getParameter("Max Step Size", 1.0);
  curStepSize = conParam.getParameter("Initial Step Size", 1.0);
  oldStepSize = 0.0;
  stepNumber = 0;                    
  agrValue = conParam.getParameter("Step Size Aggressiveness", 0.0);
  maxNonlinearIters = conParam.getParameter("Max Nonlinear Iterations", 20);
  maxConSteps = conParam.getParameter("Max Continuation Steps", 1);
  order = conParam.getParameter("Order of Continuation", ZERO);

  // Set up utilities (i.e., set print processor, etc)
  Utils::setUtils(iparams.sublist("Nonlinear Solver"));
  
  // Print out initialization parameter information
  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "\n-- User Defined Parameters Passed to LOCA --\n\n";
    iparams.print(cout,5);
  }

  // Test the initial guess
  status = testConPtr->operator()(*this);

  // Print the parameters if requested
  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n-- Status Tests Passed to LOCA --\n\n";
    testConPtr->print(cout, 5);
    cout <<"\n" << Utils::fill(72) << "\n";
  }

  //Create the solver type, throw error if invalid method is chosen
  //RPP: NOTE: need to check if this has changed for a restart!!!!
  setSolver();

  // Print out initial information
  printStart();
}

bool LOCA::reset(Abstract::Group& xgrp, Status::Test& t, 
		 const Parameter::List& p) 
{
  solnPtr = &xgrp;
  testConPtr = &t;
  iparams = p;	
  init();
  return true;
}

LOCA::~LOCA() 
{
  delete oldSolnPtr;
}


Status::StatusType LOCA::getStatus()
{
  return status;
}

Status::StatusType LOCA::iterate()
{
  // First check status
  if (status != Status::Unconverged) 
    return status;

  // Copy pointers into temporary references
  Abstract::Group& soln = *solnPtr;
  Status::Test& test = *testConPtr;

  // Copy current soln to the old soln.
  oldSoln = soln;

  // Compute RHS for new current solution.

  // Update iteration count.
  stepNumber ++;

  // Evaluate the current status.
  status = test(*this);
 
  // Return status.
  return status;
}

Status::StatusType LOCA::solve()
{
  printUpdate();

  // Iterate until converged or failed
  while (status == Status::Unconverged) {
    status = iterate();
    printUpdate();
  }

  return status;
}

const Abstract::Group& LOCA::getSolutionGroup() const
{
  return *solnPtr;
}

const Abstract::Group& LOCA::getPreviousSolutionGroup() const
{
  return oldSoln;
}

int LOCA::getNumIterations() const
{
  return stepNumber;
}

const Parameter::List& LOCA::getOutputParameters() const
{
  oparams.setParameter("Continuation Steps", stepNumber);
  oparams.setParameter("2-Norm of Residual", solnPtr->getNormRHS());
  return oparams;
}

// protected
void LOCA::printStart() 
{
  if (Utils::doPrint(Utils::LocaOuterIteration)) {
    cout << endl << Utils::repeat(cout, 80) << endl;
    cout << "Beginning LOCA Continuation: " 
	 << iparams.getParameter("LOCA Method","???") << endl;
    cout << "Initial Step Size = " << initStepSize << endl;
    cout << "Maximum Number of Continuation Steps = " << maxConSteps 
	 << endl;
    cout << "Final Parameter Value = " << finalValue << endl;
    cout << Utils::repeat(cout, 80) << endl << endl;
  }
}

// protected
void LOCA::printUpdate() 
{

  // Print results of continuation step
  if (Utils::doPrint(Utils::LocaOuterIteration)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "-- Continuation Step " << stepNumber << " -- \n";
    if (status == Status::Converged)
      cout << " (Converged!)";
    if (status == Status::Failed)
      cout << " (Failed!)";
    cout << "\n" << Utils::fill(72) << "\n" << endl;
  }
  
  if ((status != Status::Unconverged) && 
      (Utils::doPrint(Utils::LocaOuterIteration))) {
    cout << Utils::fill(72) << "\n";
    cout << "-- Final Status Test Results --\n";    
    testConPtr->print(cout);
    cout << Utils::fill(72) << "\n";
  }
}

void LOCA::printFailedStep() 
{
  //RPP: We may not need this, the failure info should be at the method level!
  cout << endl << Utils::repeat(cout, 80) << endl;
  cout << "Continuation Step Number " << stepNumber << " experienced a "
       << "convergence failure in the nonlinear or linear solver!" << endl;
  cout << "Value of continuation parameter at failed step = " << curValue
       << endl;
  cout << endl << Utils::repeat(cout, 80) << endl;
}


bool LOCA::setSolver()
{
  string newmethod = iparams.getParameter("LOCA Method", "Newton");

  if (method != newmethod) {
    
    method = newmethod;

    delete solverPtr;
    solverPtr = 0;
    
    if (method == "Zero Order Continuation") {
      //solverPtr = new ZeroOrderContinuation(grp, tests, params);
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
      cout << "ERROR: NOX::Solver::LOCA - Invalid solver choice" << endl;
      throw "LOCA Error";
    }

    if (solverPtr == NULL) {
      cerr << "NOX::Solver::LOCA::setSolver - Null pointer error for solver" 
	   << endl;
      return false;
    }

    return true;
  }
  else {

    if (solverPtr == NULL) {
      cerr << "NOX::Solver::LOCA::reset - Null pointer error for solver" 
	   << endl;
      return false;
    }

    return solverPtr->reset(*solnPtr, *testConPtr, iparams);
  }
}
