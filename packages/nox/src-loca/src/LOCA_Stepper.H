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

#ifndef LOCA_STEPPER_H
#define LOCA_STEPPER_H

#include "LOCA_Solver_Generic.H"    // class data element
#include "NOX_Parameter_List.H"	    // class data element
#include "NOX_Status_Test.H"	    // class data element

// Forward Declarations
namespace NOX {
  namespace Abstract {
    class Vector;
    class Group;
  }
}

namespace LOCA {
  
using namespace NOX;

/*!
  \brief Driver for the LOCA Library.

  This method accesses the Library of Continuation Algorithms (LOCA). 
  Zero order, first order and arc-length continuation routines are in this
  object.
  
  The following parameters are valid for the LOCA Library:

  - "Nonlinear %Solver" - Sublist containing the parameters that are 
    normally set for a standard NOX nonlinear solve.  See NOX example problems 
    for details.  This includes the "Line Search", "Direction", and 
    "Linear %Solver" sublists.

  - "LOCA Method" - String that chooses the correct algorithm for the 
    nonlinear iteration.  Valid choices are:
    <ul>
    <li> "Zero Order Continuation"
    <li> "First Order Continuation"
    <li> "Arc Length Continuation" 
    <li> "Turning Point"
    <li> "Pitchfork"
    <li> "Hopf"
    <li> "Homotopy"
    <li> "Phase Transition" </ul></li>      

  - "Initial Value" - Initial value of continuation parameter.

  - "Final Value" - Final value of continuation parameter.

  - "Initial Step Size" - Initial step size of the continuation parameter.

  - "Min Step Size" - Minimum step size of continuation parameter.

  - "Max Step Size" - Maximum step size of continuation parameter.

  - "Step Size Aggressiveness" - Controls the aggressiveness in computing the next step size.  Must be a positive number and can be greater than 1.0.  Setting this to zero yields a constant step size.

  - "Max Continuation Steps" - Maximum number of continuation steps.

  - "Order of Continuation" - Order of the continuation to use.

  - "Max Nonlinear Iterations" - Maximum number of nonlinear iterations allowed in the nonlinear solve.  This is used in computing a new step size based on the difficulty of the nonlinear solve at the last step.

  \authors Roger Pawlowski (SNL 9233), Andrew Salinger (SNL 9233)
*/

class Stepper {

public:

  //! %Stepper method constructor. 
  Stepper(const Parameter::List& p, 
	  Status::Test& test,
	  Abstract::Group& initialGuess,
	  Abstract::Vector& nullVector1,
	  Abstract::Vector& nullVector2,
	  Abstract::Group& massMatrix);

  //! Destructor
  virtual ~Stepper();

  virtual bool reset(Abstract::Group& xgrp, Status::Test& t, const Parameter::List& p);

  virtual Status::StatusType getStatus();

  //! Take one step of the continuation problem (one entire nonlinear solve)
  virtual Status::StatusType iterate();

  //! Solve the entire continuation run
  virtual Status::StatusType solve();

  virtual const Abstract::Group& getSolutionGroup() const;

  virtual const Abstract::Group& getPreviousSolutionGroup() const;

  //! Returns the number of continuation steps taken by LOCA
  virtual int getNumIterations() const;

  virtual const Parameter::List& getOutputParameters() const;

protected:
  
  //! Set all variables to initial values and calcuate the RHS.
  virtual void init();

  //! Print out the initialization information
  virtual void printStart();

  //! Print out information for the beginning of a new continuation step 
  virtual void printStartStep();

  //! Print out information for the end of a successful continuation step .
  virtual void printEndStep(NOX::Status::StatusType& solverStatus);

  //! Creates the LOCA solver algorithm to be used for each continuation step. 
  bool setSolver();

  
  //! Perturb the current vector values slightly.  
  /*! this method is used to move a converged solution at a singularity off
   * of the singularity before doing the next continuation step. 
   * This isn't pretty but has helped convergence on many turning 
   * point tracking problems:
   * 
   * \f[ this = this * (1.0+1.0E-4 * (x_{random} - 0.5))\f]
   *
   * where \f$ this \f$ is the vector to be perturbed 
   * and \f$x_{random}\f$ is a random vector ranging from 0.0 to 1.0
   * computed using the absolute values of the random() method.  
   * Returns true if successful.  
   */
  virtual bool perturb();
  
protected:
  
  //! Pointer to the LOCA::Solver:Generic object to be used in the nonlinear 
  //! solve at each continuation step.
  Solver::Generic* solverPtr;

  //! Type of continuation routine to use.
  enum Order {ZERO, FIRST, ARC_LENGTH}; 

  /*! Continuation order flag: 
   * 0 - zero-order continuation
   * 1 - first-order continuation
   * 2 - arc-length continuation
   * This flag is always 0 on the first solve, and 0 for turning point, 
   * pitchfork, or Hopf continuation. 
   */
  //Order order;  // compiler is bombing, parameterlist doesn't handle enum.
  int order;          
    
  //! Input parameters.
  Parameter::List iparams;	

  //! Stopping test to be used for nonlinear solves.
  Status::Test* testPtr;		

  //! Pointer to the current solution group.
  Abstract::Group* solnPtr;		

  //! Pointer to solution group at previous continuation step.
  Abstract::Group* oldSolnPtr;		

  //! Solution at previous continuation step.
  Abstract::Group& oldSoln;		

  //! Initial guess for null vector Y.
  Abstract::Vector* nullVectorXPtr;		

  //! Initial guess for null vector X.
  Abstract::Vector* nullVectorYPtr;		

  //! Group for the mass matrix evaluation.
  Abstract::Group* massMatrixPtr;		

  //! string containing the current solver method.
  string method;

  //! Name of the continuation parameter
  string conParamLabel;

  //! Output paramters.
  mutable Parameter::List oparams; 

  //! Starting value of continuation parameter, \f$\lambda_b \f$.
  double startValue;

  //! Value of continuation parameter at previous continuation step., \f$\lambda_{i-1} \f$.
  double oldValue;

  //! Current value of continuation parameter, \f$\lambda_i \f$.
  double curValue;

  //! Final value of continuation parameter, \f$\lambda_e \f$.
  double finalValue;

  //! Initial Step Size of the continuation parameter, \f$\Delta \lambda_0 \f$.
  double initStepSize;

  //! Minimum Step Size of the continuation parameter, \f$\Delta \lambda_{min} \f$.
  double minStepSize;

  //! Maximum Step Size of the continuation parameter, \f$\Delta \lambda_{max} \f$.
  double maxStepSize;

  //! Current step size (change in the continuation parameter), \f$\Delta \lambda_i \f$.
  double curStepSize;

  //! Step size of previous step.  Initializes to zero on first step, \f$\Delta \lambda_{i-1} \f$.
  double oldStepSize;

  //! Current step number of continuation algorithm.
  int stepNumber;                    

  //! Aggressiveness parameter used to adjust the step size when using an 
  //! adjustable step size, \f$a \f$.  Defaults to zero for a constant 
  //! step size.  Range: [0.0 to 1.0] 
  double agrValue;
  
  //! Number of nonlinear iterations to reach convergence for last 
  //! nonlinear solve -- used to pick next step size,\f$ N_{i} \f$   
  int numNonlinearSteps;    

  //! Maximum number of newton iterations per continuation step,\f$ N_{max} \f$.
  int maxNonlinearIters;

  //! Maximum number of continuation steps to take. Defaults to 1.
  int maxConSteps;

  //! %Status of nonlinear solver.
  Status::StatusType status;

};
} // namespace NOX

#endif

