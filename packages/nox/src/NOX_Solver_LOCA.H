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

#ifndef NOX_SOLVER_LOCA_H
#define NOX_SOLVER_LOCA_H

#include "NOX_Solver_Generic.H"	    // base class
#include "NOX_Parameter_List.H"	    // class data element

namespace NOX {
namespace Solver {
  
/*!
  \brief Driver for the LOCA Library.

  This method accesses the Library of Continuation Algorithms (LOCA). 
  Zero order, first order and arc-length continuation routines are in this
  object.
  
  The following parameters are valid for this solver:

  - "Continuation" - Sublist of the LOCA continuation parameters, passed to
    the NOX::Linesearch::Manager constructor. Defaults to an empty list.

  - "Nonlinear Solver" - Sublist containing the parameters that are 
    normally set for a standard NOX nonlinear solve.  See NOX example problems 
    for details.  This includes the "Line Search", "Direction", and 
    "Linear %Solver" sublists.

  - "LOCA Method" - String that chooses the correct algorithm for the 
    nonlinear iteration.  Valid choices are:
    <li> "Zero Order Continuation"
    <li> "First Order Continuation"
    <li> "Arc Length Continuation"
    <li> "Turning Point" - This is a Zero Order method
    <li> "Pitchfork" - This is a Zero Order method
    <li> "Hopf" - This is a Zero Order method
    <li> "Homotopy"
    <li> "Phase Transition" - This is a Zero Order method
 

  \authors Andrew Salinger (SNL 9233), Roger Pawlowski (SNL 9233)
*/
class LOCA : public Generic {

public:

  //! %LOCA method constructor. 
  LOCA(Abstract::Group& xgrp, 
       Status::Test& continuationTest, 
       Status::Test& nonlinearSolverTest, 
       const Parameter::List& p);

  //! Destructor
  virtual ~LOCA();

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

  //! Prints the current iteration information.
  virtual void printUpdate();

  //! Prints information about a failed continuation step.
  virtual void printFailedStep();

  //! Creates the LOCA solver algorithm to be used for each continuation step. 
  bool setSolver();

protected:
  
  //! Pointer to the NOX::Solver:Generic object to be used in the nonlinear 
  //! solve at each continuation step.
  Generic* solverPtr;

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
    
  //! Pointer to the current solution group.
  Abstract::Group* solnPtr;		

  //! Pointer to solution group at previous continuation step.
  Abstract::Group* oldSolnPtr;		

  //! Solution at previous continuation step.
  Abstract::Group& oldSoln;		

  //! Stopping test to be used for the continuation loop.
  Status::Test* testConPtr;		

  //! Stopping test to be used for nonlinear solves.
  Status::Test* testNonlinearPtr;		

  //! Input parameters.
  Parameter::List iparams;	

  //! string containing the current solver method.
  string method;

  //! Output paramters.
  mutable Parameter::List oparams; 

  //! Starting value of continuation parameter, \f$\lambda_b \f$.
  double startValue;

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
} // namespace Solver
} // namespace NOX

#endif

