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

#ifndef LOCA_STEPPER_GENERIC_H
#define LOCA_STEPPER_GENERIC_H

#include "LOCA_Parameter_Vector.H"   // class data element
#include "NOX_Parameter_List.H"	     // class data element
#include "NOX_StatusTest_Generic.H"  // for enum StatusType
#include "NOX_Common.H"              // for <string> class

// Forward Declarations
namespace LOCA {
  namespace Abstract {
    class Vector;
    class Group;
  }
  namespace Solver {
    class Generic;
  }
}

namespace LOCA {

namespace Stepper {
  
/*!
  \brief Driver for the %LOCA Library.

  This is the base class object for all stepper algorithms in LOCA.  
  
  The following parameters are valid for the LOCA::Stepper::Generic base class and are ste through the "Stepper" sublist:

  - "Stepper Method" - Choice of continuation algorithm to choose:
    <ul>
    <li> "Zero Order"
    <li> "First Order"
    <li> "Arc-Length" </ul> </li>

  - "Bifurcation Parameters" - sublist containing the names and values of the parameters used in the user's application.

  - "Continuation Parameter" - String with the name of the parameter from the "Bifurcation Parameter" sublist that will be used as the continuation parameter.

  - "Initial Value" - Initial value of continuation parameter.

  - "Final Value" - Final value of continuation parameter.

  - "Initial Step Size" - Initial step size of the continuation parameter.

  - "Min Step Size" - Minimum step size of continuation parameter.

  - "Max Step Size" - Maximum step size of continuation parameter.

  - "Step Size Aggressiveness" - Controls the aggressiveness in computing the next step size.  Must be a positive number and can be greater than 1.0.  Setting this to zero yields a constant step size.

  - "Max Continuation Steps" - Maximum number of continuation steps.

  - "Max Nonlinear Iterations" - Maximum number of nonlinear iterations allowed in the nonlinear solve.  This is used in computing a new step size based on the difficulty of the nonlinear solve at the last step.

  \authors Roger Pawlowski (SNL 9233), Andrew Salinger (SNL 9233)
*/

class Generic {

public:

  //! Constructor 
  Generic(Solver::Generic& s);

  //! Destructor
  virtual ~Generic();

  //! Reset the Stepper to start a new continuation run
  virtual bool reset(Solver::Generic& s) = 0;

  //! Reset the solver after a failed continuation step
  virtual bool resetFromFailedStep() = 0;

  //! Return the status of the continuation routine
  virtual NOX::StatusTest::StatusType getStatus() = 0;

  //! Solve for the entire continuation run
  virtual NOX::StatusTest::StatusType solve() = 0;

  //! Solve for one step of the continuation problem 
  //! (i.e. one entire nonlinear solve)
  virtual NOX::StatusTest::StatusType step() = 0;

  //! Return the current solution group.
  virtual const Abstract::Group& getSolutionGroup() const = 0;

  //! Return the solution group with the solution for the previous continuation
  //! step.
  virtual const Abstract::Group& getPreviousSolutionGroup() const = 0;

  //! Returns the number of continuation steps taken by LOCA
  virtual int getNumContinuationSteps() const;

  //! Returns the number of continuation steps that failed to converge
  virtual int getNumFailedSteps() const;

  //! Returns the total number of continuation steps attempted
  virtual int getNumTotalSteps() const;

  //! Return the output parameters from the stepper algorithm. 
  virtual const NOX::Parameter::List& getParameterList() const = 0;

protected:

  //! Print to the screen the initialization information
  virtual void printInitializationInfo();

  //! Print to the screen information for the beginning of a new continuation step 
  virtual void printStartStep();

  //! Print to the screen information for the end of a successful continuation step
  virtual void printEndStep(NOX::StatusTest::StatusType& solverStatus);

  //! Print to the screen the concluding information
  virtual void printEndInfo();

  //! Reset the objects in the base class
  virtual void resetGenericMembers(Solver::Generic& s);

  //! Compute a new step size based on results from the previous step
  virtual double computeStepSize(NOX::StatusTest::StatusType solverStatus);

  //! Checks to see if convergence of the stepper algorithm is achieved.  
  /*! Convergence is achieved if the stepper has either hit the max number of steps or has reached the final parameter value
   */
  virtual NOX::StatusTest::StatusType checkStepperStatus();

protected:

  //! Label with the requested Stepper routine 
  string stepperMethod;
  
  //! Label with the requested nonlinear solver routine
  string solverMethod;

  //! Label that identifies the continuation parameter from the list of bifurcation parameters supplied by the user.
  string conParamID;

  //! Starting value of continuation parameter, \f$\lambda_b \f$.
  double startValue;

  //! Final value of continuation parameter, \f$\lambda_e \f$.
  double finalValue;

  //! Value of continuation parameter at previous continuation step, \f$\lambda_{i-1} \f$.
  double prevValue;

  //! Current value of continuation parameter, \f$\lambda_i \f$.
  double curValue;

  //! Initial Step Size of the continuation parameter, \f$\Delta \lambda_0 \f$.
  double startStepSize;

  //! Minimum Step Size of the continuation parameter, \f$\Delta \lambda_{min} \f$.
  double minStepSize;

  //! Maximum Step Size of the continuation parameter, \f$\Delta \lambda_{max} \f$.
  double maxStepSize;

  //! Step size of previous step.  Initializes to zero on first step, \f$\Delta \lambda_{i-1} \f$.
  double prevStepSize;

  //! Current step size (change in the continuation parameter), \f$\Delta \lambda_i \f$.
  double curStepSize;

  //! Current step number of continuation algorithm (does not include failed steps).
  int stepNumber;                    

  //! Number of continuation steps that have failed.
  int numFailedSteps;

  //! Total number of steps attempeted (includes failed and successful steps).
  int numTotalSteps;

  //! Aggressiveness parameter used to adjust the step size when using an 
  //! adjustable step size, \f$ a \f$.  Defaults to zero (a constant 
  //! step size).  Range: [0.0 to 1.0] 
  double agrValue;
  
  //! Maximum number of newton iterations per continuation step,\f$ N_{max} \f$. Defaults to 15.
  int maxNonlinearSteps;

  //! Maximum number of continuation steps to take. Defaults to 100.
  int maxConSteps;

  //! Current status of the stepper.
  NOX::StatusTest::StatusType status;
  
  //! Status of the last step of the Solver.
  NOX::StatusTest::StatusType solverStatus;
  
  //! Pointer to  the parameter list owned by the solver.
  NOX::Parameter::List* paramsPtr;
  
  //! Copy of the ParameterVector of the solution group.
  ParameterVector conParams;

  //! Pointer to the solver object that was passed in through the constructor or reset methods.
  Solver::Generic* solverPtr;

private:
  //! Disallow default constructor
  Generic() {};

};
} // namespace Stepper
} // namespace LOCA

#endif

