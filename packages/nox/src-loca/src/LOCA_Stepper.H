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

#include "NOX_Parameter_List.H"	     // class data element
#include "NOX_Solver_Manager.H"	     // class data element
#include "LOCA_Continuation_Manager.H" // class data element
#include "LOCA_Predictor_Manager.H" // class data element
#include "LOCA_StepSize_Manager.H" // class data element
#include "NOX_StatusTest_Generic.H"  // for enum StatusType
#include "NOX_Common.H"              // for <string> class
#include "LOCA_Abstract_DataOutput.H"

// Forward Declarations
namespace LOCA {
  namespace Continuation {
    class Group;
  }
}

namespace LOCA {

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

class Stepper {

public:

  //! Constructor 
  Stepper(LOCA::Abstract::Group& initialGuess,
	  NOX::StatusTest::Generic& t,
	  NOX::Parameter::List& p,
	  LOCA::Abstract::DataOutput& dataOut);

  //! Destructor
  virtual ~Stepper();

  //! Reset the Stepper to start a new continuation run
  virtual bool reset(LOCA::Abstract::Group& initialGuess,
		     NOX::StatusTest::Generic& t,
		     const NOX::Parameter::List& p);

  //! Return the status of the continuation routine
  virtual NOX::StatusTest::StatusType getStatus();

  //! Solve for the entire continuation run
  virtual NOX::StatusTest::StatusType solve();

  //! Solve for one step of the continuation problem 
  //! (i.e. predictor, stepper, then call nonlinearSolve())
  virtual NOX::StatusTest::StatusType step();

  //! Call nonlinear solver, check status, print, update counters
  virtual NOX::StatusTest::StatusType nonlinearSolve();

  //@{ \name Information accessor methods.

  //! Return the current solution group.
  virtual const LOCA::Abstract::Group& getSolutionGroup() const;

  //! Returns the number of continuation steps taken by LOCA
  virtual int getNumContinuationSteps() const;

  //! Returns the number of continuation steps that failed to converge
  virtual int getNumFailedSteps() const;

  //! Returns the total number of continuation steps attempted
  virtual int getNumTotalSteps() const;

  //! Return the output parameters from the stepper algorithm. 
  virtual const NOX::Parameter::List& getParameterList() const;

  //@}

protected:

  //! Copy Constructor: Not implemented, but put here to prevent compiler form generating one
  Stepper(const Stepper& s);

  //! Reset the objects in the base class
  virtual bool init(LOCA::Abstract::Group& initialGuess);

  //! Compute a new step size based on results from the previous step
  virtual double computeStepSize(NOX::StatusTest::StatusType solverStatus);

  //! Checks to see if convergence of the stepper algorithm is achieved.  
  /*! Convergence is achieved if the stepper has either hit the max number of steps or has reached the final parameter value
   */
  virtual NOX::StatusTest::StatusType checkStepperStatus();

  //@{ \name Output methods.

  //! Print to the screen the initialization information
  virtual void printInitializationInfo();

  //! Print to the screen information for the beginning of a new continuation step 
  virtual void printStartStep();

  //! Print to the screen information for the end of a successful continuation step
  virtual void printEndStep(NOX::StatusTest::StatusType& solverStatus);

  //! Print to the screen the concluding information
  virtual void printEndInfo();

  //@}

protected:

  //! Continuation group manager
  LOCA::Continuation::Manager conGroupManager;

  //! Pointer to the Group with the initial guess for the nonlinear solver.  This is used for resetting the nonlinear solver in between successful steps.
  LOCA::Continuation::Group* curGroupPtr;

  //! Pointer to the a Group with the solution for the previous continuation step.  This is used for resetting the nonlinear solver in between failed steps.
  LOCA::Continuation::Group* prevGroupPtr;

  //! Reference to data output object.  This is used to save intermediate data generated during continuation run.
  LOCA::Abstract::DataOutput& dataOutput;

  //! Pointer to the status test.  This is used for resetting the nonlinear solver in between steps.
  NOX::StatusTest::Generic* statusTestPtr;

  //! Deep copy of the parameter list supplied initially in the constuctor. This can be changed via a call to the reset() method.
  NOX::Parameter::List paramList;

  //! Solver (a.k.a corrector)
  NOX::Solver::Manager* solverPtr;

  //! Predictor
  LOCA::Predictor::Manager predictor;

  //! Predictor direction
  LOCA::Continuation::Vector *predictorDirection;

  //! Step size manager
  LOCA::StepSize::Manager stepSizeManager;

  //! Starting value of continuation parameter, \f$\lambda_b \f$.
  double startValue;

  //! Final value of continuation parameter, \f$\lambda_e \f$.
  double finalValue;

  //! Current step size (change in the continuation parameter), \f$\Delta \lambda_i \f$.
  double stepSize;

  //! Current step number of continuation algorithm (does not include failed steps).
  int stepNumber;                    

  //! Number of continuation steps that have failed.
  int numFailedSteps;

  //! Total number of steps attempeted (includes failed and successful steps).
  int numTotalSteps;
  
  //! Maximum number of newton iterations per continuation step,\f$ N_{max} \f$. Defaults to 15.
  int maxNonlinearSteps;

  //! Maximum number of continuation steps to take. Defaults to 100.
  int maxConSteps;

  //! Current status of the stepper.
  NOX::StatusTest::StatusType stepperStatus;
  
  //! Status of the last step of the Solver.
  NOX::StatusTest::StatusType solverStatus;

  //! Flag inidicating if current step is the last step (reached final value)
  bool isLastStep;
  
};
} // end namespace LOCA

#endif

