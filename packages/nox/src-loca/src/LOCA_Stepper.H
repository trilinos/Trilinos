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

#include "LOCA_Abstract_Iterator.H"  // base class
#include "NOX_Parameter_List.H"	     // class data element
#include "NOX_Solver_Manager.H"	     // class data element
#include "LOCA_Continuation_Manager.H" // class data element
#include "LOCA_Predictor_Manager.H" // class data element
#include "LOCA_StepSize_Manager.H" // class data element
#include "NOX_StatusTest_Generic.H"  // for enum StatusType
#include "NOX_Common.H"              // for <string> class

// Forward Declarations
namespace LOCA {
  namespace Continuation {
    class Group;
  }
}

namespace LOCA {

  /*! 
    \brief Implementation of LOCA::Abstract::Iterator for computing points
    along a continuation curve.
  */
  /*!
    The %Stepper class implements the pure virtual methods of the 
    LOCA::Abstract::Iterator for iteratively computing points along a 
    continuation curve.
  */

  class Stepper : public LOCA::Abstract::Iterator {

  public:

    //! Constructor 
    Stepper(LOCA::Abstract::Group& initialGuess,
	    NOX::StatusTest::Generic& t,
	    NOX::Parameter::List& p);

    //! Destructor
    virtual ~Stepper();

    //! Reset the Stepper to start a new continuation run
    virtual bool reset(LOCA::Abstract::Group& initialGuess,
		       NOX::StatusTest::Generic& t,
		       NOX::Parameter::List& p);

    //! Return the current solution group.
    virtual const LOCA::Abstract::Group& getSolutionGroup() const;

    //! Return the output parameters from the stepper algorithm. 
    virtual const NOX::Parameter::List& getParameterList() const;

  protected:

    //! Copy Constructor: Not implemented, but put here to prevent compiler form generating one
    Stepper(const Stepper& s);

    //! Intialize iterator
    virtual LOCA::Abstract::Iterator::IteratorStatus start();

    //! Finalize iterator
    virtual LOCA::Abstract::Iterator::IteratorStatus finish();

    //! Preprocess step
    virtual LOCA::Abstract::Iterator::StepStatus 
    preprocess(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Compute step
    virtual LOCA::Abstract::Iterator::StepStatus 
    compute(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Postprocess step
    virtual LOCA::Abstract::Iterator::StepStatus 
    postprocess(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Check stopping criteria
    virtual LOCA::Abstract::Iterator::IteratorStatus 
    stop(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Compute a new step size based on results from the previous step
    LOCA::Abstract::Iterator::StepStatus
    computeStepSize(LOCA::Abstract::Iterator::StepStatus stepStatus, 
		    double& stepSize);

    //! Print to the screen the initialization information
    virtual void printInitializationInfo();

    //! Print to the screen information for the beginning of a new continuation step 
    virtual void printStartStep();

    //! Print to the screen information for the end of a successful continuation step
    virtual void printEndStep(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Print to the screen the concluding information
    virtual void printEndInfo();

  protected:

    //! Continuation group manager
    LOCA::Continuation::Manager* conGroupManagerPtr;

    //! Pointer to the Group with the initial guess for the nonlinear solver.  This is used for resetting the nonlinear solver in between successful steps.
    LOCA::Continuation::Group* curGroupPtr;

    //! Pointer to the a Group with the solution for the previous continuation step.  This is used for resetting the nonlinear solver in between failed steps.
    LOCA::Continuation::Group* prevGroupPtr;

    //! Pointer to the status test.  This is used for resetting the nonlinear solver in between steps.
    NOX::StatusTest::Generic* statusTestPtr;

    //! Pointer to parameter list passed in constructor/reset method
    NOX::Parameter::List* paramListPtr;

    //! Solver (a.k.a corrector)
    NOX::Solver::Manager* solverPtr;

    //! Predictor
    LOCA::Predictor::Manager* predictorManagerPtr;

    //! Current predictor direction
    LOCA::Continuation::Vector* curPredictorPtr;

    //! previous predictor direction
    LOCA::Continuation::Vector* prevPredictorPtr;

    //! Step size manager
    LOCA::StepSize::Manager* stepSizeManagerPtr;

    //! Starting value of continuation parameter, \f$\lambda_b \f$.
    double startValue;

    //! Largest allowable value of continuation parameter
    double maxValue;

    //! Smallest allowable value of continuation parameter
    double minValue;

    //! Current step size (change in the continuation parameter), \f$\Delta \lambda_i \f$.
    double stepSize;
  
    //! Maximum number of newton iterations per continuation step,\f$ N_{max} \f$. Defaults to 15.
    int maxNonlinearSteps;

    //! Flag indicating if this is the last step of continuation
    bool isLastStep;

    //! Target parameter value for last step (either maxValue or minValue)
    double targetValue;

    //! Flag indicating if this is an extra step to hit target value
    bool isTargetStep;

    //! Tangent factor
    double tangentFactor;

    //! Minimum tangent factor
    double minTangentFactor;

    //! Tangent factor exponent
    double tangentFactorExponent;
  };
} // end namespace LOCA

#endif

