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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_NEWSTEPPER_H
#define LOCA_NEWSTEPPER_H

#include "LOCA_Abstract_Iterator.H"    // base class
#include "NOX_Parameter_List.H"	       // class data element
#include "NOX_Solver_Manager.H"	       // class data element
#include "LOCA_Bifurcation_Manager.H"  // class data element
#include "LOCA_MultiContinuation_Manager.H" // class data element
#include "LOCA_StepSize_Manager.H"     // class data element
#include "NOX_StatusTest_Generic.H"    // for enum StatusType
#include "NOX_Common.H"                // for <string> class

// Forward Declarations
namespace LOCA {
  namespace MultiContinuation {
    class ExtendedGroup;
    class AbstractGroup;
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

  class NewStepper : public LOCA::Abstract::Iterator {

  public:

    //! Constructor 
    NewStepper(LOCA::MultiContinuation::AbstractGroup& initialGuess,
	       NOX::StatusTest::Generic& t,
	       NOX::Parameter::List& p);

    //! Destructor
    virtual ~NewStepper();

    //! Reset the Stepper to start a new continuation run
    virtual bool reset(LOCA::MultiContinuation::AbstractGroup& initialGuess,
		       NOX::StatusTest::Generic& t,
		       NOX::Parameter::List& p);

    //! Return the current solution group.
    virtual LOCA::MultiContinuation::AbstractGroup& getSolutionGroup();

    //! Return the output parameters from the stepper algorithm. 
    virtual const NOX::Parameter::List& getParameterList() const;

  protected:

    //! Copy Constructor: Not implemented, but put here to prevent compiler form generating one
    NewStepper(const NewStepper& s);

    //! Intialize iterator
    virtual LOCA::Abstract::Iterator::IteratorStatus start();

    //! Finalize iterator
    virtual LOCA::Abstract::Iterator::IteratorStatus 
    finish(LOCA::Abstract::Iterator::IteratorStatus iteratorStatus);

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

    //! Checks to see if continuation parameter is within threshold of bound
    virtual bool withinThreshold();

  protected:

    //! Bifurcation group manager
    LOCA::Bifurcation::Manager* bifGroupManagerPtr;

    //! Pointer to bifurcation group
    LOCA::MultiContinuation::AbstractGroup* bifGroupPtr;

    //! Continuation group manager
    LOCA::MultiContinuation::Manager* conGroupManagerPtr;

    //! Pointer to the Group with the initial guess for the nonlinear solver.  This is used for resetting the nonlinear solver in between successful steps.
    LOCA::MultiContinuation::ExtendedGroup* curGroupPtr;

    //! Pointer to the a Group with the solution for the previous continuation step.  This is used for resetting the nonlinear solver in between failed steps.
    LOCA::MultiContinuation::ExtendedGroup* prevGroupPtr;

    //! Pointer to the status test.  
    NOX::StatusTest::Generic* statusTestPtr;

    //! Pointer to parameter list passed in constructor/reset method
    NOX::Parameter::List* paramListPtr;

    //! Solver (a.k.a corrector)
    NOX::Solver::Manager* solverPtr;

    //! Current predictor direction
    LOCA::MultiContinuation::ExtendedVector* curPredictorPtr;

    //! previous predictor direction
    LOCA::MultiContinuation::ExtendedVector* prevPredictorPtr;

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

    //! Target parameter value for last step (either maxValue or minValue)
    double targetValue;

    //! Flag indicating if this is an extra step to hit target value
    bool isTargetStep;

    //! Flag indicating whether to do tangent factor step size scaling
    bool doTangentFactorScaling;

    //! Tangent factor
    double tangentFactor;

    //! Minimum tangent factor
    double minTangentFactor;

    //! Tangent factor exponent
    double tangentFactorExponent;

    //! Flag indicating whether to compute eigenvalues after eachs step
    bool calcEigenvalues;
  };
} // end namespace LOCA

#endif

