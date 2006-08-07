// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_STEPPER_H
#define LOCA_STEPPER_H

#include "Teuchos_RefCountPtr.hpp"

#include "LOCA_Abstract_Iterator.H"    // base class
#include "Teuchos_ParameterList.hpp"	       // class data element
#include "NOX_Solver_Manager.H"	       // class data element
#include "NOX_StatusTest_Generic.H"    // for enum StatusType
#include "NOX_Common.H"                // for <string> class

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace MultiContinuation {
    class AbstractStrategy;
    class AbstractGroup;
    class ExtendedVector;
  }
  namespace MultiPredictor {
    class AbstractStrategy;
  }
  namespace Eigensolver {
    class AbstractStrategy;
  }
  namespace SaveEigenData {
    class AbstractStrategy;
  }
  namespace StepSize {
    class AbstractStrategy;
  }
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  /*!
   * \brief Implementation of LOCA::Abstract::Iterator for computing points
   * along a continuation curve.
   */
  /*!
   * The %Stepper class implements the pure virtual methods of the
   * LOCA::Abstract::Iterator for iteratively computing points along a
   * continuation curve.
   */

  class Stepper : public LOCA::Abstract::Iterator {

  public:

    //! Constructor
    Stepper(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	    const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
	    const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t,
	    const Teuchos::RefCountPtr<Teuchos::ParameterList>& p);

    //! Destructor
    virtual ~Stepper();

    //! Reset the Stepper to start a new continuation run
    virtual bool 
    reset(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	  const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
	  const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t,
	  const Teuchos::RefCountPtr<Teuchos::ParameterList>& p);

    //! Return the current solution group.
    virtual Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup> 
    getSolutionGroup() const;

    //! Return the current bifurcation group.
    /*!
     * If the current bifurcation method is "None", then the returned group
     * is the same as getSolutionGroup(), otherwise this method returns the
     * current bifurcation group (e.g., a turning point group).
     */
    virtual Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup> 
    getBifurcationGroup() const;

    //! Return the output parameters from the stepper algorithm.
    virtual Teuchos::RefCountPtr<const Teuchos::ParameterList> 
    getList() const;

    /*! @brief Return the current nonlinear solver pointer.
     *
     * Will throw an error if the solver does not exist yet.
     */
    virtual Teuchos::RefCountPtr<const NOX::Solver::Generic> 
    getSolver() const;

  protected:

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

    //! Builds a constrained group is specified through parameter lists
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
    buildConstrainedGroup(
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp);

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

  private:
    
    //! Private to prohibit copying
    Stepper(const Stepper&);

    //! Private to prohibit copying
    Stepper& operator = (const Stepper&);

  protected:

    //! Global data
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

    //! Parsed parameter list
    Teuchos::RefCountPtr<LOCA::Parameter::SublistParser> parsedParams;

    //! Predictor strategy
    Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy> predictor;

    //! Current continuation group
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy> curGroupPtr;

    //! Continuation group from previous step
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy> prevGroupPtr;

    //! Eigensolver for computing eigenvalues after each continuation step
    Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy> eigensolver;

    //! Strategy for saving eigenvectors/values
    Teuchos::RefCountPtr<LOCA::SaveEigenData::AbstractStrategy> saveEigenData;

    //! Pointer to bifurcation group
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> bifGroupPtr;

    //! Pointer to the status test.
    Teuchos::RefCountPtr<NOX::StatusTest::Generic> statusTestPtr;

    //! Pointer to parameter list passed in constructor/reset method
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramListPtr;

    //! Pointer to stepper parameters
    Teuchos::RefCountPtr<Teuchos::ParameterList> stepperList;

    //! Solver (a.k.a corrector)
    Teuchos::RefCountPtr<NOX::Solver::Manager> solverPtr;

    //! Current predictor direction
    Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedVector> curPredictorPtr;

    //! previous predictor direction
    Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedVector> prevPredictorPtr;

    //! Step size control strategy
    Teuchos::RefCountPtr<LOCA::StepSize::AbstractStrategy> stepSizeStrategyPtr;

    //! Name of continuation parameter
    string conParamName;

    //! Index of continuation parameter
    vector<int> conParamIDs;

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

