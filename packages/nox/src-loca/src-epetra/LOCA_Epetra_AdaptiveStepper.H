// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRAADAPTIVESTEPPER_H
#define LOCA_EPETRAADAPTIVESTEPPER_H

#include "Teuchos_RCP.hpp"

#include "LOCA_Abstract_Iterator.H"    // base class
#include "Teuchos_ParameterList.hpp"           // class data element
#include "NOX_Solver_Generic.H"           // class data element
#include "NOX_StatusTest_Generic.H"    // for enum StatusType
#include "NOX_Common.H"                // for <string> class

#include "LOCA_Epetra_AdaptiveSolutionManager.H"

namespace LOCA {

namespace Epetra {

 /*!
   * \brief Implementation of LOCA::Abstract::Iterator for computing points
   * along a continuation curve.
   *
   * The %AdaptiveStepper class implements the pure virtual methods of the
   * LOCA::Abstract::Iterator for iteratively computing points along a
   * continuation curve.
   */


  class AdaptiveStepper : public LOCA::Abstract::Iterator {

  public:

    AdaptiveStepper(const Teuchos::RCP<Teuchos::ParameterList>& pList,
      const Teuchos::RCP<LOCA::Epetra::AdaptiveSolutionManager>& mgr,
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<NOX::StatusTest::Generic>& nt);

    //! Destructor
    virtual ~AdaptiveStepper();

    //! Replaces the eigensolver parameter list.
    bool
    eigensolverReset( Teuchos::RCP<Teuchos::ParameterList> & newEigensolverList );

    //! Return the current solution group.
    virtual Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
    getSolutionGroup() const;

    //! Return the current bifurcation group.
    /*!
     * If the current bifurcation method is "None", then the returned group
     * is the same as getSolutionGroup(), otherwise this method returns the
     * current bifurcation group (e.g., a turning point group).
     */
    virtual Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
    getBifurcationGroup() const;

    //! Return the output parameters from the stepper algorithm.
    virtual Teuchos::RCP<const Teuchos::ParameterList>
    getList() const;

    /*! @brief Return the current nonlinear solver pointer.
     *
     * Will throw an error if the solver does not exist yet.
     */
    virtual Teuchos::RCP<const NOX::Solver::Generic>
    getSolver() const;

    /*! @brief Return the current continuation parameter from the underlying LOCA::MultiContinuation::AbstractStrategy.
     *
     */
    virtual double
    getContinuationParameter() const;

    //! Run the iterator
    virtual LOCA::Abstract::Iterator::IteratorStatus run();

  protected:

    void setSolutionGroup(const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
      const double value);

    //! Build LOCA Factory needed to support stepping strategy
    void buildLOCAFactory();

    //! Initialize iterator
    virtual LOCA::Abstract::Iterator::IteratorStatus start();

    //! Finalize iterator
    virtual LOCA::Abstract::Iterator::IteratorStatus
      finish(LOCA::Abstract::Iterator::IteratorStatus iteratorStatus);

    //! Preprocess step
    virtual LOCA::Abstract::Iterator::StepStatus
      preprocess(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Adaptation step
    virtual LOCA::Abstract::Iterator::StepStatus
      adapt(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Compute step
    virtual LOCA::Abstract::Iterator::StepStatus
      compute(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Postprocess step
    virtual LOCA::Abstract::Iterator::StepStatus
      postprocess(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Check stopping criteria.
    LOCA::Abstract::Iterator::IteratorStatus
      stop(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Perform iteration
    virtual LOCA::Abstract::Iterator::IteratorStatus iterate();

    //! Builds a constrained group is specified through parameter lists
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
      buildConstrainedGroup(
        const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp);

    //! Compute a new step size based on results from the previous step
    LOCA::Abstract::Iterator::StepStatus
      computeStepSize(LOCA::Abstract::Iterator::StepStatus stepStatus,
            double& stepSize);

    //! Print to the screen the initialization information
    virtual void printInitializationInfo();

    //! Print to the screen information for the beginning of a new continuation step
    virtual void printStartStep();

    //! Print to the screen information for the beginning of a restart relaxation step
    virtual void printRelaxationStep();

    //! Print to the screen information for the end of a successful continuation step
    virtual void printEndStep(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Print to the screen information for the end of a relaxation step
    virtual void printRelaxationEndStep(LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Print to the screen the concluding information
    virtual void printEndInfo();

    //! Checks to see if continuation parameter is within threshold of bound
    virtual bool withinThreshold();

  private:

    //! Private to prohibit copying
    AdaptiveStepper(const AdaptiveStepper&);

    //! Private to prohibit copying
    AdaptiveStepper& operator = (const AdaptiveStepper&);

  protected:

    //! Adaptive solution manager provides interface to remesh/resize functions
    Teuchos::RCP<LOCA::Epetra::AdaptiveSolutionManager> mgr;

    //! Global data
    Teuchos::RCP<LOCA::GlobalData> globalData;

    //! Parsed parameter list
    Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

    //! Predictor strategy
    Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy> predictor;

    //! Current continuation group
    Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy> curGroupPtr;

    //! Continuation group from previous step
    Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy> prevGroupPtr;

    //! Eigensolver for computing eigenvalues after each continuation step
    Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy> eigensolver;

    //! Strategy for saving eigenvectors/values
    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> saveEigenData;

    //! Pointer to bifurcation group
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> bifGroupPtr;

    //! Pointer to the NOX status test.
    Teuchos::RCP<NOX::StatusTest::Generic> noxStatusTestPtr;

    //! Pointer to parameter list passed in constructor/reset method
    Teuchos::RCP<Teuchos::ParameterList> paramListPtr;

    //! Pointer to stepper parameters
    Teuchos::RCP<Teuchos::ParameterList> stepperList;

    //! Solver (a.k.a corrector)
    Teuchos::RCP<NOX::Solver::Generic> solverPtr;

    //! Current predictor direction
    Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> curPredictorPtr;

    //! previous predictor direction
    Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> prevPredictorPtr;

    //! Step size control strategy
    Teuchos::RCP<LOCA::StepSize::AbstractStrategy> stepSizeStrategyPtr;

    //! Other parameter lists
    Teuchos::RCP<Teuchos::ParameterList> bifurcationParams;
    Teuchos::RCP<Teuchos::ParameterList> firstStepperParams;

    //! Name of continuation parameter
    std::string conParamName;

    //! Index of continuation parameter
    std::vector<int> conParamIDs;

    //! Starting value of continuation parameter, \f$\lambda_b \f$.
    double startValue;

    //! Largest allowable value of continuation parameter
    double maxValue;

    //! Smallest allowable value of continuation parameter
    double minValue;

    //! Current step size (change in the continuation parameter), \f$\Delta \lambda_i \f$.
    double stepSize;

    //! Maximum number of Newton iterations per continuation step,\f$ N_{max} \f$. Defaults to 15.
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

    //! Flag indicating whether to compute eigenvalues after each step
    bool calcEigenvalues;

    //! Flag indicating whether to return failed upon reaching max steps
    bool return_failed_on_max_steps;

    //! Flag indicating that the failed step is only due to the maximum number of steps being exceeded
    bool max_steps_exceeded;

  };
} // end namespace Epetra
} // end namespace LOCA

#endif

