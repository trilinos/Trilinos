// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorForwardSensitivity_decl_hpp
#define Tempus_IntegratorForwardSensitivity_decl_hpp

// Tempus
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_SensitivityModelEvaluatorBase.hpp"
#include "Tempus_StepperStaggeredForwardSensitivity.hpp"

namespace Tempus {


/** \brief Time integrator implementing forward sensitivity analysis */
/**
 * This integrator implements forward parameter sensitivity analysis by
 * propagating the derivative of the solution with respect to model parameters
 * alongside the solution.  It supports sensitivity propagation methods:
 * <ul>
 *  <li> "Combined" where the sensitivity and state equations are solved
 *       simultaneously.  This is most appropriate for explicit time
 *       integration methods or implicit methods where a very-lightweight
 *       nonlinear solution strategy is used.
 *  <li> "Staggered" where the sensitivity equations are solved at each time
 *       step immediately after the state equations are solved.  This is
 *       useful for implicit methods since it saves the cost of solving the
 *       sensitivity equations while the state equations are being solved.  It
 *       generally does not work for explicit methods (and wouldn't be more
 *       efficient than the combined method even if it did).
 * </ul>
 *
 * Note that this integrator implements all of the same functions as the
 * IntegratorBasic, but is not derived from IntegratorBasic.  It also provides
 * functions for setting the sensitivity initial conditions and extracting the
 * sensitivity at the final time.  Also the vectors stored in the solution
 * history store product vectors of the state and sensitivities using
 * Thyra;:DefaultMultiVectorProductVector.
 */
template<class Scalar>
class IntegratorForwardSensitivity : virtual public Tempus::Integrator<Scalar>
{
public:

  /** \brief Constructor with ParameterList and model, and will be fully
   * initialized. */
  /*!
   * In addition to all of the regular integrator options, the supplied
   * parameter list supports the following options contained within a sublist
   * "Sensitivities" from the top-level parameter list:
   * <ul>
   *   <li> "Sensitivity Method" (default: "Combined") The sensitivity
   *        analysis method as described above.
   *   <li> "Reuse State Linear Solver" (default: false) For the staggered
   *        method, whether to reuse the model's W matrix, solver, and
   *        preconditioner when solving the sensitivity equations.  If they
   *        can be reused, substantial savings in compute time are possible.
   *   <li> "Use DfDp as Tangent" (default:  false) Reinterpret the df/dp
   *        out-arg as the tangent vector (df/dx)(x,p) * dx/dp + df/dp(x,p)
   *        as described in the Tempus::CombinedForwardSensitivityModelEvaluator
   *        documentation.
   *   <li> "Sensitivity Parameter Index" (default: 0) Model evaluator
   *        parameter index for which sensitivities will be computed.
   *   <li> "Sensitivity X Tangent Index" (default: 1) If "Use DfDp as Tangent"
   *        is true, the model evaluator parameter index for passing dx/dp
   *        as a Thyra::DefaultMultiVectorProductVector.
   *   <li> "Sensitivity X-Dot Tangent Index" (default: 2) If
   *        "Use DfDp as Tangent" is true, the model evaluator parameter index
   *        for passing dx_dot/dp as a Thyra::DefaultMultiVectorProductVector.
   *   <li> "Sensitivity X-Dot-Dot Tangent Index" (default: 3) If
   *        "Use DfDp as Tangent" is true, the model evaluator parameter index
   *        for passing dx_dot_dot/dp as a
   *        Thyra::DefaultMultiVectorProductVector (if the model supports
   *        x_dot_dot).
   * </ul>
   */
  IntegratorForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList>                pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

  /** \brief Constructor with model and "Stepper Type" and is fully initialized with default settings. */
  IntegratorForwardSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType);

  /// Destructor
  /** \brief Constructor that requires a subsequent setParameterList, setStepper, and initialize calls. */
  IntegratorForwardSensitivity();

  /// Destructor
  virtual ~IntegratorForwardSensitivity() {}

  /// \name Basic integrator methods
  //@{
    /// Advance the solution to timeMax, and return true if successful.
  virtual bool advanceTime()
    { return integrator_->advanceTime(); }
    /// Advance the solution to timeFinal, and return true if successful.
  virtual bool advanceTime(const Scalar timeFinal) override
    { return integrator_->advanceTime(timeFinal); }
  /// Perform tasks before start of integrator.
  virtual void startIntegrator()
    { integrator_->startIntegrator(); }
  /// Start time step.
  virtual void startTimeStep()
    { integrator_->startTimeStep(); }
  /// Check if time step has passed or failed.
  virtual void checkTimeStep()
    { integrator_->checkTimeStep(); }
  /// Perform tasks after end of integrator.
  virtual void endIntegrator()
    { integrator_->endIntegrator(); }
  /// Return a copy of the Tempus ParameterList
  virtual Teuchos::RCP<Teuchos::ParameterList> getTempusParameterList() override
    { return integrator_->getTempusParameterList(); }
  virtual void setTempusParameterList(Teuchos::RCP<Teuchos::ParameterList> pl) override
    { integrator_->setTempusParameterList(pl); }
  //@}

  /// \name Accessor methods
  //@{
  /// Get current time
  virtual Scalar getTime() const override
    { return integrator_->getTime(); }
  /// Get current index
  virtual Scalar getIndex() const override
    { return integrator_->getIndex(); }
  /// Get Status
  virtual Status getStatus() const override
    { return integrator_->getStatus(); }
  /// Get the Stepper
  virtual Teuchos::RCP<Stepper<Scalar> > getStepper() const override
    { return integrator_->getStepper(); }

  /// Set the Stepper
  virtual void setStepper(Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model);

  /// Set the Stepper
  virtual void setStepperWStepper(Teuchos::RCP<Stepper<Scalar> > stepper)
    { integrator_->setStepperWStepper(stepper); }
  /// Set the initial state which has the initial conditions
  virtual void setInitialState(
    Teuchos::RCP<SolutionState<Scalar> > state = Teuchos::null)
    { integrator_->setInitialState(state); }

  /// Set the initial state from Thyra::VectorBase(s)
  virtual void setInitialState(
    Scalar t0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0 = Teuchos::null,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0 = Teuchos::null,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxDp0 = Teuchos::null,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxdotDp0 = Teuchos::null,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxdotdotDp0 = Teuchos::null);

  /// Get the SolutionHistory
  virtual Teuchos::RCP<const SolutionHistory<Scalar> > getSolutionHistory() const override
    { return integrator_->getSolutionHistory(); }
  /// Set the SolutionHistory
  virtual void setSolutionHistory(
    Teuchos::RCP<SolutionHistory<Scalar> > sh = Teuchos::null)
    { integrator_->setSolutionHistory(sh); }
  /// Get the TimeStepControl
  virtual Teuchos::RCP<const TimeStepControl<Scalar> > getTimeStepControl() const override
    { return integrator_->getTimeStepControl(); }
  /// Set the TimeStepControl
  virtual void setTimeStepControl(
    Teuchos::RCP<TimeStepControl<Scalar> > tsc = Teuchos::null)
    { integrator_->setTimeStepControl(tsc); }
  /// Get the Observer
  virtual Teuchos::RCP<IntegratorObserver<Scalar> > getObserver()
    { return integrator_->getObserver(); }
  /// Set the Observer
  virtual void setObserver(
    Teuchos::RCP<IntegratorObserver<Scalar> > obs = Teuchos::null)
    { integrator_->setObserver(obs); }
  /// Initializes the Integrator after set* function calls
  virtual void initialize()
    { integrator_->initialize(); }
  virtual Teuchos::RCP<Teuchos::Time> getIntegratorTimer() const override
    { return integrator_->getIntegratorTimer(); }
  virtual Teuchos::RCP<Teuchos::Time> getStepperTimer() const override
    { return integrator_->getStepperTimer(); }

  /// Get current the solution, x
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getX() const;
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > getDxDp() const;
  /// Get current the time derivative of the solution, xdot
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getXdot() const;
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > getDxdotDp() const;
  /// Get current the second time derivative of the solution, xdotdot
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getXdotdot() const;
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > getDxdotdotDp() const;

  /// Get current state
  virtual Teuchos::RCP<SolutionState<Scalar> > getCurrentState()
    {return integrator_->getCurrentState();}
  //@}

  /// Parse when screen output should be executed
  void parseScreenOutput() { integrator_->parseScreenOutput(); }

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl)
      override;
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList() override
      { return tempus_pl_; }
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList() override;

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters()
      const override;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    std::string description() const override;
    void describe(Teuchos::FancyOStream        & out,
                  const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

protected:

  // Create sensitivity model evaluator from application model
  void
  createSensitivityModelAndStepper(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model_;
  Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar> > sens_model_;
  Teuchos::RCP<StepperStaggeredForwardSensitivity<Scalar> > sens_stepper_;
  Teuchos::RCP<IntegratorBasic<Scalar> > integrator_;
  Teuchos::RCP<Teuchos::ParameterList> tempus_pl_;
  Teuchos::RCP<Teuchos::ParameterList> sens_pl_;
  Teuchos::RCP<Teuchos::ParameterList> stepper_pl_;
  bool use_combined_method_;
};

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorForwardSensitivity<Scalar> >
integratorForwardSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorForwardSensitivity<Scalar> >
integratorForwardSensitivity(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  std::string stepperType);

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorForwardSensitivity<Scalar> >
integratorForwardSensitivity();

} // namespace Tempus

#endif // Tempus_IntegratorForwardSensitivity_decl_hpp
