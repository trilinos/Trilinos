//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperExponential_decl_hpp
#define Tempus_StepperExponential_decl_hpp

// Tempus
#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_PhiEvaluator.hpp"

/**
 * @brief Per-evaluation metadata passed from an exponential stepper to its model.
 *
 * The `Scalar` timeStepSize_ is forwarded as the Thyra step size, and the
 * `int` stageNumber_ is forwarded when the model supports a stage number.
 */
template <class Scalar>
class ExponentialODEParameters {
 public:
  /** @brief Construct metadata with a zero `Scalar` time step and stage zero. */
  ExponentialODEParameters() : timeStepSize_(Scalar(0.0)), stageNumber_(0) {}

  /**
   * @brief Construct metadata for one step or internal stage.
   *
   * @param timeStepSize `Scalar` step size supplied to the model.
   * @param stageNumber `int` internal-stage index supplied to supporting models;
   *   defaults to zero.
   */
  ExponentialODEParameters(Scalar timeStepSize, int stageNumber = 0)
    : timeStepSize_(timeStepSize), stageNumber_(stageNumber)
  {
  }

  /// `Scalar` time-step size supplied through Thyra InArgs.
  Scalar timeStepSize_;
  /// `int` internal-stage index supplied through Thyra InArgs when supported.
  int stageNumber_;
};


namespace Tempus {

/** \brief Base class for exponential steppers using an implicit Thyra model.
 *
 * The application model supplies the first-order implicit residual
 *  \f[
 *    \mathcal{F}_{\mathrm{impl}}(\dot{x},x,t) = 0,
 *  \f]
 * where \f$M = \partial\mathcal{F}_{\mathrm{impl}}/\partial\dot{x}\f$.
 * Exponential methods operate on the explicit tendency
 * \f$\mathcal{F}(x,t)=-M^{-1}\mathcal{F}_{\mathrm{impl}}(0,x,t)\f$ and
 * its Jacobian.  The PhiEvaluator performs the required mass operations
 * without explicitly forming \f$M^{-1}\f$.
 *
 * Parameter lists include the inherited stepper settings, a `PhiEvaluator`
 * sublist used to construct the evaluator, `Epsilon for RHS finite difference`
 * for the nonautonomous RHS time-derivative finite difference (positive values
 * enable it), and `Adapt PhiEvaluator Interval`.  A positive interval requests
 * evaluator adaptation on the first step and then at matching step indices;
 * a nonpositive value disables adaptation.
 *
 */
template <class Scalar>
class StepperExponential : virtual public Tempus::Stepper<Scalar> {
 public:
  /// \name Basic exponential stepper methods
  //@{

  /**
   * @brief Construct an unconfigured exponential-stepper base.
   *
   * Set a PhiEvaluator and application model, then initialize the derived
   * stepper before taking steps.
   */
  StepperExponential();

  /**
   * @brief Set the evaluator for phi-function products and mass operations.
   *
   * When a model is already set, the `Teuchos::RCP<PhiEvaluator<Scalar>>` is
   * given that model and initialized.  This invalidates the stepper setup.
   *
   * @param phiEvaluator Evaluator to associate with this stepper.
   */
  virtual void setPhiEvaluator(
    const Teuchos::RCP<Tempus::PhiEvaluator<Scalar> >& phiEvaluator);
  /** @brief Construct, set, and invalidate setup for the factory-default PhiEvaluator. */
  virtual void setDefaultPhiEvaluator();
  /** @brief Return the configured `Teuchos::RCP<PhiEvaluator<Scalar>>`, or null if unset. */
  virtual Teuchos::RCP<Tempus::PhiEvaluator<Scalar> > getPhiEvaluator() const;

  /**
   * @brief Set the implicit application model used for residual evaluations.
   *
   * The `Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>` is validated as an
   * implicit ODE/DAE model.  If a PhiEvaluator is present, it receives and
   * initializes with this model; setup is then invalidated.
   *
   * @param appModel Application ModelEvaluator defining the implicit residual.
   */
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
      override;
  /** @brief Return the configured implicit application ModelEvaluator. */
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel()
      const override;

  /**
   * @brief Set initial state data and apply the configured consistency policy.
   *
   * A null state vector is copied from the model nominal values.  For
   * `ICConsistency = "Consistent"`, the history must store xDot and this method
   * computes it from \f$M\dot{x}=-\mathcal{F}_{\mathrm{impl}}(0,x,t)\f$.
   * Other inherited policies select no synchronization, a zero derivative, or
   * the model nominal derivative.  The optional consistency check
   * reports, but does not reject, a residual mismatch.
   *
   * @param solutionHistory History containing at least one initial SolutionState.
   */
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

  /**
   * @brief Build model input arguments for an exponential residual evaluation.
   *
   * The result carries caller-supplied `x`, `xDot`, and `time`, plus
   * ExponentialODEParameters metadata where supported by the model.  It
   * sets \f$\alpha=0\f$ and \f$\beta=1\f$ so model linearizations represent
   * the residual Jacobian with respect to \f$x\f$.
   *
   * @param x `Teuchos::RCP<Thyra::VectorBase<Scalar>>` state vector.
   * @param xDot `Teuchos::RCP<Thyra::VectorBase<Scalar>>` derivative vector.
   * @param time `Scalar` evaluation time.
   * @param p Metadata containing the step size and stage number.
   * @return Thyra InArgs configured for the application ModelEvaluator.
   */
  Thyra::ModelEvaluatorBase::InArgs<Scalar>
  virtual createInArgsExponentialODE(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
      const Teuchos::RCP<ExponentialODEParameters<Scalar> >& p);

  /**
   * @brief Evaluate the implicit residual used to obtain the exponential RHS.
   *
   * Sets `xDot` to zero and evaluates
   * \f$f=\mathcal{F}_{\mathrm{impl}}(0,x,t)=-M\mathcal{F}(x,t)\f$.  The
   * PhiEvaluator converts this mass-weighted residual to the explicit tendency
   * or applies phi functions without forming the inverse mass matrix.
   *
   * @param f Output `Teuchos::RCP<Thyra::VectorBase<Scalar>>` residual vector.
   * @param x State vector passed to the application model.
   * @param xDot Derivative workspace, overwritten with zeros before evaluation.
   * @param time `Scalar` evaluation time.
   * @param p Metadata containing the step size and stage number.
   */
  virtual void evaluateExponentialODE(
      Teuchos::RCP<Thyra::VectorBase<Scalar> >& f,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
      const Teuchos::RCP<ExponentialODEParameters<Scalar> >& p);

  /**
   * @brief Accept an implicit-solver initial guess without using it.
   *
   * Exponential steppers do not perform Newton solves.
   *
   * @param initialGuess Unused `Teuchos::RCP<const Thyra::VectorBase<Scalar>>`.
   */
  virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > initialGuess) override
  {
  }

  /**
   * @brief Return a large `Scalar` initial-step estimate.
   */
  virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */)
      const override
  {
    // return a large value that should still fit into any supported scalar type
    return Scalar(Teuchos::ScalarTraits<Scalar>::rmax() / 1e2);
  }

  /** @brief Create a `StepperState<Scalar>` initialized with this stepper type. */
  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState() override;
  /**
   * @brief Set First-Same-As-Last use and invalidate setup.
   *
   *  The FSAL principle does generally not apply to exponential integrators.
   *  If it is set to true here, we will treat the method as if it was FSAL, by evaluating and
   *  storing the first stage evaluation of the future next step at the end of takeStep.
   *
   * @param a `bool` value stored as the inherited FSAL setting.
   */
  virtual void setUseFSAL(bool a) override
  {
    this->useFSAL_       = a;
    this->isInitialized_ = false;
  }

  /** @brief Return false because the stepper is backed by an implicit model. */
  virtual bool isExplicit() const override {return false;}
  /** @brief Return true because residual and linearization data come from an implicit model. */
  virtual bool isImplicit() const override {return true;}
  /** @brief Return false because this is not an explicit-implicit stepper. */
  virtual bool isExplicitImplicit() const override
    {return isExplicit() && isImplicit();}
  /** @brief Return FIRST_ORDER_ODE for the supported model equation order. */
  virtual OrderODE getOrderODE() const override {return FIRST_ORDER_ODE;}

  /**
   * @brief Check base setup and require both a model and PhiEvaluator.
   *
   * @param out Stream receiving diagnostics for missing dependencies.
   * @return True only when the base setup is valid and both dependencies are set.
   */
  virtual bool isValidSetup(Teuchos::FancyOStream& out) const override;

  /**
   * @brief Return valid parameters for the exponential-stepper base.
   *
   * Includes inherited basic settings, the `PhiEvaluator` sublist, and the
   * exponential finite-difference and adaptation settings.
   *
   * @return Valid `Teuchos::ParameterList` with current defaults.
   */
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters()
      const override;
  //@}

  /**
   * @brief Build the nonconst list of valid and default parameters shared by exponential steppers.
   *
   * `Epsilon for RHS finite difference` is a `double` finite-difference scale
   * for nonautonomous RHS corrections for Rosenbrock type integrators;
   * values at or below zero disable them.
   * `Adapt PhiEvaluator Interval` is an `int` adaptation interval; only positive
   * values enable adaptation.  The `PhiEvaluator` sublist selects and configures
   * the evaluator.  Some compatibility implicit solver and predictor entries are also accepted.
   *
   * @return Mutable valid `Teuchos::ParameterList` with current defaults.
   */
  Teuchos::RCP<Teuchos::ParameterList> getValidParametersBasicExponential() const;

  /// \name Overridden from Teuchos::Describable
  //@{
  /**
   * @brief Write base-stepper and PhiEvaluator descriptions to an output stream.
   *
   * @param out `Teuchos::FancyOStream` receiving the description.
   * @param verbLevel Requested Teuchos verbosity level, forwarded to dependencies.
   */
  virtual void describe(
      Teuchos::FancyOStream& out,
      const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  /**
   * @brief Validate and apply exponential-stepper parameter settings.
   *
   * If `PhiEvaluator` is present, its sublist is used to construct a new
   * evaluator before validation.  The method then applies inherited settings,
   * the `double` finite-difference epsilon, and the `int` adaptation interval.
   *
   * @param pl Parameter list to validate and apply; a null RCP leaves settings unchanged.
   */
  void setStepperExponentialValues(Teuchos::RCP<Teuchos::ParameterList> pl);

 protected:
  /// @brief Compute the temporal finite difference dt_Mf_deriv
  /// 
  /// dt_Mf_deriv = -dt * M * [d/dt F(x,t)]
  void computeTemporalFD(
    Teuchos::RCP<Thyra::VectorBase<Scalar>>& dt_Mf_deriv,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x,
    const Scalar t0,
    const Scalar dt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mf
  );

  /// Check if temporal derivative correction is desired for Rosenbrock integrators
  bool getTemporalDerivative()
  {
    return (temporal_finite_difference_eps_ > 0);
  }

  /// @brief Compute the nonlinear remainder R
  /// 
  ///   R = -M * (F(xr,tr) - F(x0,t0) - J_{x0} * (xr-x0) - F'(t0) * (tr-t0))
  /// including multiple of negative mass matrix (-M).
  ///
  /// @param remf is storage for the remainder R
  /// @param dt is the current time-step, not necessarily (tr-t0)
  /// @param Mf contains already evaluated -M*F(x0,t0)
  /// @param dt_Mf_deriv contains already evaluated dt*M*F'(t0)
  /// @param Mfr can optionally contain -M*F(xr,tr), if already pre-evaluated
  void computeRemf(
    Teuchos::RCP<Thyra::VectorBase<Scalar>>& remf,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& xr,
    const Scalar tr,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x0,
    const Scalar t0,
    const Scalar dt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mf,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& dt_Mf_deriv,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mfr = Teuchos::null
  );

  /// Check if PhiEvaluator adaptivity is desired and return positive interval number
  int getAdaptPhiEvaluator()
  {
    return adapt_phi_evaluator_interval_;
  }

 private:
  /// RCP to the application provided ModelEvaluator
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel_;

  /// RCP to the PhiEvaluator
  Teuchos::RCP<PhiEvaluator<Scalar> > phiEvaluator_;

  /// Finite difference step size used for RHS time derivative estimation
  /// needed for nonautonomous correction.
  double temporal_finite_difference_eps_;

  /// Number of time steps to wait between adapt PhiEvaluator calls
  int adapt_phi_evaluator_interval_;

};

}  // namespace Tempus
#endif  // Tempus_StepperExponential_decl_hpp
