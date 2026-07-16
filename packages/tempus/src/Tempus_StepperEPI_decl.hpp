//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPI_decl_hpp
#define Tempus_StepperEPI_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExponential.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperEPIAppAction.hpp"

#include "Tempus_PhiEvaluator.hpp"



namespace Tempus {
/**
 * @brief Base for exponential propagation iterative (EPI) multistep steppers.
 *
 * For an implicit model residual
 * \f[
 *   \mathcal{F}_{\mathrm{impl}}(\dot{x},x,t)=0,
 * \f]
 * this stepper obtains the explicit tendency and its current Jacobian as
 * \f[
 *   \mathcal{F}_n = \mathcal{F}(x_n,t_n)=-M^{-1}\mathcal{F}_{\mathrm{impl}}(0,x_n,t_n),
 * \f]
 * and
 * \f[
 *   J_n=-M^{-1}\frac{\partial\mathcal{F}_{\mathrm{impl}}}{\partial x}(0,x_n,t_n).
 * \f]
 * The evaluator receives the corresponding mass-weighted vectors, including
 * the explicit tendency \f$-M\mathcal{F}_n\f$,
 * the nonautonomous temporal correction \f$-\Delta t\,M\mathcal{F}'_n\f$, and the
 * mass-weighted nonlinear remainder, and performs inverse-mass operations
 * internally.
 *
 * EPI2 uses the exponential Rosenbrock-Euler update
 * \f[
 *   x_{n+1}=x_n+
 *     \Delta t\,\varphi_1(\Delta t J_n)\mathcal{F}_n+
 *     (\Delta t)^2\varphi_2(\Delta t J_n)\mathcal{F}'_n.
 * \f]
 * The final term can be omitted for autonomous systems and disabled by omitting
 * the configured finite difference estimates \f$\mathcal{F}'_n\f$ through the parameter list.
 *
 * EPI3 additionally uses the preceding state through the nonlinear remainder
 * \f[
 *   R_n(x,t)=\mathcal{F}(x,t)-\mathcal{F}_n-
 *     J_n(x-x_n)-\mathcal{F}'_n(t-t_n)
 * \f]
 * and updates by
 * \f[
 *   x_{n+1}=x_n+
 *     \Delta t\,\varphi_1(\Delta t J_n)\mathcal{F}_n+
 *     (\Delta t)^2\varphi_2(\Delta t J_n)\mathcal{F}'_n+
 *     \frac{2}{3}\Delta t\,\varphi_2(\Delta t J_n)
 *     R_n(x_{n-1},t_{n-1}).
 * \f]
 * EPI3 uses the EPI2 update on its first step because the preceding state is
 * unavailable.  In contrast to exponential Euler, EPI2 is second order for
 * non-autonomous problems (unless the temporal derivative is disabled),
 * while EPI3 uses the additional remainder term to obtain third order. Similarly to
 * EPI2, it uses the temporal derivative estimate for non-autonomous problems, and
 * will degrade to first order if this is incorectly disabled.
 *
 * The inherited `Epsilon for RHS finite difference` parameter enables a
 * finite-difference time-derivative correction when positive; a nonpositive
 * value disables it.  With `Use FSAL` and stored, synchronized xDot data, the
 * stepper reuses the current tendency and stores the next one for a later step.
 * A StepperEPIAppAction callback is invoked at BEGIN_STEP, BEFORE_EXP,
 * AFTER_EXP, and END_STEP.
 *
 * @tparam Scalar Scalar type used by the Thyra model, solution history, and
 *   phi-function evaluator.
 */
template<class Scalar>
class StepperEPI :
    virtual public Tempus::StepperExponential<Scalar>
{
public:

  /**
   * @brief Construct an EPI stepper with order two, a default PhiEvaluator, and
   * a default no-op application action.
   *
   * Set an application model and initialize the stepper before taking a step.
   */
  StepperEPI();

  /**
   * @brief Construct and initialize an EPI stepper with supplied dependencies.
   *
   * @param appModel RCP to the implicit Thyra model that defines the residual.
   * @param phiEvaluator RCP that computes phi-function products and mass
   *   operations for `appModel`.
   * @param useFSAL `bool` enabling reuse and storage of synchronized xDot data.
   * @param ICConsistency Initial-condition consistency policy.
   * @param ICConsistencyCheck `bool` requesting the inherited consistency check.
   * @param stepperEPIAppAction RCP callback; null selects the default no-op
   *   action.
   */
  StepperEPI(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>& appModel,
    const Teuchos::RCP<Tempus::PhiEvaluator<Scalar>>& phiEvaluator,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction);

  /// \name Basic stepper methods
  //@{
    /**
     * @brief Set callbacks executed at EPI step locations.
     *
     * A null `Teuchos::RCP<StepperEPIAppAction<Scalar>>` installs the default
     * no-op action and invalidates stepper setup.
     *
     * @param appAction Callback RCP, or null for the default action.
     */
    virtual void setAppAction(
      Teuchos::RCP<StepperEPIAppAction<Scalar> > appAction);

    /** @brief Return the configured EPI application-action RCP. */
    virtual Teuchos::RCP<StepperEPIAppAction<Scalar> > getAppAction() const
    { return stepperEPIAppAction_; }

    /**
     * @brief Set the `Scalar` temporal order used to select EPI2 or EPI3 logic.
     *
     * Values greater than two select the EPI3 path when sufficient history is
     * available; the value is stored internally as `double`.
     *
     * @param order Requested temporal order.
     */
    void setOrder(Scalar order) {this->order_ = order;}

    /**
     * @brief Check base exponential-stepper setup and require an application action.
     *
     * @param out `Teuchos::FancyOStream` receiving setup diagnostics.
     * @return True when the inherited setup and EPI callback are valid.
     */
    bool isValidSetup(Teuchos::FancyOStream & out) const;

    /**
     * @brief Advance the working state by its configured time step.
     *
     * The PhiEvaluator is linearized at the current state.  EPI2 requires at
     * least two SolutionStates (current=NM1 and working=N).
     * EPI3 requires history storage for at least
     * three states (NM2, current=NM1 and working=N) and uses EPI2 on its first step,
     * before the prior state NM2 is available.
     * The working-state solution status is set from the
     * PhiEvaluator result.  FSAL caching is used only when enabled and xDot is
     * present and synchronized.
     *
     * @param solutionHistory RCP containing current and working states, plus
     *   one more history entry when required.
     */
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

    /** @brief Return the configured `Scalar` temporal order. */
    virtual Scalar getOrder() const override {return Scalar(order_);}
    /** @brief Return the minimum supported `Scalar` order, two. */
    virtual Scalar getOrderMin() const override {return Scalar(2.0);}
    /** @brief Return the maximum supported `Scalar` order, three. */
    virtual Scalar getOrderMax() const override {return Scalar(3.0);}

    /** @brief Return false because EPI methods use solution-history data. */
    virtual bool isOneStepMethod() const override {return false;}
    /** @brief Return true because EPI methods are classified as multistep. */
    virtual bool isMultiStepMethod() const override {return !isOneStepMethod();}
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    /**
     * @brief Write inherited and EPI callback configuration to a stream.
     *
     * @param out `Teuchos::FancyOStream` receiving the description.
     * @param verbLevel Requested Teuchos verbosity level.
     */
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

 private:

  Teuchos::RCP<StepperEPIAppAction<Scalar> > stepperEPIAppAction_;

  /// temporal Integration order
  double order_;
};


/**
 * @brief Second-order EPI method.
 *
 * EPI2 uses the current tendency and Jacobian in the EPI2/exponential
 * Rosenbrock-Euler update.  For non-autonomous problems it is second order due to the
 * additional temporal correction and is distinct from Exponential Euler,
 * which lacks this correction and can use an inexact linear operator
 * and is generally first order.  A positive `Epsilon for RHS finite
 * difference` adds the finite-difference temporal nonautonomous correction.
 *
 * @tparam Scalar Scalar type used by the base EPI stepper.
 */
template <class Scalar>
class StepperExponential_EPI2 : virtual public StepperEPI<Scalar> {
 public:
  /** @brief Construct an unconfigured EPI2 stepper with default settings. */
  StepperExponential_EPI2()
  {
    this->setStepperName("EPI2");
    this->setStepperType("EPI2");
    this->setUseFSAL(false);
    this->setICConsistency("Consistent");
    this->setICConsistencyCheck(false);
    this->setAppAction(Teuchos::null);
    this->setOrder(2.0);
  }

  /**
   * @brief Construct and initialize EPI2 with an application model.
   *
   * @param appModel RCP to the implicit Thyra model; a null RCP leaves the
   *   stepper uninitialized.
   * @param useFSAL `bool` enabling optional xDot caching.
   * @param ICConsistency Initial-condition consistency policy.
   * @param ICConsistencyCheck `bool` requesting the inherited consistency check.
   * @param stepperEPIAppAction RCP callback; null selects the default action.
   */
  StepperExponential_EPI2(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction)
  {
    this->setStepperName("EPI2");
    this->setStepperType("EPI2");
    this->setUseFSAL(useFSAL);
    this->setICConsistency(ICConsistency);
    this->setICConsistencyCheck(ICConsistencyCheck);
    this->setOrder(2.0);

    this->setAppAction(stepperEPIAppAction);

    if (appModel != Teuchos::null) {
      this->setModel(appModel);
      this->initialize();
    }
  }
};

/**
 * @brief Create and optionally initialize an EPI2 stepper from parameters.
 *
 * The `Teuchos::ParameterList` configures inherited exponential settings,
 * including the PhiEvaluator sublist and optional nonautonomous correction.
 * A nonnull model is set and followed by initialization; a null model returns
 * an uninitialized stepper.
 *
 * @tparam Scalar Scalar type used by the model and stepper.
 * @param model RCP to the implicit Thyra application model, or null.
 * @param pl RCP to the EPI2 parameter list, or null to retain defaults.
 * @return RCP owning the configured `StepperExponential_EPI2<Scalar>`.
 */
template <class Scalar>
Teuchos::RCP<StepperExponential_EPI2<Scalar> >
createStepperExponential_EPI2(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperExponential_EPI2<Scalar>());
  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}


/**
 * @brief Third-order EPI multistep method.
 *
 * EPI3 augments the EPI2 phi-function update with the nonlinear remainder
 * formed from the preceding solution and the current Jacobian.  It therefore
 * requires storage for three SolutionStates.  On the initial step, when only
 * two states exist, it deliberately uses the EPI2 update to seed the history.
 * The optional positive `Epsilon for RHS finite difference` supplies the
 * nonautonomous correction required for non-autonomous problems;
 * optional FSAL caching reuses synchronized xDot
 * data when available.
 *
 * @tparam Scalar Scalar type used by the base EPI stepper.
 */
template <class Scalar>
class StepperExponential_EPI3 : virtual public StepperEPI<Scalar> {
 public:
  /** @brief Construct an unconfigured EPI3 stepper with default settings. */
  StepperExponential_EPI3()
  {
    this->setStepperName("EPI3");
    this->setStepperType("EPI3");
    this->setUseFSAL(false);
    this->setICConsistency("Consistent");
    this->setICConsistencyCheck(false);
    this->setAppAction(Teuchos::null);
    this->setOrder(3.0);
  }

  /**
   * @brief Construct and initialize EPI3 with an application model.
   *
   * @param appModel RCP to the implicit Thyra model; a null RCP leaves the
   *   stepper uninitialized.
   * @param useFSAL `bool` enabling optional xDot caching.
   * @param ICConsistency Initial-condition consistency policy.
   * @param ICConsistencyCheck `bool` requesting the inherited consistency check.
   * @param stepperEPIAppAction RCP callback; null selects the default action.
   */
  StepperExponential_EPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction)
  {
    this->setStepperName("EPI3");
    this->setStepperType("EPI3");
    this->setUseFSAL(useFSAL);
    this->setICConsistency(ICConsistency);
    this->setICConsistencyCheck(ICConsistencyCheck);
    this->setOrder(3.0);

    this->setAppAction(stepperEPIAppAction);

    if (appModel != Teuchos::null) {
      this->setModel(appModel);
      this->initialize();
    }
  }
};

/**
 * @brief Create and optionally initialize an EPI3 stepper from parameters.
 *
 * The `Teuchos::ParameterList` configures inherited exponential settings,
 * including the PhiEvaluator sublist and optional nonautonomous correction.
 * A nonnull model is set and followed by initialization; a null model returns
 * an uninitialized stepper.  Configure SolutionHistory storage for three
 * states before taking EPI3 steps.
 *
 * @tparam Scalar Scalar type used by the model and stepper.
 * @param model RCP to the implicit Thyra application model, or null.
 * @param pl RCP to the EPI3 parameter list, or null to retain defaults.
 * @return RCP owning the configured `StepperExponential_EPI3<Scalar>`.
 */
template <class Scalar>
Teuchos::RCP<StepperExponential_EPI3<Scalar> >
createStepperExponential_EPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperExponential_EPI3<Scalar>());
  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

} // namespace Tempus

#endif // Tempus_StepperEPI_decl_hpp
