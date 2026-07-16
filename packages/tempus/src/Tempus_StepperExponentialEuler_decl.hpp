//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperExponentialEuler_decl_hpp
#define Tempus_StepperExponentialEuler_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExponential.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperExponentialEulerAppAction.hpp"

#include "Tempus_PhiEvaluator.hpp"

namespace Tempus {

/**
 * \brief One-step exponential Euler method for an implicitly defined ODE.
 *
 * The application model supplies the residual
 * \f$\mathcal{F}_{\mathrm{impl}}(\dot{x},x,t)=0\f$.  With
 * \f$M=\partial\mathcal{F}_{\mathrm{impl}}/\partial\dot{x}\f$, the
 * associated explicit tendency is
 * \f$\mathcal{F}(x,t)=-M^{-1}\mathcal{F}_{\mathrm{impl}}(0,x,t)\f$.
 * At the current state, the stepper assembles the residual Jacobian and mass
 * matrix through the PhiEvaluator, then computes the implemented phi_1 update
 * \f[
 *   x_{n+1}=x_n+\Delta t\,\phi_1(\Delta t W_n)\mathcal{F}(x_n,t_n),
 *   \qquad \phi_1(z)=\frac{\exp(z)-1}{z},
 * \f]
 * where \f$W_n=\partial\mathcal{F}/\partial x = -M^{-1}\partial\mathcal{F}_{\mathrm{impl}}(0,x,t)/\partial x\f$.
 * The residual supplied to
 * the PhiEvaluator is mass weighted; the evaluator performs the required mass
 * solves internally rather than forming \f$M^{-1}\f$.
 *
 * The method requires a current and working SolutionState.  App actions run at
 * BEGIN_STEP, before linearization and the phi evaluation (BEFORE_EXP), after
 * the state update (AFTER_EXP), and after status and norms are recorded
 * (END_STEP).  `Use FSAL` is false by default.  When enabled and xDot storage
 * is available, the stepper evaluates and stores the new-state derivative as
 * an optional cache for the next step; it is not a First-Same-As-Last method,
 * it only behaves as one by reordering operations.
 *
 * The reported order is one.  The reported order range is one through two;
 * order two applies only to autonomous problems when the Jacobian is exact.
 *
 * @tparam Scalar Scalar type used by the model, state vectors, and time.
 */
template<class Scalar>
class StepperExponentialEuler :
    virtual public Tempus::StepperExponential<Scalar>
{
public:

  /**
   * \brief Construct an uninitialized stepper with the default PhiEvaluator.
   *
   * Sets `Use FSAL` to false, initial-condition consistency to `Consistent`,
   * disables the consistency check, and installs the default no-op app action.
   * Set an implicit ModelEvaluator and call initialize() before takeStep().
   */
  StepperExponentialEuler();

  /**
   * \brief Construct and initialize a configured exponential Euler stepper.
   *
   * The PhiEvaluator supplies mass operations and the phi_1 vector product for
   * the implicit application model.  `useFSAL` optionally caches the derivative
   * at a completed working state when that state stores xDot.  `ICConsistency`
   * selects the inherited `None`, `Zero`, `App`, or `Consistent` initial-state
   * policy; `Consistent` computes xDot from the mass-weighted residual and
   * requires xDot storage.  `ICConsistencyCheck` enables its optional residual
   * check.
   *
   * @param appModel `Teuchos::RCP` to the const implicit Thyra model defining
   *   the residual.
   * @param phiEvaluator `Teuchos::RCP` to the PhiEvaluator used for mass and
   *   phi-function operations.
   * @param useFSAL `bool` enabling the optional completed-state derivative cache.
   * @param ICConsistency `std::string` inherited initial-condition policy.
   * @param ICConsistencyCheck `bool` enabling the inherited consistency check.
   * @param stepperEEAppAction `Teuchos::RCP` callback object; null selects the
   *   default no-op action.
   */
  StepperExponentialEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Tempus::PhiEvaluator<Scalar>>& phiEvaluator,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    const Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> >& stepperEEAppAction);

  /// \name Basic stepper methods
  //@{
    /**
     * \brief Set callbacks executed at the four exponential Euler action locations.
     *
     * A null `Teuchos::RCP` installs the default no-op action.  Changing the
     * action invalidates the stepper setup.
     *
     * @param appAction `Teuchos::RCP` to the action callback, or null for the
     *   default action.
     */
    virtual void setAppAction(
      Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > appAction);

    /**
     * \brief Return the action callback used during a step.
     */
    virtual Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > getAppAction() const
    { return stepperEEAppAction_; }

    /**
     * \brief Advance the current SolutionState into the working SolutionState.
     *
     * The SolutionHistory must contain at least two states.  The working state
     * provides \f$\Delta t\f$ and is updated with the phi_1 exponential Euler
     * formula.  BEGIN_STEP runs first; BEFORE_EXP runs before residual/Jacobian
     * assembly and the phi evaluation; AFTER_EXP runs after the update; and
     * END_STEP runs after the working-state status, order, and norms are set.
     * The working SolutionStatus is the SolveStatus returned by computePhi,
     * including a nonconverged status when reported by the PhiEvaluator.
     *
     * When `Use FSAL` is enabled and xDot storage is present, the method
     * reevaluates the completed state, converts the mass-weighted residual to
     * xDot, and marks that state synchronized.  Otherwise the working state is
     * marked unsynchronized.
     *
     * @pre The stepper is initialized and `solutionHistory` has current and
     *   working SolutionStates.
     * @param solutionHistory `Teuchos::RCP` holding the current state and the
     *   working state to update.
     */
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

    /**
     * \brief Return the currently reported temporal order.
     *
     * @return `Scalar` value 1.0.
     */
    virtual Scalar getOrder() const override { return 1.0; }
    /**
     * \brief Return the minimum reported temporal order.
     *
     * @return `Scalar` value 1.0.
     */
    virtual Scalar getOrderMin() const override { return 1.0; }
    /**
     * \brief Return the maximum reported temporal order.
     *
     * Autonomous problems with an exact Jacobian can attain order two.
     *
     * @return `Scalar` value 2.0.
     */
    virtual Scalar getOrderMax() const override { return 2.0;}

    /** @brief Return true because this method advances from one current state. */
    virtual bool isOneStepMethod() const override {return true;}
    /** @brief Return false because this method does not consume older solution states. */
    virtual bool isMultiStepMethod() const override {return !isOneStepMethod();}
  //@}

  /**
   * \brief Return valid settings for this exponential Euler stepper.
   *
   * The list contains inherited stepper settings, the `PhiEvaluator` sublist,
   * and `Adapt PhiEvaluator Interval`.
   * A positive adaptation interval requests evaluator
   * adaptation on the first step and then at matching step indices;
   * a nonpositive values disable it.
   *
   * @return Const `Teuchos::RCP` to a ParameterList populated with current
   *   defaults.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /// \name Overridden from Teuchos::Describable
  //@{
    /**
     * \brief Write the base-stepper and exponential Euler configuration.
     *
     * @param out `Teuchos::FancyOStream` receiving the description.
     * @param verbLevel Requested Teuchos verbosity level.
     */
    virtual void describe(Teuchos::FancyOStream& out,
                          const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  /**
   * \brief Check that base, exponential, and app-action setup is complete.
   *
   * In addition to the inherited model and PhiEvaluator checks, requires a
   * nonnull exponential Euler app action.  Diagnostics are written to `out`
   * for each missing dependency.
   *
   * @param out `Teuchos::FancyOStream` receiving setup diagnostics.
   * @return `bool` true only when all required setup is valid.
   */
  virtual bool isValidSetup(Teuchos::FancyOStream & out) const override;

private:

  Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > stepperEEAppAction_;

};


/**
 * \brief Create an exponential Euler stepper from a model and parameter list.
 *
 * The factory constructs the default stepper, applies the inherited
 * exponential-stepper settings, and uses the optional `PhiEvaluator` sublist
 * to select and configure an evaluator.  A nonnull model is set and the
 * stepper is initialized; a null model returns a configured but uninitialized
 * stepper.
 *
 * @tparam Scalar Scalar type used by the model and stepper.
 * @param model `Teuchos::RCP` to the const implicit Thyra application model,
 *   or null to defer model setup.
 * @param pl Nonconst `Teuchos::RCP` ParameterList to validate and apply; null
 *   retains constructor defaults.
 * @return `Teuchos::RCP` to the created StepperExponentialEuler.
 */
template<class Scalar>
Teuchos::RCP<StepperExponentialEuler<Scalar> >
createStepperExponentialEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  Teuchos::RCP<Teuchos::ParameterList> pl);


} // namespace Tempus

#endif // Tempus_StepperExponentialEuler_decl_hpp
