//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorForwardSensitivity_decl_hpp
#define Tempus_IntegratorForwardSensitivity_decl_hpp

// Tempus
#include "Tempus_config.hpp"
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
 * sensitivity at the final time.  One should use the getX() and getDxDp()
 * methods for extracting the final solution and its parameter sensitivity
 * as a multi-vector.  This data can also be extracted from the solution
 * history, but is stored as a Thyra product vector which requires knowledge
 * of the internal implementation.
 */
template <class Scalar>
class IntegratorForwardSensitivity : virtual public Tempus::Integrator<Scalar> {
 public:
  /** \brief Full Constructor with model, and will be fully initialized.
   *
   * \param[in] model               The forward physics ModelEvaluator
   * \param[in] integrator          Forward state Integrator
   * \param[in] sens_model          The sensitivity ModelEvaluator
   * \param[in] sens_stepper        Tempus stepper for the sensitivity
   * integration \param[in] use_combined_method Indicates whether or not to use
   * the "Combined" sensitivity method
   *
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
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model,
      const Teuchos::RCP<IntegratorBasic<Scalar>> &integrator,
      const Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar>> &sens_model,
      const Teuchos::RCP<StepperStaggeredForwardSensitivity<Scalar>>
          &sens_stepper,
      const bool use_combined_method);

  /// Destructor
  /** \brief Constructor that requires a subsequent setStepper, and initialize
   * calls. */
  IntegratorForwardSensitivity();

  /// Destructor
  virtual ~IntegratorForwardSensitivity() {}

  /// \name Basic integrator methods
  //@{
  /// Advance the solution to timeMax, and return true if successful.
  virtual bool advanceTime() { return integrator_->advanceTime(); }
  /// Advance the solution to timeFinal, and return true if successful.
  virtual bool advanceTime(const Scalar timeFinal) override
  {
    return integrator_->advanceTime(timeFinal);
  }
  /// Perform tasks before start of integrator.
  virtual void startIntegrator() { integrator_->startIntegrator(); }
  /// Start time step.
  virtual void startTimeStep() { integrator_->startTimeStep(); }
  /// Check if time step has passed or failed.
  virtual void checkTimeStep() { integrator_->checkTimeStep(); }
  /// Perform tasks after end of integrator.
  virtual void endIntegrator() { integrator_->endIntegrator(); }
  //@}

  /// \name Accessor methods
  //@{
  /// Get current time
  virtual Scalar getTime() const override { return integrator_->getTime(); }
  /// Get current index
  virtual int getIndex() const override { return integrator_->getIndex(); }
  /// Get Status
  virtual Status getStatus() const override { return integrator_->getStatus(); }
  // Set Status
  virtual void setStatus(const Status st) override
  {
    integrator_->setStatus(st);
  }
  /// Get the Stepper
  virtual Teuchos::RCP<Stepper<Scalar>> getStepper() const override
  {
    return integrator_->getStepper();
  }

  /// Set the Stepper
  virtual void setStepper(Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> model);

  /// Set the Stepper
  virtual void setStepper(Teuchos::RCP<Stepper<Scalar>> stepper)
  {
    integrator_->setStepper(stepper);
  }
  /// Set the initial state which has the initial conditions
  virtual void initializeSolutionHistory(
      Teuchos::RCP<SolutionState<Scalar>> state = Teuchos::null)
  {
    integrator_->initializeSolutionHistory(state);
  }

  /// Set the initial state from Thyra::VectorBase(s)
  virtual void initializeSolutionHistory(
      Scalar t0, Teuchos::RCP<const Thyra::VectorBase<Scalar>> x0,
      Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdot0      = Teuchos::null,
      Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdotdot0   = Teuchos::null,
      Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> DxDp0 = Teuchos::null,
      Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> DxdotDp0 =
          Teuchos::null,
      Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> DxdotdotDp0 =
          Teuchos::null);

  /// Get the SolutionHistory
  virtual Teuchos::RCP<const SolutionHistory<Scalar>> getSolutionHistory()
      const override
  {
    return integrator_->getSolutionHistory();
  }
  /// Get the SolutionHistory
  virtual Teuchos::RCP<SolutionHistory<Scalar>> getNonConstSolutionHistory()
      override
  {
    return integrator_->getNonConstSolutionHistory();
  }
  /// Set the SolutionHistory
  virtual void setSolutionHistory(
      Teuchos::RCP<SolutionHistory<Scalar>> sh = Teuchos::null)
  {
    integrator_->setSolutionHistory(sh);
  }
  /// Get the TimeStepControl
  virtual Teuchos::RCP<const TimeStepControl<Scalar>> getTimeStepControl()
      const override
  {
    return integrator_->getTimeStepControl();
  }
  virtual Teuchos::RCP<TimeStepControl<Scalar>> getNonConstTimeStepControl()
      override
  {
    return integrator_->getNonConstTimeStepControl();
  }
  /// Set the TimeStepControl
  virtual void setTimeStepControl(
      Teuchos::RCP<TimeStepControl<Scalar>> tsc = Teuchos::null)
  {
    integrator_->setTimeStepControl(tsc);
  }
  /// Get the Observer
  virtual Teuchos::RCP<IntegratorObserver<Scalar>> getObserver()
  {
    return integrator_->getObserver();
  }
  /// Set the Observer
  virtual void setObserver(
      Teuchos::RCP<IntegratorObserver<Scalar>> obs = Teuchos::null)
  {
    integrator_->setObserver(obs);
  }
  /// Initializes the Integrator after set* function calls
  virtual void initialize() { integrator_->initialize(); }
  virtual Teuchos::RCP<Teuchos::Time> getIntegratorTimer() const override
  {
    return integrator_->getIntegratorTimer();
  }
  virtual Teuchos::RCP<Teuchos::Time> getStepperTimer() const override
  {
    return integrator_->getStepperTimer();
  }

  /**
   * @brief Get the current solution, x, only. If looking for the solution
   * vector and the sensitivities, use `SolutionState->getX()` which will return
   * a Block MultiVector with the first block containing the current solution,
   * x, and the remaining blocks are the forward sensitivities \f$dx/dp\f$.
   *
   * Use `getDxDp` to get the forward sensitivities \f$dx/dp\f$ only.
   *
   * @return The current solution, x, without the sensitivities.
   * */
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getX() const;

  /// Get the forward sensitivities \f$dx/dp\f$
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDxDp() const;
  /**
   * @brief Get current the time derivative of the solution, xdot, only.  This
   * is the first block only and not the full Block MultiVector.
   *
   * @return Get current the time derivative of the solution, xdot.
   * */
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getXDot() const;
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDXDotDp() const;
  /**
   * @brief Get current the second time derivative of the solution, xdotdot,
   * only.  This is the first block only and not the full Block MultiVector.
   *
   * Use `getDXDotDp` to get the forward sensitivities.
   *
   * @return Get current the second time derivative of the solution, xdotdot.
   * */
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getXDotDot() const;
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDXDotDotDp()
      const;

  /// Return response function g
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getG() const;
  /// Return forward sensitivity stored in Jacobian format
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDgDp() const;

  /// Get current state
  virtual Teuchos::RCP<SolutionState<Scalar>> getCurrentState()
  {
    return integrator_->getCurrentState();
  }
  //@}

  /// Parse when screen output should be executed
  void parseScreenOutput() { integrator_->parseScreenOutput(); }

  /// \name Overridden from Teuchos::Describable
  //@{
  std::string description() const override;
  void describe(Teuchos::FancyOStream &out,
                const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  //! What mode the current time integration step is in
  SensitivityStepMode getStepMode() const;

 protected:
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> model_;
  Teuchos::RCP<IntegratorBasic<Scalar>> integrator_;
  Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar>> sens_model_;
  Teuchos::RCP<StepperStaggeredForwardSensitivity<Scalar>> sens_stepper_;
  bool use_combined_method_;
};

/// Nonmember constructor
/**
 * @brief Nonmember constructor
 *
 * This nonmember constructor calls parses the `pList` provided to constructor
 * the sub-objects needed to call the full `IntegratorForwardSensitivity`
 * construtor
 *
 * @param pList  ParameterList defining the integrator options and options
 *               defining the sensitivity analysis
 * @param model  ModelEvaluator for the problem
 * @param sens_residual_model  Sensitivity residual model
 * @param sens_solve_model     Sensitivity solve model
 *
 * @return Time integrator implementing forward sensitivity
 */
template <class Scalar>
Teuchos::RCP<IntegratorForwardSensitivity<Scalar>>
createIntegratorForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &sens_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &sens_solve_model);

/// Nonmember constructor
/**
 * @brief Nonmember constructor
 *
 * This nonmember constructor calls parses the `pList` provided to constructor
 * the sub-objects needed to call the full `IntegratorForwardSensitivity`
 * construtor
 *
 * @param pList  ParameterList defining the integrator options and options
 *               defining the sensitivity analysis
 * @param model  ModelEvaluator for the problem
 * @param sens_residual_model Model evaluator for sensitivity residual
 *
 * @return Time integrator implementing forward sensitivity
 */
template <class Scalar>
Teuchos::RCP<IntegratorForwardSensitivity<Scalar>>
createIntegratorForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &sens_residual_model)
{
  return createIntegratorForwardSensitivity(pList, model, sens_residual_model,
                                            sens_residual_model);
}

/// Nonmember constructor
/**
 * @brief Nonmember constructor
 *
 * This nonmember constructor calls parses the `pList` provided to constructor
 * the sub-objects needed to call the full `IntegratorForwardSensitivity`
 * construtor
 *
 * @param pList  ParameterList defining the integrator options and options
 *               defining the sensitivity analysis
 * @param model  ModelEvaluator for the problem
 *
 * @return Time integrator implementing forward sensitivity
 */
template <class Scalar>
Teuchos::RCP<IntegratorForwardSensitivity<Scalar>>
createIntegratorForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model)
{
  return createIntegratorForwardSensitivity(pList, model, model, model);
}

/// Nonmember constructor
/**
 * @brief Default non-member constructor
 *
 * This nonmember constructor creates default state and sensitivity time
 * integrator.
 *
 * @return Time integrator implementing forward sensitivity
 */
template <class Scalar>
Teuchos::RCP<IntegratorForwardSensitivity<Scalar>>
createIntegratorForwardSensitivity();

}  // namespace Tempus

#endif  // Tempus_IntegratorForwardSensitivity_decl_hpp
