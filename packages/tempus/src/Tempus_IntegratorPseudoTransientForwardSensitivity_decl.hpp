//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorPseudoTransientForwardSensitivity_decl_hpp
#define Tempus_IntegratorPseudoTransientForwardSensitivity_decl_hpp

// Tempus
#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_SensitivityModelEvaluatorBase.hpp"

#include "Tempus_StepperStaggeredForwardSensitivity.hpp"  // For SensitivityStepMode

namespace Tempus {

/** \brief Time integrator suitable for pseudotransient forward sensitivity
 * analysis */
/**
 * For some problems, time integrators are used to compute steady-state
 * solutions (also known as pseudo-transient solvers).  When computing
 * sensitivities, it is not necessary in these cases to propagate sensitivities
 * all the way through the forward time integration.  Instead the steady-state
 * is first computed as usual, and then the sensitivities are computed using
 * a similar pseudo-transient time integration applied to the sensitivity
 * equations with the state frozen to the computed steady-state.  This
 * integrator specializes the transient sensitivity methods implemented by
 * Tempus::IntegratorForwardSensitivity to this case.
 *
 * Consider an implicit ODE f(x_dot,x,p) = 0 with a stable steady-state solution
 * x = x^s, x_dot = 0 where f(0,x^s,p) = 0 and all of the eigenvalues of
 * df/dx(0,x^s,p) are in the right half-plane (for an explicit ODE, the
 * eigenvalues must be in the left half-plane).  In the pseudo-transient method
 * a time-integrator is applied to f(x_dot,x,p) = 0 until x_dot is sufficiently
 * small.  Now consider the forward sensitivity equations:
 *       df/dx_dot*z_dot + df/dx*z + df/dp = 0
 * where z = dx/dp.  For pseudo-transient forward sensitivities, the above is
 * integrated from z(0) = 0 until z_dot is sufficiently small, in which case
 *       z^s = -(df/dx)^{-1}*(df/dp).
 * Then the final sensitivity of g is
 *       dg/dp + dg/dx*z^s.
 * One can see that z^s is the only steady-state solution of the sensitivity
 * equations, since df/dx and df/dp are constant, and must be linearly stable
 * since it has the same Jacobian matrix as the forward equations.
 *
 * One should use the getX() and getDxDp()
 * methods for extracting the final sultion and its parameter sensitivity
 * as a multi-vector.  This data can also be extracted from the solution
 * history, but is stored as a Thyra product vector which requires knowledge
 * of the internal implementation.
 */
template <class Scalar>
class IntegratorPseudoTransientForwardSensitivity
  : virtual public Tempus::Integrator<Scalar> {
 public:
  /** \brief Constructor with ParameterList and model, and will be fully
   * initialized. */
  /*!
   * In addition to all of the regular integrator options, the supplied
   * parameter list supports the following options contained within a sublist
   * "Sensitivities" from the top-level parameter list:
   * <ul>
   *   <li> "Reuse State Linear Solver" (default: false) Whether to reuse the
   *        model's W matrix, solver, and preconditioner when solving the
   *        sensitivity equations.  If they can be reused, substantial savings
   *        in compute time are possible.
   *   <li> "Force W Update" (default: false) When reusing the solver as above
   *        whether to force recomputation of W.  This can be necessary when
   *        the solver overwrites it during the solve-phase (e.g., by a
   *        factorization).
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
  IntegratorPseudoTransientForwardSensitivity(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
      const Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar>>& sens_model,
      const Teuchos::RCP<IntegratorBasic<Scalar>>& fwd_integrator,
      const Teuchos::RCP<IntegratorBasic<Scalar>>& sens_integrator,
      const bool reuse_solver, const bool force_W_update);

  /// Destructor
  /** \brief Constructor that requires a subsequent setStepper, and initialize
   * calls. */
  IntegratorPseudoTransientForwardSensitivity();

  /// Destructor
  virtual ~IntegratorPseudoTransientForwardSensitivity() {}

  /// \name Basic integrator methods
  //@{

  /// Advance the solution to timeMax, and return true if successful.
  virtual bool advanceTime();
  /// Advance the solution to timeFinal, and return true if successful.
  virtual bool advanceTime(const Scalar timeFinal) override;
  /// Get current time
  virtual Scalar getTime() const override;
  /// Get current index
  virtual int getIndex() const override;
  /// Get Status
  virtual Status getStatus() const override;
  /// Set Status
  virtual void setStatus(const Status st) override;
  /// Get the Stepper
  virtual Teuchos::RCP<Stepper<Scalar>> getStepper() const override;
  Teuchos::RCP<Stepper<Scalar>> getStateStepper() const;
  Teuchos::RCP<Stepper<Scalar>> getSensStepper() const;
  /// Get the SolutionHistory
  virtual Teuchos::RCP<const SolutionHistory<Scalar>> getSolutionHistory()
      const override;
  Teuchos::RCP<const SolutionHistory<Scalar>> getStateSolutionHistory() const;
  Teuchos::RCP<const SolutionHistory<Scalar>> getSensSolutionHistory() const;
  /// Get the SolutionHistory
  virtual Teuchos::RCP<SolutionHistory<Scalar>> getNonConstSolutionHistory()
      override;
  /// Get the TimeStepControl
  virtual Teuchos::RCP<const TimeStepControl<Scalar>> getTimeStepControl()
      const override;
  virtual Teuchos::RCP<TimeStepControl<Scalar>> getNonConstTimeStepControl()
      override;
  Teuchos::RCP<TimeStepControl<Scalar>> getStateNonConstTimeStepControl();
  Teuchos::RCP<TimeStepControl<Scalar>> getSensNonConstTimeStepControl();
  /// Get the Observer
  virtual Teuchos::RCP<IntegratorObserver<Scalar>> getObserver();
  /// Set the Observer
  virtual void setObserver(
      Teuchos::RCP<IntegratorObserver<Scalar>> obs = Teuchos::null);
  virtual Teuchos::RCP<Teuchos::Time> getIntegratorTimer() const override
  {
    return state_integrator_->getIntegratorTimer();
  }
  virtual Teuchos::RCP<Teuchos::Time> getStepperTimer() const override
  {
    return state_integrator_->getStepperTimer();
  }

  //@}

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

  /// Get current the solution, x
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getX() const;
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDxDp() const;
  /// Get current the time derivative of the solution, xdot
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getXDot() const;
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDXDotDp() const;
  /// Get current the second time derivative of the solution, xdotdot
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getXDotDot() const;
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDXDotDotDp()
      const;

  /// Return response function g
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getG() const;
  /// Return forward sensitivity stored in Jacobian format
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDgDp() const;

  /// \name Overridden from Teuchos::Describable
  //@{
  std::string description() const override;
  void describe(Teuchos::FancyOStream& out,
                const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  //! What mode the current time integration step is in
  SensitivityStepMode getStepMode() const;

 protected:
  void buildSolutionHistory();

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> model_;
  Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar>> sens_model_;
  Teuchos::RCP<IntegratorBasic<Scalar>> state_integrator_;
  Teuchos::RCP<IntegratorBasic<Scalar>> sens_integrator_;
  Teuchos::RCP<SolutionHistory<Scalar>> solutionHistory_;
  bool reuse_solver_;
  bool force_W_update_;
  SensitivityStepMode stepMode_;
};

/// Nonmember constructor
/**
 * @brief Nonmember constructor
 *
 * @param pList ParameterList to construct the Tempus state integrator, the
 *              sensitivity model evaluator, and the sensisitivity integrator
 * @param model Physics model
 * @param sens_residual_model Model evaluator for sensitivity residual
 * @param sens_solve_model Model evaluator for sensitivity solve
 *
 * @return
 */
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>>
createIntegratorPseudoTransientForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& sens_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& sens_solve_model);

/// Nonmember constructor
/**
 * @brief Nonmember constructor
 *
 * @param pList ParameterList to construct the Tempus state integrator, the
 *              sensitivity model evaluator, and the sensisitivity integrator
 * @param model Physics model
 * @param sens_residual_model Model evaluator for sensitivity residual
 *
 * @return
 */
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>>
createIntegratorPseudoTransientForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& sens_residual_model)
{
  return createIntegratorPseudoTransientForwardSensitivity(
      pList, model, sens_residual_model, sens_residual_model);
}

/// Nonmember constructor
/**
 * @brief Nonmember constructor
 *
 * @param pList ParameterList to construct the Tempus state integrator, the
 *              sensitivity model evaluator, and the sensisitivity integrator
 * @param model Physics model
 *
 * @return
 */
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>>
createIntegratorPseudoTransientForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& model)
{
  return createIntegratorPseudoTransientForwardSensitivity(pList, model, model,
                                                           model);
}

/// Nonmember constructor
/**
 * @brief Default ctor
 *
 * Instantiates a default IntegratorBasic for both the state and the sensitivity
 * integrator.
 *
 * @return IntegratorPseudoTransientForwardSensitivity
 */
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<Scalar>>
createIntegratorPseudoTransientForwardSensitivity();

}  // namespace Tempus

#endif  // Tempus_IntegratorPseudoTransientForwardSensitivity_decl_hpp
