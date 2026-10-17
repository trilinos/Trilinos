//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorAdjointSensitivity_decl_hpp
#define Tempus_IntegratorAdjointSensitivity_decl_hpp

#include <vector>

// Tempus
#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_AdjointAuxSensitivityModelEvaluator.hpp"
#include "Tempus_AdjointSensitivityUtils.hpp"
#include "Tempus_AdjointDAEJumpSolver.hpp"
#include "Tempus_DAEStructure.hpp"

#include "Tempus_StepperStaggeredForwardSensitivity.hpp"  // For SensitivityStepMode

namespace Tempus {

/** \brief Time integrator suitable for adjoint sensitivity analysis */
/**
 * This integrator implements transient adjoint sensitivities.  Given a model
 * evaluator encapsulating the equations f(x_dot,x,p) = 0, and a response
 * function g(x,p), these equations are integrated forward in time to some
 * final time T.  Then the adjoint equations
 *        F(y) = d/dt( df/dx_dot^T*y ) - df/dx^T*y = 0
 * are integrated backward in time starting at t = T where y is the adjoint
 * variable.  Nominally y is a multi-vector belonging to the vector space of f
 * and the number of columns given by the number of entries in the response g.
 * The initial conditions for y at t = T are given by
 *        y(T) = df/dx_dot(x_dot(T),x(T),p)^{-T} * dg/dx(x(T),p)^T.
 * Then the final sensitivity of g is
 *        dg/dp(T) - int_{0}^T(df/dp^T*y)dt + dx/dp(0)^T*df/dx_dot(0)^T*y(0).
 *
 * This integrator supports both implicit and explicit steppers, and provides
 * a method to compute dg/dp as described above after the reverse integration.
 * The solution history contains solutions stored as product vectors (x,y),
 * with y further stored as a product vector for each component of g.
 *
 * Because of limitations on the steppers, the implementation currently assumes
 * df/dxdot is a constant matrix.
 *
 * Besides the final-time response above, two more responses are supported
 * through the "Response Type" option of the "Sensitivities" sublist:
 * <ul>
 *   <li> "Time Integral": J = int_{0}^T g(x(t),p,t) dt.  The adjoint equations
 *        get the source term dg/dx^T, y(T) = 0, and the sensitivity is
 *        int_{0}^T(dg/dp^T - df/dp^T*y)dt + dx/dp(0)^T*df/dx_dot(0)^T*y(0).
 *   <li> "Sum Over Times": J = sum_k g(x(t_k),p,t_k), with the times t_k given
 *        by "Response Time List".  Between the times t_k the adjoint equations
 *        are homogeneous, and at each t_k the adjoint jumps by
 *        df/dx_dot^{-T}*dg/dx(t_k)^T.  The sensitivity is
 *        sum_k dg/dp(t_k)^T - int_{0}^T(df/dp^T*y)dt
 *          + dx/dp(0)^T*df/dx_dot(0)^T*y(0).
 *        The adjoint integration is restarted at every t_k, so no time step
 *        crosses a jump.  The forward time steps must end at the times t_k:
 *        with "Step Type" = "Variable" the times t_k are added as time events
 *        on which the time steps land exactly, with "Step Type" = "Constant"
 *        they must be on the time grid.
 * </ul>
 * In both cases the response must not depend on x_dot.  If the model supports
 * dg/dx_dot, it is checked to be zero.
 *
 * <b>Index-1 DAEs.</b>  If some equations do not depend on x_dot (algebraic
 * equations, e.g. constraints, or Dirichlet conditions written as
 * f_i = x_i - d_i(p,t)), df/dx_dot is singular and the jumps above are not
 * defined.  For an index-1 DAE, the "Algebraic Equations Mask" option gives
 * the algebraic equations, the rows of df/dx_dot that are zero (see
 * DAEStructure).  With
 * E = df/dx_dot, A = df/dx, and P_a, P_d the masks of the algebraic and
 * differential equations, the adjoint jump at a response time t_k is then
 * given by the solution z of
 *        C^T z = dg/dx(t_k)^T,   C = P_a A + P_d E   (at the state at t_k):
 * <ul>
 *   <li> the jump of the adjoint is P_d z (its differential components);
 *   <li> the algebraic adjoint unknowns are determined by the adjoint
 *        equations, which are a DAE as well, and the algebraic components of
 *        z add -df/dp(t_k)^T (P_a z) to the sensitivity.
 * </ul>
 * This accounts for the dependence of the response on the algebraic unknowns
 * x_a: for a semi-explicit DAE, ordering the algebraic equations and unknowns
 * first, C = [A_aa A_ad; 0 E_dd], so P_a z = A_aa^{-T} dg/dx_a^T and
 * P_d z = E_dd^{-T} (dg/dx_d^T - A_ad^T P_a z).
 * C^T z = dg/dx^T is solved with AdjointDAEJumpSolver (options in the
 * "Algebraic Jump Solver" sublist of "Sensitivities"), using only evaluations
 * of the adjoint model: the model must return the true df/dx_dot and df/dx
 * for any alpha and beta.  The "Time Integral" response has no jumps and
 * does not need this option.  The steppers must be suited to index-1 DAEs
 * (e.g. Backward Euler), and must not compute consistent initial conditions
 * ("Initial Condition Consistency" = "Consistent"), which requires a
 * nonsingular df/dx_dot.
 *
 * To extract the final solution x(T) and sensitivity dg/dp(T) one should use
 * the getX() and getDgDp() methods, which return these quantities directly.
 * One can also extract this data for all times from the solution history,
 * however the data is stored in Thyra product vectors which requires
 * knowledge of the internal implementation.
 */
template <class Scalar>
class IntegratorAdjointSensitivity : virtual public Tempus::Integrator<Scalar> {
 public:
  /** \brief Full Constructor with model, and will be fully initialized.
   *
   * \param[in] model                 The forward physics ModelEvaluator
   * \param[in] state_integrator      Forward state Integrator for the forward
   * problem \param[in] adjoint_model         ModelEvaluator for the adjoint
   * physics/problem \param[in] adjoint_aux_model     ModelEvaluator for the
   * auxiliary adjoint physics/problem \param[in] adjoint_integrator    Time
   * integrator for the adjoint problem \param[in] solution_history      The
   * forward state solution history \param[in] p_index               Sensitivity
   * parameter index \param[in] g_index               Response function index
   * \param[in] g_depends_on_p        Does response depends on parameters?
   * \param[in] f_depends_on_p        Does residual depends on parameters?
   * \param[in] ic_depends_on_p       Does the initial condition depends on
   * parameters?
   * @param adjoint_model ModelEvaluator for the adjoint problem. Optional.
   * Default value is null. \param[in] mass_matrix_is_identity Is the mass
   * matrix an identity matrix?
   *
   * In addition to all of the regular integrator options, the supplied
   * parameter list supports the following options contained within a sublist
   * "Sensitivities" from the top-level parameter list:
   * <ul>
   *    <li> "Sensitivity Parameter Index", (default: 0) The model evaluator
   *          parameter index for which sensitivities will be computed.
   *    <li> "Response Function Index", (default: 0) The model evaluator
   *         response index for which sensitivities will be computed.
   *    <li> "Response Depends on Parameters", (default: true) Whether the
   *         response function depends on the parameter vector p for which
   *         sensitivities will be computed.  If set to false, the dg/dp
   *         term in the sensitivity formula will not be computed.
   *    <li> "Residual Depends on Parameters", (default: true) Whether the
   *         model residual f depends on the parameter vector p for which
   *         sensitivities will be computed.  If set to false, its contribution
   *         to the sensitivities will not be computed.
   *    <li> "IC Depends on Parameters", (default: true) Whether the initial
   *         conditions depend on the parameter vector p for which
   *         sensitivities will be computed.  If set to false, its contribution
   *         to the sensitivities will not be computed.
   *    <li> "Mass Matrix Is Constant" (default: true) Whether the mass matrix
   *         df/dx_dot is a constant matrix.  As describe above, this is
   *         currently required to be true.
   *    <li> "Mass Matrix Is Identity" (default: false) Whether the mass matrix
   *         is the identity matrix, in which some computations can be skipped.
   *    <li> "Response Type" (default: "Final Time") "Final Time" for
   *         g(x(T),p), "Time Integral" for int_{0}^T g(x(t),p,t) dt, or
   *         "Sum Over Times" for sum_k g(x(t_k),p,t_k).
   *    <li> "Response Time List" (default: "") Comma-separated list of the
   *         times t_k of the "Sum Over Times" response.
   *    <li> "Algebraic Equations Mask" (default: null) A
   *         Teuchos::RCP<const Thyra::VectorBase<Scalar> > in the residual
   *         space, equal to 1 on the algebraic equations of an index-1 DAE
   *         and 0 elsewhere (see the class description and DAEStructure).
   *    <li> "Algebraic Jump Solver" (sublist) Options of the solver of the
   *         adjoint jumps of an index-1 DAE (see AdjointDAEJumpSolver).
   * </ul>
   *
   * To support use-cases with explicitly computed adjoint operators, the
   * constructor takes an additional model evaluator for computing the adjoint
   * W/W_op.  It is assumed the operator returned by this model evaluator is
   * the adjoint, and so will not be transposed.  It is also assumed this
   * model evaluator accepts the same inArgs as the forward model, however it
   * only requires supporting the adjoint W/W_op outArgs.
   */

  /*! \brief Version of the constructor taking a single model evaluator. */
  /*!
   * This version takes a single model evaluator for the case when the adjoint
   * is implicitly determined from the forward operator by the (conjugate)
   * transpose.
   */
  IntegratorAdjointSensitivity(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model,
      const Teuchos::RCP<IntegratorBasic<Scalar>> &state_integrator,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &adjoint_model,
      const Teuchos::RCP<AdjointAuxSensitivityModelEvaluator<Scalar>>
          &adjoint_aux_model,
      const Teuchos::RCP<IntegratorBasic<Scalar>> &adjoint_integrator,
      const Teuchos::RCP<SolutionHistory<Scalar>> &solution_history,
      const int p_index, const int g_index, const bool g_depends_on_p,
      const bool f_depends_on_p, const bool ic_depends_on_p,
      const bool mass_matrix_is_identity,
      const AdjointResponseType response_type   = AdjointResponseType::FinalTime,
      const std::vector<Scalar> &response_times = std::vector<Scalar>(),
      const Teuchos::RCP<const DAEStructure<Scalar>> &dae_structure =
          Teuchos::null,
      const Teuchos::RCP<const Teuchos::ParameterList> &dae_jump_solver_pl =
          Teuchos::null);

  /** \brief Constructor that requires a subsequent setParameterList,
   * setStepper, and initialize calls. */
  IntegratorAdjointSensitivity();

  /// Destructor
  virtual ~IntegratorAdjointSensitivity() {}

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
  /// Returns the IntegratorTimer_ for this Integrator
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

  /// Get the Observer
  virtual Teuchos::RCP<IntegratorObserver<Scalar>> getObserver();
  /// Set the Observer
  virtual void setObserver(
      Teuchos::RCP<IntegratorObserver<Scalar>> obs = Teuchos::null);
  /// Initializes the Integrator after set* function calls
  virtual void initialize();

  /// Get the current solution, x
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getX() const;
  /// Get the current time derivative of the solution, xdot
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getXDot() const;
  /// Get the current second time derivative of the solution, xdotdot
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar>> getXDotDot() const;

  /// Get the current adjoint solution, y
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getY() const;
  /// Get the current time derivative of the adjoint solution, ydot
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getYDot() const;
  /// Get the current second time derivative of the adjoint solution, ydotdot
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getYDotDot() const;

  /// Return adjoint sensitivity stored in gradient format
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> getDgDp() const;

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> getAdjointModel() const
  {
    return adjoint_model_;
  };

  /// Structure of the index-1 DAE (null for an ODE)
  Teuchos::RCP<const DAEStructure<Scalar>> getDAEStructure() const
  {
    return dae_structure_;
  }

  /// \name Overridden from Teuchos::Describable
  //@{
  std::string description() const override;
  void describe(Teuchos::FancyOStream &out,
                const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  //! What mode the current time integration step is in
  SensitivityStepMode getStepMode() const;

 protected:
  // Create sensitivity model evaluator from application model
  Teuchos::RCP<AdjointAuxSensitivityModelEvaluator<Scalar>> createAdjointModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &adjoint_model,
      const Teuchos::RCP<Teuchos::ParameterList> &inputPL);

  void buildSolutionHistory(
      const Teuchos::RCP<const SolutionHistory<Scalar>> &state_solution_history,
      const Teuchos::RCP<const SolutionHistory<Scalar>>
          &adjoint_solution_history);

  //! Solver for df/dx_dot^T, computed with the adjoint model at inArgs
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar>>
  createAdjointMassMatrixSolver(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs) const;

  //! Forward state at time t, which must be the time of a stored time step
  Teuchos::RCP<const SolutionState<Scalar>> findForwardState(
      const SolutionHistory<Scalar> &state_solution_history, const Scalar t,
      const Scalar tol) const;

  //! Integrate the adjoint problem from tau_begin to tau_end, starting from
  //! adjoint_state, which is overwritten with the adjoint state at tau_end.
  //! index is the time step index at tau_begin (negative for the initial
  //! index of the time step control), and it is overwritten with the index
  //! at tau_end.
  bool integrateAdjoint(const Scalar tau_begin, const Scalar tau_end,
                        Thyra::VectorBase<Scalar> &adjoint_state, int &index);

  //! Append the states of the adjoint integrator to adjoint_solution_history_
  //! (used when the adjoint integration is restarted at response times).
  void appendAdjointSolutionHistory(const Scalar tol);

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> model_;
  Teuchos::RCP<IntegratorBasic<Scalar>> state_integrator_;
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> adjoint_model_;
  Teuchos::RCP<AdjointAuxSensitivityModelEvaluator<Scalar>> adjoint_aux_model_;
  Teuchos::RCP<IntegratorBasic<Scalar>> adjoint_integrator_;
  Teuchos::RCP<SolutionHistory<Scalar>> solutionHistory_;
  int p_index_;
  int g_index_;
  bool g_depends_on_p_;
  bool f_depends_on_p_;
  bool ic_depends_on_p_;
  bool mass_matrix_is_identity_;
  AdjointResponseType response_type_;
  std::vector<Scalar> response_times_;
  //! Adjoint solution history combining the segments of the adjoint
  //! integration for the "Sum Over Times" response (null otherwise)
  Teuchos::RCP<SolutionHistory<Scalar>> adjoint_solution_history_;
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> dxdp_init_;
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> dgdp_;
  SensitivityStepMode stepMode_;
  //! Structure of the index-1 DAE (null for an ODE)
  Teuchos::RCP<const DAEStructure<Scalar>> dae_structure_;
  //! Solver of the adjoint jumps of the index-1 DAE (null for an ODE)
  Teuchos::RCP<const AdjointDAEJumpSolver<Scalar>> dae_jump_solver_;
  //! df/dp work objects for the algebraic contributions of the jumps
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> dfdp_op_;
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> dfdp_mv_;
};

/// Nonmember constructor
/**
 * @brief Nonmember constructor
 *
 * This nonmember constructor calls parses the `pList` provided to constructor
 * the sub-objects needed to call the full `IntegratorAdjointSensitivity`
 * construtor
 *
 * @param pList         ParameterList defining the integrator options and
 * options defining the sensitivity analysis
 * @param model         ModelEvaluator for the problem
 * @param adjoint_model ModelEvaluator for the adjoint problem. Optional.
 * Default value is null. When not provided, the adjoint_model will be
 * constructed from the forward physics model.
 *
 * @return Time integrator implementing adjoint sensitivity
 */
template <class Scalar>
Teuchos::RCP<IntegratorAdjointSensitivity<Scalar>>
createIntegratorAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &adjoint_model =
        Teuchos::null);

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorAdjointSensitivity<Scalar>>
createIntegratorAdjointSensitivity();

}  // namespace Tempus

#endif  // Tempus_IntegratorAdjointSensitivity_decl_hpp
