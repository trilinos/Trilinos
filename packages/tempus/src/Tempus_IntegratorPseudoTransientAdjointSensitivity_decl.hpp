//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorPseudoTransientAdjointSensitivity_decl_hpp
#define Tempus_IntegratorPseudoTransientAdjointSensitivity_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_AdjointSensitivityModelEvaluator.hpp"

#include "Tempus_StepperStaggeredForwardSensitivity.hpp"  // For SensitivityStepMode

namespace Tempus {

/** \brief Time integrator suitable for pseudotransient adjoint sensitivity
 * analysis */
/**
 * For some problems, time integrators are used to compute steady-state
 * solutions (also known as pseudo-transient solvers).  When computing
 * sensitivities, it is not necessary in these cases to propagate sensitivities
 * all the way through the forward time integration.  Instead the steady-state
 * is first computed as usual, and then the sensitivities are computed using
 * a similar pseudo-transient time integration applied to the adjoint
 * sensitivity equations with the state frozen to the computed steady-state.
 * This integrator specializes the transient sensitivity methods implemented by
 * Tempus::IntegratorAdjointSensitivity to this case.
 *
 * Consider an implicit ODE f(x_dot,x,p) = 0 with a stable steady-state solution
 * x = x^s, x_dot = 0 where f(0,x^s,p) = 0 and all of the eigenvalues of
 * df/dx(0,x^s,p) are in the right half-plane (for an explicit ODE, the
 * eigenvalues must be in the left half-plane).  In the pseudo-transient method
 * a time-integrator is applied to f(x_dot,x,p) = 0 until x_dot is sufficiently
 * small.  Now consider the adjoint sensitivity equations for some response
 * function g(x,p):
 *       df/dx_dot^T*y_dot + df/dx^T*y - dg/dx^T = 0
 * after the transformation tau = T - t has been applied, where T is the final
 * time.  For pseudo-transient adjoint sensitivities, the above is integrated
 * from y(0) = 0 until y_dot is sufficiently small, in which case
 *       y^s = (df/dx)^{-T}*(dg/dx)^T.
 * Then the final sensitivity of g is
 *       dg/dp^T - df/dp^T*y^s.
 * One can see that y^s is the only steady-state solution of the adjoint
 * equations, since df/dx and dg/dx are constant, and must be linearly stable
 * (since the eigenvalues of df/dx^T are the same as df/dx).
 *
 * To extract the final solution x(T) and sensitivity dg/dp one should use
 * the getX() and getDgDp() methods, which return these quantities directly.
 * One can also extract this data for all times from the solution history,
 * however the data is stored in Thyra product vectors which requires
 * knowledge of the internal implementation.
 */
template <class Scalar>
class IntegratorPseudoTransientAdjointSensitivity
  : virtual public Tempus::Integrator<Scalar>,
    virtual public Teuchos::ParameterListAcceptor {
 public:
  /** \brief Constructor with ParameterList and model, and will be fully
   * initialized. */
  /*!
   * In addition to all of the regular integrator options, the supplied
   * parameter list supports the following options contained within a sublist
   * "Sensitivities" from the top-level parameter list:
   * <ul>
   *    <li> "Sensitivity Parameter Index", (default: 0) The model evaluator
   *          parameter index for which sensitivities will be computed.
   *    <li> "Response Function Index", (default: 0) The model evaluator
   *         response index for which sensitivities will be computed.
   *    <li> "Mass Matrix Is Constant" (default: true) Whether the mass matrix
   *         df/dx_dot is a constant matrix.  As describe above, this is
   *         currently required to be true.
   *    <li> "Mass Matrix Is Identity" (default: false) Whether the mass matrix
   *         is the identity matrix, in which some computations can be skipped.
   * </ul>
   *
   * To support use-cases with explicitly computed adjoint operators, the
   * constructor takes two additional model evaluators for computing the adjoint
   * W/W_op.  It is assumed the operator returned by these model evaluators is
   * the adjoint, and so will not be transposed.  It is also assumed these
   * model evaluators accept the same inArgs as the forward model, however they
   * only requires supporting the adjoint W/W_op outArgs.  The first adjoint
   * model evaluator will only be used for computing the adjoint residual,
   * whereas the second will produce the operator used in linear solves.  The
   * passed RCP's may point to the same model evalutor object in cases where
   * these operators are the same
   */
  IntegratorPseudoTransientAdjointSensitivity(
      Teuchos::RCP<Teuchos::ParameterList> pList,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&
          adjoint_residual_model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_solve_model);

  /** \brief Constructor with model and "Stepper Type" and is fully initialized
   * with default settings. */
  IntegratorPseudoTransientAdjointSensitivity(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&
          adjoint_residual_model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_solve_model,
      std::string stepperType);

  /*! \brief Version of the constructor taking a single adjoint model evaluator.
   */
  IntegratorPseudoTransientAdjointSensitivity(
      Teuchos::RCP<Teuchos::ParameterList> pList,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_model);

  /*! \brief Version of the constructor taking a single adjoint model evaluator.
   */
  IntegratorPseudoTransientAdjointSensitivity(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_model,
      std::string stepperType);

  /*! \brief Version of the constructor taking a single model evaluator. */
  /*!
   * This version takes a single model evaluator for the case when the adjoint
   * is implicitly determined from the forward operator by the (conjugate)
   * transpose.
   */
  IntegratorPseudoTransientAdjointSensitivity(
      Teuchos::RCP<Teuchos::ParameterList> pList,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

  /*! \brief Version of the constructor taking a single model evaluator. */
  /*!
   * This version takes a single model evaluator for the case when the adjoint
   * is implicitly determined from the forward operator by the (conjugate)
   * transpose.
   */
  IntegratorPseudoTransientAdjointSensitivity(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
      std::string stepperType);

  /// Destructor
  /** \brief Constructor that requires a subsequent setParameterList,
   * setStepper, and initialize calls. */
  IntegratorPseudoTransientAdjointSensitivity();

  /// Destructor
  virtual ~IntegratorPseudoTransientAdjointSensitivity() {}

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
  virtual Teuchos::RCP<Stepper<Scalar> > getStepper() const override;
  Teuchos::RCP<Stepper<Scalar> > getStateStepper() const;
  Teuchos::RCP<Stepper<Scalar> > getSensStepper() const;
  /// Get the SolutionHistory
  virtual Teuchos::RCP<const SolutionHistory<Scalar> > getSolutionHistory()
      const override;
  Teuchos::RCP<const SolutionHistory<Scalar> > getStateSolutionHistory() const;
  Teuchos::RCP<const SolutionHistory<Scalar> > getSensSolutionHistory() const;
  /// Get the SolutionHistory
  virtual Teuchos::RCP<SolutionHistory<Scalar> > getNonConstSolutionHistory()
      override;
  /// Get the TimeStepControl
  virtual Teuchos::RCP<const TimeStepControl<Scalar> > getTimeStepControl()
      const override;
  virtual Teuchos::RCP<TimeStepControl<Scalar> > getNonConstTimeStepControl()
      override;
  Teuchos::RCP<TimeStepControl<Scalar> > getStateNonConstTimeStepControl();
  Teuchos::RCP<TimeStepControl<Scalar> > getSensNonConstTimeStepControl();
  /// Get the Observer
  virtual Teuchos::RCP<IntegratorObserver<Scalar> > getObserver();
  /// Set the Observer
  virtual void setObserver(
      Teuchos::RCP<IntegratorObserver<Scalar> > obs = Teuchos::null);
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
      Scalar t0, Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0      = Teuchos::null,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0   = Teuchos::null,
      Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > y0    = Teuchos::null,
      Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > ydot0 = Teuchos::null,
      Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > ydotdot0 =
          Teuchos::null);

  /// Get the current solution, x
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getX() const;
  /// Get the current time derivative of the solution, xdot
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getXDot() const;
  /// Get the current second time derivative of the solution, xdotdot
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getXDotDot() const;

  /// Get the current adjoint solution, y
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > getY() const;
  /// Get the current time derivative of the adjoint solution, ydot
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > getYDot() const;
  /// Get the current second time derivative of the adjoint solution, ydotdot
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > getYDotDot()
      const;

  /// Return response function g
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getG() const;
  /// Return adjoint sensitivity stored in gradient format
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > getDgDp() const;

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
  void setParameterList(
      const Teuchos::RCP<Teuchos::ParameterList>& pl) override;
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList() override;
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList() override;

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters()
      const override;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
  std::string description() const override;
  void describe(Teuchos::FancyOStream& out,
                const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  //! What mode the current time integration step is in
  SensitivityStepMode getStepMode() const;

  //! Set/get whether to do the forward integration
  void setDoForwardIntegration(const bool f) { do_forward_integration_ = f; }
  bool getDoForwardIntegration() const { return do_forward_integration_; }

  //! Set/get whether to do the adjoint integration
  void setDoAdjointIntegration(const bool f) { do_adjoint_integration_ = f; }
  bool getDoAdjointIntegration() const { return do_adjoint_integration_; }

 protected:
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  // Create sensitivity model evaluator from application model
  Teuchos::RCP<AdjointSensitivityModelEvaluator<Scalar> >
  createSensitivityModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&
          adjoint_residual_model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_solve_model,
      const Teuchos::RCP<Teuchos::ParameterList>& inputPL);

  void buildSolutionHistory();

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model_;
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > adjoint_residual_model_;
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > adjoint_solve_model_;
  Teuchos::RCP<AdjointSensitivityModelEvaluator<Scalar> > sens_model_;
  Teuchos::RCP<IntegratorBasic<Scalar> > state_integrator_;
  Teuchos::RCP<IntegratorBasic<Scalar> > sens_integrator_;
  Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > g_;
  Teuchos::RCP<DMVPV> dgdp_;
  SensitivityStepMode stepMode_;
  bool do_forward_integration_;
  bool do_adjoint_integration_;
};

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType);

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_model);

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_model,
    std::string stepperType);

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_solve_model);

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& adjoint_solve_model,
    std::string stepperType);

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity();

}  // namespace Tempus

#endif  // Tempus_IntegratorPseudoTransientAdjointSensitivity_decl_hpp
