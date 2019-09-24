// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorPseudoTransientAdjointSensitivity_decl_hpp
#define Tempus_IntegratorPseudoTransientAdjointSensitivity_decl_hpp

// Tempus
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_AdjointSensitivityModelEvaluator.hpp"

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
 */
template<class Scalar>
class IntegratorPseudoTransientAdjointSensitivity :
    virtual public Tempus::Integrator<Scalar>
{
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
   */
  IntegratorPseudoTransientAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList>                pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

  /** \brief Constructor with model and "Stepper Type" and is fully initialized with default settings. */
  IntegratorPseudoTransientAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType);

  /// Destructor
  /** \brief Constructor that requires a subsequent setParameterList, setStepper, and initialize calls. */
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
  /// Get the Stepper
  virtual Teuchos::RCP<Stepper<Scalar> > getStepper() const override;
  /// Return a copy of the Tempus ParameterList
  virtual Teuchos::RCP<Teuchos::ParameterList> getTempusParameterList() override;
  virtual void setTempusParameterList(Teuchos::RCP<Teuchos::ParameterList> pl) override;
  /// Get the SolutionHistory
  virtual Teuchos::RCP<const SolutionHistory<Scalar> > getSolutionHistory() const override;
   /// Get the TimeStepControl
  virtual Teuchos::RCP<const TimeStepControl<Scalar> > getTimeStepControl() const override;
  /// Returns the IntegratorTimer_ for this Integrator
  virtual Teuchos::RCP<Teuchos::Time> getIntegratorTimer() const override
  {return state_integrator_->getIntegratorTimer();}
  virtual Teuchos::RCP<Teuchos::Time> getStepperTimer() const override
  {return state_integrator_->getStepperTimer();}

  //@}

  /// Set the initial state from Thyra::VectorBase(s)
  virtual void initializeSolutionHistory(
    Scalar t0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0 = Teuchos::null,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0 = Teuchos::null,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > y0 = Teuchos::null,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > ydot0 = Teuchos::null,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > ydotdot0 = Teuchos::null);

  /// Get current the solution, x
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getX() const;
  /// Get current the time derivative of the solution, xdot
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getXdot() const;
  /// Get current the second time derivative of the solution, xdotdot
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getXdotdot() const;

  /// Return adjoint sensitivity stored in gradient format
  virtual Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > getDgDp() const;

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl)
      override;
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList() override;
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
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  // Create sensitivity model evaluator from application model
  Teuchos::RCP<AdjointSensitivityModelEvaluator<Scalar> >
  createSensitivityModel(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<Teuchos::ParameterList>& inputPL);

  void buildSolutionHistory();

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model_;
  Teuchos::RCP<AdjointSensitivityModelEvaluator<Scalar> > sens_model_;
  Teuchos::RCP<IntegratorBasic<Scalar> > state_integrator_;
  Teuchos::RCP<IntegratorBasic<Scalar> > sens_integrator_;
  Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory_;
  Teuchos::RCP<DMVPV> dgdp_;
};

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  std::string stepperType);

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity();

} // namespace Tempus

#endif // Tempus_IntegratorPseudoTransientAdjointSensitivity_decl_hpp
