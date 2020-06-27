// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorAdjointSensitivity_decl_hpp
#define Tempus_IntegratorAdjointSensitivity_decl_hpp

// Tempus
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_AdjointAuxSensitivityModelEvaluator.hpp"

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
 */
template<class Scalar>
class IntegratorAdjointSensitivity :
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
   * </ul>
   */
  IntegratorAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList>                pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

  /// Destructor
  /** \brief Constructor that requires a subsequent setParameterList, setStepper, and initialize calls. */
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
  /// Get the Stepper
  virtual Teuchos::RCP<Stepper<Scalar> > getStepper() const override;
  /// Return a copy of the Tempus ParameterList
  virtual Teuchos::RCP<Teuchos::ParameterList> getTempusParameterList() override;
  virtual void setTempusParameterList(Teuchos::RCP<Teuchos::ParameterList> pl) override;
  /// Get the SolutionHistory
  virtual Teuchos::RCP<const SolutionHistory<Scalar> > getSolutionHistory() const override;
   /// Get the TimeStepControl
  virtual Teuchos::RCP<const TimeStepControl<Scalar> > getTimeStepControl() const override;
  virtual Teuchos::RCP<TimeStepControl<Scalar> > getNonConstTimeStepControl() override;
  /// Returns the IntegratorTimer_ for this Integrator
  virtual Teuchos::RCP<Teuchos::Time> getIntegratorTimer() const override
  { return state_integrator_->getIntegratorTimer();}
  virtual Teuchos::RCP<Teuchos::Time> getStepperTimer() const override
  { return state_integrator_->getStepperTimer();}

  //@}

  /// Set the initial state from Thyra::VectorBase(s)
  virtual void initializeSolutionHistory(
    Scalar t0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0 = Teuchos::null,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0 = Teuchos::null,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxDp0 = Teuchos::null,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxdotDp0 = Teuchos::null,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxdotdotDp0 = Teuchos::null);

  /// Set the Observer
  virtual void setObserver(
    Teuchos::RCP<IntegratorObserver<Scalar> > obs = Teuchos::null);
  /// Initializes the Integrator after set* function calls
  virtual void initialize();

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

  // Create sensitivity model evaluator from application model
  Teuchos::RCP<AdjointAuxSensitivityModelEvaluator<Scalar> >
  createAdjointModel(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<Teuchos::ParameterList>& inputPL);

  void buildSolutionHistory(
    const Teuchos::RCP<const SolutionHistory<Scalar> >& state_solution_history,
    const Teuchos::RCP<const SolutionHistory<Scalar> >& adjoint_solution_history);

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model_;
  Teuchos::RCP<AdjointAuxSensitivityModelEvaluator<Scalar> > adjoint_model_;
  Teuchos::RCP<IntegratorBasic<Scalar> > state_integrator_;
  Teuchos::RCP<IntegratorBasic<Scalar> > adjoint_integrator_;
  Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory_;
  int p_index_;
  int g_index_;
  bool g_depends_on_p_;
  bool f_depends_on_p_;
  bool ic_depends_on_p_;
  bool mass_matrix_is_identity_;
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > dxdp_init_;
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > dgdp_;
};

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<IntegratorAdjointSensitivity<Scalar> >
integratorAdjointSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<IntegratorAdjointSensitivity<Scalar> >
integratorAdjointSensitivity();

} // namespace Tempus

#endif // Tempus_IntegratorAdjointSensitivity_decl_hpp
