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

template <class Scalar>
class ExponentialODEParameters {
 public:
  /// Constructor
  ExponentialODEParameters() : timeStepSize_(Scalar(0.0)), stageNumber_(0) {}

  /// Constructor
  ExponentialODEParameters(Scalar timeStepSize, int stageNumber = 0)
    : timeStepSize_(timeStepSize), stageNumber_(stageNumber)
  {
  }

  Scalar timeStepSize_;
  int stageNumber_;
};


namespace Tempus {

/** \brief Thyra Base interface for exponential time steppers.
 *
 * TODO: rest of the documentation is from StepperImplicit:
 *  For first-order ODEs, we can write the implicit ODE as
 *  \f[
 *    \mathcal{F}(\dot{x}_n,x_n,t_n) = 0
 *  \f]
 *  where \f$x_n\f$ is the solution vector, \f$\dot{x}\f$ is the
 *  time derivative, \f$t_n\f$ is the time and \f$n\f$ indicates
 *  the \f$n^{th}\f$ time level.  Note that \f$\dot{x}_n\f$ is
 *  different for each time stepper and is a function of other
 *  solution states, e.g., for Backward Euler,
 *  \f[
 *    \dot{x}_n(x_n) = \frac{x_n - x_{n-1}}{\Delta t}
 *  \f]
 *
 *  <b> Defining the Iteration Matrix</b>
 *
 *  Often we use Newton's method or one of its variations to solve
 *  for \f$x_n\f$, such as
 *  \f[
 *    \left[
 *    \frac{\partial}{\partial x_n}
 *    \left(
 *      \mathcal{F}(\dot{x}_n,x_n,t_n)
 *    \right)
 *    \right] \Delta x_n^\nu = - \mathcal{F}(\dot{x}_n^\nu,x_n^\nu,t_n)
 *  \f]
 *  where \f$\Delta x_n^\nu = x_n^{\nu+1} - x_n^\nu\f$ and \f$\nu\f$
 *  is the iteration index.  Using the chain rule for a function
 *  with multiple variables, we can write
 *  \f[
 *    \left[
 *    \frac{\partial \dot{x}_n(x_n) }{\partial x_n}
 *    \frac{\partial}{\partial \dot{x}_n}
 *    \left(
 *      \mathcal{F}(\dot{x}_n,x_n,t_n)
 *    \right)
 *    +
 *    \frac{\partial x_n}{\partial x_n}
 *    \frac{\partial}{\partial x_n}
 *    \left(
 *      \mathcal{F}(\dot{x}_n,x_n,t_n)
 *    \right)
 *    \right] \Delta x_n^\nu = - \mathcal{F}(\dot{x}_n^\nu,x_n^\nu,t_n)
 *  \f]
 *  Defining the iteration matrix, \f$W\f$, we have
 *  \f[
 *    W \Delta x_n^\nu = - \mathcal{F}(\dot{x}_n^\nu,x_n^\nu,t_n)
 *  \f]
 *  using \f$\mathcal{F}_n = \mathcal{F}(\dot{x}_n,x_n,t_n)\f$, where
 *  \f[
 *    W = \alpha \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *      + \beta  \frac{\partial \mathcal{F}_n}{\partial x_n}
 *  \f]
 *  and
 *  \f[
 *    W = \alpha M + \beta J
 *  \f]
 *  where
 *  \f[
 *    \alpha \equiv \frac{\partial \dot{x}_n(x_n) }{\partial x_n},
 *    \quad \quad
 *    \beta \equiv \frac{\partial x_n}{\partial x_n} = 1,
 *    \quad \quad
 *    M = \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n},
 *    \quad \quad
 *    J = \frac{\partial \mathcal{F}_n}{\partial x_n}
 *  \f]
 *  and \f$M\f$ is the mass matrix and \f$J\f$ is the Jacobian.
 *
 *  Note that sometimes it is helpful to set \f$\alpha=0\f$ and
 *  \f$\beta = 1\f$ to obtain the Jacobian, \f$J\f$, from the
 *  iteration matrix (i.e., ModelEvaluator), or set \f$\alpha=1\f$
 *  and \f$\beta = 0\f$ to obtain the mass matrix, \f$M\f$, from
 *  the iteration matrix (i.e., the ModelEvaluator).
 *
 *  As a concrete example, the time derivative for Backward Euler is
 *  \f[
 *    \dot{x}_n(x_n) = \frac{x_n - x_{n-1}}{\Delta t}
 *  \f]
 *  thus
 *  \f[
 *    \alpha \equiv \frac{\partial \dot{x}_n(x_n) }{\partial x_n}
 *           = \frac{1}{\Delta t}
 *    \quad \quad
 *    \beta \equiv \frac{\partial x_n}{\partial x_n} = 1
 *  \f]
 *  and the iteration matrix for Backward Euler is
 *  \f[
 *    W = \frac{1}{\Delta t} \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *      + \frac{\partial \mathcal{F}_n}{\partial x_n}
 *  \f]
 *
 *  <b> Dangers of multiplying through by \f$\Delta t\f$. </b>
 *  In some time-integration schemes, the application might want
 *  to multiply the governing equations by the time-step size,
 *  \f$\Delta t\f$, for scaling or other reasons.  Here we illustrate
 *  what that means and the complications that follow from it.
 *
 *  Starting with a simple implicit ODE and multiplying through by
 *  \f$\Delta t\f$, we have
 *  \f[
 *    \mathcal{F}_n = \Delta t \dot{x}_n - \Delta t \bar{f}(x_n,t_n) = 0
 *  \f]
 *  For the Backward Euler stepper, we recall from above that
 *  \f[
 *    \dot{x}_n(x_n) = \frac{x_n - x_{n-1}}{\Delta t}
 *    \quad\quad
 *    \alpha \equiv \frac{\partial \dot{x}_n(x_n) }{\partial x_n}
 *           = \frac{1}{\Delta t}
 *    \quad \quad
 *    \beta \equiv \frac{\partial x_n}{\partial x_n} = 1
 *  \f]
 *  and we can find for our simple implicit ODE, \f$\mathcal{F}_n\f$,
 *  \f[
 *    M = \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n} = \Delta t,
 *    \quad \quad
 *    J = \frac{\partial \mathcal{F}_n}{\partial x_n}
 *      = -\Delta t \frac{\partial \bar{f}_n}{\partial x_n}
 *  \f]
 *  Thus this iteration matrix, \f$W^\ast\f$, is
 *  \f[
 *    W^\ast = \alpha \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *           + \beta  \frac{\partial \mathcal{F}_n}{\partial x_n}
 *           = \frac{1}{\Delta t} \Delta t
 *           + (1) \left( - \Delta t \frac{\partial \bar{f}_n}{\partial x_n}
 *                 \right)
 *  \f]
 *  or simply
 *  \f[
 *    W^\ast = 1 - \Delta t \frac{\partial \bar{f}_n}{\partial x_n}
 *  \f]
 *  Note that \f$W^\ast\f$ is not the same as \f$W\f$ from above
 *  (i.e., \f$W = \frac{1}{\Delta t} - \frac{\partial \bar{f}_n}{\partial
 *  x_n}\f$).  But we should <b>not</b> infer from this is that
 *  \f$\alpha = 1\f$ or \f$\beta = -\Delta t\f$, as those definitions
 *  are unchanged (i.e., \f$\alpha \equiv \frac{\partial \dot{x}_n(x_n)}
 *  {\partial x_n} = \frac{1}{\Delta t}\f$ and \f$\beta \equiv
 *  \frac{\partial x_n}{\partial x_n} = 1 \f$).  However, the mass
 *  matrix, \f$M\f$, the Jacobian, \f$J\f$ and the residual,
 *  \f$-\mathcal{F}_n\f$, all need to include \f$\Delta t\f$ in
 *  their evaluations (i.e., be included in the ModelEvaluator
 *  return values for these terms).
 *
 *  <b> Dangers of explicitly including time-derivative definition.</b>
 *  If we explicitly include the time-derivative defintion for
 *  Backward Euler, we find for our simple implicit ODE,
 *  \f[
 *    \mathcal{F}_n = \frac{x_n - x_{n-1}}{\Delta t} - \bar{f}(x_n,t_n) = 0
 *  \f]
 *  that the iteration matrix is
 *  \f[
 *    W^{\ast\ast} =
 *             \alpha \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *           + \beta  \frac{\partial \mathcal{F}_n}{\partial x_n}
 *           = \frac{1}{\Delta t} (0)
 *           + (1) \left(\frac{1}{\Delta t}
 *                       - \frac{\partial \bar{f}_n}{\partial x_n}
 *                 \right)
 *  \f]
 *  or simply
 *  \f[
 *    W^{\ast\ast} =
 *      \frac{1}{\Delta t} - \frac{\partial \bar{f}_n}{\partial x_n}
 *  \f]
 *  which is the same as \f$W\f$ from above for Backward Euler, but
 *  again we should <b>not</b> infer that \f$\alpha = \frac{1}{\Delta
 *  t}\f$ or \f$\beta = -1\f$.  However the major drawback is the
 *  mass matrix, \f$M\f$, the Jacobian, \f$J\f$, and the residual,
 *  \f$-\mathcal{F}_n\f$ (i.e., the ModelEvaluator) are explicitly
 *  written only for Backward Euler.  The application would need
 *  to write separate ModelEvaluators for each stepper, thus
 *  destroying the ability to re-use the ModelEvaluator with any
 *  stepper.
 *
 */
template <class Scalar>
class StepperExponential : virtual public Tempus::Stepper<Scalar> {
 public:
  /// \name Basic exponential stepper methods
  //@{

  /** \brief Default constructor.
   *
   *  Requires subsequent setPhiEvaluator(), setModel(), and initialize()
  */
  StepperExponential();

  /// Set the PhiEvaluator
  virtual void setPhiEvaluator(
    const Teuchos::RCP<Tempus::PhiEvaluator<Scalar> >& phiEvaluator);
  /// Construct and set a default PhiEvaluator
  virtual void setDefaultPhiEvaluator();
  virtual Teuchos::RCP<Tempus::PhiEvaluator<Scalar> > getPhiEvaluator() const;

  /// Set the model
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
      override;
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel()
      const override;

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

  /// create InArgs from state and parameters
  Thyra::ModelEvaluatorBase::InArgs<Scalar>
  virtual createInArgsExponentialODE(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
      const Teuchos::RCP<ExponentialODEParameters<Scalar> >& p);

  /// evaluate an ODE residual for an exponential integrator
  virtual void evaluateExponentialODE(
      Teuchos::RCP<Thyra::VectorBase<Scalar> >& f,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
      const Teuchos::RCP<ExponentialODEParameters<Scalar> >& p);

  /// Pass initial guess to Newton solver (only relevant for implicit solvers)
  virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > initialGuess) override
  {
  }

  virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */)
      const override
  {
    // return a large value that should still fit into any supported scalar type
    return Scalar(Teuchos::ScalarTraits<Scalar>::rmax() / 1e2);
  }

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState() override;
  virtual void setUseFSAL(bool a) override
  {
    this->useFSAL_       = a;
    this->isInitialized_ = false;
  }

  virtual bool isExplicit() const override {return false;}
  /// Get the implicit/explicit type: we return true, since we rely on the implicit ModelEvaluator
  virtual bool isImplicit() const override {return true;}
  virtual bool isExplicitImplicit() const override
    {return isExplicit() && isImplicit();}
  virtual OrderODE getOrderODE() const override {return FIRST_ORDER_ODE;}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const override;

  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters()
      const override;
  //@}

  Teuchos::RCP<Teuchos::ParameterList> getValidParametersBasicExponential() const;

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(
      Teuchos::FancyOStream& out,
      const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  /// Set StepperExponential member data from the ParameterList.
  void setStepperExponentialValues(Teuchos::RCP<Teuchos::ParameterList> pl);

 protected:
  /// Compute the temporal finite difference dt_Mf_deriv
  ///   d/dt (-M * F(x,t))
  void computeTemporalFD(
    Teuchos::RCP<Thyra::VectorBase<Scalar>>& dt_Mf_deriv,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x,
    const Scalar t0,
    const Scalar dt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mf
  );

  // Check if temporal derivative is desired
  bool getTemporalDerivative()
  {
    return (temporal_finite_difference_eps_ > 0);
  }

  /// Compute the nonlinear remainder:
  ///   remf = -M * (F(xr,tr) - F(x0,t0) - J_{x0} * (xr-x0) - F'(t0) * (tr-t0))
  /// including multiple of negative mass matrix (-M).
  ///
  /// dt is the current time-step, not necessarily (tr-t0)
  /// Mf contains already evaluated -M*F(x0,t0)
  /// dt_Mf_deriv contains already evaluated dt*M*F'(t0)
  /// Mfr can optionally contain -M*F(xr,tr), if already pre-evaluated
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

  // Check if PhiEvaluator adaptivity is desired and return positive interval number
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
