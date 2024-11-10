//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperImplicit_decl_hpp
#define Tempus_StepperImplicit_decl_hpp

// Tempus
#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"

namespace Tempus {

/** \brief Thyra Base interface for implicit time steppers.
 *
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
class StepperImplicit : virtual public Tempus::Stepper<Scalar> {
 public:
  /// \name Basic implicit stepper methods
  //@{
  /// Set the model
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
      override;

  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel()
      const override
  {
    Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model;
    if (wrapperModel_ != Teuchos::null) model = wrapperModel_->getAppModel();
    return model;
  }

  virtual Teuchos::RCP<const WrapperModelEvaluator<Scalar> > getWrapperModel()
  {
    return wrapperModel_;
  }

  virtual void setDefaultSolver();

  /// Set solver.
  virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver) override;

  virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver()
      const override
  {
    return solver_;
  }

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const = 0;

  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta(const Scalar dt) const = 0;

  /// Solve implicit ODE, f(x, xDot, t, p) = 0.
  const Thyra::SolveStatus<Scalar> solveImplicitODE(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
      const Teuchos::RCP<ImplicitODEParameters<Scalar> >& p,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& y = Teuchos::null,
      const int index                                   = -1 /* index and y are for IMEX_RK_Partition */);

  /// Evaluate implicit ODE residual, f(x, xDot, t, p).
  void evaluateImplicitODE(
      Teuchos::RCP<Thyra::VectorBase<Scalar> >& f,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
      const Teuchos::RCP<ImplicitODEParameters<Scalar> >& p);

  /// Pass initial guess to Newton solver (only relevant for implicit solvers)
  virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > initialGuess) override
  {
    initialGuess_        = initialGuess;
    this->isInitialized_ = false;
  }

  /// Set parameter so that the initial guess is set to zero (=True) or use last
  /// timestep (=False).
  virtual void setZeroInitialGuess(bool zIG)
  {
    zeroInitialGuess_    = zIG;
    this->isInitialized_ = false;
  }

  virtual bool getZeroInitialGuess() const { return zeroInitialGuess_; }

  virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */)
      const override
  {
    return Scalar(1.0e+99);
  }
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(
      Teuchos::FancyOStream& out,
      const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const override;

  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters()
      const override;

  Teuchos::RCP<Teuchos::ParameterList> getValidParametersBasicImplicit() const;

  /// Set StepperImplicit member data from the ParameterList.
  void setStepperImplicitValues(Teuchos::RCP<Teuchos::ParameterList> pl);

  /// Set solver from ParameterList.
  void setStepperSolverValues(Teuchos::RCP<Teuchos::ParameterList> pl);

  /// Set the Solver Name
  void setSolverName(std::string i) { solverName_ = i; }
  /// Get the Solver Name.
  std::string getSolverName() const { return solverName_; }

 protected:
  Teuchos::RCP<WrapperModelEvaluator<Scalar> > wrapperModel_;
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > initialGuess_;
  bool zeroInitialGuess_;
  std::string solverName_;
};

}  // namespace Tempus
#endif  // Tempus_StepperImplicit_decl_hpp
