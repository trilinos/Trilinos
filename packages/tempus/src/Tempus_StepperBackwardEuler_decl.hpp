//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperBackwardEuler_decl_hpp
#define Tempus_StepperBackwardEuler_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperBackwardEulerAppAction.hpp"
#include "Tempus_StepperOptimizationInterface.hpp"

namespace Tempus {

/** \brief Backward Euler time stepper.
 *
 *  For the implicit ODE system, \f$\mathcal{F}(\dot{x},x,t) = 0\f$,
 *  the solution, \f$\dot{x}\f$ and \f$x\f$, is determined using a
 *  solver (e.g., a non-linear solver, like NOX).
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Backward Euler is
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Backward Euler \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt}
 * \setlength{\parsep}{0pt} \item {\it appAction.execute(solutionHistory,
 * stepper, BEGIN\_STEP)} \item {\bf Compute the predictor} (e.g., apply stepper
 * to $x_n$). \item {\it appAction.execute(solutionHistory, stepper,
 * BEGIN\_SOLVE)} \item {\bf Solve $\mathcal{F}_n(\dot{x}=(x_n-x_{n-1})/\Delta
 * t_n, x_n, t_n)=0$ for $x_n$} \item {\it appAction.execute(solutionHistory,
 * stepper, AFTER\_SOLVE)} \item $\dot{x}_n \leftarrow (x_n-x_{n-1})/\Delta t_n$
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *  The First-Same-As-Last (FSAL) principle is not needed with Backward Euler.
 *  The default is to set useFSAL=false, however useFSAL=true will also work
 *  but have no affect (i.e., no-op).
 *
 *  <b> Iteration Matrix, \f$W\f$.</b>
 *  Recalling that the definition of the iteration matrix, \f$W\f$, is
 *  \f[
 *    W = \alpha \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *      + \beta  \frac{\partial \mathcal{F}_n}{\partial x_n},
 *  \f]
 *  where \f$ \alpha \equiv \frac{\partial \dot{x}_n(x_n) }{\partial x_n}, \f$
 *  and \f$ \beta \equiv \frac{\partial x_n}{\partial x_n} = 1\f$, and
 *  the time derivative for Backward Euler is
 *  \f[
 *    \dot{x}_n(x_n) = \frac{x_n - x_{n-1}}{\Delta t},
 *  \f]
 *  we can determine that
 *  \f$ \alpha = \frac{1}{\Delta t} \f$
 *  and \f$ \beta = 1 \f$, and therefore write
 *  \f[
 *    W = \frac{1}{\Delta t}
 *        \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *      + \frac{\partial \mathcal{F}_n}{\partial x_n}.
 *  \f]
 */
template <class Scalar>
class StepperBackwardEuler
  : virtual public Tempus::StepperImplicit<Scalar>,
    virtual public Tempus::StepperOptimizationInterface<Scalar> {
 public:
  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
   */
  StepperBackwardEuler();

  /// Constructor
  StepperBackwardEuler(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
      const Teuchos::RCP<Stepper<Scalar> >& predictorStepper, bool useFSAL,
      std::string ICConsistency, bool ICConsistencyCheck, bool zeroInitialGuess,
      const Teuchos::RCP<StepperBackwardEulerAppAction<Scalar> >&
          stepperBEAppAction);

  /// \name Basic stepper methods
  //@{
  virtual void setAppAction(
      Teuchos::RCP<StepperBackwardEulerAppAction<Scalar> > appAction);

  virtual Teuchos::RCP<StepperBackwardEulerAppAction<Scalar> > getAppAction()
      const
  {
    return stepperBEAppAction_;
  }

  /// Set the predictor
  void setPredictor(std::string predictorType = "None");
  void setPredictor(Teuchos::RCP<Stepper<Scalar> > predictorStepper);

  /// Set the model
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
      override;

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState()
      override;
  virtual Scalar getOrder() const override { return 1.0; }
  virtual Scalar getOrderMin() const override { return 1.0; }
  virtual Scalar getOrderMax() const override { return 1.0; }

  virtual bool isExplicit() const override { return false; }
  virtual bool isImplicit() const override { return true; }
  virtual bool isExplicitImplicit() const override
  {
    return isExplicit() && isImplicit();
  }
  virtual bool isOneStepMethod() const override { return true; }
  virtual bool isMultiStepMethod() const override { return !isOneStepMethod(); }
  virtual OrderODE getOrderODE() const override { return FIRST_ORDER_ODE; }
  //@}

  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const override
  {
    return Scalar(1.0) / dt;
  }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta(const Scalar) const override { return Scalar(1.0); }

  /// Compute predictor given the supplied stepper
  virtual void computePredictor(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Return a valid ParameterList with current settings.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters()
      const override;

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(
      Teuchos::FancyOStream& out,
      const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const override;

  /// \name Implementation of StepperOptimizationInterface
  //@{
  virtual int stencilLength() const override;
  virtual void computeStepResidual(
      Thyra::VectorBase<Scalar>& residual,
      const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t, const Thyra::VectorBase<Scalar>& p,
      const int param_index) const override;
  virtual void computeStepJacobian(
      Thyra::LinearOpBase<Scalar>& jacobian,
      const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t, const Thyra::VectorBase<Scalar>& p,
      const int param_index, const int deriv_index) const override;
  virtual void computeStepParamDeriv(
      Thyra::LinearOpBase<Scalar>& deriv,
      const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t, const Thyra::VectorBase<Scalar>& p,
      const int param_index) const override;
  virtual void computeStepSolver(
      Thyra::LinearOpWithSolveBase<Scalar>& jacobian_solver,
      const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t, const Thyra::VectorBase<Scalar>& p,
      const int param_index) const override;
  //@}

 private:
  /// Implementation of computeStep*() methods
  void computeStepResidDerivImpl(
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs,
      const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t, const Thyra::VectorBase<Scalar>& p,
      const int param_index, const int deriv_index = 0) const;

 private:
  Teuchos::RCP<Stepper<Scalar> > predictorStepper_;
  Teuchos::RCP<StepperBackwardEulerAppAction<Scalar> > stepperBEAppAction_;
};

/** \brief Time-derivative interface for Backward Euler.
 *
 *  Given the state \f$x\f$, compute the Backward Euler time-derivative,
 *  \f[
 *    \dot{x}_{n} = \frac{(x_{n} - x_{n-1})}{\Delta t_{n}}.
 *  \f]
 *  \f$\ddot{x}\f$ is not used and set to null.
 */
template <typename Scalar>
class StepperBackwardEulerTimeDerivative
  : virtual public Tempus::TimeDerivative<Scalar> {
 public:
  /// Constructor
  StepperBackwardEulerTimeDerivative(
      Scalar s, Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld)
  {
    initialize(s, xOld);
  }

  /// Destructor
  virtual ~StepperBackwardEulerTimeDerivative() {}

  /// Compute the time derivative.
  virtual void compute(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xDot,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xDotDot = Teuchos::null)
  {
    xDotDot = Teuchos::null;
    // Calculate the Backward Euler x dot vector
    Thyra::V_StVpStV(xDot.ptr(), s_, *x, -s_, *xOld_);
  }

  virtual void initialize(Scalar s,
                          Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld)
  {
    s_    = s;
    xOld_ = xOld;
  }

 private:
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld_;
  Scalar s_;  // = 1.0/dt
};

/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperBackwardEuler<Scalar> > createStepperBackwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif  // Tempus_StepperBackwardEuler_decl_hpp
