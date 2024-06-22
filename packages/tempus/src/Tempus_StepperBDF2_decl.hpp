//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperBDF2_decl_hpp
#define Tempus_StepperBDF2_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperBDF2AppAction.hpp"

namespace Tempus {

/** \brief BDF2 (Backward-Difference-Formula-2) time stepper.
 *
 *  For the implicit ODE system, \f$\mathcal{F}(\dot{x},x,t) = 0\f$,
 *  the solution, \f$\dot{x}\f$ and \f$x\f$, is determined using a
 *  solver (e.g., a non-linear solver, like NOX).  This stepper allows
 *  for a variable time-step, \f$\Delta t\f$.  It is a 2-step method.
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for BDF2 is
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} BDF2 \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt}
 \setlength{\parsep}{0pt}
 *      \item {\bf if ( "Startup", i.e., number of states $<$ 3) then}
 *      \item \quad {\bf Take timestep with startup Stepper.}
 *      \item \quad {\bf return}
 *      \item {\bf endif}
 *      \item {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *      \item {\bf Get old states $x_{n-1}$ and $x_{n-2}$.}
 *      \item {\bf Set ODE parameters.}
 *      \item \quad {\bf Time derivative: }
 *            $\dot{x}_{n} = \frac{2\tau_n + \tau_{n-1}}{\tau_n + \tau_{n-1}}
 *                  \left[ \frac{x_n-x_{n-1}}{\tau_n}\right]
 *                -  \frac{\tau_n}{\tau_n + \tau_{n-1}}
 *                   \left[ \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}\right]$
 *      \item \quad {\bf Alpha: $\alpha = \frac{2\tau_n + \tau_{n-1}}{(\tau_n +
 \tau_{n-1})\tau_n}$}
 *      \item \quad {\bf Beta: $\beta = 1$}
 *      \item {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}
 *      \item {\bf Solve $\mathcal{F}(\dot{x}_n,x_n,t_n) = 0$ for $x_n$.}
 *      \item {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}
 *      \item {\bf Update} $\dot{x}_{n} = \frac{2\tau_n + \tau_{n-1}}{\tau_n +
 \tau_{n-1}}
 *                  \left[ \frac{x_n-x_{n-1}}{\tau_n}\right]
 *                -  \frac{\tau_n}{\tau_n + \tau_{n-1}}
 *                   \left[ \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}\right]$
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *  The startup stepper allows BDF2 to use user-specified Stepper for the
 *  first timestep in order to populate the SolutionHistory with past states.
 *  A one-step startup stepper is perfect for this situation, e.g., Backward
 *  Euler or RK4.  The default startup stepper is 'DIRK 1 Stage Theta Method',
 *  which is second order accurate and allows an overall second-order solution.
 *
 *  The First-Same-As-Last (FSAL) principle is not needed for BDF2.
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
 *  the time derivative for BDF2 is
 *  \f[
 *    \dot{x}_n(x_n) = \frac{2\tau_n + \tau_{n-1}}{\tau_n + \tau_{n-1}}
 *                  \left[ \frac{x_n-x_{n-1}}{\tau_n}\right]
 *                -  \frac{\tau_n}{\tau_n + \tau_{n-1}}
 *                   \left[ \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}\right],
 *  \f]
 *  where \f$\Delta t_n = \tau_n = t_n - t_{n-1}\f$,
 *  we can determine that
 *  \f$ \alpha = \frac{2\tau_n + \tau_{n-1}}{(\tau_n + \tau_{n-1})\tau_n}\f$
 *  and \f$ \beta = 1 \f$, and therefore write
 *  \f[
 *    W = \frac{2\tau_n + \tau_{n-1}}{(\tau_n + \tau_{n-1})\tau_n}
          \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *      + \frac{\partial \mathcal{F}_n}{\partial x_n}.
 *  \f]
 */
template <class Scalar>
class StepperBDF2 : virtual public Tempus::StepperImplicit<Scalar> {
 public:
  /** \brief Default constructor.
   *
   *  - Requires the following calls before takeStep():
   *    setModel() and initialize().
   */
  StepperBDF2();

  StepperBDF2(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
      const Teuchos::RCP<Stepper<Scalar> >& startUpStepper, bool useFSAL,
      std::string ICConsistency, bool ICConsistencyCheck, bool zeroInitialGuess,
      const Teuchos::RCP<StepperBDF2AppAction<Scalar> >& stepperBDF2AppAction);

  /// \name Basic stepper methods
  //@{
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

  virtual void setAppAction(
      Teuchos::RCP<StepperBDF2AppAction<Scalar> > appAction);

  virtual Teuchos::RCP<StepperBDF2AppAction<Scalar> > getAppAction() const
  {
    return stepperBDF2AppAction_;
  }

  /// Set the stepper to use in first step
  void setStartUpStepper(std::string startupStepperType);
  void setStartUpStepper(Teuchos::RCP<Stepper<Scalar> > startupStepper);

  /// Initialize during construction and after changing input parameters.
  virtual void initialize();

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
  virtual Scalar getOrder() const { return order_; }
  virtual Scalar getOrderMin() const { return 1.0; }
  virtual Scalar getOrderMax() const { return 2.0; }

  virtual bool isExplicit() const { return false; }
  virtual bool isImplicit() const { return true; }
  virtual bool isExplicitImplicit() const
  {
    return isExplicit() && isImplicit();
  }
  virtual bool isOneStepMethod() const { return false; }
  virtual bool isMultiStepMethod() const { return !isOneStepMethod(); }
  virtual OrderODE getOrderODE() const { return FIRST_ORDER_ODE; }
  //@}

  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const { return getAlpha(dt, dt); }
  virtual Scalar getAlpha(const Scalar dt, const Scalar dtOld) const
  {
    return (Scalar(2.0) * dt + dtOld) / (dt * (dt + dtOld));
  }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta(const Scalar) const { return Scalar(1.0); }

  /// Compute the first time step given the supplied startup stepper
  virtual void computeStartUp(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const;

 private:
  Teuchos::RCP<Stepper<Scalar> > startUpStepper_;
  Teuchos::RCP<StepperBDF2AppAction<Scalar> > stepperBDF2AppAction_;
  Scalar order_ = Scalar(2.0);
};

/** \brief Time-derivative interface for BDF2.
 *
 *  Given the state \f$x_n\f$, compute the BDF2 time-derivative,
 *  \f[
 *    \dot{x}_{n} = \frac{2\tau_n + \tau_{n-1}}{\tau_n + \tau_{n-1}}
 *                  \left[ \frac{x_n-x_{n-1}}{\tau_n}\right]
 *                -  \frac{\tau_n}{\tau_n + \tau_{n-1}}
 *                   \left[ \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}\right]
 *  \f]
 *  where
 *  \f[
 *   \tau_n = t_n - t_{n-1}.
 *   \f]
 *  \f$\ddot{x}\f$ is not used and set to null.
 */
template <typename Scalar>
class StepperBDF2TimeDerivative
  : virtual public Tempus::TimeDerivative<Scalar> {
 public:
  /// Constructor
  StepperBDF2TimeDerivative(
      Scalar dt, Scalar dtOld,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOldOld)
  {
    initialize(dt, dtOld, xOld, xOldOld);
  }

  /// Destructor
  virtual ~StepperBDF2TimeDerivative() {}

  /// Compute the time derivative.
  virtual void compute(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xDot,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xDotDot = Teuchos::null)
  {
    xDotDot = Teuchos::null;
    // Calculate the BDF2 x dot vector
    const Scalar a = ((Scalar(2.0) * dt_ + dtOld_) / (dt_ + dtOld_)) / dt_;
    const Scalar b = (dt_ / (dt_ + dtOld_)) / dtOld_;
    // xDot = a*(x_n - x_{n-1}) - b*(x_{n-1} - x_{n-2})
    Thyra::V_StVpStV(xDot.ptr(), a, *x, -(a + b), *xOld_);
    Thyra::Vp_StV(xDot.ptr(), b, *xOldOld_);
  }

  virtual void initialize(
      Scalar dt, Scalar dtOld,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOldOld)
  {
    dt_      = dt;
    dtOld_   = dtOld;
    xOld_    = xOld;
    xOldOld_ = xOldOld;
  }

 private:
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOldOld_;
  Scalar dt_;     // = t_n - t_{n-1}
  Scalar dtOld_;  // = t_{n-1} - t_{n-2}
};

/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperBDF2<Scalar> > createStepperBDF2(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif  // Tempus_StepperBDF2_decl_hpp
