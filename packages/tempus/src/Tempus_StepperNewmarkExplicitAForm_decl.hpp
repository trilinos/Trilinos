//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkExplicitAForm_decl_hpp
#define Tempus_StepperNewmarkExplicitAForm_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExplicit.hpp"
#include "Tempus_StepperNewmarkExplicitAFormAppAction.hpp"

namespace Tempus {

/** \brief Newmark Explicit time stepper.
 *
 *  This is the specific case of the more general Newmark time stepper
 *  where this stepper is explicit (\f$\beta = 0\f$) (i.e., no solver used).
 *
 *  The governing equation is solved by this stepper is
 *  \f[
 *    \mathbf{M}\, \ddot{\mathbf{x}} + \mathbf{C}\, \dot{\mathbf{x}}
 *    + \mathbf{K}\, \mathbf{x} = \mathbf{F}(t)
 *  \f]
 *  For the A-form (i.e., solving for the acceleration,
 *  \f$\mathbf{a} = \ddot{\mathbf{x}}\f$), we have the following explicit ODE
 *  \f[
 *     \mathbf{a} = -\mathbf{M}^{-1}\left[ \mathbf{C}\, \mathbf{v}
 *    + \mathbf{K}\, \mathbf{d} - \mathbf{F}(t) \right]
 *    = \bar{\mathbf{f}}(\mathbf{d}, \mathbf{v}, t)
 *  \f]
 *  where \f$\mathbf{v} = \dot{\mathbf{x}}\f$ and \f$\mathbf{d} = \mathbf{x}\f$.
 *
 *  <b> Algorithm </b>
 *  The algorithm for the Newmark Explicit A-form is
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Newmark Explicit A-form \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt}
 * \setlength{\parsep}{0pt} \item {\it appAction.execute(solutionHistory,
 * stepper, BEGIN\_STEP)} \item {\bf if (Not ``Use FSAL'' or (previous step
 * failed)) then} \item \quad  {\it appAction.execute(solutionHistory, stepper,
 * BEFORE\_EXPLICIT\_EVAL)} \item \quad  $\mathbf{a}^{n-1} =
 *                   \bar{\mathbf{f}}(\mathbf{d}^{n-1}, \mathbf{v}^{n-1},
 * t^{n-1})$ \item {\bf endif} \item $\mathbf{d}^{\ast} = \mathbf{d}^{n-1} +
 * \Delta t \mathbf{v}^{n-1}
 *                               + \Delta t^2 \mathbf{a}^{n-1} / 2$
 *      \item $\mathbf{v}^{\ast} =
 *            \mathbf{v}^{n-1} + \Delta t (1-\gamma) \mathbf{a}^{n-1}$
 *      \item $\mathbf{a}^{\ast} =
 *            \bar{\mathbf{f}}(\mathbf{d}^{\ast}, \mathbf{v}^{\ast}, t^{n-1})$
 *      \item $\mathbf{d}^n = \mathbf{d}^{\ast}$
 *      \item $\mathbf{v}^n =
 *            \mathbf{v}^{\ast} + \Delta t \gamma \mathbf{a}^{\ast}$
 *      \item {\bf if (``Use FSAL'') then}
 *      \item \quad  {\it appAction.execute(solutionHistory, stepper,
 * BEFORE\_EXPLICIT\_EVAL)} \item \quad  $\mathbf{a}^n =
 *                   \bar{\mathbf{f}}(\mathbf{d}^n, \mathbf{v}^n, t^n)$
 *      \item {\bf endif}
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *  Note that with useFSAL=false \f$x_n\f$ and \f$\dot{x}_{n-1}\f$ are not
 *  at the same time level at the end of the time step (i.e., they are not
 *  sync'ed).
 *
 *  To have them at the same time level, we can use the First-Same-As-Last
 *  (FSAL) principle where the function evaulation from the last time step
 *  can be used as the first function evalulation of the current step.
 *
 *  The default is to use FSAL (useFSAL=true), but will also work
 *  with useFSAL=false.  Using useFSAL=true does assume that the
 *  solution, \f$x\f$, its time derivative, \f$\dot{x}\f$, and its
 *  second time derivative, \f$\ddot{x}\f$, are consistent at the
 *  initial conditions (ICs), i.e.,
 *  \f$\ddot{x}_{0} = \bar{f}(x_{0},\dot{x}_{0},t_{0})\f$.
 *  This can be ensured by setting setICConsistency("Consistent"),
 *  and checked with setICConsistencyCheck(true).
 *
 */
template <class Scalar>
class StepperNewmarkExplicitAForm
  : virtual public Tempus::StepperExplicit<Scalar> {
 public:
  /** \brief Default constructor.
   *
   *  - Requires subsequent setModel() and initialize() calls before calling
   *    takeStep().
   */
  StepperNewmarkExplicitAForm();

  /// Constructor
  StepperNewmarkExplicitAForm(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
      Scalar gamma,
      const Teuchos::RCP<StepperNewmarkExplicitAFormAppAction<Scalar> >&
          stepperAppAction);

  virtual Teuchos::RCP<StepperNewmarkExplicitAFormAppAction<Scalar> >
  getAppAction() const
  {
    return stepperNewmarkExpAppAction_;
  }

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
  virtual Scalar getOrder() const
  {
    if (gamma_ == 0.5)
      return 2.0;
    else
      return 1.0;
  }
  virtual Scalar getOrderMin() const { return 1.0; }
  virtual Scalar getOrderMax() const { return 2.0; }
  virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */) const
  {
    return Scalar(1.0e+99);
  }

  virtual bool isExplicit() const { return true; }
  virtual bool isImplicit() const { return false; }
  virtual bool isExplicitImplicit() const
  {
    return isExplicit() && isImplicit();
  }
  virtual bool isOneStepMethod() const { return true; }
  virtual bool isMultiStepMethod() const { return !isOneStepMethod(); }
  virtual void setUseFSAL(bool a)
  {
    this->useFSAL_       = a;
    this->isInitialized_ = false;
  }
  virtual OrderODE getOrderODE() const { return SECOND_ORDER_ODE; }
  //@}

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const;

  void predictVelocity(Thyra::VectorBase<Scalar>& vPred,
                       const Thyra::VectorBase<Scalar>& v,
                       const Thyra::VectorBase<Scalar>& a,
                       const Scalar dt) const;

  void predictDisplacement(Thyra::VectorBase<Scalar>& dPred,
                           const Thyra::VectorBase<Scalar>& d,
                           const Thyra::VectorBase<Scalar>& v,
                           const Thyra::VectorBase<Scalar>& a,
                           const Scalar dt) const;

  void correctVelocity(Thyra::VectorBase<Scalar>& v,
                       const Thyra::VectorBase<Scalar>& vPred,
                       const Thyra::VectorBase<Scalar>& a,
                       const Scalar dt) const;

  void setGamma(Scalar gamma)
  {
    gamma_ = gamma;

    TEUCHOS_TEST_FOR_EXCEPTION(
        (gamma_ > 1.0) || (gamma_ < 0.0), std::logic_error,
        "Error in 'Newmark Explicit a-Form' stepper: invalid value of Gamma = "
            << gamma_ << ".  Please select 0 <= Gamma <= 1. \n");

    this->isInitialized_ = false;
  }

  virtual void setAppAction(
      Teuchos::RCP<StepperNewmarkExplicitAFormAppAction<Scalar> > appAction);

 protected:
  Scalar gammaDefault_;
  Scalar gamma_;
  Teuchos::RCP<StepperNewmarkExplicitAFormAppAction<Scalar> >
      stepperNewmarkExpAppAction_;
};

/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperNewmarkExplicitAForm<Scalar> >
createStepperNewmarkExplicitAForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif  // Tempus_StepperNewmarkExplicitAForm_decl_hpp
