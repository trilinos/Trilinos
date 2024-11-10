//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperForwardEuler_decl_hpp
#define Tempus_StepperForwardEuler_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExplicit.hpp"
#include "Tempus_StepperForwardEulerAppAction.hpp"

namespace Tempus {

/** \brief Forward Euler time stepper.
 *
 *  For the explicit ODE system,
 *  \f[
 *    \dot{x} = \bar{f}(x,t),
 *  \f]
 *  the Forward Euler stepper can be written as
 *  \f[
 *    x_{n} = x_{n-1} + \Delta t\, \bar{f}(x_{n-1},t_{n-1})
 *  \f]
 *  Forward Euler is an explicit time stepper (i.e., no solver used).
 *  Note that the time derivative by definition is
 *  \f[
 *    \dot{x}_{n} = \bar{f}(x_{n},t_{n}),
 *  \f]
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Forward Euler is
 *
 *  \f{center}
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Forward Euler \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt}
 * \setlength{\parsep}{0pt} \item {\it appAction.execute(solutionHistory,
 * stepper, BEGIN\_STEP)} \item {\bf if (Not ``Use FSAL'' or (previous step
 * failed)) then} \item \quad  {\it appAction.execute(solutionHistory, stepper,
 * BEFORE\_EXPLICIT\_EVAL)} \item \quad  $\dot{x}_{n-1} \leftarrow
 * \bar{f}(x_{n-1},t_{n-1})$ \item {\bf endif} \item $x_{n} \leftarrow x_{n-1} +
 * \Delta t\, \dot{x}_{n-1}$ \hfill {\it * Forward Euler update.} \item {\it
 * appAction.execute(solutionHistory, stepper, END\_STEP)} \item {\bf if (``Use
 * FSAL'') then} \item \quad  {\it appAction.execute(solutionHistory, stepper,
 * BEFORE\_EXPLICIT\_EVAL)} \item \quad  $\dot{x}_n \leftarrow
 * \bar{f}(x_{n},t_{n})$ \item {\bf endif} \item {\it
 * appAction.execute(solutionHistory, stepper, END\_STEP)} \end{enumerate}
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
 *  The default for Forward Euler is to use FSAL (useFSAL=true), but will
 *  also work with useFSAL=false.  Using useFSAL=true does assume that the
 *  solution, \f$x\f$, and its time derivative, \f$\dot{x}\f$, are consistent
 *  at the initial conditions (ICs), i.e.,
 *  \f$\dot{x}_{0} = \bar{f}(x_{0},t_{0})\f$.  This can be ensured by setting
 *  setICConsistency("Consistent"), and checked with
 *  setICConsistencyCheck(true).
 */
template <class Scalar>
class StepperForwardEuler : virtual public Tempus::StepperExplicit<Scalar> {
 public:
  /** \brief Default constructor.
   *
   *  - Requires subsequent setModel() and initialize() calls before calling
   *    takeStep().
   */
  StepperForwardEuler();

  /// Constructor
  StepperForwardEuler(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
      const Teuchos::RCP<StepperForwardEulerAppAction<Scalar> >&
          stepperFEAppAction);

  virtual void setAppAction(
      Teuchos::RCP<StepperForwardEulerAppAction<Scalar> > appAction);

  virtual Teuchos::RCP<StepperForwardEulerAppAction<Scalar> > getAppAction()
      const
  {
    return stepperFEAppAction_;
  }

  /// Set the initial conditions, make them consistent, and set needed memory.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
  virtual Scalar getOrder() const { return 1.0; }
  virtual Scalar getOrderMin() const { return 1.0; }
  virtual Scalar getOrderMax() const { return 1.0; }
  virtual void setUseFSAL(bool a)
  {
    this->useFSAL_       = a;
    this->isInitialized_ = false;
  }
  virtual OrderODE getOrderODE() const { return FIRST_ORDER_ODE; }
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const;

 protected:
  Teuchos::RCP<StepperForwardEulerAppAction<Scalar> > stepperFEAppAction_;
};

/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperForwardEuler<Scalar> > createStepperForwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif  // Tempus_StepperForwardEuler_decl_hpp
