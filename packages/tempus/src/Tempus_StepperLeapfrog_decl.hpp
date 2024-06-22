//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperLeapfrog_decl_hpp
#define Tempus_StepperLeapfrog_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExplicit.hpp"
#include "Tempus_StepperLeapfrogAppAction.hpp"
#include "Tempus_StepperLeapfrogAppActionComposite.hpp"

namespace Tempus {

/** \brief Leapfrog time stepper.
 *
 *  For the governing equation,
 *  \f[
 *    M(t) \ddot{x} + K(t) x = F(t),
 *  \f]
 *  one can write the explicit ODE system,
 *  \f[
 *    \ddot{x} = f(x,t),
 *  \f]
 *  where
 *  \f[
 *    f(x,t) = \left(M(t)\right)^{-1} \left( F(t) - K(t) x \right).
 *  \f]
 *  The Leapfrog stepper can be written as
 *  \f{eqnarray*}{
 *    x_{n+1}         & = & x_{n} + \Delta t\, \dot{x}_{n+1/2} \\
 *    \ddot{x}_{n+1}  & = & f(x_{n+1},t_{n+1}) \\
 *    \dot{x}_{n+3/2} & = & \dot{x}_{n+1/2} + \Delta t\, \ddot{x}_{n+1}
 *  \f}
 *  where the position and velocity are leapfrogged over each other.
 *  On startup the velocity half-step can be obtained with
 *  \f{eqnarray*}{
 *    \dot{x}_{n+1/2} & = & \dot{x}_{n} + \frac{1}{2} \Delta t\, \ddot{x}_{n} \\
 *    \dot{x}_{n+1/2} & = & \dot{x}_{n} + \frac{1}{2} \Delta t\, f(x_{n},t_{n})
 *  \f}
 *  and to complete the time step, the final velocity half-step is obtained
 *  with
 *  \f{eqnarray*}{
 *    \dot{x}_{n+1}&=&\dot{x}_{n+1/2} +\frac{1}{2} \Delta t\, \ddot{x}_{n+1} \\
 *    \dot{x}_{n+1}&=&\dot{x}_{n+1/2} +\frac{1}{2} \Delta t\, f(x_{n+1},t_{n+1})
 *  \f}
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Leapfrog is
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Leapfrog \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt}
 * \setlength{\parsep}{0pt} \item {\it appAction.execute(solutionHistory,
 * stepper, BEGIN\_STEP)} \item {\bf if (``Startup") then} \hfill {\it * Take
 * half-step startup.} \item \quad  $\dot{x}_{n+1/2} = \dot{x}_n +
 * \frac{1}{2}\Delta t \ddot{x}_n$ \item {\bf endif} \item {\it
 * appAction.execute(solutionHistory, stepper, BEFORE\_X\_UPDATE)} \item
 * $x_{n+1} = x_n + \Delta t \dot{x}_{n+1/2}$ \item {\it
 * appAction.execute(solutionHistory, stepper, BEFORE\_EXPLICIT\_EVAL)} \item
 * $\ddot{x}_{n+1} = f(x_{n+1},t_{n+1})$ \item {\it
 * appAction.execute(solutionHistory, stepper, BEFORE\_XDOT\_UPDATE)} \item {\bf
 * if (``Ending") then} \hfill {\it * Take half-step to get solution at the same
 * time level.} \item \quad $\dot{x}_{n+1} \leftarrow \dot{x}_{n+1/2}
 * +\frac{1}{2} \Delta t\, \ddot{x}_{n+1}$ \item {\bf else} \item \quad
 * $\dot{x}_{n+3/2} \leftarrow \dot{x}_{n+1/2} + \Delta t\, \ddot{x}_{n+1}$
 *      \item {\bf endif}
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *  where one can begin with \f$(x_n,\dot{x}_{n+1/2})\f$ or
 * \f$(x_n,\dot{x}_{n})\f$ and/or end with \f$(x_{n+1},\dot{x}_{n+3/2})\f$ or
 *  \f$(x_{n+1},\dot{x}_{n+1})\f$.
 *
 *  The First-Same-As-Last (FSAL) principle is not used with Leapfrog
 *  because of the algorithm's prescribed order of solution update.
 *  The default is to set useFSAL=false, and useFSAL=true will
 *  issue a warning that it will have no affect.
 */
template <class Scalar>
class StepperLeapfrog : virtual public Tempus::StepperExplicit<Scalar> {
 public:
  /** \brief Default constructor.
   *
   *  - Requires subsequent setModel() and initialize() calls before calling
   *    takeStep().
   */
  StepperLeapfrog();

  /// Constructor
  StepperLeapfrog(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
      const Teuchos::RCP<StepperLeapfrogAppAction<Scalar> >&
          stepperLFAppAction);

  /// \name Basic stepper methods
  //@{
  virtual void setAppAction(
      Teuchos::RCP<StepperLeapfrogAppAction<Scalar> > appAction);

  virtual Teuchos::RCP<StepperLeapfrogAppAction<Scalar> > getAppAction() const
  {
    return stepperLFAppAction_;
  }

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
  virtual Scalar getOrder() const { return 2.0; }
  virtual Scalar getOrderMin() const { return 2.0; }
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
  virtual OrderODE getOrderODE() const { return SECOND_ORDER_ODE; }
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const;

 protected:
  Teuchos::RCP<StepperLeapfrogAppAction<Scalar> > stepperLFAppAction_;
};

/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperLeapfrog<Scalar> > createStepperLeapfrog(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif  // Tempus_StepperLeapfrog_decl_hpp
