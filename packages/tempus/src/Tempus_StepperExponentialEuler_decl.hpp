//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperExponentialEuler_decl_hpp
#define Tempus_StepperExponentialEuler_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExponential.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperExponentialEulerAppAction.hpp"

#include "Tempus_PhiEvaluator.hpp"

namespace Tempus {

/** \brief Exponential Euler time stepper.
 *
 *  For the explicit ODE system, \f$\dot{x} = \mathcal{F}(x,t)\f$,
 *  the solution, \f$x\f$, is determined using explicit evaluations of \f$\mathcal{F}\f$
 *  and exponentials of a linear operator \f$W\f$.
 *  TODO: We rely on the implicit model evaluator, this should be mentioned right here,
 *        and discussed in more detail later.
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Exponential Euler is:
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Exponential Euler \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt} \setlength{\parsep}{0pt}
 *      \item {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *      \item {\bf Compute the explicit evaluation}  $f_{n-1} = \mathcal{F}(x_{n-1}, t_{n-1})$.
 *      \item {\it appAction.execute(solutionHistory, stepper, BEFORE\_EXP)}
 *      \item {\bf Compute $d_{n-1} = \phi_1(\Delta{t_{n-1}} W_{n-1}) f_{n-1}$.}
 *      \item {\it appAction.execute(solutionHistory, stepper, AFTER\_EXP)}
 *      \item $\dot{x}_n \leftarrow x_{n-1} + \Delta{t_n}d_n$
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *  The First-Same-As-Last (FSAL) principle is not applicable for Exponential Euler.
 *  The default is to set useFSAL=false. However, setting useFSAL=true will instruct the
 *  Stepper to evaluate the forcing term \f$f_{n}\f$ at the new time step after completing the step
 *  sucessfully. The result will be saved in the xDot value of the SolutionHistory, and used
 *  in the beginning of the next step. As a side-effect, it will set tim-dependent Dirichlet
 *  conditions for \f$x_n\f$ (if those are part of the problem).
 *
 *  <b> Iteration Matrix, \f$W\f$.</b>
 *  Recalling that the definition of the iteration matrix, \f$W\f$, is
 *  \f[
 *    W = \frac{\partial \mathcal{F}}{\partial x},
 *  \f]
 *  to obtain that from the implicit model evaluator, we set
 *  \f$ \alpha = 0 \f$ and \f$ \beta = 1 \f$.
 *
 *  TODO: Need to document our use of implicit model evaluator:
 *        We rely on the PhiEvaluator to assumble weak Jacobian (including Mass) and mass matrix,
 *        and convert between integrated vectors Mf and f on the fly.
 */
template<class Scalar>
class StepperExponentialEuler :
    virtual public Tempus::StepperExponential<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  Requires subsequent setModel() and initialize()
   *  calls before calling takeStep().
  */
  StepperExponentialEuler();

  /// Constructor
  StepperExponentialEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    const Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> >& stepperEEAppAction);

  /// \name Basic stepper methods
  //@{
    virtual void setAppAction(
      Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > appAction);

    virtual Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > getAppAction() const
    { return stepperEEAppAction_; }

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

    virtual Scalar getOrder() const override { return 1.0; }
    // The exponential Euler is order 1 in general, and order 2 for autonomous problems
    // (and assuming exact Jacobian)
    virtual Scalar getOrderMin() const override { return 1.0; }
    virtual Scalar getOrderMax() const override { return 2.0;}

    virtual bool isOneStepMethod() const override {return true;}
    virtual bool isMultiStepMethod() const override {return !isOneStepMethod();}
  //@}

  /// Return a valid ParameterList with current settings.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream& out,
                          const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream & out) const override;

private:

  Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > stepperEEAppAction_;

};


/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template<class Scalar>
Teuchos::RCP<StepperExponentialEuler<Scalar> >
createStepperExponentialEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  Teuchos::RCP<Teuchos::ParameterList> pl);


} // namespace Tempus

#endif // Tempus_StepperExponentialEuler_decl_hpp
