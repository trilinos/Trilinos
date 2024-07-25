//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitAForm_decl_hpp
#define Tempus_StepperNewmarkImplicitAForm_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluatorSecondOrder.hpp"
#include "Tempus_StepperNewmarkImplicitAFormAppAction.hpp"

namespace Tempus {

/** \brief Newmark time stepper in acceleration form (a-form).
 *
 *  Here, we implement the Newmark scheme in predictor/corrector form;
 *  see equations (34)-(35) in: A. Mota, W. Klug, M. Ortiz,  "Finite element
 *  simulation of firearm injury to the human cranium",  Computational
 *  Mechanics 31(1) 115-121 (2003).
 *
 *  Newmark has two parameters: \f$\beta\f$
 *  and \f$\gamma\f$, both of which need to be in the range \f$[0,1]\f$.
 *  Newmark can be an explicit or implicit method, depending on
 *  the value of the \f$\beta\f$ parameter. If \f$\beta = 0\f$, the method
 *  is explicit.  Regardless of whether the method is implicit
 *  or explicit, a linear solve is required.  This linear solve can be
 *  optimized, however, for the explicit case by lumping the mass matrix.
 *  This optimization can be invoked by running "Newmark Explicit d-Form"
 *  Stepper through the Piro::TempusSolver class.
 *
 *  Newmark is second order accurate if \f$\gamma =  0.5\f$; otherwise it
 *  is first order accurate.  Some additional properties about the Newmark
 *  scheme can be found
 *<a href="http://opensees.berkeley.edu/wiki/index.php/Newmark_Method">here</a>.
 *
 *  The governing equation solved by this stepper is
 *  \f[
 *    \mathbf{M}\, \ddot{\mathbf{x}} + \mathbf{C}\, \dot{\mathbf{x}}
 *    + \mathbf{K}\, \mathbf{x} + \mathbf{F}(t) = 0
 *  \f]
 *  For the A-form (i.e., solving for the acceleration,
 *  \f$\mathbf{a} = \ddot{\mathbf{x}}\f$), we have the following implicit ODE
 *  \f[
 *       \mathbf{M}\, \mathbf{a} + \mathbf{C}\, \mathbf{v}
 *     + \mathbf{K}\, \mathbf{d} + \mathbf{F}(t) =
 *       \mathbf{f}(\mathbf{d}, \mathbf{v}, \mathbf{a}, t) = 0
 *  \f]
 *  where \f$\mathbf{v} = \dot{\mathbf{x}}\f$ and \f$\mathbf{d} = \mathbf{x}\f$.
 *
 *  <b> Algorithm </b>
 *  The algorithm for the Newmark implicit A-form with predictors and
 *  correctors is
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Newmark Implicit A-form \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt} \setlength{\parsep}{0pt}
 *      \item {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *      \item $\mathbf{d}^{\ast} = \mathbf{d}^{n-1} + \Delta t \mathbf{v}^{n-1}
 *                                     + \Delta t^2 (1-2 \beta) \mathbf{a}^{n-1} / 2$
 *      \item $\mathbf{v}^{\ast} = \mathbf{v}^{n-1} + \Delta t (1-\gamma) \mathbf{a}^{n-1}$
 *      \item {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}
 *      \item {\bf Solve
 *                 $\mathbf{f}(\mathbf{d}^n, \mathbf{v}^n, \mathbf{a}^n, t^n) = 0$
 *                 for $\mathbf{a}^n$ where} \\
 *                 $\mathbf{d}^n = \mathbf{d}^{\ast} + \beta \Delta t^2 \mathbf{a}^n$ \\
 *                 $\mathbf{v}^n = \mathbf{v}^{\ast} + \gamma \Delta t \mathbf{a}^n$
 *      \item {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}
 *      \item $\mathbf{d}^n = \mathbf{d}^{\ast} + \beta \Delta t^2 \mathbf{a}^n$
 *      \item $\mathbf{v}^n = \mathbf{v}^{\ast} + \gamma \Delta t \mathbf{a}^n$
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *  The First-Same-As-Last (FSAL) principle is not needed with Newmark
 *  Implicit A-Form.  The default is to set useFSAL=false, however
 *  useFSAL=true will also work but have no affect (i.e., no-op).
 *
 */
template <class Scalar>
class StepperNewmarkImplicitAForm
  : virtual public Tempus::StepperImplicit<Scalar> {
 public:
  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
   */
  StepperNewmarkImplicitAForm();

  /// Constructor
  StepperNewmarkImplicitAForm(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
      bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
      bool zeroInitialGuess, std::string schemeName, Scalar beta, Scalar gamma,
      const Teuchos::RCP<StepperNewmarkImplicitAFormAppAction<Scalar> >&
          stepperAppAction);

  /// \name Basic stepper methods
  //@{
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

  virtual void setAppAction(
      Teuchos::RCP<StepperNewmarkImplicitAFormAppAction<Scalar> > appAction);

  virtual Teuchos::RCP<StepperNewmarkImplicitAFormAppAction<Scalar> >
  getAppAction() const
  {
    return stepperNewmarkImpAppAction_;
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

  virtual bool isExplicit() const { return false; }
  virtual bool isImplicit() const { return true; }
  virtual bool isExplicitImplicit() const
  {
    return isExplicit() && isImplicit();
  }
  virtual bool isOneStepMethod() const { return true; }
  virtual bool isMultiStepMethod() const { return !isOneStepMethod(); }
  virtual void setUseFSAL(bool a) { this->setUseFSALTrueOnly(a); }
  virtual OrderODE getOrderODE() const { return SECOND_ORDER_ODE; }
  //@}

  /// Return W_xDotxDot_coeff = d(xDotDot)/d(xDotDot).
  virtual Scalar getW_xDotDot_coeff(const Scalar) const { return Scalar(1.0); }
  /// Return alpha = d(xDot)/d(xDotDot).
  virtual Scalar getAlpha(const Scalar dt) const { return gamma_ * dt; }
  /// Return beta  = d(x)/d(xDotDot).
  virtual Scalar getBeta(const Scalar dt) const { return beta_ * dt * dt; }

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

  void correctDisplacement(Thyra::VectorBase<Scalar>& d,
                           const Thyra::VectorBase<Scalar>& dPred,
                           const Thyra::VectorBase<Scalar>& a,
                           const Scalar dt) const;

  void setSchemeName(std::string schemeName);
  void setBeta(Scalar beta);
  void setGamma(Scalar gamma);

 private:
  std::string schemeName_;
  Scalar beta_;
  Scalar gamma_;

  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::RCP<StepperNewmarkImplicitAFormAppAction<Scalar> >
      stepperNewmarkImpAppAction_;
};

/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperNewmarkImplicitAForm<Scalar> >
createStepperNewmarkImplicitAForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif  // Tempus_StepperNewmarkImplicitAForm_decl_hpp
