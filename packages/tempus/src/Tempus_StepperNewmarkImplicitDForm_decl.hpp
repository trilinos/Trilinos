//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitDForm_decl_hpp
#define Tempus_StepperNewmarkImplicitDForm_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluatorSecondOrder.hpp"
#include "Tempus_StepperNewmarkImplicitDFormAppAction.hpp"

namespace Tempus {

/** \brief Newmark time stepper.
 *
 *  Here, we implement the Newmark scheme in displacement predictor/corrector
 *  form; see equations (34)-(35) in: A. Mota, W. Klug, M. Ortiz,
 *  "Finite element simulation of firearm injury to the human cranium",
 *  Computational Mechanics 31(1) 115-121 (2003).
 *
 *  Newmark has two parameters: \f$\beta\f$
 *  and \f$\gamma\f$, both of which need to be in the range \f$[0,1]\f$.
 *  Newmark can be an explicit or implicit method, depending on
 *  the value of the \f$\beta\f$ parameter. If \f$\beta = 0\f$, the method
 *  is explicit; but note that the d-form of the Newmark scheme is not defined
 *  for the explicit case.
 *
 *  Newmark is second order accurate if \f$\gamma =  0.5\f$; otherwise it
 *  is first order accurate.  Some additional properties about the Newmark
 *  Beta scheme can be found
 *  <a
 * href="http://opensees.berkeley.edu/wiki/index.php/Newmark_Method">here</a>.
 *
 *  <b> Algorithm </b>
 *  The algorithm for the Newmark implicit D-form with predictors and
 *  correctors is
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Newmark Implicit D-form \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt}
 * \setlength{\parsep}{0pt} \item {\it appAction.execute(solutionHistory,
 * stepper, BEGIN\_STEP)} \item $\mathbf{d}^{\ast} = \mathbf{d}^{n-1} + \Delta t
 * \mathbf{v}^{n-1}
 *                               + \Delta t^2 (1-2 \beta) \mathbf{a}^{n-1} / 2$
 *      \item $\mathbf{v}^{\ast} = \mathbf{v}^{n-1} + \Delta t (1-\gamma)
 * \mathbf{a}^{n-1}$ \item {\it appAction.execute(solutionHistory, stepper,
 * BEFORE\_SOLVE)} \item {\bf Solve
 *            $\mathbf{f}(\mathbf{d}^n, \mathbf{v}^n, \mathbf{a}^n, t^n) = 0$
 *            for $\mathbf{d}^n$ where} \\
 *            $\mathbf{a}^n = (\mathbf{d}^n - \mathbf{d}^{\ast})/(\beta \Delta
 * t^2)$ \\
 *            $\mathbf{v}^n = \mathbf{v}^{\ast} + \gamma \Delta t \mathbf{a}^n$
 *      \item {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}
 *      \item $\mathbf{a}^n = (\mathbf{d}^n - \mathbf{d}^{\ast})/(\beta \Delta
 * t^2)$ \item $\mathbf{v}^n = \mathbf{v}^{\ast} + \gamma \Delta t \mathbf{a}^n$
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *  The First-Same-As-Last (FSAL) principle is not used with the
 *  Newmark implicit D-Form method.  The default is to set useFSAL=false,
 *  however useFSAL=true will also work but have no affect (i.e., no-op).
 *
 */
template <class Scalar>
class StepperNewmarkImplicitDForm
  : virtual public Tempus::StepperImplicit<Scalar> {
 public:
  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
   */
  StepperNewmarkImplicitDForm();

  /// Constructor
  StepperNewmarkImplicitDForm(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>& appModel,
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar>>& solver,
      bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
      bool zeroInitialGuess, std::string schemeName, Scalar beta, Scalar gamma,
      const Teuchos::RCP<StepperNewmarkImplicitDFormAppAction<Scalar>>&
          stepperAppAction);

  /// \name Basic stepper methods
  //@{
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>& appModel);

  virtual Teuchos::RCP<StepperNewmarkImplicitDFormAppAction<Scalar>>
  getAppAction() const
  {
    return stepperNewmarkImpAppAction_;
  }

  virtual void setAppAction(
      Teuchos::RCP<StepperNewmarkImplicitDFormAppAction<Scalar>> appAction);

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar>>& /* solutionHistory */)
  {
  }

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar>>& solutionHistory);

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar>> getDefaultStepperState();
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
  virtual OrderODE getOrderODE() const { return SECOND_ORDER_ODE; }
  //@}

  /// Return W_xDotxDot_coeff = d(xDotDot)/d(x).
  virtual Scalar getW_xDotDot_coeff(const Scalar dt) const
  {
    return Scalar(1.0) / (beta_ * dt * dt);
  }
  /// Return alpha = d(xDot)/d(x).
  virtual Scalar getAlpha(const Scalar dt) const
  {
    return gamma_ / (beta_ * dt);
  }
  /// Return beta  = d(x)/d(x).
  virtual Scalar getBeta(const Scalar) const { return Scalar(1.0); }

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

  void correctAcceleration(Thyra::VectorBase<Scalar>& a,
                           const Thyra::VectorBase<Scalar>& dPred,
                           const Thyra::VectorBase<Scalar>& d,
                           const Scalar dt) const;

  void setSchemeName(std::string schemeName);
  void setBeta(Scalar beta);
  void setGamma(Scalar gamma);

 protected:
  std::string schemeName_;
  Scalar beta_;
  Scalar gamma_;

  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::RCP<StepperNewmarkImplicitDFormAppAction<Scalar>>
      stepperNewmarkImpAppAction_;
};

/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperNewmarkImplicitDForm<Scalar>>
createStepperNewmarkImplicitDForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif  // Tempus_StepperNewmarkImplicitDForm_decl_hpp
