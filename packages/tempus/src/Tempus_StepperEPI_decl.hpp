//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPI_decl_hpp
#define Tempus_StepperEPI_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExponential.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperEPIAppAction.hpp"

#include "Tempus_PhiEvaluator.hpp"



namespace Tempus {
// TODO: FURKAN: FIX THE DESCRIPTION
/** \brief Exponential Euler time stepper.
 *
 *  For the explicit ODE system, \f$\dot{x} = \mathcal{F}(x,t)\f$,
 *  the solution, \f$x\f$, is determined using explicit evaluations of \f$\mathcal{F}\f$
 *  and exponentials of a linear operator \f$W\f$.
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
 *      \item {\bf Compute $d_{n-1} = \exp(\Delta{t_{n-1}} W_{n-1}) f_{n-1}$.}
 *      \item {\it appAction.execute(solutionHistory, stepper, AFTER\_EXP)}
 *      \item $\dot{x}_n \leftarrow x_{n-1} + \Delta{t_n}d_n$
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *  The First-Same-As-Last (FSAL) principle is not needed with Exponential Euler.
 *  The default is to set useFSAL=false, however useFSAL=true will also work
 *  but have no affect (i.e., no-op).
 *
 *  <b> Iteration Matrix, \f$W\f$.</b>
 *  Recalling that the definition of the iteration matrix, \f$W\f$, is
 *  \f[
 *    W = \frac{\partial \mathcal{F}}{\partial x},
 *  \f]
 *  to obtain that from the implicit model evaluator, we set
 *  \f$ \alpha = 0 \f$ and \f$ \beta = 1 \f$.
 *  TODO: Need to use exponential model evaluator:
 *        explicit model evaluator plus W, or implicit model evaluator wrapper.
 */
template<class Scalar>
class StepperEPI :
    virtual public Tempus::StepperExponential<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
  */
  StepperEPI();

  /// Constructor
  StepperEPI(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction);

  /// \name Basic stepper methods
  //@{
    virtual void setAppAction(
      Teuchos::RCP<StepperEPIAppAction<Scalar> > appAction);

    virtual Teuchos::RCP<StepperEPIAppAction<Scalar> > getAppAction() const
    { return stepperEPIAppAction_; }

    /// Set the order
    void setOrder(Scalar order) {this->order_ = order;}

    /// Determine if we have a valid setup
    bool isValidSetup(Teuchos::FancyOStream & out) const;

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

    virtual Scalar getOrder() const override {return Scalar(order_);}
    virtual Scalar getOrderMin() const override {return Scalar(2.0);}
    virtual Scalar getOrderMax() const override {return Scalar(3.0);}

    virtual bool isOneStepMethod() const override {return false;}
    virtual bool isMultiStepMethod() const override {return !isOneStepMethod();}
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

 private:

  Teuchos::RCP<StepperEPIAppAction<Scalar> > stepperEPIAppAction_;

  /// temporal Integration order
  double order_;
};


/// Nonmember constructors - ModelEvaluator and ParameterList
// ----------------------------------------------------------------------------
/** \brief EPI2 Definition.
 *
 *  See Tempus_StepperEPI for additional details.
 */
template <class Scalar>
class StepperExponential_EPI2 : virtual public StepperEPI<Scalar> {
 public:
  StepperExponential_EPI2()
  {
    this->setStepperName("EPI2");
    this->setStepperType("EPI2");
    this->setUseFSAL(false);
    this->setICConsistency("Consistent");
    this->setICConsistencyCheck(false);
    this->setAppAction(Teuchos::null);
    this->setOrder(2.0);
  }

  StepperExponential_EPI2(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction)
  {
    this->setStepperName("EPI2");
    this->setStepperType("EPI2");
    this->setUseFSAL(useFSAL);
    this->setICConsistency(ICConsistency);
    this->setICConsistencyCheck(ICConsistencyCheck);
    this->setOrder(2.0);

    this->setAppAction(stepperEPIAppAction);

    if (appModel != Teuchos::null) {
      this->setModel(appModel);
      this->initialize();
    }
  }
};

/// Nonmember constructor for EPI2
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperExponential_EPI2<Scalar> >
createStepperExponential_EPI2(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperExponential_EPI2<Scalar>());
  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}


// ----------------------------------------------------------------------------
/** \brief EPI3 Definition.
 *
 *  See Tempus_StepperEPI for additional details.
 */
template <class Scalar>
class StepperExponential_EPI3 : virtual public StepperEPI<Scalar> {
 public:
  StepperExponential_EPI3()
  {
    this->setStepperName("EPI3");
    this->setStepperType("EPI3");
    this->setUseFSAL(false);
    this->setICConsistency("Consistent");
    this->setICConsistencyCheck(false);
    this->setAppAction(Teuchos::null);
    this->setOrder(3.0);
  }

  StepperExponential_EPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction)
  {
    this->setStepperName("EPI3");
    this->setStepperType("EPI3");
    this->setUseFSAL(useFSAL);
    this->setICConsistency(ICConsistency);
    this->setICConsistencyCheck(ICConsistencyCheck);
    this->setOrder(3.0);

    this->setAppAction(stepperEPIAppAction);

    if (appModel != Teuchos::null) {
      this->setModel(appModel);
      this->initialize();
    }
  }
};

/// Nonmember constructor for EPI3
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperExponential_EPI3<Scalar> >
createStepperExponential_EPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperExponential_EPI3<Scalar>());
  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

} // namespace Tempus

#endif // Tempus_StepperEPI_decl_hpp
