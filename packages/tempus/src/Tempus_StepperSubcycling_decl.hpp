//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperSubcycling_decl_hpp
#define Tempus_StepperSubcycling_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExplicit.hpp"
#include "Tempus_StepperSubcyclingAppAction.hpp"
#include "Tempus_IntegratorBasic.hpp"

namespace Tempus {

/** \brief Subcycling time stepper.
 *
 *  This stepper wraps an IntegratorBasic object to perform the
 *  subcycling, thus it has all the capbilities of an IntegratorBasic
 *  with the following specializations and defaults:
 *   - Main Integrator operates independently from subcycling Stepper,
 *     and can be driven by another physics.
 *   - Need separate TimeStepControl for the subcycling so it can be
 *     driven by the physics (e.g., variable time stepping).
 *   - Finish the subcycling exactly on the full timestep.
 *   - No solution I/O within the subcycling.
 *   - No restart capability within subcycling, but still have restart
 *     capability from the full timestep.
 *   - Do not need to keep a solution history of the subcycling.
 *
 *  <b> Algorithm </b>
 *  The algorithm for Subcycling stepper is
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Subcycling \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt}
 * \setlength{\parsep}{0pt} \item {\it appAction.execute(solutionHistory,
 * stepper, BEGIN\_STEP)} \item {\bf Advance solution, $x_{n}$ from $x_{n-1}$,
 * by cycling substeppers.} \item {\it appAction.execute(solutionHistory,
 * stepper, END\_STEP)} \end{enumerate} \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 */
template <class Scalar>
class StepperSubcycling : virtual public Tempus::Stepper<Scalar> {
 public:
  /** \brief Default constructor.
   *
   *  - Requires subsequent setModel() and initialize() call before
   *    calling takeStep().
   */
  StepperSubcycling();

  /// Constructor
  StepperSubcycling(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      const Teuchos::RCP<IntegratorBasic<Scalar> >& integrator, bool useFSAL,
      std::string ICConsistency, bool ICConsistencyCheck,
      const Teuchos::RCP<StepperSubcyclingAppAction<Scalar> >&
          stepperSCAppAction);

  /// \name Basic stepper methods
  //@{
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

  virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel);

  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const
  {
    return scIntegrator_->getStepper()->getModel();
  }

  virtual void setAppAction(Teuchos::RCP<StepperSubcyclingAppAction<Scalar> >
                                appAction = Teuchos::null);

  virtual Teuchos::RCP<StepperSubcyclingAppAction<Scalar> > getAppAction() const
  {
    return stepperSCAppAction_;
  }

  /// Initialize during construction and after changing input parameters.
  virtual void initialize();

  virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const;

  /// Set the initial conditions, make them consistent, and set needed memory.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Pass initial guess to Newton solver (only relevant for implicit solvers)
  //  thus a no-op for explicit steppers.
  virtual void setInitialGuess(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess);

  virtual void setSolver(
      Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver = Teuchos::null);

  virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const;

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();

  virtual bool isExplicit() const;
  virtual bool isImplicit() const;
  virtual bool isExplicitImplicit() const;
  virtual bool isOneStepMethod() const;
  virtual bool isMultiStepMethod() const;

  virtual Scalar getOrder() const;
  virtual Scalar getOrderMin() const;
  virtual Scalar getOrderMax() const;
  virtual OrderODE getOrderODE() const;
  //@}

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Functions to set the subcycling stepper values.
  //@{
  virtual void setSubcyclingStepper(Teuchos::RCP<Stepper<Scalar> > stepper);
  virtual void setSubcyclingMinTimeStep(Scalar MinTimeStep);
  virtual void setSubcyclingInitTimeStep(Scalar InitTimeStep);
  virtual void setSubcyclingMaxTimeStep(Scalar MaxTimeStep);
  virtual void setSubcyclingMaxFailures(int MaxFailures);
  virtual void setSubcyclingMaxConsecFailures(int MaxConsecFailures);
  virtual void setSubcyclingScreenOutputIndexInterval(int i);
  virtual void setSubcyclingScreenOutputIndexList(std::string s);
  virtual void setSubcyclingTimeStepControlStrategy(
      Teuchos::RCP<TimeStepControlStrategy<Scalar> > tscs);
  virtual void setSubcyclingIntegratorObserver(
      Teuchos::RCP<IntegratorObserver<Scalar> > obs);
  virtual void setSubcyclingPrintDtChanges(bool printDtChanges);
  //@}

  /// \name Functions to get the subcycling stepper values.
  //@{
  virtual Teuchos::RCP<const Stepper<Scalar> > getSubcyclingStepper() const;
  virtual Scalar getSubcyclingMinTimeStep() const;
  virtual Scalar getSubcyclingInitTimeStep() const;
  virtual Scalar getSubcyclingMaxTimeStep() const;
  virtual std::string getSubcyclingStepType() const;
  virtual int getSubcyclingMaxFailures() const;
  virtual int getSubcyclingMaxConsecFailures() const;
  virtual int getSubcyclingScreenOutputIndexInterval() const;
  virtual std::string getSubcyclingScreenOutputIndexList() const;
  virtual Teuchos::RCP<TimeStepControlStrategy<Scalar> >
  getSubcyclingTimeStepControlStrategy() const;
  virtual Teuchos::RCP<IntegratorObserver<Scalar> >
  getSubcyclingIntegratorObserver() const;
  virtual bool getSubcyclingPrintDtChanges() const;
  //@}

 protected:
  Teuchos::RCP<StepperSubcyclingAppAction<Scalar> > stepperSCAppAction_;
  Teuchos::RCP<IntegratorBasic<Scalar> > scIntegrator_;
};

/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperSubcycling<Scalar> > createStepperSubcycling(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif  // Tempus_StepperSubcycling_decl_hpp
