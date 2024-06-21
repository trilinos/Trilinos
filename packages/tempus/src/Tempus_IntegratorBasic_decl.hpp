//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorBasic_decl_hpp
#define Tempus_IntegratorBasic_decl_hpp

#include "Teuchos_Time.hpp"

#include "Thyra_ModelEvaluator.hpp"

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_Integrator.hpp"
#include "Tempus_TimeStepControl.hpp"
#include "Tempus_IntegratorObserverBasic.hpp"

namespace Tempus {

/** \brief Basic time integrator
 */
template <class Scalar>
class IntegratorBasic : virtual public Tempus::Integrator<Scalar> {
 public:
  /// Default constructor (requires calls to setModel and setSolutionHistory for
  /// initial conditions before calling initialize() to be fully constructed).
  IntegratorBasic();

  /// Full constructor
  IntegratorBasic(Teuchos::RCP<Stepper<Scalar> > stepper,
                  Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory,
                  Teuchos::RCP<TimeStepControl<Scalar> > timeStepControl,
                  Teuchos::RCP<IntegratorObserver<Scalar> > integratorObserver,
                  std::vector<int> outputScreenIndices,
                  int outputScreenInterval);

  /// Copy (a shallow copy)
  virtual void copy(Teuchos::RCP<IntegratorBasic<Scalar> > iB);

  /// Destructor
  virtual ~IntegratorBasic() {}

  /// \name Basic integrator methods
  //@{
  /// Advance the solution to timeMax, and return true if successful.
  virtual bool advanceTime();
  /// Advance the solution to timeFinal, and return true if successful.
  virtual bool advanceTime(const Scalar timeFinal) override;
  /// Perform tasks before start of integrator.
  virtual void startIntegrator();
  /// Start time step.
  virtual void startTimeStep();
  /// Check if time step has passed or failed.
  virtual void checkTimeStep();
  /// Perform tasks after end of integrator.
  virtual void endIntegrator();
  //@}

  /// \name Accessor methods
  //@{
  /// Get current time
  virtual Scalar getTime() const override
  {
    return solutionHistory_->getCurrentTime();
  }
  /// Get current index
  virtual int getIndex() const override
  {
    return solutionHistory_->getCurrentIndex();
  }
  /// Get Status
  virtual Status getStatus() const override { return integratorStatus_; }
  /// Set Status
  virtual void setStatus(const Status st) override { integratorStatus_ = st; }
  /// Get the Stepper
  virtual Teuchos::RCP<Stepper<Scalar> > getStepper() const override
  {
    return stepper_;
  }
  /// Set the model on the stepper.
  virtual void setModel(
      Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model);
  /// Set the Stepper
  virtual void setStepper(Teuchos::RCP<Stepper<Scalar> > stepper);
  /// Set the initial state which has the initial conditions
  virtual void initializeSolutionHistory(
      Teuchos::RCP<SolutionState<Scalar> > state = Teuchos::null);
  /// Set the initial state from Thyra::VectorBase(s)
  virtual void initializeSolutionHistory(
      Scalar t0, Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0    = Teuchos::null,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0 = Teuchos::null);
  /// Get the SolutionHistory
  virtual Teuchos::RCP<const SolutionHistory<Scalar> > getSolutionHistory()
      const override
  {
    return solutionHistory_;
  }
  /// Get the SolutionHistory
  virtual Teuchos::RCP<SolutionHistory<Scalar> > getNonConstSolutionHistory()
      override
  {
    return solutionHistory_;
  }
  /// Set the SolutionHistory
  virtual void setSolutionHistory(
      Teuchos::RCP<SolutionHistory<Scalar> > sh = Teuchos::null);
  /// Get the TimeStepControl
  virtual Teuchos::RCP<const TimeStepControl<Scalar> > getTimeStepControl()
      const override
  {
    return timeStepControl_;
  }
  virtual Teuchos::RCP<TimeStepControl<Scalar> > getNonConstTimeStepControl()
      override
  {
    return timeStepControl_;
  }
  /// Set the TimeStepControl
  virtual void setTimeStepControl(
      Teuchos::RCP<TimeStepControl<Scalar> > tsc = Teuchos::null);
  /// Get the Observer
  virtual Teuchos::RCP<IntegratorObserver<Scalar> > getObserver()
  {
    return integratorObserver_;
  }
  /// Set the Observer
  virtual void setObserver(
      Teuchos::RCP<IntegratorObserver<Scalar> > obs = Teuchos::null);
  /// Initializes the Integrator after set* function calls
  virtual void initialize();
  /// Return true if IntegratorBasic is initialized.
  bool isInitialized() { return isInitialized_; }

  // TODO: finish this
  /// Returns the IntegratorTimer_ for this Integrator
  virtual Teuchos::RCP<Teuchos::Time> getIntegratorTimer() const override
  {
    return integratorTimer_;
  }
  virtual Teuchos::RCP<Teuchos::Time> getStepperTimer() const override
  {
    return stepperTimer_;
  }

  /// Get current the solution, x
  virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getX() const
  {
    return solutionHistory_->getCurrentState()->getX();
  }
  /// Get current the time derivative of the solution, xdot
  virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getXDot() const
  {
    return solutionHistory_->getCurrentState()->getXDot();
  }
  /// Get current the second time derivative of the solution, xdotdot
  virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getXDotDot() const
  {
    return solutionHistory_->getCurrentState()->getXDotDot();
  }

  /// Get current state
  virtual Teuchos::RCP<SolutionState<Scalar> > getCurrentState()
  {
    return solutionHistory_->getCurrentState();
  }

  virtual void setScreenOutputIndexInterval(int i)
  {
    outputScreenInterval_ = i;
  }

  virtual int getScreenOutputIndexInterval() const
  {
    return outputScreenInterval_;
  }

  virtual void setScreenOutputIndexList(std::vector<int> indices)
  {
    outputScreenIndices_ = indices;
  }

  /// Parse when screen output should be executed
  virtual void setScreenOutputIndexList(std::string str);

  virtual std::vector<int> getScreenOutputIndexList() const
  {
    return outputScreenIndices_;
  }

  virtual std::string getScreenOutputIndexListString() const;
  //@}

  void parseScreenOutput()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "  IntegratorBasic::parseScreenOutput() --  "
                               "Should call setScreenOutputIndexList()\n");
  }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
  std::string description() const override;
  void describe(Teuchos::FancyOStream& out,
                const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  /// Set the Integrator Name
  void setIntegratorName(std::string i) { integratorName_ = i; }
  /// Get the Integrator Name.
  std::string getIntegratorName() const { return integratorName_; }

  /// Get the Integrator Type.
  std::string getIntegratorType() const { return integratorType_; }

 protected:
  /// Set the Integrator Type
  void setIntegratorType(std::string i);
  std::string integratorName_;  ///< integrator name used for I/O.
  std::string integratorType_;  ///< the integrator type.

  Teuchos::RCP<Stepper<Scalar> > stepper_;
  Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory_;
  Teuchos::RCP<TimeStepControl<Scalar> > timeStepControl_;
  Teuchos::RCP<IntegratorObserver<Scalar> > integratorObserver_;

  std::vector<int> outputScreenIndices_;  ///< Vector of screen output indices.
  int outputScreenInterval_;              ///< screen output interval.

  /** The integratorStatus is primarily in the WORKING Status, and
   *  PASSED/FAILED are noted at the end of the run.  A FAILED value
   *  is used to jump out of the time-integration loop.
   */
  Status integratorStatus_;
  bool isInitialized_;

  Teuchos::RCP<Teuchos::Time> integratorTimer_;
  Teuchos::RCP<Teuchos::Time> stepperTimer_;
};

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorBasic<Scalar> > createIntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList> pList, bool runInitialize = true);

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorBasic<Scalar> > createIntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    bool runInitialize = true);

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorBasic<Scalar> > createIntegratorBasic(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType);

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorBasic<Scalar> > createIntegratorBasic();

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorBasic<Scalar> > createIntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models,
    bool runInitialize = true);

}  // namespace Tempus

#endif  // Tempus_IntegratorBasic_decl_hpp
