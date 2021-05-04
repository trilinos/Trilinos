// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorBasicOld_decl_hpp
#define Tempus_IntegratorBasicOld_decl_hpp

#include "Thyra_ModelEvaluator.hpp"

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_Integrator.hpp"
#include "Tempus_TimeStepControl.hpp"
#include "Tempus_IntegratorObserverBasic.hpp"
#include "Tempus_IntegratorObserverComposite.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace Tempus {


/** \brief Basic time integrator
 */
template<class Scalar>
class IntegratorBasicOld
  : virtual public Tempus::Integrator<Scalar>,
    virtual public Teuchos::ParameterListAcceptor
{
public:

  /// Constructor with ParameterList and model, and will be fully initialized.
  IntegratorBasicOld(
    Teuchos::RCP<Teuchos::ParameterList>                pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

  /// Constructor with model and "Stepper Type" and is fully initialized with default settings.
  IntegratorBasicOld(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType);

  /// Constructor that requires a subsequent setParameterList, setStepper, and initialize calls.
  IntegratorBasicOld();

  /// Constructor with ParameterList and models, and will be fully initialized.
  IntegratorBasicOld(
    Teuchos::RCP<Teuchos::ParameterList>                pList,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models);

  /// Destructor
  virtual ~IntegratorBasicOld() {}

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
    /// Return a copy of the Tempus ParameterList
    virtual Teuchos::RCP<Teuchos::ParameterList> getTempusParameterList()
      override { return tempusPL_; }
    virtual void setTempusParameterList(
      Teuchos::RCP<Teuchos::ParameterList> pl) override
    {
      if (tempusPL_==Teuchos::null) tempusPL_=Teuchos::parameterList("Tempus");
      if (pl != Teuchos::null) *tempusPL_ = *pl;
      this->setParameterList(Teuchos::null);
    }
  //@}

  /// \name Accessor methods
  //@{
    /// Get current time
    virtual Scalar getTime() const override
    {return solutionHistory_->getCurrentTime();}
    /// Get current index
    virtual int getIndex() const override
    {return solutionHistory_->getCurrentIndex();}
    /// Get Status
    virtual Status getStatus() const override
    {return integratorStatus_;}
    /// Get the Stepper
    virtual Teuchos::RCP<Stepper<Scalar> > getStepper() const override
    {return stepper_;}
    /// Set the Stepper
    virtual void setStepper(Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model);
    /// Set the Stepper
    virtual void setStepper(
      std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models);
    /// Set the Stepper
    virtual void setStepperWStepper(Teuchos::RCP<Stepper<Scalar> > stepper);
    /// Set the initial state which has the initial conditions
    virtual void initializeSolutionHistory(
      Teuchos::RCP<SolutionState<Scalar> > state = Teuchos::null);
    /// Set the initial state from Thyra::VectorBase(s)
    virtual void initializeSolutionHistory(Scalar t0,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0 = Teuchos::null,
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0 = Teuchos::null);
    /// Get the SolutionHistory
    virtual Teuchos::RCP<const SolutionHistory<Scalar> > getSolutionHistory() const override
      {return solutionHistory_;}
    /// Set the SolutionHistory
    virtual void setSolutionHistory(
      Teuchos::RCP<SolutionHistory<Scalar> > sh = Teuchos::null);
    /// Get the TimeStepControl
    virtual Teuchos::RCP<const TimeStepControl<Scalar> > getTimeStepControl() const override
      {return timeStepControl_;}
    virtual Teuchos::RCP<TimeStepControl<Scalar> > getNonConstTimeStepControl() override
      {return timeStepControl_;}
    /// Set the TimeStepControl
    virtual void setTimeStepControl(
      Teuchos::RCP<TimeStepControl<Scalar> > tsc = Teuchos::null);
    /// Get the Observer
    virtual Teuchos::RCP<IntegratorObserverComposite<Scalar> > getObserver()
      {return integratorObserver_;}
    /// Set the Observer
    virtual void setObserver(
      Teuchos::RCP<IntegratorObserver<Scalar> > obs = Teuchos::null);
    /// Initializes the Integrator after set* function calls
    virtual void initialize();
    //TODO: finish this
    /// Returns the IntegratorTimer_ for this Integrator
    virtual Teuchos::RCP<Teuchos::Time> getIntegratorTimer() const override
    { return integratorTimer_;}
    virtual Teuchos::RCP<Teuchos::Time> getStepperTimer() const override
    { return stepperTimer_;}

    /// Get current the solution, x
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getX() const
      {return solutionHistory_->getCurrentState()->getX();}
    /// Get current the time derivative of the solution, xdot
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getXDot() const
      {return solutionHistory_->getCurrentState()->getXDot();}
    /// Get current the second time derivative of the solution, xdotdot
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getXDotDot() const
      {return solutionHistory_->getCurrentState()->getXDotDot();}

    /// Get current state
    virtual Teuchos::RCP<SolutionState<Scalar> > getCurrentState()
      {return solutionHistory_->getCurrentState();}

    Teuchos::RCP<Teuchos::ParameterList> getIntegratorParameterList()
      { return integratorPL_; }

    //virtual Teuchos::RCP<Teuchos::Time> getIntegratorTimer() const
      //{return integratorTimer_;}

    virtual void setScreenOutputIndexInterval(int i)
    { integratorPL_->set("Screen Output Index Interval", i); }

    virtual int getScreenOutputIndexInterval() const
    { return integratorPL_->get<int>("Screen Output Index Interval"); }

    virtual void setScreenOutputIndexList(std::string s)
    { integratorPL_->set("Screen Output Index List", s); }

    virtual std::string getScreenOutputIndexList() const
    { return integratorPL_->get<std::string>("Screen Output Index List", ""); }
  //@}

  /// Parse when screen output should be executed
  void parseScreenOutput();

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl)
      override;
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList() override;
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList() override;
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters()
      const override;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    std::string description() const override;
    void describe(Teuchos::FancyOStream        & out,
                  const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

protected:

  Teuchos::RCP<Teuchos::ParameterList>      tempusPL_;
  Teuchos::RCP<Teuchos::ParameterList>      integratorPL_;
  Teuchos::RCP<SolutionHistory<Scalar> >    solutionHistory_;
  Teuchos::RCP<TimeStepControl<Scalar> >    timeStepControl_;
  Teuchos::RCP<IntegratorObserverComposite<Scalar> > integratorObserver_;
  Teuchos::RCP<Stepper<Scalar> >            stepper_;

  Teuchos::RCP<Teuchos::Time> integratorTimer_;
  Teuchos::RCP<Teuchos::Time> stepperTimer_;
  Scalar runtime_;

  std::vector<int> outputScreenIndices_;  ///< Vector of screen output indices.

  /** The integratorStatus is primarily in the WORKING Status, and
   *  PASSED/FAILED are noted at the end of the run.  A FAILED value
   *  is used to jump out of the time-integration loop.
   */
  Status integratorStatus_;
  bool isInitialized_;
};

/// Nonmember constructor
template<class Scalar>
Teuchos::RCP<IntegratorBasicOld<Scalar> > integratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

/// Nonmember constructor
template<class Scalar>
Teuchos::RCP<IntegratorBasicOld<Scalar> > integratorBasic(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  std::string stepperType);

/// Nonmember constructor
template<class Scalar>
Teuchos::RCP<IntegratorBasicOld<Scalar> > integratorBasic();

/// Nonmember constructor
template<class Scalar>
Teuchos::RCP<IntegratorBasicOld<Scalar> > integratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                pList,
  std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models);


} // namespace Tempus

#endif // Tempus_IntegratorBasicOld_decl_hpp
