// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorBasic_decl_hpp
#define Tempus_IntegratorBasic_decl_hpp

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Time.hpp"
// Tempus
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
// Tempus
#include "Tempus_Integrator.hpp"
#include "Tempus_TimeStepControl.hpp"
#include "Tempus_IntegratorObserverBasic.hpp"

#include <string>

namespace Tempus {


/** \brief Basic time integrator
 */
template<class Scalar>
class IntegratorBasic : virtual public Tempus::Integrator<Scalar>
{
public:

  /** \brief Constructor with ParameterList and model, and will be fully initialized. */
  IntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList>                pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

  /** \brief Constructor with model and "Stepper Type" and is fully initialized with default settings. */
  IntegratorBasic(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType);

  /// Destructor
  /** \brief Constructor that requires a subsequent setParameterList, setStepper, and initialize calls. */
  IntegratorBasic();

  /// Destructor
  virtual ~IntegratorBasic() {}

  /// \name Basic integrator methods
  //@{
    /// Advance the solution to timeMax, and return true if successful.
    virtual bool advanceTime();
    /// Advance the solution to timeFinal, and return true if successful.
    virtual bool advanceTime(const Scalar timeFinal);
    /// Perform tasks before start of integrator.
    virtual void startIntegrator();
    /// Start time step.
    virtual void startTimeStep();
    /// Only accept step after meeting time step criteria.
    virtual void acceptTimeStep();
    /// Perform tasks after end of integrator.
    virtual void endIntegrator();
    /// Return a copy of the Tempus ParameterList
    virtual Teuchos::RCP<Teuchos::ParameterList> getTempusParameterList()
    { return tempusPL_; }
    virtual void setTempusParameterList(Teuchos::RCP<Teuchos::ParameterList> pl)
    {
      if (tempusPL_==Teuchos::null) tempusPL_=Teuchos::parameterList("Tempus");
      if (pl != Teuchos::null) *tempusPL_ = *pl;
      this->setParameterList(Teuchos::null);
    }
  //@}

  /// \name Accessor methods
  //@{
    /// Get current time
    virtual Scalar getTime() const {return solutionHistory_->getCurrentTime();}
    /// Get current index
    virtual Scalar getIndex()const {return solutionHistory_->getCurrentIndex();}
    /// Get the Stepper
    virtual Teuchos::RCP<Stepper<Scalar> > getStepper() const
      {return stepper_;}
    /// Set the Stepper
    virtual void setStepper(Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model);
    /// Set the Stepper
    virtual void setStepperWStepper(Teuchos::RCP<Stepper<Scalar> > stepper);
    /// Set the initial state which has the initial conditions
    virtual void setInitialState(
      Teuchos::RCP<SolutionState<Scalar> > state = Teuchos::null);
    /// Set the initial state from Thyra::VectorBase(s)
    virtual void setInitialState(Scalar t0,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > x0,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xdot0 = Teuchos::null,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xdotdot0 = Teuchos::null);
    /// Get the SolutionHistory
    virtual Teuchos::RCP<SolutionHistory<Scalar> > getSolutionHistory()
      {return solutionHistory_;}
    /// Set the SolutionHistory
    virtual void setSolutionHistory(
      Teuchos::RCP<SolutionHistory<Scalar> > sh = Teuchos::null);
    /// Get the TimeStepControl
    virtual Teuchos::RCP<TimeStepControl<Scalar> > getTimeStepControl()
      {return timeStepControl_;}
    /// Set the TimeStepControl
    virtual void setTimeStepControl(
      Teuchos::RCP<TimeStepControl<Scalar> > tsc = Teuchos::null);
    /// Get the Observer
    virtual Teuchos::RCP<IntegratorObserver<Scalar> > getObserver()
      {return integratorObserver_;}
    /// Set the Observer
    virtual void setObserver(
      Teuchos::RCP<IntegratorObserver<Scalar> > obs = Teuchos::null);
    /// Initializes the Integrator after set* function calls
    virtual void initialize();


    /// Get current the solution, x
    virtual Teuchos::RCP<Thyra::VectorBase<double> > getX() const
      {return solutionHistory_->getCurrentState()->getX();}
    /// Get current the time derivative of the solution, xdot
    virtual Teuchos::RCP<Thyra::VectorBase<double> > getXdot() const
      {return solutionHistory_->getCurrentState()->getXDot();}
    /// Get current the second time derivative of the solution, xdotdot
    virtual Teuchos::RCP<Thyra::VectorBase<double> > getXdotdot() const
      {return solutionHistory_->getCurrentState()->getXDotDot();}

    /// Get current state
    virtual Teuchos::RCP<SolutionState<Scalar> > getCurrentState()
      {return solutionHistory_->getCurrentState();}
  //@}

  /// Parse when screen output should be executed
  void parseScreenOutput();

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    std::string description() const;
    void describe(Teuchos::FancyOStream        & out,
                  const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

protected:

  Teuchos::RCP<Teuchos::ParameterList>      tempusPL_;
  Teuchos::RCP<Teuchos::ParameterList>      integratorPL_;
  Teuchos::RCP<SolutionHistory<Scalar> >    solutionHistory_;
  Teuchos::RCP<TimeStepControl<Scalar> >    timeStepControl_;
  Teuchos::RCP<IntegratorObserver<Scalar> > integratorObserver_;
  Teuchos::RCP<Stepper<Scalar> >            stepper_;

  Teuchos::RCP<Teuchos::Time> integratorTimer_;
  Teuchos::RCP<Teuchos::Time> stepperTimer_;

  std::vector<int> outputScreenIndices_;  ///< Vector of screen output indices.

  /** The integratorStatus is primarily in the WORKING Status, and
   *  PASSED/FAILED are noted at the end of the run.  A FAILED value
   *  is used to jump out of the time-integration loop.
   */
  Status integratorStatus_;
  bool isInitialized_;
};

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integratorBasic(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  std::string stepperType);

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integratorBasic();

} // namespace Tempus

#endif // Tempus_IntegratorBasic_decl_hpp
