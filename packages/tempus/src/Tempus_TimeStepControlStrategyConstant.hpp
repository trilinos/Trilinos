// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeStepControlStrategy_Constant_hpp
#define Tempus_TimeStepControlStrategy_Constant_hpp

#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControlStrategy.hpp"
#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperState.hpp"


namespace Tempus {

/** \brief StepControlStrategy class for TimeStepControl
 *
 */
template<class Scalar>
class TimeStepControlStrategyConstant
  : virtual public TimeStepControlStrategy<Scalar>
{
public:

  /// Default Constructor
  TimeStepControlStrategyConstant()
    : constantTimeStep_(0.0)
  {
    this->setStrategyType("Constant");
    this->setStepType("Constant");
    this->setName("Constant");
    this->initialize();
  }

  /// Full Constructor
  TimeStepControlStrategyConstant(Scalar constantTimeStep,
    std::string name = "Constant")
    : constantTimeStep_(constantTimeStep)
  {
    this->setStrategyType("Constant");
    this->setStepType("Constant");
    this->setName(name);
    this->initialize();
  }

  /// Destructor
  virtual ~TimeStepControlStrategyConstant(){}

  /** \brief Determine the time step size.*/
  virtual void setNextTimeStep(const TimeStepControl<Scalar> & tsc,
    Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory,
    Status & integratorStatus) override
  {
    using Teuchos::RCP;

    this->checkInitialized();

    RCP<SolutionState<Scalar> >workingState=solutionHistory->getWorkingState();
    const Scalar errorAbs = workingState->getErrorAbs();
    const Scalar errorRel = workingState->getErrorRel();
    Scalar dt = workingState->getTimeStep();

    RCP<Teuchos::FancyOStream> out = tsc.getOStream();
    Teuchos::OSTab ostab(out,1,"setNextTimeStep");


    // Check constant time step
    if ( dt != tsc.getInitTimeStep() ) {
      tsc.printDtChanges(workingState->getIndex(), dt, tsc.getInitTimeStep(),
                         "Resetting constant dt.");
      dt = tsc.getInitTimeStep();
    }

    // Stepper failure
    if (workingState->getSolutionStatus() == Status::FAILED) {
      *out << "Failure - Stepper failed and can not change time step size!\n"
           << "    Time step type == CONSTANT_STEP_SIZE\n" << std::endl;
      integratorStatus = FAILED;
      return;
    }

    // Absolute error failure
    if (errorAbs > tsc.getMaxAbsError()) {
      *out << "Failure - Absolute error failed and can not change time step!\n"
           << "  Time step type == CONSTANT_STEP_SIZE\n"
           << "  (errorAbs ="<<errorAbs<<") > (errorMaxAbs ="
           << tsc.getMaxAbsError() << ")" << std::endl;
      integratorStatus = FAILED;
      return;
    }

    // Relative error failure
    if (errorRel > tsc.getMaxRelError()) {
      *out << "Failure - Relative error failed and can not change time step!\n"
         << "  Time step type == CONSTANT_STEP_SIZE\n"
         << "  (errorRel ="<<errorRel<<") > (errorMaxRel ="
         << tsc.getMaxRelError() << ")" << std::endl;
      integratorStatus = FAILED;
      return;
    }

    // update dt
    workingState->setTimeStep(dt);

    // Set time from initial time, dt, and index to avoid numerical roundoff.
    const Scalar initTime = tsc.getInitTime();
    const int initIndex   = tsc.getInitIndex();
    const int index       = workingState->getIndex();
    const Scalar time     = (index-initIndex)*dt + initTime;
    workingState->setTime(time);
  }


  /// \name Overridden from Teuchos::Describable
  //@{
    std::string description() const override
    { return "Tempus::TimeStepControlStrategyConstant"; }

    void describe(Teuchos::FancyOStream          &out,
                  const Teuchos::EVerbosityLevel verbLevel) const override
    {
      Teuchos::OSTab ostab(out,2,"describe");
      out << description() << std::endl
          << "Strategy Type = " << this->getStrategyType() << std::endl
          << "Step Type     = " << this->getStepType() << std::endl
          << "Time Step     = " << getConstantTimeStep() << std::endl;
    }
  //@}

  /// Return ParameterList with current values.
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override
  {
    Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Time Step Control Strategy");

    pl->set<std::string>("Strategy Type", this->getStrategyType(), "Constant");
    pl->set<double>("Time Step", getConstantTimeStep());

    return pl;
  }


  virtual void initialize() const override
  {
    this->isInitialized_ = true;   // Only place where this is set to true!
  }

  virtual Scalar getConstantTimeStep() const { return constantTimeStep_; }

  virtual void setConstantTimeStep(Scalar dt)
  { constantTimeStep_ = dt; this->isInitialized_ = false; }


private:

  Scalar constantTimeStep_;   ///< Constant time step size.

};


/// Nonmember constructor.
template <class Scalar>
Teuchos::RCP<TimeStepControlStrategyConstant<Scalar> >
createTimeStepControlStrategyConstant(
  const Teuchos::RCP<Teuchos::ParameterList> & pList,
  std::string name = "Constant")
{
  auto tscs = Teuchos::rcp(new TimeStepControlStrategyConstant<Scalar>());
  if (pList == Teuchos::null) return tscs;

  TEUCHOS_TEST_FOR_EXCEPTION(
    pList->get<std::string>("Strategy Type", "Constant") != "Constant",
    std::logic_error,
    "Error - Strategy Type != 'Constant'.  (='"
    +pList->get<std::string>("Strategy Type")+"')\n");

  pList->validateParametersAndSetDefaults(*tscs->getValidParameters(), 0);

  tscs->setConstantTimeStep(pList->get<double>("Time Step"));

  tscs->setName(name);
  tscs->initialize();

  return tscs;
}


/// Nonmember function to return ParameterList with default values.
template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList> getTimeStepControlStrategyConstantPL()
{
  auto t = rcp(new Tempus::TimeStepControlStrategyConstant<Scalar>());
  return Teuchos::rcp_const_cast<Teuchos::ParameterList> (t->getValidParameters());
}


} // namespace Tempus
#endif // Tempus_TimeStepControlStrategy_Constant_hpp
