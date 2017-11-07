// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeStepControl_impl_hpp
#define Tempus_TimeStepControl_impl_hpp

// Teuchos
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"


namespace Tempus {

// TimeStepControl definitions:
template<class Scalar>
TimeStepControl<Scalar>::TimeStepControl()
  : outputAdjustedDt_(false), dtAfterOutput_(0.0)
{
  tscPL_->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setParameterList(tscPL_);
}

template<class Scalar>
TimeStepControl<Scalar>::TimeStepControl(
  Teuchos::RCP<Teuchos::ParameterList> pList)
  : outputAdjustedDt_(false), dtAfterOutput_(0.0)
{
  this->setParameterList(pList);
}

template<class Scalar>
TimeStepControl<Scalar>::TimeStepControl(const TimeStepControl<Scalar>& tsc_)
  : tscPL_           (tsc_.tscPL_           ),
    outputIndices_   (tsc_.outputIndices_   ),
    outputTimes_     (tsc_.outputTimes_     ),
    outputAdjustedDt_(tsc_.outputAdjustedDt_),
    dtAfterOutput_   (tsc_.dtAfterOutput_   )
{}


template<class Scalar>
void TimeStepControl<Scalar>::getNextTimeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> > & solutionHistory,
  Status & integratorStatus)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::TimeStepControl::getNextTimeStep()");
  {
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionStateMetaData<Scalar> > metaData_ = workingState->getMetaData();
    const Scalar time = metaData_->getTime();
    const int iStep = metaData_->getIStep();
    const Scalar errorAbs = metaData_->getErrorAbs();
    const Scalar errorRel = metaData_->getErrorRel();
    int order = metaData_->getOrder();
    Scalar dt = metaData_->getDt();
    bool output = metaData_->getOutput();

    RCP<StepperState<Scalar> > stepperState = workingState->getStepperState();

    output = false;

    // If last time step was adjusted for output, reinstate previous dt.
    if (outputAdjustedDt_ == true) {
      dt = dtAfterOutput_;
      outputAdjustedDt_ = false;
      dtAfterOutput_ = 0.0;
    }

    if (dt < getMinTimeStep()) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"getNextTimeStep");
      *out << "Warning - Time step size (=" << dt << ") is less than\n"
           << "  minimum time step size (=" << getMinTimeStep() << ").\n"
           << "  Resetting to minimum time step size." << std::endl;
      dt = getMinTimeStep();
    }

    if (getStepType() == "Constant") {

      dt = getInitTimeStep();

      // Stepper failure
      if (stepperState->stepperStatus_ == Status::FAILED) {
        if (order+1 <= getMaxOrder()) {
          order++;
          RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out,1,"getNextTimeStep");
          *out << "Warning - Stepper failure with constant time step.\n"
               << "  Try increasing order.  order = " << order << std::endl;
        } else {
          RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out,1,"getNextTimeStep");
          *out << "Failure - Stepper failed and can not change time step size "
               << "or order!\n"
               << "    Time step type == CONSTANT_STEP_SIZE\n"
               << "    order = " << order << std::endl;
          integratorStatus = FAILED;
          return;
        }
      }

      // Absolute error failure
      if (errorAbs > getMaxAbsError()) {
        if (order+1 <= getMaxOrder()) {
          order++;
          RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out,1,"getNextTimeStep");
          *out
            <<"Warning - Absolute error is too large with constant time step.\n"
            <<"  (errorAbs ="<<errorAbs<<") > (errorMaxAbs ="
            <<getMaxAbsError()<<")"
            <<"  Try increasing order.  order = " << order << std::endl;
        } else {
          RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out,1,"getNextTimeStep");
          *out
            <<"Failure - Absolute error failed and can not change time step "
            <<"size or order!\n"
            <<"  Time step type == CONSTANT_STEP_SIZE\n"
            <<"  order = " << order
            <<"  (errorAbs ="<<errorAbs<<") > (errorMaxAbs ="
            <<getMaxAbsError()<<")"
            << std::endl;
          integratorStatus = FAILED;
          return;
        }
      }

      // Relative error failure
      if (errorRel > getMaxRelError()) {
        if (order+1 <= getMaxOrder()) {
          order++;
          RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out,1,"getNextTimeStep");
          *out
            <<"Warning - Relative error is too large with constant time step.\n"
            <<"  (errorRel ="<<errorRel<<") > (errorMaxRel ="
            <<getMaxRelError()<<")"
            <<"  Try increasing order.  order = " << order << std::endl;
        } else {
          RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out,1,"getNextTimeStep");
          *out
            <<"Failure - Relative error failed and can not change time step "
            <<"size or order!\n"
            <<"  Time step type == CONSTANT_STEP_SIZE\n"
            <<"  order = " << order
            <<"  (errorRel ="<<errorRel<<") > (errorMaxRel ="
            <<getMaxRelError()<<")"
            << std::endl;
          integratorStatus = FAILED;
          return;
        }
      }

      // Consistency checks
      TEUCHOS_TEST_FOR_EXCEPTION(
        (order < getMinOrder() || order > getMaxOrder()), std::out_of_range,
        "Error - Solution order is out of range and can not change "
        "time step size!\n"
        "    Time step type == CONSTANT_STEP_SIZE\n"
        "    [order_min, order_max] = [" <<getMinOrder()<< ", "
        <<getMaxOrder()<< "]\n"
        "    order = " << order << "\n");

    } else { // VARIABLE_STEP_SIZE

      // \todo The following controls should be generalized to plugable options.
      if (stepperState->stepperStatus_ == Status::FAILED) dt /=2;
      if (errorAbs > getMaxAbsError()) dt /= 2;
      if (errorRel > getMaxRelError()) dt /= 2;
      if (order < getMinOrder()) dt *= 2;
      if (order > getMaxOrder()) dt /= 2;

      if (dt < getMinTimeStep()) dt = getMinTimeStep();
      if (dt > getMaxTimeStep()) dt = getMaxTimeStep();
    }

    // Adjust time step to hit final time or correct for small
    // numerical differences.
    Scalar reltol = 1.0e-6;
    if ((time + dt > getFinalTime() ) ||
        (std::abs((time+dt-getFinalTime())/(time+dt)) < reltol))
      dt = getFinalTime() - time;

    // Check if we need to output this step index
    std::vector<int>::const_iterator it =
      std::find(outputIndices_.begin(), outputIndices_.end(), iStep+1);
    if (it != outputIndices_.end()) output = true;

    // Adjust time step to hit output times.
    for (size_t i=0; i < outputTimes_.size(); ++i) {
      const Scalar oTime = outputTimes_[i];
      if (time < oTime && oTime <= time+dt+getMinTimeStep()) {
        output = true;
        outputAdjustedDt_ = true;
        dtAfterOutput_ = dt;
        if (time < oTime && oTime <= time+dt-getMinTimeStep()) {
          // Next output time is not near next time
          // (>getMinTimeStep() away from it).
          // Take time step to hit output time.
          dt = oTime - time;
        } else if (std::abs((time+dt-oTime)/(time+dt)) < reltol) {
          // Next output time IS VERY near next time (<reltol away from it),
          // e.g., adjust for numerical roundoff.
          dt = oTime - time;
        } else if (time+dt-getMinTimeStep() < oTime &&
                   oTime <= time+dt+getMinTimeStep()) {
          // Next output time IS near next time (<getMinTimeStep() away from it)
          // Take two time steps to get to next output time.
          dt = (oTime - time)/2.0;
        }
        break;
      }
    }

    // Time step always needs to keep time within range.
    TEUCHOS_TEST_FOR_EXCEPTION(
      (time + dt < getInitTime()), std::out_of_range,
      "Error - Time step does not move time INTO time range.\n"
      "    [timeMin, timeMax] = [" << getInitTime() << ", "
      << getFinalTime() << "]\n"
      "    T + dt = " << time <<" + "<< dt <<" = " << time + dt << "\n");

    TEUCHOS_TEST_FOR_EXCEPTION(
      (time + dt > getFinalTime()), std::out_of_range,
      "Error - Time step move time OUT OF time range.\n"
      "    [timeMin, timeMax] = [" << getInitTime() << ", "
      << getFinalTime() << "]\n"
      "    T + dt = " << time <<" + "<< dt <<" = " << time + dt << "\n");

    metaData_->setOrder(order);
    metaData_->setDt(dt);
    metaData_->setOutput(output);
  }
  return;
}


/// Test if time is within range: include timeMin and exclude timeMax.
template<class Scalar>
bool TimeStepControl<Scalar>::timeInRange(const Scalar time) const{
  const Scalar relTol = 1.0e-14;
  bool tir = (getInitTime()*(1.0-relTol) <= time and
              time < getFinalTime()*(1.0-relTol));
  return tir;
}


template<class Scalar>
bool TimeStepControl<Scalar>::indexInRange(const int iStep) const{
  bool iir = (getInitIndex() <= iStep and iStep < getFinalIndex());
  return iir;
}

template<class Scalar>
void TimeStepControl<Scalar>::setNumTimeSteps(int numTimeSteps) {
  if (numTimeSteps > 0) {
    tscPL_->set<int>        ("Number of Time Steps", numTimeSteps);
    const int initIndex = getInitIndex();
    tscPL_->set<int>        ("Final Time Index", initIndex + numTimeSteps);
    const double initTime = tscPL_->get<double>("Initial Time");
    const double finalTime = tscPL_->get<double>("Final Time");
    const double initTimeStep = (finalTime - initTime)/numTimeSteps;
    tscPL_->set<double>     ("Initial Time Step", initTimeStep);
    tscPL_->set<std::string>("Integrator Step Type", "Constant");

    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"setParameterList");
    *out << "Warning - Found 'Number of Time Steps' = " << getNumTimeSteps()
         << "  Set the following parameters: \n"
         << "  'Final Time Index'     = " << getFinalIndex() << "\n"
         << "  'Initial Time Step'    = " << getInitTimeStep() << "\n"
         << "  'Integrator Step Type' = " << getStepType() << std::endl;
  }
}


template<class Scalar>
std::string TimeStepControl<Scalar>::description() const
{
  std::string name = "Tempus::TimeStepControl";
  return(name);
}


template<class Scalar>
void TimeStepControl<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  if (verbLevel == Teuchos::VERB_EXTREME) {
    out << description() << "::describe:" << std::endl
        << "pList        = " << tscPL_    << std::endl;
  }
}


template <class Scalar>
void TimeStepControl<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pList));
  pList->validateParameters(*this->getValidParameters());
  pList->validateParametersAndSetDefaults(*this->getValidParameters());
  tscPL_ = pList;

  // Override parameters
  setNumTimeSteps(getNumTimeSteps());

  TEUCHOS_TEST_FOR_EXCEPTION(
    (getInitTime() > getFinalTime() ), std::logic_error,
    "Error - Inconsistent time range.\n"
    "    (timeMin = "<<getInitTime()<<") > (timeMax = "<<getFinalTime()<<")\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    (getMinTimeStep() < Teuchos::ScalarTraits<Scalar>::zero() ),
    std::logic_error,
    "Error - Negative minimum time step.  dtMin = "<<getMinTimeStep()<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (getMaxTimeStep() < Teuchos::ScalarTraits<Scalar>::zero() ),
    std::logic_error,
    "Error - Negative maximum time step.  dtMax = "<<getMaxTimeStep()<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (getMinTimeStep() > getMaxTimeStep() ), std::logic_error,
    "Error - Inconsistent time step range.\n"
    "  (dtMin = "<<getMinTimeStep()<<") > (dtMax = "<<getMaxTimeStep()<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (getInitTimeStep() < Teuchos::ScalarTraits<Scalar>::zero() ),
    std::logic_error,
    "Error - Negative initial time step.  dtInit = "<<getInitTimeStep()<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (getInitTimeStep() < getMinTimeStep() ||
     getInitTimeStep() > getMaxTimeStep() ),
    std::out_of_range,
    "Error - Initial time step is out of range.\n"
    << "    [dtMin, dtMax] = [" << getMinTimeStep() << ", "
                                << getMaxTimeStep() << "]\n"
    << "    dtInit = " << getInitTimeStep() << "\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    (getInitIndex() > getFinalIndex() ), std::logic_error,
    "Error - Inconsistent time index range.\n"
    "  (iStepMin = "<<getInitIndex()<<") > (iStepMax = "
    <<getFinalIndex()<<")\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    (getMaxAbsError() < Teuchos::ScalarTraits<Scalar>::zero() ),
    std::logic_error,
    "Error - Negative maximum time step.  errorMaxAbs = "
    <<getMaxAbsError()<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (getMaxRelError() < Teuchos::ScalarTraits<Scalar>::zero() ),
    std::logic_error,
    "Error - Negative maximum time step.  errorMaxRel = "
    <<getMaxRelError()<<")\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    (getMinOrder() < Teuchos::ScalarTraits<Scalar>::zero() ),
    std::logic_error,
    "Error - Negative minimum order.  orderMin = "<<getMinOrder()<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (getMaxOrder() < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative maximum order.  orderMax = "<<getMaxOrder()<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (getMinOrder() > getMaxOrder() ), std::logic_error,
    "Error - Inconsistent order range.\n"
    "    (orderMin = "<<getMinOrder()<<") > (orderMax = "
    <<getMaxOrder()<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (getInitOrder() < getMinOrder() || getInitOrder() > getMaxOrder()),
    std::out_of_range,
    "Error - Initial order is out of range.\n"
    << "    [orderMin, orderMax] = [" << getMinOrder() << ", "
                                      << getMaxOrder() << "]\n"
    << "    order = " << getInitOrder()  << "\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    (getStepType() != "Constant" and getStepType() != "Variable"),
    std::out_of_range,
      "Error - 'Integrator Step Type' does not equal none of these:\n"
    << "  'Constant' - Integrator will take constant time step sizes.\n"
    << "  'Variable' - Integrator will allow changes to the time step size.\n"
    << "  stepType = " << getStepType()  << "\n");

  // Parse output times
  {
    outputTimes_.clear();
    std::string str = tscPL_->get<std::string>("Output Time List");
    std::string delimiters(",");
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
      std::string token = str.substr(lastPos,pos-lastPos);
      outputTimes_.push_back(Scalar(std::stod(token)));
      if(pos==std::string::npos)
        break;

      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }

    Scalar outputTimeInterval = tscPL_->get<double>("Output Time Interval");
    Scalar output_t = getInitTime();
    while (output_t <= getFinalTime()) {
      outputTimes_.push_back(output_t);
      output_t += outputTimeInterval;
    }

    // order output times
    std::sort(outputTimes_.begin(),outputTimes_.end());
  }

  // Parse output indices
  {
    outputIndices_.clear();
    std::string str = tscPL_->get<std::string>("Output Index List");
    std::string delimiters(",");
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
      std::string token = str.substr(lastPos,pos-lastPos);
      outputIndices_.push_back(int(std::stoi(token)));
      if(pos==std::string::npos)
        break;

      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }

    Scalar outputIndexInterval = tscPL_->get<int>("Output Index Interval");
    Scalar output_i = getInitIndex();
    while (output_i <= getFinalIndex()) {
      outputIndices_.push_back(output_i);
      output_i += outputIndexInterval;
    }

    // order output indices
    std::sort(outputIndices_.begin(),outputIndices_.end());
  }

  return;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
TimeStepControl<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

  const double stdMin = std::numeric_limits<double>::epsilon();
  const double stdMax = std::numeric_limits<double>::max();
  pl->set<double>("Initial Time"      , 0.0    , "Initial time");
  pl->set<double>("Final Time"        , stdMax , "Final time");
  pl->set<int>   ("Initial Time Index", 0      , "Initial time index");
  pl->set<int>   ("Final Time Index"  , 1000000, "Final time index");
  pl->set<double>("Minimum Time Step" , stdMin , "Minimum time step size");
  pl->set<double>("Initial Time Step" , stdMin , "Initial time step size");
  pl->set<double>("Maximum Time Step" , stdMax , "Maximum time step size");
  pl->set<int>   ("Minimum Order", 0,
    "Minimum time-integration order.  If set to zero (default), the\n"
    "Stepper minimum order is used.");
  pl->set<int>   ("Initial Order", 0,
    "Initial time-integration order.  If set to zero (default), the\n"
    "Stepper minimum order is used.");
  pl->set<int>   ("Maximum Order", 0,
    "Maximum time-integration order.  If set to zero (default), the\n"
    "Stepper maximum order is used.");
  pl->set<double>("Maximum Absolute Error", 1.0e-08, "Maximum absolute error");
  pl->set<double>("Maximum Relative Error", 1.0e-08, "Maximum relative error");

  pl->set<std::string>("Integrator Step Type", "Variable",
    "'Integrator Step Type' indicates whether the Integrator will allow "
    "the time step to be modified.\n"
    "  'Constant' - Integrator will take constant time step sizes.\n"
    "  'Variable' - Integrator will allow changes to the time step size.\n");

  pl->set<std::string>("Output Time List", "",
    "Comma deliminated list of output times");
  pl->set<std::string>("Output Index List","",
    "Comma deliminated list of output indices");
  pl->set<double>("Output Time Interval", stdMax, "Output time interval");
  pl->set<int>   ("Output Index Interval", 1000000, "Output index interval");

  pl->set<int>   ("Maximum Number of Stepper Failures", 10,
    "Maximum number of Stepper failures");
  pl->set<int>   ("Maximum Number of Consecutive Stepper Failures", 5,
    "Maximum number of consecutive Stepper failures");
  pl->set<int>   ("Number of Time Steps", -1,
    "The number of constant time steps.  The actual step size gets computed\n"
    "on the fly given the size of the time domain.  Overides and resets\n"
    "  'Final Time Index'     = 'Initial Time Index' + 'Number of Time Steps'\n"
    "  'Initial Time Step'    = "
    "('Final Time' - 'Initial Time')/'Number of Time Steps'\n"
    "  'Integrator Step Type' = 'Constant'\n");

  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
TimeStepControl<Scalar>::getNonconstParameterList()
{
  return(tscPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
TimeStepControl<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = tscPL_;
  tscPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_TimeStepControl_impl_hpp
