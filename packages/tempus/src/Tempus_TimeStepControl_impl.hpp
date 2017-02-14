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
  pList_->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setParameterList(pList_);
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
  : timeMin_                (tsc_.timeMin_                ),
    timeMax_                (tsc_.timeMax_                ),
    dtMin_                  (tsc_.dtMin_                  ),
    dtInit_                 (tsc_.dtInit_                 ),
    dtMax_                  (tsc_.dtMax_                  ),
    iStepMin_               (tsc_.iStepMin_               ),
    iStepMax_               (tsc_.iStepMax_               ),
    errorMaxAbs_            (tsc_.errorMaxAbs_            ),
    errorMaxRel_            (tsc_.errorMaxRel_            ),
    orderMin_               (tsc_.orderMin_               ),
    orderInit_              (tsc_.orderInit_              ),
    orderMax_               (tsc_.orderMax_               ),
    stepType_               (tsc_.stepType_               ),
    outputIndices_          (tsc_.outputIndices_          ),
    outputTimes_            (tsc_.outputTimes_            ),
    nFailuresMax_           (tsc_.nFailuresMax_           ),
    nConsecutiveFailuresMax_(tsc_.nConsecutiveFailuresMax_),
    pList_                  (tsc_.pList_                  ),
    outputAdjustedDt_       (tsc_.outputAdjustedDt_       ),
    dtAfterOutput_          (tsc_.dtAfterOutput_          )
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
    RCP<SolutionStateMetaData<Scalar> > metaData_ = workingState->metaData_;
    const Scalar time = metaData_->getTime();
    const int iStep = metaData_->getIStep();
    const Scalar errorAbs = metaData_->getErrorAbs();
    const Scalar errorRel = metaData_->getErrorRel();
    int order = metaData_->getOrder();
    Scalar dt = metaData_->getDt();
    bool output = metaData_->getOutput();

    RCP<StepperState<Scalar> > stepperState = workingState->stepperState_;

    output = false;

    // If last time step was adjusted for output, reinstate previous dt.
    if (outputAdjustedDt_ == true) {
      dt = dtAfterOutput_;
      outputAdjustedDt_ = false;
      dtAfterOutput_ = 0.0;
    }

    if (dt < dtMin_) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"getNextTimeStep");
      *out << "Warning - Time step size (=" << dt << ") is less than\n"
           << "  minimum time step size (=" << dtMin_ << ").\n"
           << "  Resetting to minimum time step size." << std::endl;
      dt = dtMin_;
    }

    if (stepType_ == "Constant") {

      dt = dtInit_;

      // Stepper failure
      if (stepperState->stepperStatus_ == Status::FAILED) {
        if (order+1 <= orderMax_) {
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
      if (errorAbs > errorMaxAbs_) {
        if (order+1 <= orderMax_) {
          order++;
          RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out,1,"getNextTimeStep");
          *out
            <<"Warning - Absolute error is too large with constant time step.\n"
            <<"  (errorAbs ="<<errorAbs<<") > (errorMaxAbs ="<<errorMaxAbs_<<")"
            <<"  Try increasing order.  order = " << order << std::endl;
        } else {
          RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out,1,"getNextTimeStep");
          *out
            <<"Failure - Absolute error failed and can not change time step "
            <<"size or order!\n"
            <<"  Time step type == CONSTANT_STEP_SIZE\n"
            <<"  order = " << order
            <<"  (errorAbs ="<<errorAbs<<") > (errorMaxAbs ="<<errorMaxAbs_<<")"
            << std::endl;
          integratorStatus = FAILED;
          return;
        }
      }

      // Relative error failure
      if (errorRel > errorMaxRel_) {
        if (order+1 <= orderMax_) {
          order++;
          RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out,1,"getNextTimeStep");
          *out
            <<"Warning - Relative error is too large with constant time step.\n"
            <<"  (errorRel ="<<errorRel<<") > (errorMaxRel ="<<errorMaxRel_<<")"
            <<"  Try increasing order.  order = " << order << std::endl;
        } else {
          RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out,1,"getNextTimeStep");
          *out
            <<"Failure - Relative error failed and can not change time step "
            <<"size or order!\n"
            <<"  Time step type == CONSTANT_STEP_SIZE\n"
            <<"  order = " << order
            <<"  (errorRel ="<<errorRel<<") > (errorMaxRel ="<<errorMaxRel_<<")"
            << std::endl;
          integratorStatus = FAILED;
          return;
        }
      }

      // Consistency checks
      TEUCHOS_TEST_FOR_EXCEPTION(
        (order < orderMin_ || order > orderMax_), std::out_of_range,
        "Error - Solution order is out of range and can not change "
        "time step size!\n"
        "    Time step type == CONSTANT_STEP_SIZE\n"
        "    [order_min, order_max] = [" <<orderMin_<< ", " <<orderMax_<< "]\n"
        "    order = " << order << "\n");

    } else { // VARIABLE_STEP_SIZE

      // \todo The following controls should be generalized to plugable options.
      if (stepperState->stepperStatus_ == Status::FAILED) dt /=2;
      if (errorAbs > errorMaxAbs_) dt /= 2;
      if (errorRel > errorMaxRel_) dt /= 2;
      if (order < orderMin_) dt *= 2;
      if (order > orderMax_) dt /= 2;

      if (dt < dtMin_) dt = dtMin_;
      if (dt > dtMax_) dt = dtMax_;
    }

    // Adjust time step to hit final time or correct for small
    // numerical differences.
    Scalar reltol = 1.0e-6;
    if ((time + dt > timeMax_ ) ||
        (std::abs((time+dt-timeMax_)/(time+dt)) < reltol))
      dt = timeMax_ - time;

    // Check if we need to output this step index
    std::vector<int>::const_iterator it =
      std::find(outputIndices_.begin(), outputIndices_.end(), iStep+1);
    if (it != outputIndices_.end()) output = true;

    // Adjust time step to hit output times.
    for (size_t i=0; i < outputTimes_.size(); ++i) {
      const Scalar oTime = outputTimes_[i];
      if (time < oTime && oTime <= time+dt+dtMin_) {
        output = true;
        outputAdjustedDt_ = true;
        dtAfterOutput_ = dt;
        if (time < oTime && oTime <= time+dt-dtMin_) {
          // Next output time is not near next time (>dtMin_ away from it).
          // Take time step to hit output time.
          dt = oTime - time;
        } else if (std::abs((time+dt-oTime)/(time+dt)) < reltol) {
          // Next output time IS VERY near next time (<reltol away from it),
          // e.g., adjust for numerical roundoff.
          dt = oTime - time;
        } else if (time+dt-dtMin_ < oTime && oTime <= time+dt+dtMin_) {
          // Next output time IS near next time (<dtMin_ away from it).
          // Take two time steps to get to next output time.
          dt = (oTime - time)/2.0;
        }
        break;
      }
    }

    // Time step always needs to keep time within range.
    TEUCHOS_TEST_FOR_EXCEPTION(
      (time + dt < timeMin_), std::out_of_range,
      "Error - Time step does not move time INTO time range.\n"
      "    [timeMin, timeMax] = [" << timeMin_ << ", " << timeMax_ << "]\n"
      "    T + dt = " << time <<" + "<< dt <<" = " << time + dt << "\n");

    TEUCHOS_TEST_FOR_EXCEPTION(
      (time + dt > timeMax_), std::out_of_range,
      "Error - Time step move time OUT OF time range.\n"
      "    [timeMin, timeMax] = [" << timeMin_ << ", " << timeMax_ << "]\n"
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
  bool tir = (timeMin_*(1.0-relTol) <= time and time < timeMax_*(1.0-relTol));
  return tir;
}


template<class Scalar>
bool TimeStepControl<Scalar>::indexInRange(const int iStep) const{
  bool iir = (iStepMin_ <= iStep and iStep < iStepMax_);
  return iir;
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
        << "timeMin      = " << timeMin_      << std::endl
        << "timeMax      = " << timeMax_      << std::endl
        << "dtMin        = " << dtMin_        << std::endl
        << "dtInit       = " << dtInit_       << std::endl
        << "dtMax        = " << dtMax_        << std::endl
        << "iStepMin     = " << iStepMin_     << std::endl
        << "iStepMax     = " << iStepMax_     << std::endl
        << "orderMin     = " << orderMin_     << std::endl
        << "orderInit    = " << orderInit_    << std::endl
        << "orderMax     = " << orderMax_     << std::endl
        << "errorMaxAbs  = " << errorMaxAbs_  << std::endl
        << "errorMaxRel  = " << errorMaxRel_  << std::endl
        << "stepType     = " << stepType_     << std::endl
        << "nFailuresMax = " << nFailuresMax_ << std::endl
        << "nConsecutiveFailuresMax = " << nConsecutiveFailuresMax_ << std::endl
        << "pList        = " << pList_        << std::endl;
  }
}


template <class Scalar>
void TimeStepControl<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pList));
  pList->validateParameters(*this->getValidParameters());
  pList->validateParametersAndSetDefaults(*this->getValidParameters());
  pList_ = pList;

  Teuchos::readVerboseObjectSublist(&*pList_,this);

  timeMin_     = pList_->get<double>("Initial Time");
  timeMax_     = pList_->get<double>("Final Time");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (timeMin_ > timeMax_ ), std::logic_error,
    "Error - Inconsistent time range.\n"
    "    (timeMin = "<<timeMin_<<") > (timeMax = "<<timeMax_<<")\n");

  dtMin_       = pList_->get<double>("Minimum Time Step");
  dtInit_      = pList_->get<double>("Initial Time Step");
  dtMax_       = pList_->get<double>("Maximum Time Step");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dtMin_ < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative minimum time step.  dtMin = "<<dtMin_<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dtMax_ < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative maximum time step.  dtMax = "<<dtMax_<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dtMin_ > dtMax_ ), std::logic_error,
    "Error - Inconsistent time step range.\n"
    "    (dtMin = "<<dtMin_<<") > (dtMax = "<<dtMax_<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dtInit_ < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative initial time step.  dtInit = "<<dtInit_<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dtInit_ < dtMin_ || dtInit_ > dtMax_ ),
    std::out_of_range,
    "Error - Initial time step is out of range.\n"
    << "    [dtMin, dtMax] = [" << dtMin_ << ", " << dtMax_ << "]\n"
    << "    dtInit = " << dtInit_ << "\n");

  iStepMin_    = pList_->get<int>("Initial Time Index");
  iStepMax_    = pList_->get<int>("Final Time Index");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (iStepMin_ > iStepMax_ ), std::logic_error,
    "Error - Inconsistent time index range.\n"
    "    (iStepMin = "<<iStepMin_<<") > (iStepMax = "<<iStepMax_<<")\n");

  errorMaxAbs_ = pList_->get<double>("Maximum Absolute Error");
  errorMaxRel_ = pList_->get<double>("Maximum Relative Error");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (errorMaxAbs_ < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative maximum time step.  errorMaxAbs = "<<errorMaxAbs_<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (errorMaxRel_ < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative maximum time step.  errorMaxRel = "<<errorMaxRel_<<")\n");

  orderMin_    = pList_->get<int>("Minimum Order");
  orderInit_   = pList_->get<int>("Initial Order");
  orderMax_    = pList_->get<int>("Maximum Order");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (orderMin_ < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative minimum order.  orderMin = "<<orderMin_<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (orderMax_ < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative maximum order.  orderMax = "<<orderMax_<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (orderMin_ > orderMax_ ), std::logic_error,
    "Error - Inconsistent order range.\n"
    "    (orderMin = "<<orderMin_<<") > (orderMax = "<<orderMax_<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (orderInit_ < orderMin_ || orderInit_ > orderMax_), std::out_of_range,
    "Error - Initial order is out of range.\n"
    << "    [orderMin, orderMax] = [" << orderMin_ << ", " << orderMax_ << "]\n"
    << "    order = " << orderInit_  << "\n");

  stepType_ = pList_->get<std::string>("Integrator Step Type", "Variable");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (stepType_ != "Constant" and stepType_ != "Variable"), std::out_of_range,
      "Error - 'Integrator Step Type' does not equal none of these:\n"
    << "  'Constant' - Integrator will take constant time step sizes.\n"
    << "  'Variable' - Integrator will allow changes to the time step size.\n"
    << "  stepType = " << stepType_  << "\n");

  // Parse output times
  {
    outputTimes_.clear();
    std::string str = pList_->get<std::string>("Output Time List");
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

    Scalar outputTimeInterval = pList_->get<double>("Output Time Interval");
    Scalar output_t = timeMin_;
    while (output_t <= timeMax_) {
      outputTimes_.push_back(output_t);
      output_t += outputTimeInterval;
    }

    // order output times
    std::sort(outputTimes_.begin(),outputTimes_.end());
  }

  // Parse output indices
  {
    outputIndices_.clear();
    std::string str = pList_->get<std::string>("Output Index List");
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

    Scalar outputIndexInterval = pList_->get<int>("Output Index Interval");
    Scalar output_i = iStepMin_;
    while (output_i <= iStepMax_) {
      outputIndices_.push_back(output_i);
      output_i += outputIndexInterval;
    }

    // order output indices
    std::sort(outputIndices_.begin(),outputIndices_.end());
  }

  nFailuresMax_ = pList_->get<int>("Maximum Number of Stepper Failures");
  nConsecutiveFailuresMax_ = pList_->get<int>(
    "Maximum Number of Consecutive Stepper Failures");
  return;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
TimeStepControl<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    const double stdMin = std::numeric_limits<double>::epsilon();
    const double stdMax = std::numeric_limits<double>::max();
    pl->set("Initial Time"      , 0.0    , "Initial time");
    pl->set("Final Time"        , stdMax , "Final time");
    pl->set("Initial Time Index", 0      , "Initial time index");
    pl->set("Final Time Index"  , 1000000, "Final time index");
    pl->set("Minimum Time Step" , stdMin , "Minimum time step size");
    pl->set("Initial Time Step" , stdMin , "Initial time step size");
    pl->set("Maximum Time Step" , stdMax , "Maximum time step size");
    pl->set("Maximum Absolute Error", 1.0e-08, "Maximum absolute error");
    pl->set("Maximum Relative Error", 1.0e-08, "Maximum relative error");
    pl->set("Minimum Order", 0,
      "Minimum time-integration order.  If set to zero (default), the\n"
      "Stepper minimum order is used.");
    pl->set("Initial Order", 0,
      "Initial time-integration order.  If set to zero (default), the\n"
      "Stepper minimum order is used.");
    pl->set("Maximum Order", 0,
      "Maximum time-integration order.  If set to zero (default), the\n"
      "Stepper maximum order is used.");

    pl->set("Integrator Step Type", "Variable",
      "'Integrator Step Type' indicates whether the Integrator will allow "
      "the time step to be modified.\n"
      "  'Constant' - Integrator will take constant time step sizes.\n"
      "  'Variable' - Integrator will allow changes to the time step size.\n");

    pl->set("Output Time List", "", "Comma deliminated list of output times");
    pl->set("Output Index List","", "Comma deliminated list of output indices");
    pl->set("Output Time Interval", stdMax, "Output time interval");
    pl->set("Output Index Interval", 1000000, "Output index interval");

    pl->set("Maximum Number of Stepper Failures", 10,
      "Maximum number of Stepper failures");
    pl->set("Maximum Number of Consecutive Stepper Failures", 5,
      "Maximum number of consecutive Stepper failures");

    validPL = pl;

  }
  return validPL;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
TimeStepControl<Scalar>::getNonconstParameterList()
{
  return(pList_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
TimeStepControl<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = pList_;
  pList_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_TimeStepControl_impl_hpp
