// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeStepControl_impl_hpp
#define Tempus_TimeStepControl_impl_hpp

#include "Teuchos_TimeMonitor.hpp"

#include "Tempus_TimeStepControlStrategyConstant.hpp"
#include "Tempus_TimeStepControlStrategyComposite.hpp"
#include "Tempus_TimeStepControlStrategyBasicVS.hpp"
#include "Tempus_TimeStepControlStrategyIntegralController.hpp"


namespace Tempus {

template<class Scalar>
TimeStepControl<Scalar>::TimeStepControl()
  : isInitialized_      (false),
    initTime_           (0.0),
    finalTime_          (1.0e+99),
    minTimeStep_        (0.0),
    initTimeStep_       (1.0),
    maxTimeStep_        (1.0e+99),
    initIndex_          (0),
    finalIndex_         (1000000),
    maxAbsError_        (1.0e-08),
    maxRelError_        (1.0e-08),
    maxFailures_        (10),
    maxConsecFailures_  (5),
    numTimeSteps_       (-1),
    printDtChanges_     (true),
    outputExactly_      (true),
    //outputIndices_      (),
    //outputTimes_        (),
    outputIndexInterval_(1000000),
    outputTimeInterval_ (1.0e+99),
    outputAdjustedDt_(false),
    dtAfterOutput_(0.0)
{
  setTimeStepControlStrategy();
  this->initialize();
}


template<class Scalar>
TimeStepControl<Scalar>::TimeStepControl(
  Scalar              initTime,
  Scalar              finalTime,
  Scalar              minTimeStep,
  Scalar              initTimeStep,
  Scalar              maxTimeStep,
  int                 initIndex,
  int                 finalIndex,
  Scalar              maxAbsError,
  Scalar              maxRelError,
  int                 maxFailures,
  int                 maxConsecFailures,
  int                 numTimeSteps,
  bool                printDtChanges,
  bool                outputExactly,
  std::vector<int>    outputIndices,
  std::vector<Scalar> outputTimes,
  int                 outputIndexInterval,
  Scalar              outputTimeInterval,
  Teuchos::RCP<TimeStepControlStrategy<Scalar>> stepControlStrategy)
  : isInitialized_      (false              ),
    initTime_           (initTime           ),
    finalTime_          (finalTime          ),
    minTimeStep_        (minTimeStep        ),
    initTimeStep_       (initTimeStep       ),
    maxTimeStep_        (maxTimeStep        ),
    initIndex_          (initIndex          ),
    finalIndex_         (finalIndex         ),
    maxAbsError_        (maxAbsError        ),
    maxRelError_        (maxRelError        ),
    maxFailures_        (maxFailures        ),
    maxConsecFailures_  (maxConsecFailures  ),
    numTimeSteps_       (numTimeSteps       ),
    printDtChanges_     (printDtChanges     ),
    outputExactly_      (outputExactly      ),
    outputIndices_      (outputIndices      ),
    outputTimes_        (outputTimes        ),
    outputIndexInterval_(outputIndexInterval),
    outputTimeInterval_ (outputTimeInterval ),
    outputAdjustedDt_   (false              ),
    dtAfterOutput_      (0.0                ),
    stepControlStrategy_(stepControlStrategy)
{
  setNumTimeSteps(getNumTimeSteps());
  this->initialize();
}


template<class Scalar>
void TimeStepControl<Scalar>::initialize() const
{
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
    (getStepType() != "Constant" && getStepType() != "Variable"),
    std::out_of_range,
      "Error - 'Step Type' does not equal one of these:\n"
    << "  'Constant' - Integrator will take constant time step sizes.\n"
    << "  'Variable' - Integrator will allow changes to the time step size.\n"
    << "  stepType = " << getStepType() << "\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    (stepControlStrategy_ == Teuchos::null), std::logic_error,
    "Error - Strategy is unset!\n");

  stepControlStrategy_->initialize();

  isInitialized_ = true;   // Only place where this is set to true!
}


template<class Scalar>
void TimeStepControl<Scalar>::
printDtChanges(int istep, Scalar dt_old, Scalar dt_new, std::string reason) const
{
  if (!getPrintDtChanges()) return;

  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  out->setOutputToRootOnly(0);
  Teuchos::OSTab ostab(out,0,"printDtChanges");

  std::stringstream message;
  message << std::scientific
                   <<std::setw(6)<<std::setprecision(3)<<istep
    << " *  (dt = "<<std::setw(9)<<std::setprecision(3)<<dt_old
    <<   ", new = "<<std::setw(9)<<std::setprecision(3)<<dt_new
    << ")  " << reason << std::endl;
  *out << message.str();
}


template<class Scalar>
void TimeStepControl<Scalar>::checkInitialized()
{
  if ( !isInitialized_ ) {
    this->describe( *(this->getOStream()), Teuchos::VERB_MEDIUM);
    TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error - " << this->description() << " is not initialized!");
  }
}


template<class Scalar>
void TimeStepControl<Scalar>::setNextTimeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> > & solutionHistory,
  Status & integratorStatus)
{
  using Teuchos::RCP;

  checkInitialized();

  TEMPUS_FUNC_TIME_MONITOR("Tempus::TimeStepControl::setNextTimeStep()");
  {
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar lastTime = solutionHistory->getCurrentState()->getTime();
    const int iStep = workingState->getIndex();
    Scalar dt   = workingState->getTimeStep();
    Scalar time = workingState->getTime();
    bool output = false;

    RCP<StepperState<Scalar> > stepperState = workingState->getStepperState();

    // If last time step was adjusted for output, reinstate previous dt.
    if (getStepType() == "Variable") {
      if (outputAdjustedDt_ == true) {
        printDtChanges(iStep, dt, dtAfterOutput_, "Reset dt after output.");
        dt = dtAfterOutput_;
        time = lastTime + dt;
        outputAdjustedDt_ = false;
        dtAfterOutput_ = 0.0;
      }

      if (dt <= 0.0) {
        printDtChanges(iStep, dt, getInitTimeStep(), "Reset dt to initial dt.");
        dt = getInitTimeStep();
        time = lastTime + dt;
      }

      if (dt < getMinTimeStep()) {
        printDtChanges(iStep, dt, getMinTimeStep(), "Reset dt to minimum dt.");
        dt = getMinTimeStep();
        time = lastTime + dt;
      }
    }

    // Update dt for the step control strategy to be informed
    workingState->setTimeStep(dt);
    workingState->setTime(time);

    // Call the step control strategy (to update dt if needed)
    stepControlStrategy_->setNextTimeStep(*this, solutionHistory,
                                          integratorStatus);

    // Get the dt (probably have changed by stepControlStrategy_)
    dt = workingState->getTimeStep();
    time = workingState->getTime();

    if (getStepType() == "Variable") {
      if (dt < getMinTimeStep()) { // decreased below minimum dt
        printDtChanges(iStep, dt, getMinTimeStep(),
          "dt is too small.  Resetting to minimum dt.");
        dt = getMinTimeStep();
        time = lastTime + dt;
      }
      if (dt > getMaxTimeStep()) { // increased above maximum dt
        printDtChanges(iStep, dt, getMaxTimeStep(),
          "dt is too large.  Resetting to maximum dt.");
        dt = getMaxTimeStep();
        time = lastTime + dt;
      }
    }


    // Check if we need to output this step index
    std::vector<int>::const_iterator it =
      std::find(getOutputIndices().begin(), getOutputIndices().end(), iStep);
    if (it != getOutputIndices().end()) output = true;

    const int iInterval = getOutputIndexInterval();
    if ( (iStep - getInitIndex()) % iInterval == 0) output = true;

    // Check if we need to output in the next timestep based on
    // outputTimes_ or "Output Time Interval".
    Scalar reltol = 1.0e-6;
    Scalar endTime = lastTime+dt+getMinTimeStep();
    // getMinTimeStep() = dt for constant time step
    // so we can't add it on here
    if (getStepType() == "Constant") endTime = lastTime+dt;
    bool checkOutput = false;
    Scalar oTime = getInitTime();
    for (size_t i=0; i < outputTimes_.size(); ++i) {
      oTime = outputTimes_[i];
      if (lastTime < oTime && oTime <= endTime) {
        checkOutput = true;
        break;
      }
    }
    const Scalar tInterval = getOutputTimeInterval();
    Scalar oTime2 =  ceil((lastTime-getInitTime())/tInterval)*tInterval
                   + getInitTime();
    if (lastTime < oTime2 && oTime2 <= endTime) {
      if (checkOutput == true) {
        if (oTime2 < oTime) oTime = oTime2;  // Use the first output time.
      } else {
        checkOutput = true;
        oTime = oTime2;
      }
    }

    if (checkOutput == true) {
      const bool outputExactly = getOutputExactly();
      if (getStepType() == "Variable" && outputExactly == true) {
        // Adjust time step to hit output times.
        if ( time > oTime ) {
          output = true;
          printDtChanges(iStep, dt, oTime - lastTime,
            "Adjusting dt to hit the next output time.");
          // Next output time is not near next time
          // (>getMinTimeStep() away from it).
          // Take time step to hit output time.
          outputAdjustedDt_ = true;
          dtAfterOutput_ = dt;
          dt = oTime - lastTime;
          time = lastTime + dt;
        } else if (std::fabs((time-oTime)/(time)) < reltol) {
          output = true;
          printDtChanges(iStep, dt, oTime - lastTime,
            "Adjusting dt for numerical roundoff to hit the next output time.");
          // Next output time IS VERY near next time (<reltol away from it),
          // e.g., adjust for numerical roundoff.
          outputAdjustedDt_ = true;
          dtAfterOutput_ = dt;
          dt = oTime - lastTime;
          time = lastTime + dt;
        } else  if (std::fabs((time + getMinTimeStep() - oTime)/oTime) < reltol ) {
          printDtChanges(iStep, dt, (oTime - lastTime)/2.0,
            "The next output time is within the minimum dt of the next time. "
            "Adjusting dt to take two steps.");
          // Next output time IS near next time
          // (<getMinTimeStep() away from it).
          // Take two time steps to get to next output time.
          dt = (oTime - lastTime)/2.0;
          time = lastTime + dt;
        }
      } else {
        // Stepping over output time and want this time step for output,
        // but do not want to change dt. Either because of 'Constant' time
        // step or user specification, "Output Exactly On Output Times"=false.
        output = true;
      }
    }

    // Adjust time step to hit final time or correct for small
    // numerical differences.
    if (getStepType() == "Variable") {
      if ((lastTime + dt > getFinalTime() ) ||
          (std::fabs((lastTime+dt-getFinalTime())/(lastTime+dt)) < reltol)) {
        printDtChanges(iStep, dt, getFinalTime() - lastTime,
          "Adjusting dt to hit final time.");
        dt = getFinalTime() - lastTime;
        time = lastTime + dt;
      }
    }

    // Check for negative time step.
    TEUCHOS_TEST_FOR_EXCEPTION( dt <= Scalar(0.0), std::out_of_range,
      "Error - Time step is not positive.  dt = " << dt <<"\n");

    // Time step always needs to keep time within range.
    TEUCHOS_TEST_FOR_EXCEPTION(
      (lastTime + dt < getInitTime()), std::out_of_range,
      "Error - Time step does not move time INTO time range.\n"
      "    [timeMin, timeMax] = [" << getInitTime() << ", "
      << getFinalTime() << "]\n"
      "    T + dt = " << lastTime <<" + "<< dt <<" = " << lastTime + dt <<"\n");

    if (getStepType() == "Variable") {
      TEUCHOS_TEST_FOR_EXCEPTION(
        (lastTime + dt > getFinalTime()), std::out_of_range,
        "Error - Time step move time OUT OF time range.\n"
        "    [timeMin, timeMax] = [" << getInitTime() << ", "
        << getFinalTime() << "]\n"
        "    T + dt = " << lastTime <<" + "<< dt <<" = " << lastTime + dt <<"\n");
    }

    workingState->setTimeStep(dt);
    workingState->setTime(time);
    workingState->setOutput(output);
  }
  return;
}


/// Test if time is within range: include initTime and exclude finalTime.
template<class Scalar>
bool TimeStepControl<Scalar>::timeInRange(const Scalar time) const
{
  // Get absolute tolerance 1.0e-(i+14), i.e., 14 digits of accuracy.
  const int relTol = 14;
  const int i =
    (getInitTime() == 0) ? 0 : 1 + (int)std::floor(std::log10(std::fabs(getInitTime()) ) );
  const Scalar absTolInit = std::pow(10, i-relTol);
  const int j =
    (getFinalTime() == 0) ? 0 : 1 + (int)std::floor(std::log10(std::fabs(getFinalTime()) ) );
  const Scalar absTolFinal = std::pow(10, j-relTol);

  const bool test1 = getInitTime() - absTolInit <= time;
  const bool test2 = time < getFinalTime() - absTolFinal;

  return (test1 && test2);
}


/// Test if index is within range: include initIndex and exclude finalIndex.
template<class Scalar>
bool TimeStepControl<Scalar>::indexInRange(const int iStep) const{
  bool iir = (getInitIndex() <= iStep && iStep < getFinalIndex());
  return iir;
}


template<class Scalar>
void TimeStepControl<Scalar>::setNumTimeSteps(int numTimeSteps)
{
  TEUCHOS_TEST_FOR_EXCEPTION( getStepType() != "Constant", std::out_of_range,
      "Error - Can only use setNumTimeSteps() when 'Step Type' == 'Constant'.\n");

  if (numTimeSteps >= 0) {
    numTimeSteps_ = numTimeSteps;
    setFinalIndex(getInitIndex() + numTimeSteps_);
    Scalar initTimeStep;
    if (numTimeSteps_ == 0)
      initTimeStep = Scalar(0.0);
    else
      initTimeStep = (getFinalTime() - getInitTime())/numTimeSteps_;
    setInitTimeStep(initTimeStep);
    setMinTimeStep (initTimeStep);
    setMaxTimeStep (initTimeStep);

    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out,1,"setNumTimeSteps");
    *out << "Warning - setNumTimeSteps() Setting 'Number of Time Steps' = " << getNumTimeSteps()
         << "  Set the following parameters: \n"
         << "  'Final Time Index'     = " << getFinalIndex() << "\n"
         << "  'Initial Time Step'    = " << getInitTimeStep() << "\n"
         << "  'Step Type'            = " << getStepType() << std::endl;

    isInitialized_ = false;
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

    std::vector<int> idx = getOutputIndices();
    std::ostringstream listIdx;
    if (!idx.empty()) {
      for(std::size_t i = 0; i < idx.size()-1; ++i) listIdx << idx[i] << ", ";
      listIdx << idx[idx.size()-1];
    }

    std::vector<Scalar> times = getOutputTimes();
    std::ostringstream listTimes;
    if (!times.empty()) {
      for(std::size_t i = 0; i < times.size()-1; ++i) listTimes << times[i] << ", ";
      listTimes << times[times.size()-1];
    }

    auto l_out = Teuchos::fancyOStream( out.getOStream() );
    l_out->setOutputToRootOnly(0);
    *l_out << description() << "::describe:" << std::endl
           << "stepType           = " << getStepType()            << std::endl
           << "initTime           = " << getInitTime()            << std::endl
           << "finalTime          = " << getFinalTime()           << std::endl
           << "minTimeStep        = " << getMinTimeStep()         << std::endl
           << "initTimeStep       = " << getInitTimeStep()        << std::endl
           << "maxTimeStep        = " << getMaxTimeStep()         << std::endl
           << "initIndex          = " << getInitIndex()           << std::endl
           << "finalIndex         = " << getFinalIndex()          << std::endl
           << "maxAbsError        = " << getMaxAbsError()         << std::endl
           << "maxRelError        = " << getMaxRelError()         << std::endl
           << "maxFailures        = " << getMaxFailures()         << std::endl
           << "maxConsecFailures  = " << getMaxConsecFailures()   << std::endl
           << "numTimeSteps       = " << getNumTimeSteps()        << std::endl
           << "printDtChanges     = " << getPrintDtChanges()      << std::endl
           << "outputExactly      = " << getOutputExactly()       << std::endl
           << "outputIndices      = " << listIdx.str()            << std::endl
           << "outputTimes        = " << listTimes.str()          << std::endl
           << "outputIndexInterval= " << getOutputIndexInterval() << std::endl
           << "outputTimeInterval = " << getOutputTimeInterval()  << std::endl
           << "outputAdjustedDt   = " << outputAdjustedDt_        << std::endl
           << "dtAfterOutput      = " << dtAfterOutput_           << std::endl
           << "stepControlSrategy = " << std::endl;
           stepControlStrategy_->describe(out, verbLevel);
  }
}


template<class Scalar>
void TimeStepControl<Scalar>::setTimeStepControlStrategy(
  Teuchos::RCP<TimeStepControlStrategy<Scalar> > tscs)
{
  using Teuchos::rcp;

  if ( tscs != Teuchos::null ) {
    stepControlStrategy_ = tscs;
  } else {
    stepControlStrategy_ =
      rcp(new TimeStepControlStrategyConstant<Scalar>(getInitTimeStep()));
  }
  isInitialized_ = false;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
TimeStepControl<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList("Time Step Control");

  pl->set<double>("Initial Time"          , getInitTime()    , "Initial time");
  pl->set<double>("Final Time"            , getFinalTime()   , "Final time");
  pl->set<double>("Minimum Time Step"     , getMinTimeStep() , "Minimum time step size");
  pl->set<double>("Initial Time Step"     , getInitTimeStep(), "Initial time step size");
  pl->set<double>("Maximum Time Step"     , getMaxTimeStep() , "Maximum time step size");
  pl->set<int>   ("Initial Time Index"    , getInitIndex()   , "Initial time index");
  pl->set<int>   ("Final Time Index"      , getFinalIndex()  , "Final time index");
  pl->set<int>   ("Number of Time Steps", getNumTimeSteps(),
    "The number of constant time steps.  The actual step size gets computed\n"
    "on the fly given the size of the time domain.  Overides and resets\n"
    "  'Final Time Index'     = 'Initial Time Index' + 'Number of Time Steps'\n"
    "  'Initial Time Step'    = "
    "('Final Time' - 'Initial Time')/'Number of Time Steps'\n");
  pl->set<double>("Maximum Absolute Error", getMaxAbsError() , "Maximum absolute error");
  pl->set<double>("Maximum Relative Error", getMaxRelError() , "Maximum relative error");

  pl->set<bool>  ("Print Time Step Changes", getPrintDtChanges(),
    "Print timestep size when it changes");

  pl->set<bool>("Output Exactly On Output Times", getOutputExactly(),
    "This determines if the timestep size will be adjusted to exactly land\n"
    "on the output times for 'Variable' timestepping (default=true).\n"
    "When set to 'false' or for 'Constant' time stepping, the timestep\n"
    "following the output time will be flagged for output.\n");

  pl->set<int>   ("Output Index Interval", getOutputIndexInterval(), "Output index interval");
  pl->set<double>("Output Time Interval", getOutputTimeInterval(), "Output time interval");

  {
    std::vector<int> idx = getOutputIndices();
    std::ostringstream list;
    if (!idx.empty()) {
      for(std::size_t i = 0; i < idx.size()-1; ++i) list << idx[i] << ", ";
      list << idx[idx.size()-1];
    }
    pl->set<std::string>("Output Index List", list.str(),
      "Comma deliminated list of output indices");
  }
  {
    std::vector<Scalar> times = getOutputTimes();
    std::ostringstream list;
    if (!times.empty()) {
      for(std::size_t i = 0; i < times.size()-1; ++i) list << times[i] << ", ";
      list << times[times.size()-1];
    }
    pl->set<std::string>("Output Time List", list.str(),
      "Comma deliminated list of output times");
  }

  pl->set<int>   ("Maximum Number of Stepper Failures", getMaxFailures(),
    "Maximum number of Stepper failures");
  pl->set<int>   ("Maximum Number of Consecutive Stepper Failures", getMaxConsecFailures(),
    "Maximum number of consecutive Stepper failures");

  pl->set("Time Step Control Strategy", *stepControlStrategy_->getValidParameters());

  return pl;
}


// Nonmember constructor - ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<TimeStepControl<Scalar> > createTimeStepControl(
  Teuchos::RCP<Teuchos::ParameterList> const& pList, bool runInitialize)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  auto tsc = Teuchos::rcp(new TimeStepControl<Scalar>());
  if (pList == Teuchos::null) return tsc;

  pList->validateParametersAndSetDefaults(*tsc->getValidParameters(), 0);

  tsc->setInitTime(         pList->get<double>("Initial Time"));
  tsc->setFinalTime(        pList->get<double>("Final Time"));
  tsc->setMinTimeStep(      pList->get<double>("Minimum Time Step"));
  tsc->setInitTimeStep(     pList->get<double>("Initial Time Step"));
  tsc->setMaxTimeStep(      pList->get<double>("Maximum Time Step"));
  tsc->setInitIndex(        pList->get<int>   ("Initial Time Index"));
  tsc->setFinalIndex(       pList->get<int>   ("Final Time Index"));
  tsc->setMaxAbsError(      pList->get<double>("Maximum Absolute Error"));
  tsc->setMaxRelError(      pList->get<double>("Maximum Relative Error"));
  tsc->setMaxFailures(      pList->get<int>   ("Maximum Number of Stepper Failures"));
  tsc->setMaxConsecFailures(pList->get<int>   ("Maximum Number of Consecutive Stepper Failures"));
  tsc->setPrintDtChanges(   pList->get<bool>  ("Print Time Step Changes"));
  tsc->setNumTimeSteps(     pList->get<int>   ("Number of Time Steps"));

  tsc->setOutputExactly(    pList->get<bool>  ("Output Exactly On Output Times"));

  // Parse output indices
  {
    std::vector<int> outputIndices;
    outputIndices.clear();
    std::string str = pList->get<std::string>("Output Index List");
    std::string delimiters(",");
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
      std::string token = str.substr(lastPos,pos-lastPos);
      outputIndices.push_back(int(std::stoi(token)));
      if(pos==std::string::npos) break;

      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }

    // order output indices
    std::sort(outputIndices.begin(),outputIndices.end());
    tsc->setOutputIndices(outputIndices);
  }

  tsc->setOutputIndexInterval(pList->get<int>("Output Index Interval"));

  // Parse output times
  {
    std::vector<Scalar> outputTimes;
    outputTimes.clear();
    std::string str = pList->get<std::string>("Output Time List");
    std::string delimiters(",");
    // Skip delimiters at the beginning
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find the first delimiter
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
      // Found a token, add it to the vector
      std::string token = str.substr(lastPos,pos-lastPos);
      outputTimes.push_back(Scalar(std::stod(token)));
      if(pos==std::string::npos) break;

      lastPos = str.find_first_not_of(delimiters, pos); // Skip delimiters
      pos = str.find_first_of(delimiters, lastPos);     // Find next delimiter
    }

    // order output times
    std::sort(outputTimes.begin(),outputTimes.end());
    outputTimes.erase(std::unique(outputTimes.begin(),
                                  outputTimes.end()   ),
                                  outputTimes.end()     );
    tsc->setOutputTimes(outputTimes);
  }

  tsc->setOutputTimeInterval(pList->get<double>("Output Time Interval"));


  if ( !pList->isParameter("Time Step Control Strategy") ) {

    tsc->setTimeStepControlStrategy();  // i.e, set default Constant timestep strategy.

  } else {

    RCP<ParameterList> tscsPL =
      Teuchos::sublist(pList, "Time Step Control Strategy", true);

    auto strategyType = tscsPL->get<std::string>("Strategy Type");
    if (strategyType == "Constant") {
      tsc->setTimeStepControlStrategy(
        createTimeStepControlStrategyConstant<Scalar>(tscsPL));
    } else if (strategyType == "Basic VS") {
      tsc->setTimeStepControlStrategy(
        createTimeStepControlStrategyBasicVS<Scalar>(tscsPL));
    } else if (strategyType == "Integral Controller") {
      tsc->setTimeStepControlStrategy(
        createTimeStepControlStrategyIntegralController<Scalar>(tscsPL));
    } else if (strategyType == "Composite") {
      tsc->setTimeStepControlStrategy(
        createTimeStepControlStrategyComposite<Scalar>(tscsPL));
    } else {
      RCP<Teuchos::FancyOStream> out =
        Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setOutputToRootOnly(0);
      Teuchos::OSTab ostab(out,1, "createTimeStepControl()");
      *out << "Warning -- Did not find a Tempus strategy to create!\n"
           << "'Strategy Type' = '" << strategyType << "'\n"
           << "Should call setTimeStepControlStrategy() with this\n"
           << "(app-specific?) strategy, and initialize().\n" << std::endl;
    }
  }

  if (runInitialize) tsc->initialize();
  return tsc;
}


} // namespace Tempus
#endif // Tempus_TimeStepControl_impl_hpp
