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

//Step control strategy
#include "Tempus_TimeStepControlStrategyConstant.hpp"
#include "Tempus_TimeStepControlStrategyComposite.hpp"
#include "Tempus_TimeStepControlStrategyBasicVS.hpp"
#include "Tempus_TimeStepControlStrategyIntegralController.hpp"

//Thyra
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

template<class Scalar>
TimeStepControl<Scalar>::TimeStepControl(
  Teuchos::RCP<Teuchos::ParameterList> pList)
  : outputAdjustedDt_(false), dtAfterOutput_(0.0)
{
  this->setParameterList(pList);
}

template<class Scalar>
TimeStepControl<Scalar>::TimeStepControl(const TimeStepControl<Scalar>& tsc_)
  : tscPL_              (tsc_.tscPL_              ),
    outputIndices_      (tsc_.outputIndices_      ),
    outputTimes_        (tsc_.outputTimes_        ),
    outputAdjustedDt_   (tsc_.outputAdjustedDt_   ),
    dtAfterOutput_      (tsc_.dtAfterOutput_      ),
    stepControlStategy_ (tsc_.stepControlStategy_ )
{}


template<class Scalar>
void TimeStepControl<Scalar>::getNextTimeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> > & solutionHistory,
  Status & integratorStatus)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::TimeStepControl::getNextTimeStep()");
  {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"getNextTimeStep");
    bool printChanges = solutionHistory->getVerbLevel() !=
                        Teuchos::as<int>(Teuchos::VERB_NONE);
    bool printChangesHi = solutionHistory->getVerbLevel() >=
                          Teuchos::as<int>(Teuchos::VERB_HIGH);

    auto changeDT = [] (Scalar dt_old, Scalar dt_new, std::string reason) {
      std::stringstream message;
      message <<
    "     (dt = "<<std::scientific<<std::setw(9)<<std::setprecision(3)<<dt_old
    << ", new = "<<std::scientific<<std::setw(9)<<std::setprecision(3)<<dt_new
    << ")  " << reason << std::endl;
      return message.str();
    };

    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionStateMetaData<Scalar> > metaData = workingState->getMetaData();
    const Scalar lastTime = solutionHistory->getCurrentState()->getTime();
    const int iStep = metaData->getIStep();
    volatile int order = metaData->getOrder();
    Scalar dt = metaData->getDt();
    bool output = metaData->getOutput();

    RCP<StepperState<Scalar> > stepperState = workingState->getStepperState();

    output = false;

    // If last time step was adjusted for output, reinstate previous dt.
    if (outputAdjustedDt_ == true) {
      if (printChanges) *out << changeDT(dt, dtAfterOutput_,
        "Reset dt after output.");
      dt = dtAfterOutput_;
      outputAdjustedDt_ = false;
      dtAfterOutput_ = 0.0;
    }

    if (dt <= 0.0) {
      if (printChanges) *out << changeDT(dt, getInitTimeStep(),
        "Reset dt to initial dt.");
      dt = getInitTimeStep();
    }

    if (dt < getMinTimeStep()) {
      if (printChanges) *out << changeDT(dt, getMinTimeStep(),
        "Reset dt to minimum dt.");
      dt = getMinTimeStep();
    }

    // update dt in metaData for the step control strategy to be informed
    metaData->setDt(dt);

    // call the step control strategy (to update order/dt if needed)
    stepControlStategy_->getNextTimeStep(*this, solutionHistory, integratorStatus);

    // get the order and dt (probably have changed by stepControlStategy_)
    order = metaData->getOrder();
    dt = metaData->getDt();

      if (dt < getMinTimeStep()) { // decreased below minimum dt
        if (printChanges) *out << changeDT(dt, getMinTimeStep(),
          "dt is too small.  Resetting to minimum dt.");
        dt = getMinTimeStep();
      }
      if (dt > getMaxTimeStep()) { // increased above maximum dt
        if (printChanges) *out << changeDT(dt, getMaxTimeStep(),
          "dt is too large.  Resetting to maximum dt.");
        dt = getMaxTimeStep();
      }


    // Check if we need to output this step index
    std::vector<int>::const_iterator it =
      std::find(outputIndices_.begin(), outputIndices_.end(), iStep);
    if (it != outputIndices_.end()) output = true;

    // Adjust time step to hit output times.
    Scalar reltol = 1.0e-6;
    for (size_t i=0; i < outputTimes_.size(); ++i) {
      const Scalar oTime = outputTimes_[i];
      if (lastTime < oTime && oTime <= lastTime+dt+getMinTimeStep()) {
        if (std::abs((lastTime+dt-oTime)/(lastTime+dt)) < reltol) {
          if (printChangesHi) *out << changeDT(dt, oTime - lastTime,
            "Adjusting dt for numerical roundoff to hit the next output time.");
          // Next output time IS VERY near next time (<reltol away from it),
          // e.g., adjust for numerical roundoff.
          output = true;
          outputAdjustedDt_ = true;
          dtAfterOutput_ = dt;
          dt = oTime - lastTime;
        } else if (lastTime*(1.0+reltol) < oTime &&
                   oTime < (lastTime+dt-getMinTimeStep())*(1.0+reltol)) {
          if (printChanges) *out << changeDT(dt, oTime - lastTime,
            "Adjusting dt to hit the next output time.");
          // Next output time is not near next time
          // (>getMinTimeStep() away from it).
          // Take time step to hit output time.
          output = true;
          outputAdjustedDt_ = true;
          dtAfterOutput_ = dt;
          dt = oTime - lastTime;
        } else {
          if (printChanges) *out << changeDT(dt, (oTime - lastTime)/2.0,
            "The next output time is within the minimum dt of the next time.  "
            "Adjusting dt to take two steps.");
          // Next output time IS near next time (<getMinTimeStep() away from it)
          // Take two time steps to get to next output time.
          dt = (oTime - lastTime)/2.0;
        }
        break;
      }
    }

    // Adjust time step to hit final time or correct for small
    // numerical differences.
    if ((lastTime + dt > getFinalTime() ) ||
        (std::abs((lastTime+dt-getFinalTime())/(lastTime+dt)) < reltol)) {
      if (printChangesHi) *out << changeDT(dt, getFinalTime() - lastTime,
        "Adjusting dt to hit the final time.");
      dt = getFinalTime() - lastTime;
    }

    // Time step always needs to keep time within range.
    TEUCHOS_TEST_FOR_EXCEPTION(
      (lastTime + dt < getInitTime()), std::out_of_range,
      "Error - Time step does not move time INTO time range.\n"
      "    [timeMin, timeMax] = [" << getInitTime() << ", "
      << getFinalTime() << "]\n"
      "    T + dt = " << lastTime <<" + "<< dt <<" = " << lastTime + dt <<"\n");

    TEUCHOS_TEST_FOR_EXCEPTION(
      (lastTime + dt > getFinalTime()), std::out_of_range,
      "Error - Time step move time OUT OF time range.\n"
      "    [timeMin, timeMax] = [" << getInitTime() << ", "
      << getFinalTime() << "]\n"
      "    T + dt = " << lastTime <<" + "<< dt <<" = " << lastTime + dt <<"\n");

    metaData->setOrder(order);
    metaData->setDt(dt);
    metaData->setTime(lastTime + dt);
    metaData->setOutput(output);
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
void TimeStepControl<Scalar>::setNumTimeSteps(int numTimeSteps)
{
  if (numTimeSteps >= 0) {
    tscPL_->set<int>        ("Number of Time Steps", numTimeSteps);
    const int initIndex = getInitIndex();
    tscPL_->set<int>        ("Final Time Index", initIndex + numTimeSteps);
    const double initTime = tscPL_->get<double>("Initial Time");
    const double finalTime = tscPL_->get<double>("Final Time");
    double initTimeStep = (finalTime - initTime)/numTimeSteps;
    if (numTimeSteps == 0) initTimeStep = Scalar(0.0);
    tscPL_->set<double>     ("Initial Time Step", initTimeStep);
    tscPL_->set<std::string>("Integrator Step Type", "Constant");

    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"setNumTimeSteps");
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
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (tscPL_ == Teuchos::null) {
      tscPL_ = Teuchos::parameterList("TimeStepControl");
      *tscPL_ = *(this->getValidParameters());
    }
  } else {
    tscPL_ = pList;
  }
  tscPL_->validateParametersAndSetDefaults(*this->getValidParameters(), 0);

  // Override parameters
  setNumTimeSteps(getNumTimeSteps());

  // set the time step control strategy
  setTimeStepControlStrategy();

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
    // Skip delimiters at the beginning
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find the first delimiter
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
      // Found a token, add it to the vector
      std::string token = str.substr(lastPos,pos-lastPos);
      outputTimes_.push_back(Scalar(std::stod(token)));
      if(pos==std::string::npos) break;

      lastPos = str.find_first_not_of(delimiters, pos); // Skip delimiters
      pos = str.find_first_of(delimiters, lastPos);     // Find next delimiter
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
      if(pos==std::string::npos) break;

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
void TimeStepControl<Scalar>::setTimeStepControlStrategy(
  Teuchos::RCP<TimeStepControlStrategy<Scalar> > tscs)
{
   using Teuchos::RCP;
   using Teuchos::ParameterList;

   if (stepControlStategy_ == Teuchos::null){
      stepControlStategy_ =
         Teuchos::rcp(new TimeStepControlStrategyComposite<Scalar>());
   }

   if (tscs == Teuchos::null) {
      // Create stepControlStategy_ if null, otherwise keep current parameters.

      if (getStepType() == "Constant"){
         stepControlStategy_->addStrategy(
               Teuchos::rcp(new TimeStepControlStrategyConstant<Scalar>()));
      } else if (getStepType() == "Variable") {
         // add TSCS from "Time Step Control Strategy List"

         RCP<ParameterList> tscsPL =
            Teuchos::sublist(tscPL_,"Time Step Control Strategy",true);
         // Construct from TSCS sublist
         std::vector<std::string> tscsLists;

         // string tokenizer
         tscsLists.clear();
         std::string str = tscsPL->get<std::string>("Time Step Control Strategy List");
         std::string delimiters(",");
         // Skip delimiters at the beginning
         std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
         // Find the first delimiter
         std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
         while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
            // Found a token, add it to the vector
            std::string token = str.substr(lastPos,pos-lastPos);
            tscsLists.push_back(token);
            if(pos==std::string::npos) break;

            lastPos = str.find_first_not_of(delimiters, pos); // Skip delimiters
            pos = str.find_first_of(delimiters, lastPos);     // Find next delimiter
         }

         // For each sublist name tokenized, add the TSCS
         for( auto el: tscsLists){

            RCP<Teuchos::ParameterList> pl =
               Teuchos::rcp(new ParameterList(tscsPL->sublist(el)));

            RCP<TimeStepControlStrategy<Scalar>> ts;

            // construct appropriate TSCS
            if(pl->get<std::string>("Name") == "Integral Controller")
               ts = Teuchos::rcp(new TimeStepControlStrategyIntegralController<Scalar>(pl));
            else if(pl->get<std::string>("Name") == "Basic VS")
               ts = Teuchos::rcp(new TimeStepControlStrategyBasicVS<Scalar>(pl));

            stepControlStategy_->addStrategy(ts);
         }
      }

   } else {
      // just add the new tscs to the vector of strategies
      stepControlStategy_->addStrategy(tscs);
   }

}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
TimeStepControl<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

  const double stdMin = 1.0e+04*std::numeric_limits<double>::epsilon();
  const double stdMax =         std::numeric_limits<double>::max();
  pl->set<double>("Initial Time"         , 0.0    , "Initial time");
  pl->set<double>("Final Time"           , stdMax , "Final time");
  pl->set<int>   ("Initial Time Index"   , 0      , "Initial time index");
  pl->set<int>   ("Final Time Index"     , 1000000, "Final time index");
  pl->set<double>("Minimum Time Step"    , stdMin , "Minimum time step size");
  pl->set<double>("Initial Time Step"    , stdMin , "Initial time step size");
  pl->set<double>("Maximum Time Step"    , stdMax , "Maximum time step size");
  //From (Denner, 2014), amplification factor can be at most 1.91 for stability.
/*
  pl->set<double>("Amplification Factor" , 1.75   , "Amplification factor");
  pl->set<double>("Reduction Factor"     , 0.5    , "Reduction factor");
  //FIXME? may need to modify default values of monitoring function
  //IKT, 1/5/17: from (Denner, 2014), it seems a reasonable choice for
  //eta_min is 0.1*eta_max.  Numerical tests confirm this.
  //TODO: Change default value of eta_min to 1.0e-2?
  pl->set<double>("Minimum Value Monitoring Function", 1.0e-6, "Min value eta");
  pl->set<double>("Maximum Value Monitoring Function", 1.0e-1, "Max value eta");
*/
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

   Teuchos::RCP<Teuchos::ParameterList> tscsPL = Teuchos::parameterList("Time Step Control Strategy");
   tscsPL->set<std::string>("Time Step Control Strategy List","");
   pl->set("Time Step Control Strategy", *tscsPL);
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
