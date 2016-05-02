#ifndef TEMPUS_SOLUTIONHISTORY_IMPL_HPP
#define TEMPUS_SOLUTIONHISTORY_IMPL_HPP

// Teuchos
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
// Tempus
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_SolutionStateMetaData.hpp"

//#include "Thyra_VectorStdOps.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;

namespace {

  static std::string Invalid_name    = "Invalid";
  static std::string KeepNewest_name = "Keep Newest";
  static std::string Undo_name       = "Undo";
  static std::string Static_name     = "Static";
  static std::string Unlimited_name  = "Unlimited";
  static std::string Storage_name    = "Storage Type";
  static std::string Storage_default = KeepNewest_name;

  static std::string StorageLimit_name    = "Storage Limit";
  static int         StorageLimit_default = 1;

  Teuchos::Array<std::string>
    HistoryPolicies = Teuchos::tuple<std::string>(
      Invalid_name,
      KeepNewest_name,
      Undo_name,
      Static_name,
      Unlimited_name);

  const RCP<Teuchos::StringToIntegralParameterEntryValidator<Tempus::HistoryPolicy> >
    PolicyValidator = Teuchos::rcp(
        new Teuchos::StringToIntegralParameterEntryValidator<Tempus::HistoryPolicy>(
          HistoryPolicies,
          Teuchos::tuple<Tempus::HistoryPolicy>(
            Tempus::HISTORY_POLICY_INVALID,
            Tempus::HISTORY_POLICY_KEEP_NEWEST,
            Tempus::HISTORY_POLICY_UNDO,
            Tempus::HISTORY_POLICY_STATIC,
            Tempus::HISTORY_POLICY_UNLIMITED),
          Storage_name));

} // namespace


namespace Tempus {

template<class Scalar>
SolutionHistory<Scalar>::SolutionHistory(RCP<ParameterList> pList_)
{
  // Create history, an array of solution states.
  history = rcp(new Teuchos::Array<SolutionState<Scalar> >);

  if (pList_ == Teuchos::null) {
    pList     = Teuchos::null;
    //interpolator  = Teuchos::null;
    storage_limit = StorageLimit_default;
    historyPolicy = HISTORY_POLICY_KEEP_NEWEST;
  } else {
    this->setParameterList(pList_);
  }

  if (Teuchos::as<int>(this->getVerbLevel()) >=
      Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"SolutionHistory::SolutionHistory");
    *out << this->description() << std::endl;
  }
}


template<class Scalar>
void SolutionHistory<Scalar>::addState(
  const RCP<SolutionState<Scalar> >& state_)
{
  // Check that we're not going to exceed our storage limit:
  if (Teuchos::as<int>(history->size()+1) > storage_limit) {
    switch (historyPolicy) {
    case HISTORY_POLICY_INVALID: {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - history policy is HISTORY_POLICY_INVALID.\n");
      break;
    }
    case HISTORY_POLICY_STATIC: {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - history will overflow and policy is HISTORY_POLICY_STATIC.  "
        "This state can not be added\n");
      break;
    }
    case HISTORY_POLICY_KEEP_NEWEST:
    case HISTORY_POLICY_UNDO:
    case HISTORY_POLICY_STATIC: {
      if (state_->getTime() > history->front()->getTime()) {
        // Case:  State is newer than the oldest state in history.
        // Remove state from the beginning of history, then add new state.
        history->erase(history->begin());
      } else {
        // Case:  State is older than the oldest state in history.
        RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,"SolutionHistory::addState");
        *out << "Warning, state is older than oldest state in history.  "
             << "State not added!\n" << std::endl;
        return;
      }
      break;
    }
    case HISTORY_POLICY_UNLIMITED:
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown history policy.\n");
    }
  }

  // Add new state in chronological order.
  if (history->size() == 0) {
    history->push_back(state_);
  } else {
    typename Teuchos::Array<SolutionState<Scalar> >::iterator
      state_it = history->begin();
    for (state_it; state_it < history->end(); ++state_it) {
      if (state_ < state_it) break;
    }
    history->insert(state_it, state_);
  }
}

template<class Scalar>
void SolutionHistory<Scalar>::removeState(
  const RCP<SolutionState<Scalar> >& state_)
{
  if (history->size() != 0) {
    typename Teuchos::Array<SolutionState<Scalar> >::reverse_iterator
      state_it = history->rbegin();
    for (state_it; state_it < history->rend(); ++state_it) {
      if (state_ == state_it) break;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(state_it == history->rend(), std::logic_error,
      "Error - removeState() Could not remove state = "
      // << state_it->describe()
      );

    history->erase(state_it);
  }
  return;
}


template<class Scalar>
void SolutionHistory<Scalar>::removeState(const Scalar time)
{
  RCP<SolutionState<Scalar> > tmpState = findState(time);
  removeState(tmpState);
}


template<class Scalar>
RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::findState(const Scalar time) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(minTime() <= time and time <= maxTime()), std::logic_error,
    "Error - SolutionHistory::findState() Requested time out of range!\n"
    "        [Min, Max] = [" << minTime() << ", " << maxTime() << "]\n"
    "        time = "<< time <<"\n");

  Scalar reltol = 1.0e-14;
  typename Teuchos::Array<SolutionState<Scalar> >::iterator
    state_it = history->begin();
  // Linear search
  for (state_it; state_it < history->end(); ++state_it) {
    if (abs((state_it.getTime()-time)/state_it.getTime()) < reltol) break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(state_it == history->end(), std::logic_error,
    "Error - SolutionHistory::findState()!\n"
    "        Did not find a SolutionState with time = " <<time<< std::endl);

  return history[state_it];
}


template<class Scalar>
RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::interpolateState(const Scalar time) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(minTime() <= time and time <= maxTime()), std::logic_error,
    "Error - SolutionHistory::getTime() Requested time out of range!\n"
    "        [Min, Max] = [" << minTime() << ", " << maxTime() << "]\n"
    "        time = "<< time <<"\n");

  // Interpolate the state.
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "Error - Implement interpolation.\n");
  //SolutionState<Scalar> state_out;
  //interpolate<Scalar>(*interpolator, history, time_vec, &state_out);
}


template<class Scalar>
RCP<SolutionState<Scalar> > SolutionHistory<Scalar>::getCurrentState() const
{
  return currentState;
}


template<class Scalar>
RCP<SolutionState<Scalar> > SolutionHistory<Scalar>::getWorkingState() const
{
  return workingState;
}

/** Initialize the working state
 *
 *  <b>Tasks:</b>
 *    -# Get the working state by creating a new SolutionState or
 *       use the last SolutionState in the SolutionHistory.
 *    -# Set working state to current state (deepcopy).
 *       This is making the initial guess for the next time step
 *       equal to the current time step.
 */
template<class Scalar>
RCP<SolutionState<Scalar> > SolutionHistory<Scalar>::initWorkingState()
{
  TEUCHOS_TEST_FOR_EXCEPTION(history->size() == 0, std::logic_error,
    "Error - SolutionHistory::initWorkingState()\n"
    "Can not initialize working state without at least one solution state!\n");

  if (history->size() < storage_limit) {
    workingState->clone(currentState);
  } else {
    workingState = history->back();
    workingState->deepcopy(currentState);
  }
}


template<class Scalar>
void SolutionHistory<Scalar>::promoteWorkingState()
{
  RCP<SolutionStateMetaData<Scalar> > md = workingState->metaData;
  md->time += md->dt;
  md->iStep++;
  md->nConsecutiveFailures = std::max(0,md->nConsecutiveFailures-1);
  md->status = SolutionStatus::PASSING;
  md->isAccepted = true;
  md->isRestartable = true;
  md->isInterpolated = false;
  currentState = workingState;
}


template<class Scalar>
void SolutionHistory<Scalar>::setStorage(int storage)
{
  int storage_limit = std::max(1,storage);

  TEUCHOS_TEST_FOR_EXCEPTION(
    (Teuchos::as<int>(history->size()) > storage_limit), std::logic_error,
    "Error - requested storage limit = " << storage_limit
    << " is smaller than the current number of states stored = "
    << history->size() << "!\n");
}


template<class Scalar>
int SolutionHistory<Scalar>::getStorage() const
{
  return(storage_limit);
}


template <class Scalar>
HistoryPolicy SolutionHistory<Scalar>::getHistoryPolicy()
{
  return historyPolicy;
}


template<class Scalar>
Scalar SolutionHistory<Scalar>::minTime() const
{
  return (history->front())->getTime();
}


template<class Scalar>
Scalar SolutionHistory<Scalar>::maxTime() const
{
  return (history->back())->getTime();
}


template<class Scalar>
std::string SolutionHistory<Scalar>::description() const
{
  std::string name = "Tempus::SolutionHistory";
  return(name);
}


template<class Scalar>
void SolutionHistory<Scalar>::describe(
  Teuchos::FancyOStream          &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  if ((Teuchos::as<int>(verbLevel)==Teuchos::as<int>(Teuchos::VERB_DEFAULT)) ||
      (Teuchos::as<int>(verbLevel)>=Teuchos::as<int>(Teuchos::VERB_LOW)    )  ){
    out << description() << "::describe" << std::endl;
    //out << "interpolator     = " << interpolator->description() << std::endl;
    out << "storage_limit    = " << storage_limit << std::endl;
    out << "historyPolicy    = " << historyPolicy << std::endl;
    out << "number of states = " << history->size() << std::endl;
    out << "time range       = (" << history->front()->getTime() << ", "
                                  << history->back()->getTime() << ")" << std::endl;
  } else if (Teuchos::as<int>(verbLevel) >=
             Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    out << "SolutionStates: " << std::endl;
    for (Teuchos::Ordinal i=0; i<history->size() ; ++i) {
      out << "SolutionState[" << i << "] = " << std::endl;
      (*history)[i].describe(out,this->getVerbLevel());
    }
  }
}


template <class Scalar>
void SolutionHistory<Scalar>::setParameterList(
  RCP<ParameterList> const& pList_)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pList_));
  pList_->validateParameters(*this->getValidParameters());
  pList = pList_;

  Teuchos::readVerboseObjectSublist(&*pList,this);

  //setInterpolator(interpolator);

  HistoryPolicy policy_ = PolicyValidator->getIntegralValue(
      *pList, Storage_name, Storage_default);

  if (policy_ != HISTORY_POLICY_INVALID)
    historyPolicy = policy_;
  else
    historyPolicy = HISTORY_POLICY_KEEP_NEWEST;

  int storage_limit = pList->get(StorageLimit_name, StorageLimit_default);
  setStorage(storage_limit);
}


template<class Scalar>
RCP<const ParameterList> SolutionHistory<Scalar>::getValidParameters() const
{
  static RCP<ParameterList> validPL;

  if (is_null(validPL)) {

    RCP<ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    pl->set(Storage_name, Storage_default,
      "History Policy sets the memory storage.  "
      "'Keep Newest' - will retain the single newest solution state.  "
      "'Undo' - will retain two solution states in order to do a single undo.  "
      "'Static' - will retain 'Storage Limit' number of solution states.  "
      "'Unlimited' - will not remove any solution states!",
      PolicyValidator);

    pl->set(StorageLimit_name, StorageLimit_default,
      "Storage limit for the solution history.");

    validPL = pl;

  }
  return validPL;
}


template <class Scalar>
RCP<ParameterList>
SolutionHistory<Scalar>::getNonconstParameterList()
{
  return(pList);
}


template <class Scalar>
RCP<ParameterList> SolutionHistory<Scalar>::unsetParameterList()
{
  RCP<ParameterList> temp_plist = pList;
  pList = Teuchos::null;
  return(temp_plist);
}


//template<class Scalar>
//void SolutionHistory<Scalar>::setInterpolator(
//  const RCP<InterpolatorBase<Scalar> >& interpolator_)
//{
//  if (interpolator_ == Teuchos::null) {
//    interpolator = linearInterpolator<Scalar>();
//  } else {
//    interpolator = interpolator_;
//  }
//  if (Teuchos::as<int>(this->getVerbLevel()) >=
//      Teuchos::as<int>(Teuchos::VERB_HIGH)) {
//    RCP<Teuchos::FancyOStream> out = this->getOStream();
//    Teuchos::OSTab ostab(out,1,"SolutionHistory::setInterpolator");
//    *out << "interpolator = " << interpolator_->description() << std::endl;
//  }
//}
//
//template<class Scalar>
//RCP<InterpolatorBase<Scalar> >
//  SolutionHistory<Scalar>::getNonconstInterpolator()
//{
//  return interpolator;
//}
//
//template<class Scalar>
//RCP<const InterpolatorBase<Scalar> >
//  SolutionHistory<Scalar>::getInterpolator() const
//{
//  return interpolator;
//}
//
//template<class Scalar>
//RCP<InterpolatorBase<Scalar> > SolutionHistory<Scalar>::unSetInterpolator()
//{
//  RCP<InterpolatorBase<Scalar> > old_interpolator = interpolator;
//  interpolator = linearInterpolator<Scalar>();
//  return old_interpolator;
//}


} // namespace Tempus
#endif // TEMPUS_SOLUTIONHISTORY_IMPL_HPP
