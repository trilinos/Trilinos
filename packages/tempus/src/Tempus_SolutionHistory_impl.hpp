#ifndef TEMPUS_SOLUTIONHISTORY_IMPL_HPP
#define TEMPUS_SOLUTIONHISTORY_IMPL_HPP

#include "Tempus_SolutionHistory.hpp"

//#include "Thyra_VectorStdOps.hpp"

#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace {

  static std::string Invalid_name      = "Invalid";
  static std::string KeepNewest_name   = "Keep Newest";
  static std::string Undo_name         = "Undo";
  static std::string Static_name       = "Static";
  static std::string Unlimited_name    = "Unlimited";
  static std::string Selection_name    = "History Policy";
  static std::string Selection_default = KeepNewest_name;

  static std::string StorageLimit_name    = "Storage Limit";
  static int         StorageLimit_default = 0;

  Teuchos::Array<std::string>
    HistoryPolicies = Teuchos::tuple<std::string>(
      Invalid_name,
      KeepNewest_name,
      Undo_name,
      Static_name,
      Unlimited_name);

  const Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<tempus::HistoryPolicy> >
    PolicyValidator = Teuchos::rcp(
        new Teuchos::StringToIntegralParameterEntryValidator<tempus::HistoryPolicy>(
          HistoryPolicies,
          Teuchos::tuple<tempus::HistoryPolicy>(
            tempus::HISTORY_POLICY_INVALID,
            tempus::HISTORY_POLICY_KEEP_NEWEST,
            tempus::HISTORY_POLICY_UNDO,
            tempus::HISTORY_POLICY_STATIC,
            tempus::HISTORY_POLICY_UNLIMITED),
          Selection_name));

} // namespace


namespace tempus {


template<class Scalar>
SolutionHistory<Scalar>::SolutionHistory(
  RCP<Teuchos::ParameterList> paramList_ = Teuchos::null )
{
  // Create history, an array of solution states.
  history = rcp(new Array<SolutionState<Scalar> >);

  if ( paramList_ == Teuchos::null ) {
    paramList     = Teuchos::null;
    interpolator  = Teuchos::null;
    storage_limit = 1;
    historyPolicy = HISTORY_POLICY_KEEP_NEWEST;
  } else {
    this->setParameterList(paramList_);
  }

  if ( Teuchos::as<int>(this->getVerbLevel()) >=
       Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"SolutionHistory::SolutionHistory");
    *out << this->description() << std::endl;
  }
}


template<class Scalar>
void SolutionHistory<Scalar>::setStorage( int storage )
{
  int storage_limit = std::max(1,storage);

  TEUCHOS_TEST_FOR_EXCEPTION(
    (Teuchos::as<int>(history->size()) > storage_limit), std::logic_error,
    "Error, requested storage limit = " << storage_limit
    << " is smaller than the current number of states stored = "
    << history->size() << "!\n");
  }
}


template<class Scalar>
int SolutionHistory<Scalar>::getStorage() const
{
  return(storage_limit);
}


template<class Scalar>
void SolutionHistory<Scalar>::setInterpolator(
  const RCP<InterpolatorBase<Scalar> >& interpolator_)
{
  if (interpolator_ == Teuchos::null) {
    interpolator = linearInterpolator<Scalar>();
  } else {
    interpolator = interpolator_;
  }
  if ( Teuchos::as<int>(this->getVerbLevel()) >=
       Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"SolutionHistory::setInterpolator");
    *out << "interpolator = " << interpolator_->description() << std::endl;
  }
}

template<class Scalar>
RCP<InterpolatorBase<Scalar> >
  SolutionHistory<Scalar>::getNonconstInterpolator()
{
  return interpolator;
}

template<class Scalar>
RCP<const InterpolatorBase<Scalar> >
  SolutionHistory<Scalar>::getInterpolator() const
{
  return interpolator;
}

template<class Scalar>
RCP<InterpolatorBase<Scalar> > SolutionHistory<Scalar>::unSetInterpolator()
{
  RCP<InterpolatorBase<Scalar> > old_interpolator = interpolator;
  interpolator = linearInterpolator<Scalar>();
  return old_interpolator;
}


template<class Scalar>
void SolutionHistory<Scalar>::addState( const RCP<SolutionState<Scalar> >& state_ )
{
  // Check that we're not going to exceed our storage limit:
  if (Teuchos::as<int>(history->size()+1) > storage_limit) {
    switch (historyPolicy) {
    case HISTORY_POLICY_INVALID: {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "Error, history policy is HISTORY_POLICY_INVALID.\n");
      break;
    }
    case HISTORY_POLICY_STATIC: {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "Error, history will overflow and policy is HISTORY_POLICY_STATIC.  "
        "This state can not be added\n");
      break;
    }
    case HISTORY_POLICY_KEEP_NEWEST:
    case HISTORY_POLICY_UNDO:
    case HISTORY_POLICY_STATIC: {
      if (state_.time > history->front().time) {
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
    default: {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "Error, unknown history policy.\n");
    }
  }

  // Add new state in chronological order.
  Array<SolutionState<Scalar> >::iterator state_it = history->begin();
  for (state_it ; state_it < history->end() ; ++state_it) {
    if (state_.time < state_it.time ) break;
  }
  history->insert(state_it,state_);
}


template<class Scalar>
RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::getState( const Scalar time ) const
{
  reltol = 1.0e-14;
  Array<SolutionState<Scalar> >::iterator state_it = history->begin();
  // Linear search
  for (state_it ; state_it < history->end() ; ++state_it) {
    if ( abs((state_it.time-time)/state_it.time) < reltol ) break;
  }

  if ( state_it != history->end() )
    // Return original state.
    return state_it;
  else {
    // Interpolate the state.
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
      "Error, Implement interpolation.\n");
    SolutionState<Scalar> state_out;
    interpolate<Scalar>(*interpolator, history, time_vec, &state_out);
  }
}


template<class Scalar>
TimeRange<Scalar> SolutionHistory<Scalar>::getTimeRange() const
{
  TimeRange<Scalar> timerange;
  if (history->size() > 0)
    timerange = TimeRange<Scalar>(history->front().time,history->back().time);
  return(timerange);
}


template<class Scalar>
void SolutionHistory<Scalar>::getTimes( Array<Scalar>* time_vec ) const
{
  int N = history->size();
  time_vec->clear();
  time_vec->reserve(N);
  for (int i=0 ; i<N ; ++i)
    time_vec->push_back((*history)[i].time);
}


template<class Scalar>
void SolutionHistory<Scalar>::removeStates( Array<Scalar>& time_vec )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  int N = time_vec.size();

  SolutionState<Scalar> ss_tmp();
  for (int i=0; i<N ; ++i) {
    ss_tmp.time = time_vec[i];
    Array<SolutionState<Scalar> >::iterator state_it =
      std::find(history->begin(),history->end(), ss_tmp);
    TEUCHOS_TEST_FOR_EXCEPTION( state_it == history->end(), std::logic_error,
      "Error, time_vec[" << i << "] = " << time_vec[i]
       << "is not a time in the history!\n");
    history->erase(state_it);
  }
}


template<class Scalar>
int SolutionHistory<Scalar>::getInterplatorOrder() const
{
  return(interpolator->order());
}


template<class Scalar>
std::string SolutionHistory<Scalar>::description() const
{
  std::string name = "tempus::SolutionHistory";
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
    out << "interpolator     = " << interpolator->description() << std::endl;
    out << "storage_limit    = " << storage_limit << std::endl;
    out << "historyPolicy    = " << historyPolicy << std::endl;
    out << "number of states = " << history->size() << std::endl;
    out << "time range       = (" << history->front().time << ", "
                                  << history->back().time << ")" << std::endl;
  } else if ( Teuchos::as<int>(verbLevel) >=
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
  RCP<Teuchos::ParameterList> const& paramList_)
{
  TEUCHOS_TEST_FOR_EXCEPT( is_null(paramList_) );
  paramList_->validateParameters(*this->getValidParameters());
  paramList = paramList_;

  Teuchos::readVerboseObjectSublist(&*paramList,this);

  setInterpolator(interpolator);

  HistoryPolicy policy_ = PolicyValidator->getIntegralValue(
      *paramList, Selection_name, Selection_default);

  if (policy_ != HISTORY_POLICY_INVALID)
    historyPolicy = policy_;
  else
    historyPolicy = HISTORY_POLICY_KEEP_NEWEST;

  int storage_limit = paramList->get( StorageLimit_name, StorageLimit_default);
  setStorage(storage_limit);
}


template<class Scalar>
RCP<const Teuchos::ParameterList> SolutionHistory<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    pl->set( Selection_name, Selection_default,
      "History Policy sets the memory storage.  "
      "'Keep Newest' - will retain the single newest solution state.  "
      "'Undo' - will retain two solution states in order to do a single undo.  "
      "'Static' - will retain 'Storage Limit' number of solution states.  "
      "'Unlimited' - will not remove any solution states!"
      PolicyValidator);

    pl->set( StorageLimit_name, StorageLimit_default,
      "Storage limit for the solution history.");

    validPL = pl;

  }
  return validPL;
}


template <class Scalar>
RCP<Teuchos::ParameterList>
SolutionHistory<Scalar>::getNonconstParameterList()
{
  return(paramList);
}


template <class Scalar>
RCP<Teuchos::ParameterList> SolutionHistory<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = paramList;
  paramList = Teuchos::null;
  return(temp_param_list);
}

template <class Scalar>
HistoryPolicy SolutionHistory<Scalar>::getHistoryPolicy()
{
  return historyPolicy;
}

} // namespace tempus
#endif // TEMPUS_SOLUTIONHISTORY_IMPL_HPP
