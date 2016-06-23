#ifndef TEMPUS_SOLUTIONSTATEMETADATA_IMPL_HPP
#define TEMPUS_SOLUTIONSTATEMETADATA_IMPL_HPP


namespace Tempus {

// SolutionStateMetaData definitions:
template<class Scalar>
SolutionStateMetaData<Scalar>::SolutionStateMetaData()
  :time_          (0.0),
   iStep_         (0),
   dt_            (0.0),
   errorAbs_      (0.0),
   errorRel_      (0.0),
   order_         (1),
   nFailures_     (0),
   nConsecutiveFailures_(0),
   solutionStatus_(WORKING),
   output_        (false),
   isRestartable_ (true),
   isInterpolated_(false),
   accuracy_      (0.0)
{}

template<class Scalar>
SolutionStateMetaData<Scalar>::SolutionStateMetaData(
  const Scalar time,
  const int    iStep,
  const Scalar dt,
  const Scalar errorAbs,
  const Scalar errorRel,
  const int    order,
  const int    nFailures,
  const int    nConsecutiveFailures,
  const Status solutionStatus,
  const bool   output,
  const bool   isRestartable,
  const bool   isInterpolated,
  const Scalar accuracy)
  :time_          (time),
   iStep_         (iStep),
   dt_            (dt),
   errorAbs_      (errorAbs),
   errorRel_      (errorRel),
   order_         (order),
   nFailures_     (nFailures),
   nConsecutiveFailures_(nConsecutiveFailures),
   solutionStatus_(solutionStatus),
   output_        (output),
   isRestartable_ (isRestartable),
   isInterpolated_(isInterpolated),
   accuracy_      (accuracy)
{}

template<class Scalar>
SolutionStateMetaData<Scalar>::SolutionStateMetaData(const SolutionStateMetaData<Scalar>& ssmd_)
  :time_          (ssmd_.time_),
   iStep_         (ssmd_.iStep_),
   dt_            (ssmd_.dt_),
   errorAbs_      (ssmd_.errorAbs_),
   errorRel_      (ssmd_.errorRel_),
   order_         (ssmd_.order_),
   nFailures_     (ssmd_.nFailures_),
   nConsecutiveFailures_(ssmd_.nConsecutiveFailures_),
   solutionStatus_(ssmd_.solutionStatus_),
   output_        (ssmd_.output_),
   isRestartable_ (ssmd_.isRestartable_),
   isInterpolated_(ssmd_.isInterpolated_),
   accuracy_      (ssmd_.accuracy_)
{}


template<class Scalar>
Teuchos::RCP<SolutionStateMetaData<Scalar> > SolutionStateMetaData<Scalar>::clone()
{
  Teuchos::RCP<SolutionStateMetaData<Scalar> > md =
    rcp(new SolutionStateMetaData<Scalar> (
      time_,
      iStep_,
      dt_,
      errorAbs_,
      errorRel_,
      order_,
      nFailures_,
      nConsecutiveFailures_,
      solutionStatus_,
      output_,
      isRestartable_,
      isInterpolated_,
      accuracy_));

  return md;
}


template<class Scalar>
std::string SolutionStateMetaData<Scalar>::description() const
{
  std::string name = "Tempus::SolutionStateMetaData";
  return(name);
}


template<class Scalar>
void SolutionStateMetaData<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  if (verbLevel == Teuchos::VERB_EXTREME) {
    out << description() << "::describe:" << std::endl
        << "time           = " << time_ << std::endl
        << "iStep          = " << iStep_ << std::endl
        << "dt             = " << dt_ << std::endl
        << "errorAbs       = " << errorAbs_ << std::endl
        << "errorRel       = " << errorRel_ << std::endl
        << "order          = " << order_ << std::endl
        << "nFailures      = " << nFailures_ << std::endl
        << "nConsecutiveFailures = " << nConsecutiveFailures_ << std::endl
        << "solutionStatus = " << toString(solutionStatus_) << std::endl
        << "output         = " << output_ << std::endl
        << "isRestartable  = " << isRestartable_ << std::endl
        << "isInterpolated = " << isInterpolated_ << std::endl
        << "accuracy       = " << accuracy_ << std::endl;
  }
}

} // namespace Tempus
#endif // TEMPUS_SOLUTIONSTATEMETADATA_IMPL_HPP
