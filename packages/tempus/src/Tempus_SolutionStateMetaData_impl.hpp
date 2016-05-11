#ifndef TEMPUS_SOLUTIONSTATEMETADATA_IMPL_HPP
#define TEMPUS_SOLUTIONSTATEMETADATA_IMPL_HPP


namespace Tempus {

// SolutionStateMetaData definitions:
template<class Scalar>
SolutionStateMetaData<Scalar>::SolutionStateMetaData()
  :time          (0.0),
   iStep         (0),
   dt            (0.0),
   suggestedDt   (0.0),
   errorAbs      (0.0),
   errorRel      (0.0),
   order         (1),
   nFailures     (0),
   nConsecutiveFailures(0),
   status        (PASSED),
   output        (false),
   isAccepted    (false),
   isRestartable (true),
   isInterpolated(false),
   accuracy      (0.0)
{}

template<class Scalar>
SolutionStateMetaData<Scalar>::SolutionStateMetaData(
  const Scalar time_,
  const int    iStep_,
  const Scalar dt_,
  const Scalar suggestedDt_,
  const Scalar errorAbs_,
  const Scalar errorRel_,
  const int    order_,
  const int    nFailures_,
  const int    nConsecutiveFailures_,
  const SolutionStatus status_,
  const bool   output_,
  const bool   isAccepted_,
  const bool   isInterpolated_,
  const bool   isRestartable_,
  const Scalar accuracy_)
  :time          (time_),
   iStep         (iStep_),
   dt            (dt_),
   suggestedDt   (suggestedDt_),
   errorAbs      (errorAbs_),
   errorRel      (errorRel_),
   order         (order_),
   nFailures     (nFailures_),
   nConsecutiveFailures(nConsecutiveFailures_),
   status        (status_),
   output        (output_),
   isAccepted    (isAccepted_),
   isRestartable (isRestartable_),
   isInterpolated(isInterpolated_),
   accuracy      (accuracy_)
{}

template<class Scalar>
SolutionStateMetaData<Scalar>::SolutionStateMetaData(const SolutionStateMetaData<Scalar>& ssmd_)
  :time          (ssmd_.time),
   iStep         (ssmd_.iStep),
   dt            (ssmd_.dt),
   suggestedDt   (ssmd_.suggestedDt),
   errorAbs      (ssmd_.errorAbs),
   errorRel      (ssmd_.errorRel),
   order         (ssmd_.order),
   nFailures     (ssmd_.nFailures),
   nConsecutiveFailures(ssmd_.nConsecutiveFailures),
   status        (ssmd_.status),
   output        (ssmd_.output),
   isAccepted    (ssmd_.isAccepted),
   isRestartable (ssmd_.isRestartable),
   isInterpolated(ssmd_.isInterpolated),
   accuracy      (ssmd_.accuracy)
{}


template<class Scalar>
Teuchos::RCP<SolutionStateMetaData<Scalar> > SolutionStateMetaData<Scalar>::clone()
{
  Teuchos::RCP<SolutionStateMetaData<Scalar> > md =
    rcp(new SolutionStateMetaData<Scalar> (
      time,
      iStep,
      dt,
      suggestedDt,
      errorAbs,
      errorRel,
      order,
      nFailures,
      nConsecutiveFailures,
      status,
      output,
      isAccepted,
      isRestartable,
      isInterpolated,
      accuracy));

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
        << "time           = " << time << std::endl
        << "iStep          = " << iStep << std::endl
        << "dt             = " << dt << std::endl
        << "suggestedDt    = " << suggestedDt << std::endl
        << "errorAbs       = " << errorAbs << std::endl
        << "errorRel       = " << errorRel << std::endl
        << "order          = " << order << std::endl
        << "nFailures      = " << nFailures << std::endl
        << "nConsecutiveFailures = " << nConsecutiveFailures << std::endl
        << "status         = " << toString(status) << std::endl
        << "output         = " << output << std::endl
        << "isAccepted     = " << isAccepted    << std::endl
        << "isRestartable  = " << isRestartable << std::endl
        << "isInterpolated = " << isInterpolated << std::endl
        << "accuracy       = " << accuracy << std::endl;
  }
}

} // namespace Tempus
#endif // TEMPUS_SOLUTIONSTATEMETADATA_IMPL_HPP
