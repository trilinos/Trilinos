// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_SolutionStateMetaData_impl_hpp
#define Tempus_SolutionStateMetaData_impl_hpp

#include "Tempus_config.hpp"

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
   nRunningFailures_(0),
   nConsecutiveFailures_(0),
   solutionStatus_(WORKING),
   output_        (false),
   outputScreen_  (false),
   isSynced_      (true),
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
  const int    nRunningFailures,
  const int    nConsecutiveFailures,
  const Status solutionStatus,
  const bool   output,
  const bool   outputScreen,
  const bool   isSynced,
  const bool   isInterpolated,
  const Scalar accuracy)
  :time_          (time),
   iStep_         (iStep),
   dt_            (dt),
   errorAbs_      (errorAbs),
   errorRel_      (errorRel),
   order_         (order),
   nFailures_     (nFailures),
   nRunningFailures_(nRunningFailures),
   nConsecutiveFailures_(nConsecutiveFailures),
   solutionStatus_(solutionStatus),
   output_        (output),
   outputScreen_  (outputScreen),
   isSynced_      (isSynced),
   isInterpolated_(isInterpolated),
   accuracy_      (accuracy)
{}

template<class Scalar>
SolutionStateMetaData<Scalar>::SolutionStateMetaData(const SolutionStateMetaData<Scalar>& ssmd)
  :time_          (ssmd.time_),
   iStep_         (ssmd.iStep_),
   dt_            (ssmd.dt_),
   errorAbs_      (ssmd.errorAbs_),
   errorRel_      (ssmd.errorRel_),
   order_         (ssmd.order_),
   nFailures_     (ssmd.nFailures_),
   nRunningFailures_(ssmd.nRunningFailures_),
   nConsecutiveFailures_(ssmd.nConsecutiveFailures_),
   solutionStatus_(ssmd.solutionStatus_),
   output_        (ssmd.output_),
   outputScreen_  (ssmd.outputScreen_),
   isSynced_      (ssmd.isSynced_),
   isInterpolated_(ssmd.isInterpolated_),
   accuracy_      (ssmd.accuracy_)
{}


template<class Scalar>
Teuchos::RCP<SolutionStateMetaData<Scalar> > SolutionStateMetaData<Scalar>::clone() const
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
      nRunningFailures_,
      nConsecutiveFailures_,
      solutionStatus_,
      output_,
      outputScreen_,
      isSynced_,
      isInterpolated_,
      accuracy_));

  return md;
}


template<class Scalar>
void SolutionStateMetaData<Scalar>::
copy(const Teuchos::RCP<const SolutionStateMetaData<Scalar> >& ssmd)
{
  time_           = ssmd->time_;
  iStep_          = ssmd->iStep_;
  dt_             = ssmd->dt_;
  errorAbs_       = ssmd->errorAbs_;
  errorRel_       = ssmd->errorRel_;
  order_          = ssmd->order_;
  nFailures_      = ssmd->nFailures_;
  nRunningFailures_= ssmd->nRunningFailures_;
  nConsecutiveFailures_ = ssmd->nConsecutiveFailures_;
  solutionStatus_ = ssmd->solutionStatus_;
  output_         = ssmd->output_;
  outputScreen_   = ssmd->outputScreen_;
  isSynced_       = ssmd->isSynced_;
  isInterpolated_ = ssmd->isInterpolated_;
  accuracy_       = ssmd->accuracy_;
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
        << "nRunningFailures= " << nRunningFailures_<< std::endl
        << "nConsecutiveFailures = " << nConsecutiveFailures_ << std::endl
        << "solutionStatus = " << toString(solutionStatus_) << std::endl
        << "output         = " << output_ << std::endl
        << "outputScreen   = " << outputScreen_ << std::endl
        << "isSynced       = " << isSynced_ << std::endl
        << "isInterpolated = " << isInterpolated_ << std::endl
        << "accuracy       = " << accuracy_ << std::endl;
  }
}

} // namespace Tempus
#endif // Tempus_SolutionStateMetaData_impl_hpp
