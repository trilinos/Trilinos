// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorObserverLogging_impl_hpp
#define Tempus_IntegratorObserverLogging_impl_hpp

#include "Tempus_IntegratorObserver.hpp"
#include "Tempus_TimeStepControl.hpp"

namespace Tempus {

template<class Scalar>
IntegratorObserverLogging<Scalar>::IntegratorObserverLogging(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory,
  const Teuchos::RCP<TimeStepControl<Scalar> >& timeStepControl)
  : nameObserveStartIntegrator_ ("observeStartIntegrator" ),
    nameObserveStartTimeStep_   ("observeStartTimeStep"   ),
    nameObserveNextTimeStep_    ("observeNextTimeStep"    ),
    nameObserveBeforeTakeStep_  ("observeBeforeTakeStep"  ),
    nameObserveAfterTakeStep_   ("observeAfterTakeStep"   ),
    nameObserveAcceptedTimeStep_("observeAcceptedTimeStep"),
    nameObserveEndIntegrator_   ("observeEndIntegrator"   ),
    solutionHistory_(solutionHistory), timeStepControl_(timeStepControl)
{
  counters_ = Teuchos::rcp(new std::map<std::string,int>);
  order_ = Teuchos::rcp(new std::list<std::string>);
  this->resetLogCounters();
}

template<class Scalar>
IntegratorObserverLogging<Scalar>::~IntegratorObserverLogging(){}

template<class Scalar>
void IntegratorObserverLogging<Scalar>::observeStartIntegrator()
{ logCall(nameObserveStartIntegrator_); }

template<class Scalar>
void IntegratorObserverLogging<Scalar>::observeStartTimeStep()
{ logCall(nameObserveStartTimeStep_); }

template<class Scalar>
void IntegratorObserverLogging<Scalar>::
observeNextTimeStep(Status & integratorStatus)
{ logCall(nameObserveNextTimeStep_); }

template<class Scalar>
void IntegratorObserverLogging<Scalar>::observeBeforeTakeStep()
{ logCall(nameObserveBeforeTakeStep_); }

template<class Scalar>
void IntegratorObserverLogging<Scalar>::observeAfterTakeStep()
{ logCall(nameObserveAfterTakeStep_); }

template<class Scalar>
void IntegratorObserverLogging<Scalar>::
observeAcceptedTimeStep(Status & integratorStatus)
{ logCall(nameObserveAcceptedTimeStep_); }

template<class Scalar>
void IntegratorObserverLogging<Scalar>::
observeEndIntegrator(const Status integratorStatus)
{ logCall(nameObserveEndIntegrator_); }

template<class Scalar>
void IntegratorObserverLogging<Scalar>::resetLogCounters()
{
  (*counters_)[nameObserveStartIntegrator_ ] = 0;
  (*counters_)[nameObserveStartTimeStep_   ] = 0;
  (*counters_)[nameObserveNextTimeStep_    ] = 0;
  (*counters_)[nameObserveBeforeTakeStep_  ] = 0;
  (*counters_)[nameObserveAfterTakeStep_   ] = 0;
  (*counters_)[nameObserveAcceptedTimeStep_] = 0;
  (*counters_)[nameObserveEndIntegrator_   ] = 0;
  order_->clear();
}

template<class Scalar>
Teuchos::RCP<const std::map<std::string,int> >
IntegratorObserverLogging<Scalar>::getCounters()
{ return counters_; }

template<class Scalar>
Teuchos::RCP<const std::list<std::string> >
IntegratorObserverLogging<Scalar>::getOrder()
{ return order_; }

template<class Scalar>
void IntegratorObserverLogging<Scalar>::logCall(const std::string call) const
{
  (*counters_)[call] += 1;
  order_->push_back(call);
}

} // namespace Tempus
#endif // Tempus_IntegratorObserverLogging_impl_hpp
