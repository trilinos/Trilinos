// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorObserverComposite_impl_hpp
#define Tempus_IntegratorObserverComposite_impl_hpp

#include "Tempus_IntegratorObserver.hpp"
#include "Tempus_TimeStepControl.hpp"

namespace Tempus {

template<class Scalar>
IntegratorObserverComposite<Scalar>::IntegratorObserverComposite(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory,
  const Teuchos::RCP<TimeStepControl<Scalar> >& timeStepControl)
  : solutionHistory_(solutionHistory), timeStepControl_(timeStepControl){}

template<class Scalar>
IntegratorObserverComposite<Scalar>::~IntegratorObserverComposite(){}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::observeStartIntegrator()
{ 
  for(auto& o : observers_)
    o->observeStartIntegrator();
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::observeStartTimeStep()
{
  for(auto& o : observers_)
    o->observeStartTimeStep();
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeNextTimeStep(Status & integratorStatus)
{
  for(auto& o : observers_)
    o->observeNextTimeStep(integratorStatus);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::observeBeforeTakeStep()
{
  for(auto& o : observers_)
    o->observeBeforeTakeStep();
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::observeAfterTakeStep()
{
  for(auto& o : observers_)
    o->observeAfterTakeStep();
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeAcceptedTimeStep(Status & integratorStatus)
{
  for(auto& o : observers_)
    o->observeAcceptedTimeStep(integratorStatus);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeEndIntegrator(const Status integratorStatus)
{
  for(auto& o : observers_)  
    o->observeEndIntegrator(integratorStatus);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
setSolutionHistory(Teuchos::RCP<SolutionHistory<Scalar> > sh)
{
  solutionHistory_ = sh;
  for (auto& o : observers_)
    o->setSolutionHistory(sh);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
setTimeStepControl(Teuchos::RCP<TimeStepControl<Scalar> > tsc)
{
  timeStepControl_ = tsc;
  for (auto& o : observers_)
    o->setTimeStepControl(tsc);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
addObserver(const Teuchos::RCP<IntegratorObserver<Scalar> > &observer)
{
    observer->setSolutionHistory(solutionHistory_);
    observer->setTimeStepControl(timeStepControl_);
    observers_.push_back(observer);
}

} // namespace Tempus
#endif // Tempus_IntegratorObserverComposite_impl_hpp
