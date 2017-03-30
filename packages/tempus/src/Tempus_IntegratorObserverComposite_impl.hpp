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
    using Teuchos::as;

    for( int i = 0; i < as<int>(observers_.size()); i++)
        observers_[i]->observeStartIntegrator();
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::observeStartTimeStep()
{
    using Teuchos::as;

    for( int i = 0; i < as<int>(observers_.size()); i++)
        observers_[i]->observeStartTimeStep();
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeNextTimeStep(Status & integratorStatus)
{ 
    using Teuchos::as;

    for( int i = 0; i < as<int>(observers_.size()); i++)
        observers_[i]->observeNextTimeStep( integratorStatus);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::observeBeforeTakeStep()
{ 
    using Teuchos::as;

    for( int i = 0; i < as<int>(observers_.size()); i++)
        observers_[i]->observeBeforeTakeStep();
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::observeAfterTakeStep()
{ 
    using Teuchos::as;

    for( int i = 0; i < as<int>(observers_.size()); i++)
        observers_[i]->observeAfterTakeStep();
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeAcceptedTimeStep(Status & integratorStatus)
{ 
    using Teuchos::as;

    for( int i = 0; i < as<int>(observers_.size()); i++)
        observers_[i]->observeAcceptedTimeStep( integratorStatus);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeEndIntegrator(const Status integratorStatus)
{ 
    using Teuchos::as;

    for( int i = 0; i < as<int>(observers_.size()); i++)
        observers_[i]->observeEndIntegrator( integratorStatus);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
setSolutionHistory(Teuchos::RCP<SolutionHistory<Scalar> > sh)
{ solutionHistory_ = sh; return; }

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
setTimeStepControl(Teuchos::RCP<TimeStepControl<Scalar> > tsc)
{ timeStepControl_ = tsc; return; }

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
addObserver(const Teuchos::RCP<IntegratorObserver<Scalar> > &observer)
{
    observers_.push_back(observer);
}

} // namespace Tempus
#endif // Tempus_IntegratorObserverComposite_impl_hpp
