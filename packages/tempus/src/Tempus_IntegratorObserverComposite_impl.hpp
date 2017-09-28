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
IntegratorObserverComposite<Scalar>::IntegratorObserverComposite(){}

template<class Scalar>
IntegratorObserverComposite<Scalar>::~IntegratorObserverComposite(){}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeStartIntegrator(const Integrator<Scalar>& integrator)
{ 
  for(auto& o : observers_)
    o->observeStartIntegrator(integrator);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeStartTimeStep(const Integrator<Scalar>& integrator)
{
  for(auto& o : observers_)
    o->observeStartTimeStep(integrator);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeNextTimeStep(const Integrator<Scalar>& integrator)
{
  for(auto& o : observers_)
    o->observeNextTimeStep(integrator);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeBeforeTakeStep(const Integrator<Scalar>& integrator)
{
  for(auto& o : observers_)
    o->observeBeforeTakeStep(integrator);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeAfterTakeStep(const Integrator<Scalar>& integrator)
{
  for(auto& o : observers_)
    o->observeAfterTakeStep(integrator);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeAcceptedTimeStep(const Integrator<Scalar>& integrator)
{
  for(auto& o : observers_)
    o->observeAcceptedTimeStep(integrator);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
observeEndIntegrator(const Integrator<Scalar>& integrator)
{
  for(auto& o : observers_)  
    o->observeEndIntegrator(integrator);
}

template<class Scalar>
void IntegratorObserverComposite<Scalar>::
addObserver(const Teuchos::RCP<IntegratorObserver<Scalar> > &observer)
{
  observers_.push_back(observer);
}

} // namespace Tempus
#endif // Tempus_IntegratorObserverComposite_impl_hpp
