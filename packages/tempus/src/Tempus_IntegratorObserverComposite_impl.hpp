//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorObserverComposite_impl_hpp
#define Tempus_IntegratorObserverComposite_impl_hpp

#include "Tempus_IntegratorObserver.hpp"
#include "Tempus_TimeStepControl.hpp"

namespace Tempus {

template <class Scalar>
IntegratorObserverComposite<Scalar>::IntegratorObserverComposite()
{
}

template <class Scalar>
IntegratorObserverComposite<Scalar>::~IntegratorObserverComposite()
{
}

template <class Scalar>
void IntegratorObserverComposite<Scalar>::observeStartIntegrator(
    const Integrator<Scalar>& integrator)
{
  for (auto& o : observers_) o->observeStartIntegrator(integrator);
}

template <class Scalar>
void IntegratorObserverComposite<Scalar>::observeStartTimeStep(
    const Integrator<Scalar>& integrator)
{
  for (auto& o : observers_) o->observeStartTimeStep(integrator);
}

template <class Scalar>
void IntegratorObserverComposite<Scalar>::observeNextTimeStep(
    const Integrator<Scalar>& integrator)
{
  for (auto& o : observers_) o->observeNextTimeStep(integrator);
}

template <class Scalar>
void IntegratorObserverComposite<Scalar>::observeBeforeTakeStep(
    const Integrator<Scalar>& integrator)
{
  for (auto& o : observers_) o->observeBeforeTakeStep(integrator);
}

template <class Scalar>
void IntegratorObserverComposite<Scalar>::observeAfterTakeStep(
    const Integrator<Scalar>& integrator)
{
  for (auto& o : observers_) o->observeAfterTakeStep(integrator);
}

template <class Scalar>
void IntegratorObserverComposite<Scalar>::observeAfterCheckTimeStep(
    const Integrator<Scalar>& integrator)
{
  for (auto& o : observers_) o->observeAfterCheckTimeStep(integrator);
}

template <class Scalar>
void IntegratorObserverComposite<Scalar>::observeEndTimeStep(
    const Integrator<Scalar>& integrator)
{
  for (auto& o : observers_) o->observeEndTimeStep(integrator);
}

template <class Scalar>
void IntegratorObserverComposite<Scalar>::observeEndIntegrator(
    const Integrator<Scalar>& integrator)
{
  for (auto& o : observers_) o->observeEndIntegrator(integrator);
}

template <class Scalar>
void IntegratorObserverComposite<Scalar>::addObserver(
    const Teuchos::RCP<IntegratorObserver<Scalar> >& observer)
{
  observers_.push_back(observer);
}

template <class Scalar>
void IntegratorObserverComposite<Scalar>::clearObservers()
{
  observers_.clear();
}

}  // namespace Tempus
#endif  // Tempus_IntegratorObserverComposite_impl_hpp
