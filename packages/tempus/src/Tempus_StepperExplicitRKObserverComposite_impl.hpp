// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicitRKObserverComposite_impl_hpp
#define Tempus_StepperExplicitRKObserverComposite_impl_hpp

#include "Tempus_StepperExplicitRKObserver.hpp"
#include "Tempus_TimeStepControl.hpp"

namespace Tempus {

template<class Scalar>
StepperExplicitRKObserverComposite<Scalar>::StepperExplicitRKObserverComposite(){}

template<class Scalar>
StepperExplicitRKObserverComposite<Scalar>::~StepperExplicitRKObserverComposite(){}

template<class Scalar>
void StepperExplicitRKObserverComposite<Scalar>::
observeBeginTakeStep(Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeBeginTakeStep(sh,stepper);
}

template<class Scalar>
void StepperExplicitRKObserverComposite<Scalar>::observeBeginStage(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperExplicitRK<Scalar> & stepperExplicitRK)
{
  for(auto& o : observers_)
    o->observeBeginStage(sh,stepperExplicitRK);
}

template<class Scalar>
void StepperExplicitRKObserverComposite<Scalar>::observeBeforeExplicit(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperExplicitRK<Scalar> & stepperExplicitRK)
{
  for(auto& o : observers_)
    o->observeBeforeExplicit(sh,stepperExplicitRK);
}

template<class Scalar>
void StepperExplicitRKObserverComposite<Scalar>::observeEndStage(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperExplicitRK<Scalar> & stepperExplicitRK)
{
  for(auto& o : observers_)
    o->observeEndStage(sh,stepperExplicitRK);
}

template<class Scalar>
void StepperExplicitRKObserverComposite<Scalar>::
observeEndTakeStep(Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeEndTakeStep(sh,stepper);
}

template<class Scalar>
void StepperExplicitRKObserverComposite<Scalar>::
addObserver(const Teuchos::RCP<StepperExplicitRKObserver<Scalar> > &observer)
{
  observers_.push_back(observer);
}

template<class Scalar>
void StepperExplicitRKObserverComposite<Scalar>::
clearObservers() { observers_.clear();}

} // namespace Tempus
#endif // Tempus_StepperExplicitRKObserverComposite_impl_hpp
