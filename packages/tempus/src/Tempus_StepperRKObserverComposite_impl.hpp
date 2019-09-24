// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKObserverComposite_impl_hpp
#define Tempus_StepperRKObserverComposite_impl_hpp

#include "Tempus_StepperRKObserver.hpp"
#include "Tempus_TimeStepControl.hpp"

namespace Tempus {

template<class Scalar>
StepperRKObserverComposite<Scalar>::StepperRKObserverComposite(){}

template<class Scalar>
StepperRKObserverComposite<Scalar>::~StepperRKObserverComposite(){}

template<class Scalar>
void StepperRKObserverComposite<Scalar>::
observeBeginTakeStep(Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeBeginTakeStep(sh,stepper);
}

template<class Scalar>
void StepperRKObserverComposite<Scalar>::observeBeginStage(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeBeginStage(sh,stepper);
}

template<class Scalar>
void StepperRKObserverComposite<Scalar>::observeBeforeImplicitExplicitly(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeBeforeImplicitExplicitly(sh,stepper);
}

template<class Scalar>
void StepperRKObserverComposite<Scalar>::observeBeforeSolve(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeBeforeSolve(sh,stepper);
}

template<class Scalar>
void StepperRKObserverComposite<Scalar>::observeAfterSolve(
        Teuchos::RCP<SolutionHistory<Scalar> >  sh,
        Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeAfterSolve(sh,stepper);
}

template<class Scalar>
void StepperRKObserverComposite<Scalar>::observeBeforeExplicit(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeBeforeExplicit(sh,stepper);
}
template<class Scalar>
void StepperRKObserverComposite<Scalar>::observeEndStage(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeEndStage(sh,stepper);
}

template<class Scalar>
void StepperRKObserverComposite<Scalar>::
observeEndTakeStep(Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeEndTakeStep(sh,stepper);
}

template<class Scalar>
void StepperRKObserverComposite<Scalar>::
addObserver(const Teuchos::RCP<StepperRKObserver<Scalar> > &observer)
{
  observers_.push_back(observer);
}

template<class Scalar>
void StepperRKObserverComposite<Scalar>::
clearObservers() { observers_.clear();}

} // namespace Tempus
#endif // Tempus_StepperRKObserverComposite_impl_hpp
