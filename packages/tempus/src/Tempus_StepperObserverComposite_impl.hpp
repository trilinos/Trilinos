// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperObserverComposite_impl_hpp
#define Tempus_StepperObserverComposite_impl_hpp

#include "Tempus_StepperObserver.hpp"
#include "Tempus_TimeStepControl.hpp"

namespace Tempus {

template<class Scalar>
StepperObserverComposite<Scalar>::StepperObserverComposite(){}

template<class Scalar>
StepperObserverComposite<Scalar>::~StepperObserverComposite(){}

template<class Scalar>
void StepperObserverComposite<Scalar>::
observeBeginTakeStep(Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeBeginTakeStep(sh,stepper);
}

template<class Scalar>
void StepperObserverComposite<Scalar>::
observeEndTakeStep(Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper)
{
  for(auto& o : observers_)
    o->observeEndTakeStep(sh,stepper);
}

template<class Scalar>
void StepperObserverComposite<Scalar>::
addObserver(const Teuchos::RCP<StepperObserver<Scalar> > &observer)
{
  observers_.push_back(observer);
}

template<class Scalar>
void StepperObserverComposite<Scalar>::
clearObservers() { observers_.clear();}

} // namespace Tempus
#endif // Tempus_StepperObserverComposite_impl_hpp
