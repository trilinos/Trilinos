// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKObserverLogging_impl_hpp
#define Tempus_StepperRKObserverLogging_impl_hpp

#include "Tempus_StepperRKObserver.hpp"
#include "Tempus_TimeStepControl.hpp"

namespace Tempus {

template<class Scalar>
StepperRKObserverLogging<Scalar>::StepperRKObserverLogging()
  : nameObserveBeginTakeStep_            ( "observeBeginTakeStep"),
    nameObserveBeginStage_               ( "observeBeginStage"),
    nameObserveBeforeImplicitExplicitly_ ( "observeBeforeImplicitExplicitly"),
    nameObserveBeforeSolve_              ( "observeBeforeSolve"),
    nameObserveAfterSolve_               ( "observeAfterSolve"),
    nameObserveBeforeExplicit_           ( "observeBeforeExplicit"),
    nameObserveEndStage_                 ( "observeEndStage"),
    nameObserveEndTakeStep_              ( "observeEndTakeStep")
{
  counters_ = Teuchos::rcp(new std::map<std::string,int>);
  order_ = Teuchos::rcp(new std::list<std::string>);
  this->resetLogCounters();
}

template<class Scalar>
StepperRKObserverLogging<Scalar>::~StepperRKObserverLogging(){}

template<class Scalar>
void StepperRKObserverLogging<Scalar>::
observeBeginTakeStep(
        Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
        Stepper<Scalar> & /* stepper */)
{ logCall(nameObserveBeginTakeStep_); }

template<class Scalar>
void StepperRKObserverLogging<Scalar>::
observeBeginStage(
        Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
        Stepper<Scalar> & /* stepper */)
{ logCall(nameObserveBeginStage_); }

template<class Scalar>
void StepperRKObserverLogging<Scalar>::
observeBeforeImplicitExplicitly(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */)
{ logCall(nameObserveBeforeImplicitExplicitly_); }

template<class Scalar>
void StepperRKObserverLogging<Scalar>::
observeBeforeSolve(
        Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
        Stepper<Scalar> & /* stepper */)
{ logCall(nameObserveBeforeSolve_); }

template<class Scalar>
void StepperRKObserverLogging<Scalar>::
observeAfterSolve(
        Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
        Stepper<Scalar> & /* stepper */)
{ logCall(nameObserveAfterSolve_); }

template<class Scalar>
void StepperRKObserverLogging<Scalar>::
observeBeforeExplicit(
        Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
        Stepper<Scalar> & /* stepper */)
{ logCall(nameObserveBeforeExplicit_); }

template<class Scalar>
void StepperRKObserverLogging<Scalar>::
observeEndStage(
        Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
        Stepper<Scalar> & /* stepper */)
{ logCall(nameObserveEndStage_); }

template<class Scalar>
void StepperRKObserverLogging<Scalar>::
observeEndTakeStep(
        Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
        Stepper<Scalar> & /* stepper */)
{ logCall(nameObserveEndTakeStep_); }

template<class Scalar>
void StepperRKObserverLogging<Scalar>::resetLogCounters()
{
  (*counters_)[nameObserveBeginTakeStep_            ] = 0;
  (*counters_)[nameObserveBeginStage_               ] = 0;
  (*counters_)[nameObserveBeforeImplicitExplicitly_ ] = 0;
  (*counters_)[nameObserveBeforeSolve_              ] = 0;
  (*counters_)[nameObserveAfterSolve_               ] = 0;
  (*counters_)[nameObserveBeforeExplicit_           ] = 0;
  (*counters_)[nameObserveEndStage_                 ] = 0;
  (*counters_)[nameObserveEndTakeStep_              ] = 0;
  order_->clear();
}

template<class Scalar>
Teuchos::RCP<const std::map<std::string,int> >
StepperRKObserverLogging<Scalar>::getCounters()
{ return counters_; }

template<class Scalar>
Teuchos::RCP<const std::list<std::string> >
StepperRKObserverLogging<Scalar>::getOrder()
{ return order_; }

template<class Scalar>
void StepperRKObserverLogging<Scalar>::logCall(const std::string call) const
{
  (*counters_)[call] += 1;
  order_->push_back(call);
}

} // namespace Tempus
#endif // Tempus_StepperRKObserverLogging_impl_hpp
