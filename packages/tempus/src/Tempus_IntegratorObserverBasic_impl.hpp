// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorObserverBasic_impl_hpp
#define Tempus_IntegratorObserverBasic_impl_hpp


namespace Tempus {

template<class Scalar>
IntegratorObserverBasic<Scalar>::IntegratorObserverBasic(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory,
  const Teuchos::RCP<TimeStepControl<Scalar> >& timeStepControl)
  : solutionHistory_(solutionHistory), timeStepControl_(timeStepControl)
{}

template<class Scalar>
IntegratorObserverBasic<Scalar>::~IntegratorObserverBasic(){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::observeStartIntegrator(){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::observeStartTimeStep(){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeNextTimeStep(Status & integratorStatus){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::observeBeforeTakeStep(){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::observeAfterTakeStep(){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeAcceptedTimeStep(Status & integratorStatus){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeEndIntegrator(const Status integratorStatus){}

} // namespace Tempus
#endif // Tempus_IntegratorObserverBasic_impl_hpp
