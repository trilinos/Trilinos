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
IntegratorObserverBasic<Scalar>::IntegratorObserverBasic(){}

template<class Scalar>
IntegratorObserverBasic<Scalar>::~IntegratorObserverBasic(){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeStartIntegrator(const Integrator<Scalar>& ){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeStartTimeStep(const Integrator<Scalar>& ){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeNextTimeStep(const Integrator<Scalar>& ){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeBeforeTakeStep(const Integrator<Scalar>& ){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeAfterTakeStep(const Integrator<Scalar>& ){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeAcceptedTimeStep(const Integrator<Scalar>& ){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeEndIntegrator(const Integrator<Scalar>& ){}

} // namespace Tempus
#endif // Tempus_IntegratorObserverBasic_impl_hpp
