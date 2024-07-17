//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorObserverNoOp_impl_hpp
#define Tempus_IntegratorObserverNoOp_impl_hpp

#include "Tempus_Stepper.hpp"

namespace Tempus {

template <class Scalar>
IntegratorObserverNoOp<Scalar>::IntegratorObserverNoOp()
{
}

template <class Scalar>
IntegratorObserverNoOp<Scalar>::~IntegratorObserverNoOp()
{
}

template <class Scalar>
void IntegratorObserverNoOp<Scalar>::observeStartIntegrator(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverNoOp<Scalar>::observeStartTimeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverNoOp<Scalar>::observeNextTimeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverNoOp<Scalar>::observeBeforeTakeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverNoOp<Scalar>::observeAfterTakeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverNoOp<Scalar>::observeAfterCheckTimeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverNoOp<Scalar>::observeEndTimeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverNoOp<Scalar>::observeEndIntegrator(
    const Integrator<Scalar>& /* integrator */)
{
}

}  // namespace Tempus
#endif  // Tempus_IntegratorObserverNoOp_impl_hpp
