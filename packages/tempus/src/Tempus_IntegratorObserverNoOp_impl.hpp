// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

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
