//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorObserverSubcycling_impl_hpp
#define Tempus_IntegratorObserverSubcycling_impl_hpp

#include "Tempus_Stepper.hpp"

namespace Tempus {

template <class Scalar>
IntegratorObserverSubcycling<Scalar>::IntegratorObserverSubcycling()
{
}

template <class Scalar>
IntegratorObserverSubcycling<Scalar>::~IntegratorObserverSubcycling()
{
}

template <class Scalar>
void IntegratorObserverSubcycling<Scalar>::observeStartIntegrator(
    const Integrator<Scalar>& integrator)
{
  const Teuchos::RCP<Teuchos::FancyOStream> out = integrator.getOStream();
  out->setOutputToRootOnly(0);
  Teuchos::OSTab ostab(out, 0, "ScreenOutput");
  *out << "\n    Begin Subcycling "
          "-------------------------------------------------------\n";
  // << "  Step       Time         dt  Abs Error  Rel Error  Order  nFail
  // dCompTime"
  // << std::endl;
}

template <class Scalar>
void IntegratorObserverSubcycling<Scalar>::observeStartTimeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverSubcycling<Scalar>::observeNextTimeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverSubcycling<Scalar>::observeBeforeTakeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverSubcycling<Scalar>::observeAfterTakeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverSubcycling<Scalar>::observeAfterCheckTimeStep(
    const Integrator<Scalar>& /* integrator */)
{
}

template <class Scalar>
void IntegratorObserverSubcycling<Scalar>::observeEndTimeStep(
    const Integrator<Scalar>& integrator)
{
  using Teuchos::RCP;
  auto cs = integrator.getSolutionHistory()->getCurrentState();

  if ((cs->getOutputScreen() == true) || (cs->getOutput() == true) ||
      (cs->getTime() == integrator.getTimeStepControl()->getFinalTime())) {
    const Scalar steppertime = integrator.getStepperTimer()->totalElapsedTime();
    // reset the stepper timer
    integrator.getStepperTimer()->reset();

    const Teuchos::RCP<Teuchos::FancyOStream> out = integrator.getOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 0, "ScreenOutput");
    *out << std::scientific << std::setw(6) << std::setprecision(3)
         << cs->getIndex() << std::setw(11) << std::setprecision(3)
         << cs->getTime() << std::setw(11) << std::setprecision(3)
         << cs->getTimeStep() << std::setw(11) << std::setprecision(3)
         << cs->getErrorAbs() << std::setw(11) << std::setprecision(3)
         << cs->getErrorRel() << std::fixed << std::setw(7)
         << std::setprecision(1) << cs->getOrder() << std::scientific
         << std::setw(7) << std::setprecision(3) << cs->getNFailures()
         << std::setw(11) << std::setprecision(3) << steppertime << std::endl;
  }
}

template <class Scalar>
void IntegratorObserverSubcycling<Scalar>::observeEndIntegrator(
    const Integrator<Scalar>& integrator)
{
  const Teuchos::RCP<Teuchos::FancyOStream> out = integrator.getOStream();
  out->setOutputToRootOnly(0);
  Teuchos::OSTab ostab(out, 0, "ScreenOutput");
  *out << "    End Subcycling "
          "---------------------------------------------------------\n\n";
}

}  // namespace Tempus
#endif  // Tempus_IntegratorObserverSubcycling_impl_hpp
