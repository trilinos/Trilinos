// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorObserverBasic_impl_hpp
#define Tempus_IntegratorObserverBasic_impl_hpp

#include "Tempus_Stepper.hpp"

namespace Tempus {

template<class Scalar>
IntegratorObserverBasic<Scalar>::IntegratorObserverBasic(){}

template<class Scalar>
IntegratorObserverBasic<Scalar>::~IntegratorObserverBasic(){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeStartIntegrator(const Integrator<Scalar>& integrator){

  std::time_t begin = std::time(nullptr);
  const Teuchos::RCP<Teuchos::FancyOStream> out = integrator.getOStream();
  out->setOutputToRootOnly(0);
  Teuchos::OSTab ostab(out,0,"ScreenOutput");
  *out << "\nTempus - IntegratorBasic\n"
       << std::asctime(std::localtime(&begin)) << "\n"
       << "  Stepper = " << integrator.getStepper()->description() << "\n"
       << "  Simulation Time Range  [" << integrator.getTimeStepControl()->getInitTime()
       << ", " << integrator.getTimeStepControl()->getFinalTime() << "]\n"
       << "  Simulation Index Range [" << integrator.getTimeStepControl()->getInitIndex()
       << ", " << integrator.getTimeStepControl()->getFinalIndex() << "]\n"
       << "============================================================================\n"
       << "  Step       Time         dt  Abs Error  Rel Error  Order  nFail  dCompTime"
       << std::endl;
}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeStartTimeStep(const Integrator<Scalar>& /* integrator */){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeNextTimeStep(const Integrator<Scalar>& /* integrator */){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeBeforeTakeStep(const Integrator<Scalar>& /* integrator */){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeAfterTakeStep(const Integrator<Scalar>& /* integrator */){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeAfterCheckTimeStep(const Integrator<Scalar>& /* integrator */){}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeEndTimeStep(const Integrator<Scalar>& integrator){

  using Teuchos::RCP;
  auto cs = integrator.getSolutionHistory()->getCurrentState();

  if ((cs->getOutputScreen() == true) or
      (cs->getOutput() == true) or
      (cs->getTime() == integrator.getTimeStepControl()->getFinalTime())) {

     const Scalar steppertime = integrator.getStepperTimer()->totalElapsedTime();
     // reset the stepper timer
     integrator.getStepperTimer()->reset();

     const Teuchos::RCP<Teuchos::FancyOStream> out = integrator.getOStream();
     out->setOutputToRootOnly(0);
     Teuchos::OSTab ostab(out,0,"ScreenOutput");
     *out<<std::scientific
        <<std::setw( 6)<<std::setprecision(3)<<cs->getIndex()
        <<std::setw(11)<<std::setprecision(3)<<cs->getTime()
        <<std::setw(11)<<std::setprecision(3)<<cs->getTimeStep()
        <<std::setw(11)<<std::setprecision(3)<<cs->getErrorAbs()
        <<std::setw(11)<<std::setprecision(3)<<cs->getErrorRel()
        <<std::fixed     <<std::setw( 7)<<std::setprecision(1)<<cs->getOrder()
        <<std::scientific<<std::setw( 7)<<std::setprecision(3)<<cs->getNFailures()
        <<std::setw(11)<<std::setprecision(3)<<steppertime
        <<std::endl;
  }

}

template<class Scalar>
void IntegratorObserverBasic<Scalar>::
observeEndIntegrator(const Integrator<Scalar>& integrator){

  std::string exitStatus;
  //const Scalar runtime = integrator.getIntegratorTimer()->totalElapsedTime();
  if (integrator.getSolutionHistory()->getCurrentState()->getSolutionStatus() ==
      Status::FAILED or integrator.getStatus() == Status::FAILED) {
    exitStatus = "Time integration FAILURE!";
  } else {
    exitStatus = "Time integration complete.";
  }
  std::time_t end = std::time(nullptr);
  const Scalar runtime = integrator.getIntegratorTimer()->totalElapsedTime();
  const Teuchos::RCP<Teuchos::FancyOStream> out = integrator.getOStream();
  out->setOutputToRootOnly(0);
  Teuchos::OSTab ostab(out,0,"ScreenOutput");
  *out << "============================================================================\n"
       << "  Total runtime = " << runtime << " sec = "
       << runtime/60.0 << " min\n"
       << std::asctime(std::localtime(&end))
       << "\nNumber of Accepted Steps = " << integrator.getSolutionHistory()->getCurrentState()->getIndex()
       << "\nNumber of Failures = " << integrator.getSolutionHistory()->getCurrentState()->getMetaData()->getNRunningFailures()
       << "\n"
       << exitStatus << "\n"
       << std::endl;
}

} // namespace Tempus
#endif // Tempus_IntegratorObserverBasic_impl_hpp
