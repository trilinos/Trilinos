//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperFactory_impl_hpp
#define Tempus_StepperFactory_impl_hpp

#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"
#include "Tempus_StepperTrapezoidal.hpp"
#include "Tempus_StepperBDF2.hpp"
#include "Tempus_StepperNewmarkImplicitAForm.hpp"
#include "Tempus_StepperNewmarkImplicitDForm.hpp"
#include "Tempus_StepperNewmarkExplicitAForm.hpp"
#include "Tempus_StepperHHTAlpha.hpp"
#include "Tempus_StepperRKButcherTableau.hpp"
#include "Tempus_StepperIMEX_RK.hpp"
#include "Tempus_StepperIMEX_RK_Partition.hpp"
#include "Tempus_StepperOperatorSplit.hpp"
#include "Tempus_StepperSubcycling.hpp"
#include "Tempus_StepperLeapfrog.hpp"

namespace Tempus {

template <class Scalar>
Teuchos::RCP<Stepper<Scalar> > StepperFactory<Scalar>::createStepper(
    std::string stepperType,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  if (stepperType == "") stepperType = "Forward Euler";
  return this->createStepper(stepperType, Teuchos::null, model);
}

template <class Scalar>
Teuchos::RCP<Stepper<Scalar> > StepperFactory<Scalar>::createStepper(
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  std::string stepperType = "Forward Euler";
  if (stepperPL != Teuchos::null)
    stepperType = stepperPL->get<std::string>("Stepper Type", "Forward Euler");
  return this->createStepper(stepperType, stepperPL, model);
}

template <class Scalar>
Teuchos::RCP<Stepper<Scalar> > StepperFactory<Scalar>::createStepper(
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models)
{
  std::string stepperType = "Operator Split";
  if (stepperPL != Teuchos::null)
    stepperType = stepperPL->get<std::string>("Stepper Type");

  if (stepperType == "Operator Split")
    return createStepperOperatorSplit(models, stepperPL);
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Unknown 'Stepper Type' = " << stepperType);
  }
}

template <class Scalar>
Teuchos::RCP<Stepper<Scalar> > StepperFactory<Scalar>::createStepper(
    std::string stepperType, Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  using Teuchos::rcp;
  if (stepperType == "Forward Euler")
    return createStepperForwardEuler(model, stepperPL);
  else if (stepperType == "Backward Euler")
    return createStepperBackwardEuler(model, stepperPL);
  else if (stepperType == "Trapezoidal Method")
    return createStepperTrapezoidal(model, stepperPL);
  else if (stepperType == "BDF2")
    return createStepperBDF2(model, stepperPL);
  else if (stepperType == "Newmark Implicit a-Form")
    return createStepperNewmarkImplicitAForm(model, stepperPL);
  else if (stepperType == "Newmark Implicit d-Form")
    return createStepperNewmarkImplicitDForm(model, stepperPL);
  else if (stepperType == "Newmark Explicit a-Form")
    return createStepperNewmarkExplicitAForm(model, stepperPL);
  else if (stepperType == "HHT-Alpha")
    return createStepperHHTAlpha(model, stepperPL);
  else if (stepperType == "General ERK")
    return createStepperERK_General(model, stepperPL);
  else if (stepperType == "RK Forward Euler" || stepperType == "RK1")
    return createStepperERK_ForwardEuler(model, stepperPL);
  else if (stepperType == "RK Explicit 4 Stage")
    return createStepperERK_4Stage4thOrder(model, stepperPL);
  else if (stepperType == "RK Explicit 3/8 Rule")
    return createStepperERK_3_8Rule(model, stepperPL);
  else if (stepperType == "RK Explicit 4 Stage 3rd order by Runge")
    return createStepperERK_4Stage3rdOrderRunge(model, stepperPL);
  else if (stepperType == "RK Explicit 5 Stage 3rd order by Kinnmark and Gray")
    return createStepperERK_5Stage3rdOrderKandG(model, stepperPL);
  else if (stepperType == "RK Explicit 3 Stage 3rd order")
    return createStepperERK_3Stage3rdOrder(model, stepperPL);
  else if (stepperType == "RK Explicit 3 Stage 3rd order TVD" ||
           stepperType == "SSPERK33" || stepperType == "SSPRK3")
    return createStepperERK_3Stage3rdOrderTVD(model, stepperPL);
  else if (stepperType == "RK Explicit 3 Stage 3rd order by Heun")
    return createStepperERK_3Stage3rdOrderHeun(model, stepperPL);
  else if (stepperType == "RK Explicit Midpoint")
    return createStepperERK_Midpoint(model, stepperPL);
  else if (stepperType == "RK Explicit Trapezoidal" ||
           stepperType == "Heuns Method" || stepperType == "SSPERK22" ||
           stepperType == "SSPRK2")
    return createStepperERK_Trapezoidal(model, stepperPL);
  else if (stepperType == "RK Explicit Ralston" || stepperType == "RK2")
    return createStepperERK_Ralston(model, stepperPL);
  else if (stepperType == "SSPERK54")
    return createStepperERK_SSPERK54(model, stepperPL);
  else if (stepperType == "Bogacki-Shampine 3(2) Pair")
    return createStepperERK_BogackiShampine32(model, stepperPL);
  else if (stepperType == "Merson 4(5) Pair")
    return createStepperERK_Merson45(model, stepperPL);
  else if (stepperType == "General DIRK")
    return createStepperDIRK_General(model, stepperPL);
  else if (stepperType == "RK Backward Euler")
    return createStepperDIRK_BackwardEuler(model, stepperPL);
  else if (stepperType == "SDIRK 2 Stage 2nd order")
    return createStepperSDIRK_2Stage2ndOrder(model, stepperPL);
  else if (stepperType == "SSPDIRK22")
    return createStepperSDIRK_SSPDIRK22(model, stepperPL);
  else if (stepperType == "SDIRK 3 Stage 2nd order")
    return createStepperSDIRK_3Stage2ndOrder(model, stepperPL);
  else if (stepperType == "SSPDIRK32")
    return createStepperSDIRK_SSPDIRK32(model, stepperPL);
  else if (stepperType == "SSPDIRK23")
    return createStepperSDIRK_SSPDIRK23(model, stepperPL);
  else if (stepperType == "SSPDIRK33")
    return createStepperSDIRK_SSPDIRK33(model, stepperPL);
  else if (stepperType == "SDIRK 2 Stage 3rd order")
    return createStepperSDIRK_2Stage3rdOrder(model, stepperPL);
  else if (stepperType == "EDIRK 2 Stage 3rd order")
    return createStepperEDIRK_2Stage3rdOrder(model, stepperPL);
  else if (stepperType == "DIRK 1 Stage Theta Method")
    return createStepperDIRK_1StageTheta(model, stepperPL);
  else if (stepperType == "EDIRK 2 Stage Theta Method")
    return createStepperEDIRK_2StageTheta(model, stepperPL);
  else if (stepperType == "RK Trapezoidal Rule" ||
           stepperType == "RK Crank-Nicolson")
    return createStepperEDIRK_TrapezoidalRule(model, stepperPL);
  else if (stepperType == "RK Implicit Midpoint")
    return createStepperSDIRK_ImplicitMidpoint(model, stepperPL);
  else if (stepperType == "RK Implicit 1 Stage 1st order Radau IA")
    return createStepperDIRK_1Stage1stOrderRadauIA(model, stepperPL);
  else if (stepperType == "RK Implicit 2 Stage 2nd order Lobatto IIIB")
    return createStepperDIRK_2Stage2ndOrderLobattoIIIB(model, stepperPL);
  else if (stepperType == "SDIRK 5 Stage 4th order")
    return createStepperSDIRK_5Stage4thOrder(model, stepperPL);
  else if (stepperType == "SDIRK 3 Stage 4th order")
    return createStepperSDIRK_3Stage4thOrder(model, stepperPL);
  else if (stepperType == "SDIRK 5 Stage 5th order")
    return createStepperSDIRK_5Stage5thOrder(model, stepperPL);
  else if (stepperType == "SDIRK 2(1) Pair")
    return createStepperSDIRK_21Pair(model, stepperPL);
  else if (stepperType == "IMEX RK 1st order" || stepperType == "SSP1_111" ||
           stepperType == "IMEX RK SSP2" || stepperType == "IMEX RK SSP3" ||
           stepperType == "SSP3_332" || stepperType == "SSP2_222" ||
           stepperType == "SSP2_222_L" || stepperType == "SSP2_222_A" ||
           stepperType == "IMEX RK ARS 233" || stepperType == "ARS 233" ||
           stepperType == "General IMEX RK")
    return createStepperIMEX_RK(model, stepperType, stepperPL);
  else if (stepperType == "Partitioned IMEX RK 1st order" ||
           stepperType == "Partitioned IMEX RK SSP2" ||
           stepperType == "Partitioned IMEX RK ARS 233" ||
           stepperType == "General Partitioned IMEX RK")
    return createStepperIMEX_RK_Partition(model, stepperType, stepperPL);
  else if (stepperType == "Leapfrog")
    return createStepperLeapfrog(model, stepperPL);
  else if (stepperType == "Subcycling")
    return createStepperSubcycling(model, stepperPL);
  else {
    Teuchos::RCP<Teuchos::FancyOStream> out =
        Teuchos::VerboseObjectBase::getDefaultOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 1, "StepperFactory::createStepper");
    *out << "Unknown Stepper Type!  ('" + stepperType + "').\n"
         << "Here is a list of available Steppers.\n"
         << "  One-Step Methods:\n"
         << "    'Forward Euler'\n"
         << "    'Backward Euler'\n"
         << "    'Trapezoidal Method'\n"
         << "  Multi-Step Methods:\n"
         << "    'BDF2'\n"
         << "  Second-order PDE Methods:\n"
         << "    'Leapfrog'\n"
         << "    'Newmark Implicit a-Form'\n"
         << "    'Newmark Implicit d-Form'\n"
         << "    'Newmark Explicit a-Form'\n"
         << "    'HHT-Alpha'\n"
         << "  Explicit Runge-Kutta Methods:\n"
         << "    'RK Forward Euler (RK1)'\n"
         << "    'RK Explicit 4 Stage'\n"
         << "    'RK Explicit 3/8 Rule'\n"
         << "    'RK Explicit 4 Stage 3rd order by Runge'\n"
         << "    'RK Explicit 5 Stage 3rd order by Kinnmark and Gray'\n"
         << "    'RK Explicit 3 Stage 3rd order'\n"
         << "    'RK Explicit 3 Stage 3rd order TVD'\n"
         << "    'RK Explicit 3 Stage 3rd order by Heun'\n"
         << "    'RK Explicit Midpoint'\n"
         << "    'RK Explicit Trapezoidal' or 'Heuns Method'\n"
         << "    'Bogacki-Shampine 3(2) Pair'\n"
         << "    'SSPERK22 (SSPRK2)'\n"
         << "    'SSPERK33 (SSPRK3)'\n"
         << "    'SSPERK54'\n"
         << "    'General ERK'\n"
         << "  Implicit Runge-Kutta Methods:\n"
         << "    'RK Backward Euler'\n"
         << "    'DIRK 1 Stage Theta Method'\n"
         << "    'RK Implicit Midpoint'\n"
         << "    'SDIRK 1 Stage 1st order'\n"
         << "    'SDIRK 2 Stage 2nd order'\n"
         << "    'SDIRK 2 Stage 3rd order'\n"
         << "    'EDIRK 2 Stage 3rd order'\n"
         << "    'EDIRK 2 Stage Theta Method'\n"
         << "    'SDIRK 3 Stage 4th order'\n"
         << "    'SDIRK 5 Stage 4th order'\n"
         << "    'SDIRK 5 Stage 5th order'\n"
         << "    'SSPDIRK22'\n"
         << "    'SSPDIRK32'\n"
         << "    'SSPDIRK23'\n"
         << "    'SSPDIRK33'\n"
         << "    'SDIRK 2(1) Pair'\n"
         << "    'SDIRK 3 Stage 2nd order'\n"
         << "    'RK Trapezoidal Rule' or 'RK Crank-Nicolson'\n"
         << "    'General DIRK'\n"
         << "  Implicit-Explicit (IMEX) Methods:\n"
         << "    'IMEX RK 1st order'\n"
         << "    'SSP1_111'\n"
         << "    'IMEX RK SSP2 (SSP2_222)'\n"
         << "    'IMEX RK SSP3 (SSP3_332)'\n"
         << "    'IMEX RK ARS 233'\n"
         << "    'General IMEX RK'\n"
         << "    'Partitioned IMEX RK 1st order'\n"
         << "    'Partitioned IMEX RK SSP2'\n"
         << "    'Partitioned IMEX RK ARS 233'\n"
         << "    'General Partitioned IMEX RK'\n"
         << "  Steppers with subSteppers:\n"
         << "    'Operator Split'\n"
         << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Unknown 'Stepper Type' = " << stepperType);
  }
}

}  // namespace Tempus
#endif  // Tempus_StepperFactory_impl_hpp
