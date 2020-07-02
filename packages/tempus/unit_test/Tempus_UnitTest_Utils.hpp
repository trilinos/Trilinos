// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_UnitTest_Utils_hpp
#define Tempus_UnitTest_Utils_hpp

#include "Tempus_StepperFactory.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperRKModifierBase.hpp"
#include "Tempus_StepperRKModifierXBase.hpp"
#include "Tempus_StepperRKObserverBase.hpp"
#include "Tempus_StepperRKModifierDefault.hpp"
#include "Tempus_StepperRKModifierXDefault.hpp"
#include "Tempus_StepperRKObserverDefault.hpp"
#include "Tempus_StepperRKAppActionComposite.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ParameterList;

using Tempus::StepperFactory;

/** \brief Unit test utility for Stepper construction through StepperFactory.
 */
void testFactoryConstruction(std::string stepperType,
  const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model)
{
  RCP<StepperFactory<double> > sf = Teuchos::rcp(new StepperFactory<double>());

  // Test using stepperType
  // Passing in model.
  auto stepper = sf->createStepper(stepperType, model);
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  // With setting model.
  stepper = sf->createStepper(stepperType);
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test using ParameterList.
  // Passing in model.
  auto stepperPL = rcp_const_cast<ParameterList>(stepper->getValidParameters());
  stepper = sf->createStepper(stepperPL, model);
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  // With setting model.
  stepper = sf->createStepper(stepperPL);
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
}


/** \brief Unit test utility for ExplicitRK Stepper construction and accessors.
 */
void testExplicitRKAccessorsFullConstruction(
  const RCP<Tempus::StepperExplicitRK<double> >& stepper)
{
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto modifier  = rcp(new Tempus::StepperRKModifierDefault<double>());
  auto modifierX = rcp(new Tempus::StepperRKModifierXDefault<double>());
  auto observer  = rcp(new Tempus::StepperRKObserverDefault<double>());
  bool useFSAL              = stepper->getUseFSALDefault();
  std::string ICConsistency = stepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheckDefault();
  bool useEmbedded          = stepper->getUseEmbeddedDefault();

  // Test the set functions.
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  auto obs    = rcp(new Tempus::StepperRKObserverComposite<double>());
  stepper->setObserver(obs);                           stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
  stepper->setAppAction(modifier);                     stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAppAction(modifierX);                    stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAppAction(observer);                     stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseEmbedded(useEmbedded);                stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  std::string stepperType = stepper->getStepperType();
  // Full argument list construction.
  if        (stepperType == "RK Explicit 3 Stage 3rd order") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_3Stage3rdOrder<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_3Stage3rdOrder<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Explicit 3 Stage 3rd order by Heun") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_3Stage3rdOrderHeun<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_3Stage3rdOrderHeun<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Explicit 3 Stage 3rd order TVD") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_3Stage3rdOrderTVD<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_3Stage3rdOrderTVD<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Explicit 3/8 Rule") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_3_8Rule<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_3_8Rule<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Explicit 4 Stage 3rd order by Runge") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_4Stage3rdOrderRunge<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_4Stage3rdOrderRunge<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Explicit 4 Stage") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_4Stage4thOrder<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_4Stage4thOrder<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Explicit 5 Stage 3rd order by Kinnmark and Gray") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_5Stage3rdOrderKandG<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_5Stage3rdOrderKandG<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "Bogacki-Shampine 3(2) Pair") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_BogackiShampine32<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_BogackiShampine32<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Forward Euler") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_ForwardEuler<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_ForwardEuler<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "Merson 4(5) Pair") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_Merson45<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_Merson45<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Explicit Midpoint") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_Midpoint<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_Midpoint<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Explicit Trapezoidal") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperERK_Trapezoidal<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperERK_Trapezoidal<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  }
  else TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error - unknown stepperType = "+stepperType);
}


/** \brief Unit test utility for ExplicitRK Stepper construction and accessors.
 */
void testDIRKAccessorsFullConstruction(
  const RCP<Tempus::StepperDIRK<double> >& stepper)
{
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto modifier  = rcp(new Tempus::StepperRKModifierDefault<double>());
  auto modifierX = rcp(new Tempus::StepperRKModifierXDefault<double>());
  auto observer  = rcp(new Tempus::StepperRKObserverDefault<double>());
  auto solver    = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSALDefault();
  std::string ICConsistency = stepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheckDefault();
  bool useEmbedded          = stepper->getUseEmbeddedDefault();
  bool zeroInitialGuess     = stepper->getZeroInitialGuess();

  // Test the set functions.
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  auto obs    = rcp(new Tempus::StepperRKObserverComposite<double>());
  stepper->setObserver(obs);                           stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
  stepper->setAppAction(modifier);                     stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAppAction(modifierX);                    stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAppAction(observer);                     stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setSolver(solver);                          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseEmbedded(useEmbedded);                stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setZeroInitialGuess(zeroInitialGuess);      stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  std::string stepperType = stepper->getStepperType();
  // Full argument list construction.
  if      (stepperType == "RK Backward Euler") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperDIRK_BackwardEuler<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperDIRK_BackwardEuler<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SDIRK 2 Stage 2nd order") {
    double gamma = 0.2928932188134524;
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_2Stage2ndOrder<double>(
      model, obs, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, gamma));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_2Stage2ndOrder<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier, gamma));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SDIRK 3 Stage 2nd order") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_3Stage2ndOrder<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_3Stage2ndOrder<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SDIRK 2 Stage 3rd order") {
    std::string gammaType = "3rd Order A-stable";
    double gamma          = 0.7886751345948128;
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_2Stage3rdOrder<double>(
      model, obs, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, gammaType, gamma));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_2Stage3rdOrder<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier, gammaType, gamma));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "EDIRK 2 Stage 3rd order") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperEDIRK_2Stage3rdOrder<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperEDIRK_2Stage3rdOrder<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "DIRK 1 Stage Theta Method") {
    double theta = 0.5;
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperDIRK_1StageTheta<double>(
      model, obs, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess), theta);
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperDIRK_1StageTheta<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier, theta));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "EDIRK 2 Stage Theta Method") {
    double theta = 0.5;
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperEDIRK_2StageTheta<double>(
      model, obs, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, theta));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperEDIRK_2StageTheta<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier, theta));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
    s->setTheta(theta);
    s->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Trapezoidal Rule") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperEDIRK_TrapezoidalRule<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperEDIRK_TrapezoidalRule<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Implicit Midpoint") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_ImplicitMidpoint<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_ImplicitMidpoint<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SSPDIRK22") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_SSPDIRK22<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_SSPDIRK22<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SSPDIRK32") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_SSPDIRK32<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_SSPDIRK32<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SSPDIRK23") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_SSPDIRK23<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_SSPDIRK23<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SSPDIRK33") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_SSPDIRK33<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_SSPDIRK33<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Implicit 1 Stage 1st order Radau IA") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperDIRK_1Stage1stOrderRadauIA<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperDIRK_1Stage1stOrderRadauIA<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "RK Implicit 2 Stage 2nd order Lobatto IIIB") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperDIRK_2Stage2ndOrderLobattoIIIB<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperDIRK_2Stage2ndOrderLobattoIIIB<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SDIRK 5 Stage 4th order") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_5Stage4thOrder<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_5Stage4thOrder<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SDIRK 3 Stage 4th order") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_3Stage4thOrder<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_3Stage4thOrder<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SDIRK 5 Stage 5th order") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_5Stage5thOrder<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_5Stage5thOrder<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else if (stepperType == "SDIRK 2(1) Pair") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    { auto s = rcp(new Tempus::StepperSDIRK_21Pair<double>(
      model, obs, solver, useFSAL,
      ICConsistency, ICConsistencyCheck, useEmbedded, zeroInitialGuess));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized()); }
#endif
    auto s = rcp(new Tempus::StepperSDIRK_21Pair<double>(
      model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
      useEmbedded, zeroInitialGuess, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!s->isInitialized());
  } else TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error - unknown stepperType = "+stepperType);

}


/** \brief Unit test class for RK Stepper Modifier AppAction.
 */
class StepperRKModifierTest
  : virtual public Tempus::StepperRKModifierBase<double>
{
public:

  /// Constructor
  StepperRKModifierTest()
    : testBEGIN_STEP(false),
      testBEGIN_STAGE(false),
      testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),
      testBEFORE_EXPLICIT_EVAL(false),
      testEND_STAGE(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(-1.5),
      testType("")
  {}

  /// Destructor
  virtual ~StepperRKModifierTest(){}

  /// Modify RK Stepper at action location.
  virtual void modify(
    Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
    Teuchos::RCP<Tempus::StepperRKBase<double> > stepper,
    const typename Tempus::StepperRKAppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperRKAppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        auto x = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        testType = stepper->getStepperType() + " - Modifier";
        stepper->setStepperType(testType);
        break;
      }
      case StepperRKAppAction<double>::BEGIN_STAGE:
      {
        testBEGIN_STAGE = true;
        break;
      }
      case StepperRKAppAction<double>::BEFORE_SOLVE:
      {
        testBEFORE_SOLVE = true;
        testDt = sh->getWorkingState()->getTimeStep()/10.0;
        break;
      }
      case StepperRKAppAction<double>::AFTER_SOLVE:
      {
        testAFTER_SOLVE = true;
        break;
      }
      case StepperRKAppAction<double>::BEFORE_EXPLICIT_EVAL:
      {
        testBEFORE_EXPLICIT_EVAL = true;
        break;
      }
      case StepperRKAppAction<double>::END_STAGE:
      {
        testEND_STAGE = true;
        break;
      }
      case StepperRKAppAction<double>::END_STEP:
      {
        testEND_STEP = true;
        auto x = sh->getWorkingState()->getX();
        testWorkingValue = get_ele(*(x), 0);
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }
  }

  bool testBEGIN_STEP;
  bool testBEGIN_STAGE;
  bool testBEFORE_SOLVE;
  bool testAFTER_SOLVE;
  bool testBEFORE_EXPLICIT_EVAL;
  bool testEND_STAGE;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testType;
};


/** \brief Unit test class for RK Stepper Observer AppAction.
 */
class StepperRKObserverTest
  : virtual public Tempus::StepperRKObserverBase<double>
{
public:

  /// Constructor
  StepperRKObserverTest()
    : testBEGIN_STEP(false),
      testBEGIN_STAGE(false),
      testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),
      testBEFORE_EXPLICIT_EVAL(false),
      testEND_STAGE(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(-1.5),
      testType("")
  {}

  /// Destructor
  virtual ~StepperRKObserverTest(){}

  /// Observe RK Stepper at action location.
  virtual void observe(
    Teuchos::RCP<const Tempus::SolutionHistory<double> > sh,
    Teuchos::RCP<const Tempus::StepperRKBase<double> > stepper,
    const typename Tempus::StepperRKAppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperRKAppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        auto x = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperRKAppAction<double>::BEGIN_STAGE:
      {
        testBEGIN_STAGE = true;
        break;
      }
      case StepperRKAppAction<double>::BEFORE_SOLVE:
      {
        testBEFORE_SOLVE = true;
        testDt = sh->getWorkingState()->getTimeStep()/10.0;
        break;
      }
      case StepperRKAppAction<double>::AFTER_SOLVE:
      {
        testAFTER_SOLVE = true;
        testType = stepper->getStepperType() + " - Observer";
        break;
      }
      case StepperRKAppAction<double>::BEFORE_EXPLICIT_EVAL:
      {
        testBEFORE_EXPLICIT_EVAL = true;
        break;
      }
      case StepperRKAppAction<double>::END_STAGE:
      {
        testEND_STAGE = true;
        break;
      }
      case StepperRKAppAction<double>::END_STEP:
      {
        testEND_STEP = true;
        auto x = sh->getWorkingState()->getX();
        testWorkingValue = get_ele(*(x), 0);
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }
  }

  bool testBEGIN_STEP;
  bool testBEGIN_STAGE;
  bool testBEFORE_SOLVE;
  bool testAFTER_SOLVE;
  bool testBEFORE_EXPLICIT_EVAL;
  bool testEND_STAGE;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testType;
};


class StepperRKModifierXTest
  : virtual public Tempus::StepperRKModifierXBase<double>
{
public:

  /// Constructor
  StepperRKModifierXTest()
    : testX_BEGIN_STEP(false),
      testX_BEGIN_STAGE(false),
      testX_BEFORE_SOLVE(false),
      testX_AFTER_SOLVE(false),
      testX_BEFORE_EXPLICIT_EVAL(false),
      testXDOT_END_STAGE(false),
      testX_END_STEP(false),
      testX(-0.99),
      testEndStageX(-0.99),
      testDt(-1.5),
      testTime(-1.5),
      testStageNumber(-1),
      testStageX(-0.99)
  {}

  /// Destructor
  virtual ~StepperRKModifierXTest(){}

  /// Observe RK Stepper at end of takeStep.
  virtual void modify(
    Teuchos::RCP<Thyra::VectorBase<double> > x,
    const double time, const double dt,
    const int stageNumber,
    const typename Tempus::StepperRKModifierXBase<double>::MODIFIER_TYPE modType)
  {
    switch(modType) {
      case StepperRKModifierXBase<double>::X_BEGIN_STEP:
      {
        testX_BEGIN_STEP = true;
        testX = get_ele(*(x), 0);
        break;
      }
      case StepperRKModifierXBase<double>::X_BEGIN_STAGE:
      {
        testX_BEGIN_STAGE = true;
        break;
      }
      case StepperRKModifierXBase<double>::X_BEFORE_SOLVE:
      {
        testX_BEFORE_SOLVE = true;
        testDt = dt;
        break;
      }
      case StepperRKModifierXBase<double>::X_AFTER_SOLVE:
      {
        testX_AFTER_SOLVE = true;
        testTime = time;
        break;
      }
      case StepperRKModifierXBase<double>::X_BEFORE_EXPLICIT_EVAL:
      {
        testX_BEFORE_EXPLICIT_EVAL = true;
        testStageNumber = stageNumber;
        testStageX = get_ele(*(x), 0);  // x should be the stage value.
        break;
      }
      case StepperRKModifierXBase<double>::X_END_STAGE:
      {
        testXDOT_END_STAGE = true;
        testEndStageX = get_ele(*(x), 0);
        break;
      }
      case StepperRKModifierXBase<double>::X_END_STEP:
      {
        testX_END_STEP = true;
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }
  }

  bool testX_BEGIN_STEP;
  bool testX_BEGIN_STAGE;
  bool testX_BEFORE_SOLVE;
  bool testX_AFTER_SOLVE;
  bool testX_BEFORE_EXPLICIT_EVAL;
  bool testXDOT_END_STAGE;
  bool testX_END_STEP;
  double testX;
  double testEndStageX;
  double testDt;
  double testTime;
  int testStageNumber;
  double testStageX;
};


/** \brief Unit test utility for Stepper RK AppAction.
 */
void testRKAppAction(
  const Teuchos::RCP<Tempus::StepperRKBase<double> >& stepper,
  const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model,
  Teuchos::FancyOStream &out, bool &success)
{
  auto testTypeOrig = stepper->getStepperType();

  // Test Modifier.
  {
    stepper->setModel(model);
    auto modifier = rcp(new StepperRKModifierTest());
    stepper->setAppAction(modifier);
    stepper->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
    auto testType = testTypeOrig + " - Modifier";

    // Create a SolutionHistory.
    auto solutionHistory = Tempus::createSolutionHistoryME(model);

    // Take one time step.
    stepper->setInitialConditions(solutionHistory);
    solutionHistory->initWorkingState();
    double dt = 0.1;
    solutionHistory->getWorkingState()->setTimeStep(dt);
    stepper->takeStep(solutionHistory);

    // Testing that each ACTION_LOCATION has been called.
    TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
    TEST_COMPARE(modifier->testBEGIN_STAGE, ==, true);
    TEST_COMPARE(modifier->testBEFORE_SOLVE, ==, true);
    TEST_COMPARE(modifier->testAFTER_SOLVE, ==, true);
    TEST_COMPARE(modifier->testBEFORE_EXPLICIT_EVAL, ==, true);
    TEST_COMPARE(modifier->testEND_STAGE, ==, true);
    TEST_COMPARE(modifier->testEND_STEP, ==, true);

    // Testing that values can be set through the modifier.
    auto x = solutionHistory->getCurrentState()->getX();
    TEST_FLOATING_EQUALITY(modifier->testCurrentValue,get_ele(*(x), 0),1.0e-15);
    x = solutionHistory->getWorkingState()->getX();
    TEST_FLOATING_EQUALITY(modifier->testWorkingValue,get_ele(*(x), 0),1.0e-15);
    auto Dt = solutionHistory->getWorkingState()->getTimeStep();
    TEST_FLOATING_EQUALITY(modifier->testDt, Dt/10.0, 1.0e-15);

    TEST_COMPARE(modifier->testType, ==, testType);
  }

  // Test Observer.
  {
    stepper->setModel(model);
    auto observer = rcp(new StepperRKObserverTest());
    stepper->setAppAction(observer);
    stepper->setStepperType(testTypeOrig);
    stepper->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

    // Create a SolutionHistory.
    auto solutionHistory = Tempus::createSolutionHistoryME(model);

    // Take one time step.
    stepper->setInitialConditions(solutionHistory);
    solutionHistory->initWorkingState();
    double dt = 0.1;
    solutionHistory->getWorkingState()->setTimeStep(dt);
    stepper->takeStep(solutionHistory);

    // Testing that each ACTION_LOCATION has been called.
    TEST_COMPARE(observer->testBEGIN_STEP, ==, true);
    TEST_COMPARE(observer->testBEGIN_STAGE, ==, true);
    TEST_COMPARE(observer->testBEFORE_SOLVE, ==, true);
    TEST_COMPARE(observer->testAFTER_SOLVE, ==, true);
    TEST_COMPARE(observer->testBEFORE_EXPLICIT_EVAL, ==, true);
    TEST_COMPARE(observer->testEND_STAGE, ==, true);
    TEST_COMPARE(observer->testEND_STEP, ==, true);

    // Testing that values can be observed through the observer.
    auto x = solutionHistory->getCurrentState()->getX();
    TEST_FLOATING_EQUALITY(observer->testCurrentValue,get_ele(*(x), 0),1.0e-15);
    x = solutionHistory->getWorkingState()->getX();
    TEST_FLOATING_EQUALITY(observer->testWorkingValue,get_ele(*(x), 0),1.0e-15);
    auto Dt = solutionHistory->getWorkingState()->getTimeStep();
    TEST_FLOATING_EQUALITY(observer->testDt, Dt/10.0, 1.0e-15);

    auto testType = testTypeOrig + " - Observer";
    TEST_COMPARE(observer->testType, ==, testType);
  }

  // Test ModifierX.
  {
    stepper->setModel(model);
    auto modifierX = rcp(new StepperRKModifierXTest());
    stepper->setAppAction(modifierX);
    stepper->setStepperType(testTypeOrig);
    stepper->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

    // Create a SolutionHistory.
    auto solutionHistory = Tempus::createSolutionHistoryME(model);

    // Take one time step.
    stepper->setInitialConditions(solutionHistory);
    solutionHistory->initWorkingState();
    double dt = 0.1;
    solutionHistory->getWorkingState()->setTimeStep(dt);
    stepper->takeStep(solutionHistory);

    // Testing that each ACTION_LOCATION has been called.
    TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
    TEST_COMPARE(modifierX->testX_BEGIN_STAGE, ==, true);
    TEST_COMPARE(modifierX->testX_BEFORE_SOLVE, ==, true);
    TEST_COMPARE(modifierX->testX_AFTER_SOLVE, ==, true);
    TEST_COMPARE(modifierX->testX_BEFORE_EXPLICIT_EVAL, ==, true);
    TEST_COMPARE(modifierX->testXDOT_END_STAGE, ==, true);
    TEST_COMPARE(modifierX->testX_END_STEP, ==, true);

    // Testing that values can be set through the modifierX.
    auto x = solutionHistory->getCurrentState()->getX();
    TEST_FLOATING_EQUALITY(modifierX->testX, get_ele(*(x), 0), 1.0e-15);
    auto Dt = solutionHistory->getWorkingState()->getTimeStep();
    TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-15);

    auto time = solutionHistory->getWorkingState()->getTime();
    TEST_FLOATING_EQUALITY(modifierX->testTime, time, 1.0e-15);

    // Stage Number should be -1 outside stage loop.
    TEST_COMPARE(stepper->getStageNumber(), ==, -1);
    // The last stage number through X_BEFORE_EXPLICIT_EVAL should be
    // the number of stages minus one.
    TEST_COMPARE(modifierX->testStageNumber,==,stepper->getNumberOfStages()-1);

    if (rcp_dynamic_cast<Tempus::StepperIMEX_RK<double>>(stepper) != Teuchos::null ||
        rcp_dynamic_cast<Tempus::StepperIMEX_RK_Partition<double>>(stepper) != Teuchos::null) {
      TEST_FLOATING_EQUALITY(modifierX->testStageX, 2.0, 1.0e-15);
      TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 2.0,1.0e-15);
    } else if (stepper->isImplicit()) {
      // For implicit steppers, stageX is separate memory which can be tested.
      auto stageX = stepper->getStageX();
      TEST_FLOATING_EQUALITY(modifierX->testStageX, get_ele(*(stageX), 0),1.0e-15);
      TEST_FLOATING_EQUALITY(modifierX->testEndStageX, get_ele(*(stageX), 0),1.0e-15);
    } else {
      // For explicit steppers, stageX is under written and not available
      // outside takeStep, so direct comparisons are needed.
      if (rcp_dynamic_cast<Tempus::StepperERK_3_8Rule<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09966666666666668, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09966666666666668, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_3Stage3rdOrder<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.1, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.1, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_3Stage3rdOrderHeun<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.06666666666666667, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.06666666666666667, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_3Stage3rdOrderTVD<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.05, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.05, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_4Stage3rdOrderRunge<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09950000000000001, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09950000000000001, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_4Stage4thOrder<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09975000000000001, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09975000000000001, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_5Stage3rdOrderKandG<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.06662222222222222, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.06662222222222222, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_BogackiShampine32<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09983333333333333, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09983333333333333, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_ForwardEuler<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.0, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.0, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_General<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09975000000000001, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09975000000000001, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_Merson45<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09983333333333332, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09983333333333332, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_Midpoint<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.05, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.05, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_Trapezoidal<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.1, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.1, 1.0e-15);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error - unknown stepperType = "+stepper->getStepperType());
      }
    }
  }


  // Test Composite.
  {
    stepper->setModel(model);

    auto modifier = rcp(new StepperRKModifierTest());
    auto observer = rcp(new StepperRKObserverTest());
    auto modifierX = rcp(new StepperRKModifierXTest());
    auto composite = rcp(new Tempus::StepperRKAppActionComposite<double>());

    composite->addRKAppAction(modifier);
    composite->addRKAppAction(observer);
    composite->addRKAppAction(modifierX);
    stepper->setAppAction(composite);

    stepper->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

    // Create a SolutionHistory.
    auto solutionHistory = Tempus::createSolutionHistoryME(model);

    // Take one time step.
    stepper->setInitialConditions(solutionHistory);
    solutionHistory->initWorkingState();
    double dt = 0.1;
    solutionHistory->getWorkingState()->setTimeStep(dt);
    stepper->takeStep(solutionHistory);

    auto xCS = solutionHistory->getCurrentState()->getX();
    auto xWS = solutionHistory->getWorkingState()->getX();
    auto Dt = solutionHistory->getWorkingState()->getTimeStep();

    // Test Modifier.
    // Testing that each ACTION_LOCATION has been called.
    TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
    TEST_COMPARE(modifier->testBEGIN_STAGE, ==, true);
    TEST_COMPARE(modifier->testBEFORE_SOLVE, ==, true);
    TEST_COMPARE(modifier->testAFTER_SOLVE, ==, true);
    TEST_COMPARE(modifier->testBEFORE_EXPLICIT_EVAL, ==, true);
    TEST_COMPARE(modifier->testEND_STAGE, ==, true);
    TEST_COMPARE(modifier->testEND_STEP, ==, true);

    // Testing that values can be set through the modifier.
    TEST_FLOATING_EQUALITY(modifier->testCurrentValue, get_ele(*(xCS), 0),1.0e-15);
    TEST_FLOATING_EQUALITY(modifier->testWorkingValue, get_ele(*(xWS), 0),1.0e-15);
    TEST_FLOATING_EQUALITY(modifier->testDt, Dt/10.0, 1.0e-15);

    auto testType = testTypeOrig + " - Modifier";
    TEST_COMPARE(modifier->testType, ==, testType);


    // Test Observer.
    // Testing that each ACTION_LOCATION has been called.
    TEST_COMPARE(observer->testBEGIN_STEP, ==, true);
    TEST_COMPARE(observer->testBEGIN_STAGE, ==, true);
    TEST_COMPARE(observer->testBEFORE_SOLVE, ==, true);
    TEST_COMPARE(observer->testAFTER_SOLVE, ==, true);
    TEST_COMPARE(observer->testBEFORE_EXPLICIT_EVAL, ==, true);
    TEST_COMPARE(observer->testEND_STAGE, ==, true);
    TEST_COMPARE(observer->testEND_STEP, ==, true);

    // Testing that values can be observed through the observer.
    TEST_FLOATING_EQUALITY(observer->testCurrentValue, get_ele(*(xCS), 0),1.0e-15);
    TEST_FLOATING_EQUALITY(observer->testWorkingValue, get_ele(*(xWS), 0),1.0e-15);
    TEST_FLOATING_EQUALITY(observer->testDt, Dt/10.0, 1.0e-15);

    testType = testType + " - Observer";
    TEST_COMPARE(observer->testType, ==, testType);


    // Test ModifierX.
    // Testing that each ACTION_LOCATION has been called.
    TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
    TEST_COMPARE(modifierX->testX_BEGIN_STAGE, ==, true);
    TEST_COMPARE(modifierX->testX_BEFORE_SOLVE, ==, true);
    TEST_COMPARE(modifierX->testX_AFTER_SOLVE, ==, true);
    TEST_COMPARE(modifierX->testX_BEFORE_EXPLICIT_EVAL, ==, true);
    TEST_COMPARE(modifierX->testXDOT_END_STAGE, ==, true);
    TEST_COMPARE(modifierX->testX_END_STEP, ==, true);

    // Testing that values can be set through the modifierX.
    TEST_FLOATING_EQUALITY(modifierX->testX, get_ele(*(xCS), 0), 1.0e-15);
    TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-15);

    auto time = solutionHistory->getWorkingState()->getTime();
    TEST_FLOATING_EQUALITY(modifierX->testTime, time, 1.0e-15);

    // Stage Number should be -1 outside stage loop.
    TEST_COMPARE(stepper->getStageNumber(), ==, -1);
    // The last stage number through X_BEFORE_EXPLICIT_EVAL should be
    // the number of stages minus one.
    TEST_COMPARE(modifierX->testStageNumber,==,stepper->getNumberOfStages()-1);

    if (rcp_dynamic_cast<Tempus::StepperIMEX_RK<double>>(stepper) != Teuchos::null ||
        rcp_dynamic_cast<Tempus::StepperIMEX_RK_Partition<double>>(stepper) != Teuchos::null) {
      TEST_FLOATING_EQUALITY(modifierX->testStageX, 2.0, 1.0e-15);
      TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 2.0,1.0e-15);
    } else if (stepper->isImplicit()) {
      // For implicit steppers, stageX is separate memory which can be tested.
      auto stageX = stepper->getStageX();
      TEST_FLOATING_EQUALITY(modifierX->testStageX, get_ele(*(stageX), 0),1.0e-15);
      TEST_FLOATING_EQUALITY(modifierX->testEndStageX, get_ele(*(stageX), 0),1.0e-15);
    } else {
      // For explicit steppers, stageX is under written and not available
      // outside takeStep, so direct comparisons are needed.
      if (rcp_dynamic_cast<Tempus::StepperERK_3_8Rule<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09966666666666668, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09966666666666668, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_3Stage3rdOrder<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.1, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.1, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_3Stage3rdOrderHeun<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.06666666666666667, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.06666666666666667, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_3Stage3rdOrderTVD<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.05, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.05, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_4Stage3rdOrderRunge<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09950000000000001, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09950000000000001, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_4Stage4thOrder<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09975000000000001, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09975000000000001, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_5Stage3rdOrderKandG<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.06662222222222222, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.06662222222222222, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_BogackiShampine32<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09983333333333333, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09983333333333333, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_ForwardEuler<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.0, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.0, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_General<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09975000000000001, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09975000000000001, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_Merson45<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.09983333333333332, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.09983333333333332, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_Midpoint<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.05, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.05, 1.0e-15);
      } else if (rcp_dynamic_cast<Tempus::StepperERK_Trapezoidal<double>>(stepper) != Teuchos::null) {
        TEST_FLOATING_EQUALITY(modifierX->testStageX, 0.1, 1.0e-15);
        TEST_FLOATING_EQUALITY(modifierX->testEndStageX, 0.1, 1.0e-15);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error - unknown stepperType = "+stepper->getStepperType());
      }
    }
  }
}


} // namespace Tempus_Test
#endif // Tempus_UnitTest_Utils_hpp
