// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_UnitTest_Utils.hpp"
#include "Tempus_StepperRKButcherTableau.hpp"

#include "../TestModels/VanDerPol_IMEX_ExplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEXPart_ImplicitModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <vector>

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::StepperFactory;
using Tempus::StepperExplicitRK;


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK_Partition, Default_Construction)
{
  // Setup the explicit VanDerPol ModelEvaluator
  const bool useProductVector = true;
  auto explicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>(Teuchos::null, useProductVector));
  auto implicitModel = rcp(new Tempus_Test::VanDerPol_IMEXPart_ImplicitModel<double>());

  // Setup the IMEX Pair ModelEvaluator
  const int numExplicitBlocks = 1;
  const int parameterIndex = 4;
  auto model = rcp(new Tempus::WrapperModelEvaluatorPairPartIMEX_Basic<double>(
                         explicitModel, implicitModel,
                         numExplicitBlocks, parameterIndex));


  // Default construction.
  auto stepper = rcp(new Tempus::StepperIMEX_RK_Partition<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto modifier  = rcp(new Tempus::StepperRKModifierDefault<double>());
  auto modifierX = rcp(new Tempus::StepperRKModifierXDefault<double>());
  auto observer  = rcp(new Tempus::StepperRKObserverDefault<double>());
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSALDefault();
  std::string ICConsistency = stepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheckDefault();
  bool zeroInitialGuess     = stepper->getZeroInitialGuess();
  std::string stepperType   = "IMEX RK SSP2";
  auto stepperERK = Teuchos::rcp(new Tempus::StepperERK_Trapezoidal<double>());
  auto explicitTableau      = stepperERK->getTableau();
  auto stepperSDIRK = Teuchos::rcp(new Tempus::StepperSDIRK_2Stage3rdOrder<double>());
  auto implicitTableau      = stepperSDIRK->getTableau();
  int order                 = 2;


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
  stepper->setZeroInitialGuess(zeroInitialGuess);      stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  stepper->setStepperType(stepperType);                stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setExplicitTableau(explicitTableau);        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setImplicitTableau(implicitTableau);        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setOrder(order);                            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  TEUCHOS_TEST_FOR_EXCEPT(explicitTableau != stepper->getTableau());
  TEUCHOS_TEST_FOR_EXCEPT(explicitTableau != stepper->getExplicitTableau());
  TEUCHOS_TEST_FOR_EXCEPT(implicitTableau != stepper->getImplicitTableau());

  // Full argument list construction.
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  stepper = rcp(new Tempus::StepperIMEX_RK_Partition<double>(
    model, obs, solver, useFSAL,
    ICConsistency, ICConsistencyCheck, zeroInitialGuess,
    stepperType, explicitTableau, implicitTableau, order));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
  stepper = rcp(new Tempus::StepperIMEX_RK_Partition<double>(
    model, solver, useFSAL, ICConsistency, ICConsistencyCheck,
    zeroInitialGuess, modifier, stepperType, explicitTableau,
    implicitTableau, order));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  std::cout << "ordero = " << stepper->getOrder() << std::endl;
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK_Partition, StepperFactory_Construction)
{
  // Setup the explicit VanDerPol ModelEvaluator
  const bool useProductVector = true;
  auto explicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>(Teuchos::null, useProductVector));
  auto implicitModel = rcp(new Tempus_Test::VanDerPol_IMEXPart_ImplicitModel<double>());

  // Setup the IMEX Pair ModelEvaluator
  const int numExplicitBlocks = 1;
  const int parameterIndex = 4;
  auto model = rcp(new Tempus::WrapperModelEvaluatorPairPartIMEX_Basic<double>(
                         explicitModel, implicitModel,
                         numExplicitBlocks, parameterIndex));

  testFactoryConstruction("Partitioned IMEX RK SSP2", model);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK_Partition, StepperFactory_Construction_General_wo_Parameterlist)
{
  // Setup the explicit VanDerPol ModelEvaluator
  const bool useProductVector = true;
  auto explicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>(Teuchos::null, useProductVector));
  auto implicitModel = rcp(new Tempus_Test::VanDerPol_IMEXPart_ImplicitModel<double>());

  // Setup the IMEX Pair ModelEvaluator
  const int numExplicitBlocks = 1;
  const int parameterIndex = 4;
  auto model = rcp(new Tempus::WrapperModelEvaluatorPairPartIMEX_Basic<double>(
                         explicitModel, implicitModel,
                         numExplicitBlocks, parameterIndex));

  RCP<StepperFactory<double> > sf = Teuchos::rcp(new StepperFactory<double>());

  auto stepper = sf->createStepper("General Partitioned IMEX RK", model);
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK_Partition, StepperFactory_Construction_General_wo_Parameterlist_Model)
{
  // Setup the explicit VanDerPol ModelEvaluator
  const bool useProductVector = true;
  auto explicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>(Teuchos::null, useProductVector));
  auto implicitModel = rcp(new Tempus_Test::VanDerPol_IMEXPart_ImplicitModel<double>());

  // Setup the IMEX Pair ModelEvaluator
  const int numExplicitBlocks = 1;
  const int parameterIndex = 4;
  auto model = rcp(new Tempus::WrapperModelEvaluatorPairPartIMEX_Basic<double>(
                         explicitModel, implicitModel,
                         numExplicitBlocks, parameterIndex));

  RCP<StepperFactory<double> > sf = Teuchos::rcp(new StepperFactory<double>());

  auto stepper = sf->createStepper("General Partitioned IMEX RK");
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK_Partition, AppAction)
{
  auto stepper = rcp(new Tempus::StepperIMEX_RK_Partition<double>());

  // Setup the explicit VanDerPol ModelEvaluator
  const bool useProductVector = true;
  auto explicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>(Teuchos::null, useProductVector));
  auto implicitModel = rcp(new Tempus_Test::VanDerPol_IMEXPart_ImplicitModel<double>());

  // Setup the IMEX Pair ModelEvaluator
  const int numExplicitBlocks = 1;
  const int parameterIndex = 4;
  auto model = rcp(new Tempus::WrapperModelEvaluatorPairPartIMEX_Basic<double>(
                         explicitModel, implicitModel,
                         numExplicitBlocks, parameterIndex));

  testRKAppAction(stepper, model, out, success);
}


} // namespace Tempus_Test
