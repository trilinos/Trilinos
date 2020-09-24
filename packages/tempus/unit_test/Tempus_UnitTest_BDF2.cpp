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
#include "Tempus_StepperHHTAlphaModifierBase.hpp"
#include "Tempus_StepperHHTAlphaModifierXBase.hpp"
#include "Tempus_StepperHHTAlphaObserverBase.hpp"
#include "Tempus_StepperHHTAlphaModifierDefault.hpp"
#include "Tempus_StepperHHTAlphaModifierXDefault.hpp"
#include "Tempus_StepperHHTAlphaObserverDefault.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
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

// Comment out any of the following tests to exclude from build/run.


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, Default_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperBDF2<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Default values for construction.
  auto obs    = rcp(new Tempus::StepperBDF2Observer<double>());
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  auto startUpStepper = rcp(new Tempus::StepperDIRK_1StageTheta<double>());
  startUpStepper->setModel(model);  // Can use the same model since both steppers are implicit ODEs.
  startUpStepper->initialize();

  auto defaultStepper = rcp(new Tempus::StepperBDF2<double>());
  bool useFSAL              = defaultStepper->getUseFSALDefault();
  std::string ICConsistency = defaultStepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = defaultStepper->getICConsistencyCheckDefault();
  bool zeroInitialGuess     = defaultStepper->getZeroInitialGuess();

  // Test the set functions.
  stepper->setObserver(obs);                           stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setSolver(solver);                          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setStartUpStepper(startUpStepper);          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setZeroInitialGuess(zeroInitialGuess);      stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Full argument list construction.
  stepper = rcp(new Tempus::StepperBDF2<double>(
    model, obs, solver, startUpStepper, useFSAL,
    ICConsistency, ICConsistencyCheck, zeroInitialGuess));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 2);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("BDF2", model);
}


} // namespace Tempus_Test
