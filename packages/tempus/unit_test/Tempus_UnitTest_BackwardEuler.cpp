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
#define CONSTRUCTION
#define STEPPERFACTORY_CONSTRUCTION


#ifdef CONSTRUCTION
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, Construction)
{
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperBackwardEuler<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Default values for construction.
  auto obs    = rcp(new Tempus::StepperBackwardEulerObserver<double>());
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  auto predictorStepper = rcp(new Tempus::StepperForwardEuler<double>());
  predictorStepper->setModel(model);  // Can use the same model since both steppers are implicit ODEs.
  predictorStepper->initialize();

  auto defaultStepper = rcp(new Tempus::StepperBackwardEuler<double>());
  bool useFSAL              = defaultStepper->getUseFSALDefault();
  std::string ICConsistency = defaultStepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = defaultStepper->getICConsistencyCheckDefault();
  bool zeroInitialGuess     = defaultStepper->getZeroInitialGuess();

  // Test the set functions.
  stepper->setObserver(obs);                           stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setSolver(solver);                          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setPredictor(predictorStepper);             stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setZeroInitialGuess(zeroInitialGuess);      stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Full argument list construction.
  stepper = rcp(new Tempus::StepperBackwardEuler<double>(
    model, obs, solver, predictorStepper, useFSAL,
    ICConsistency, ICConsistencyCheck, zeroInitialGuess));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

}
#endif // CONSTRUCTION


#ifdef STEPPERFACTORY_CONSTRUCTION
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("Backward Euler", model);
}
#endif // STEPPERFACTORY_CONSTRUCTION


} // namespace Tempus_Test
