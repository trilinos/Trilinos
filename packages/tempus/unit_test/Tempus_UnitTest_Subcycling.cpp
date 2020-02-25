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
//#include "Tempus_UnitTest_Utils.hpp"

#include "../TestModels/SinCosModel.hpp"
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

// Comment out any of the following tests to exclude from build/run.
#define CONSTRUCTION
#define MAXTIMESTEPDOESNOTCHANGEDURING_TAKESTEP


#ifdef CONSTRUCTION
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Subcycling, Default_Construction)
{
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperSubcycling<double>());
  auto sf = Teuchos::rcp(new Tempus::StepperFactory<double>());
  auto stepperBE = sf->createStepperBackwardEuler(model, Teuchos::null);
  stepper->setSubcyclingStepper(stepperBE);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto obs    = rcp(new Tempus::StepperSubcyclingObserver<double>());
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSALDefault();
  std::string ICConsistency = stepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheckDefault();

  // Test the set functions
  stepper->setSolver(solver);                          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setObserver(obs);                           stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Full argument list construction.
  auto scIntegrator = Teuchos::rcp(new Tempus::IntegratorBasic<double>());
  auto stepperFE = sf->createStepperForwardEuler(model, Teuchos::null);
  scIntegrator->setStepperWStepper(stepperFE);
  scIntegrator->initialize();

  stepper = rcp(new Tempus::StepperSubcycling<double>(
    model, obs, scIntegrator, useFSAL, ICConsistency, ICConsistencyCheck));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

}
#endif // CONSTRUCTION


#ifdef MAXTIMESTEPDOESNOTCHANGEDURING_TAKESTEP
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Subcycling, MaxTimeStepDoesNotChangeDuring_takeStep)
{
  // Setup the stepper ----------------------------------------
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  auto stepper = rcp(new Tempus::StepperSubcycling<double>());
  auto sf = Teuchos::rcp(new Tempus::StepperFactory<double>());
  auto stepperBE = sf->createStepperBackwardEuler(model, Teuchos::null);
  stepper->setSubcyclingStepper(stepperBE);
  stepper->initialize();

  // Setup SolutionHistory ------------------------------------
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
    stepper->getModel()->getNominalValues();
  auto icSolution =rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = rcp(new Tempus::SolutionState<double>(icSolution));
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->addState(icState);
  solutionHistory->initWorkingState();

  // Test
  stepper->setSubcyclingMaxTimeStep(0.5);
  double maxTimeStep_Set = stepper->getSubcyclingMaxTimeStep();
  stepper->takeStep(solutionHistory);
  double maxTimeStep_After = stepper->getSubcyclingMaxTimeStep();

  TEST_FLOATING_EQUALITY(maxTimeStep_Set, maxTimeStep_After, 1.0e-14 );
}
#endif // MAXTIMESTEPDOESNOTCHANGEDURING_TAKESTEP


} // namespace Tempus_Test
