//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tempus_StepperBackwardEuler.hpp"
#include "Tempus_TimeStepControlStrategyConstant.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

using Tempus::StepperFactory;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IntegratorBasic, Default_Construction)
{
  // Default Construction.
  auto integrator = Teuchos::rcp(new Tempus::IntegratorBasic<double>());
  TEST_ASSERT(integrator->isInitialized() == false);

  TEST_COMPARE(integrator->getIntegratorName(), ==, "Integrator Basic");
  TEST_COMPARE(integrator->getIntegratorType(), ==, "Integrator Basic");
  TEST_COMPARE(integrator->getStepper()->getStepperName(), ==, "Forward Euler");
  TEST_ASSERT(integrator->getStepper()->getModel() == Teuchos::null);
  TEST_ASSERT(integrator->getSolutionHistory() != Teuchos::null);
  TEST_COMPARE(integrator->getSolutionHistory()->getNumStates(), ==, 0);
  TEST_ASSERT(integrator->getTimeStepControl() != Teuchos::null);
  TEST_ASSERT(integrator->getTimeStepControl()->getStepType() == "Constant");
  TEST_ASSERT(integrator->getObserver() != Teuchos::null);

  // Setup ModelEvaluator -------------------------------------
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  integrator->setModel(model);

  // Setup SolutionHistory ------------------------------------
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState         = Tempus::createSolutionStateX(icSolution);
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->addState(icState);

  integrator->setSolutionHistory(solutionHistory);
  integrator->initialize();

  TEST_ASSERT(integrator->isInitialized() == true);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IntegratorBasic, Full_Construction)
{
  auto stepper = rcp(new Tempus::StepperBackwardEuler<double>());
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  stepper->setModel(model);

  // Setup SolutionHistory ------------------------------------
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState         = Tempus::createSolutionStateX(icSolution);
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->addState(icState);

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());

  // Setup IntegratorObserver ---------------------------------
  auto integratorObserver = rcp(new Tempus::IntegratorObserverBasic<double>());

  std::vector<int> outputScreenIndices{10, 20, 30};
  int outputScreenInterval = 72;

  // Full argument list construction.
  auto integrator = Teuchos::rcp(new Tempus::IntegratorBasic<double>(
      stepper, solutionHistory, timeStepControl, integratorObserver,
      outputScreenIndices, outputScreenInterval));

  TEST_ASSERT(integrator->isInitialized() == true);

  TEST_COMPARE(integrator->getIntegratorName(), ==, "Integrator Basic");
  TEST_COMPARE(integrator->getIntegratorType(), ==, "Integrator Basic");
  TEST_COMPARE(integrator->getStepper()->getStepperName(), ==,
               "Backward Euler");
  TEST_ASSERT(integrator->getStepper()->getModel() != Teuchos::null);
  TEST_ASSERT(integrator->getSolutionHistory() != Teuchos::null);
  TEST_COMPARE(integrator->getSolutionHistory()->getNumStates(), ==, 1);
  TEST_ASSERT(integrator->getTimeStepControl() != Teuchos::null);
  TEST_ASSERT(integrator->getTimeStepControl()->getStepType() == "Constant");
  TEST_ASSERT(integrator->getObserver() != Teuchos::null);
  TEST_ASSERT(integrator->getScreenOutputIndexList() == outputScreenIndices);
  TEST_ASSERT(integrator->getScreenOutputIndexInterval() ==
              outputScreenInterval);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IntegratorBasic, Describe)
{
  // 1) Setup the ParameterList (here we start with params from .xml file)
  RCP<ParameterList> pl =
      Teuchos::getParametersFromXmlFile("Tempus_IntegratorBasic_default.xml");

  // 2) Setup the ModelEvaluator
  auto model = Teuchos::rcp(new Tempus_Test::SinCosModel<double>());

  // 3) Setup the Integrator
  RCP<ParameterList> tempusPL = sublist(pl, "Tempus", true);
  RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::createIntegratorBasic<double>(tempusPL, model);

  std::ostringstream ss;
  Teuchos::RCP<Teuchos::FancyOStream> myOut =
      Teuchos::fancyOStream(Teuchos::rcpFromRef(ss));

  integrator->describe(*myOut, Teuchos::VERB_EXTREME);

  auto testS = ss.str();

  // Find major headers.
  auto npos = std::string::npos;
  TEST_ASSERT(npos != testS.find("--- Tempus::IntegratorBasic ---"));
  TEST_ASSERT(npos != testS.find("--- Tempus::SolutionHistory"));
  TEST_ASSERT(npos != testS.find("--- SolutionState (index =     0; time =     "
                                 "    0; dt =     1e+99) ---"));
  TEST_ASSERT(npos != testS.find("--- Tempus::SolutionStateMetaData ---"));
  TEST_ASSERT(npos != testS.find("--- Tempus::StepperState"));
  TEST_ASSERT(npos != testS.find("--- Tempus::PhysicsState"));
  TEST_ASSERT(npos != testS.find("--- Tempus::TimeStepControl ---"));
  TEST_ASSERT(npos !=
              testS.find("--- Tempus::TimeStepControlStrategyConstant ---"));
  TEST_ASSERT(npos != testS.find("--- Stepper ---"));
  TEST_ASSERT(npos != testS.find("stepperType_        = Forward Euler"));
  TEST_ASSERT(npos != testS.find("--- StepperExplicit ---"));

  integrator->setStatus(Tempus::Status::FAILED);
  TEST_ASSERT(integrator->getStatus() == Tempus::Status::FAILED);
  integrator->setStatus(Tempus::Status::WORKING);
  TEST_ASSERT(integrator->getStatus() == Tempus::Status::WORKING);
  integrator->setStatus(Tempus::Status::PASSED);
  TEST_ASSERT(integrator->getStatus() == Tempus::Status::PASSED);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IntegratorBasic, checkTimeStep)
{
  auto model      = Teuchos::rcp(new Tempus_Test::SinCosModel<double>());
  auto integrator = Tempus::createIntegratorBasic<double>(
      model, std::string("Backward Euler"));

  // Ensure initial status is working and unchanged by checkTimeStep.
  integrator->getNonConstSolutionHistory()->initWorkingState();
  TEST_ASSERT(integrator->getStatus() == Tempus::Status::WORKING);

  auto tsc = integrator->getNonConstTimeStepControl();
  auto ws  = integrator->getSolutionHistory()->getWorkingState();

  // integrator->checkTimeStep();

  // Test "Too many TimeStep failures"
  ws->setNFailures(11);
  integrator->checkTimeStep();
  TEST_ASSERT(integrator->getStatus() == Tempus::Status::FAILED);
  // Reset test.
  ws->setNFailures(0);
  integrator->setStatus(Tempus::Status::WORKING);

  // Test "Too many consecutive TimeStep failures"
  ws->setNConsecutiveFailures(6);
  integrator->checkTimeStep();
  TEST_ASSERT(integrator->getStatus() == Tempus::Status::FAILED);
  // Reset test.
  ws->setNConsecutiveFailures(0);
  integrator->setStatus(Tempus::Status::WORKING);

  // Test "Timestep size is at the minimum timestep size and the step failed."
  ws->setTimeStep(1.0);
  ws->setSolutionStatus(Tempus::Status::FAILED);
  tsc->setMinTimeStep(1.0);
  integrator->checkTimeStep();
  TEST_ASSERT(integrator->getStatus() == Tempus::Status::FAILED);
  // Reset test.
  ws->setSolutionStatus(Tempus::Status::PASSED);
  tsc->setMinTimeStep(0.1);
  integrator->setStatus(Tempus::Status::WORKING);

  // Test "Stepper failure."
  ws->setSolutionStatus(Tempus::Status::FAILED);
  integrator->checkTimeStep();
  TEST_ASSERT(ws->getNFailures() == 1);
  TEST_ASSERT(ws->getNRunningFailures() == 1);
  TEST_ASSERT(ws->getNConsecutiveFailures() == 1);
  TEST_ASSERT(ws->getSolutionStatus() == Tempus::Status::FAILED);
  // Reset test.
  ws->setNFailures(0);
  ws->setNRunningFailures(0);
  ws->setNConsecutiveFailures(0);
  ws->setSolutionStatus(Tempus::Status::PASSED);

  // Test "Constant time step failure."
  auto tscs = rcp(new Tempus::TimeStepControlStrategyConstant<double>());
  tsc->setTimeStepControlStrategy(tscs);
  ws->setTimeStep(0.1);
  tsc->setInitTimeStep(1.0);
  integrator->checkTimeStep();
  TEST_ASSERT(ws->getNFailures() == 1);
  TEST_ASSERT(ws->getNRunningFailures() == 1);
  TEST_ASSERT(ws->getNConsecutiveFailures() == 1);
  TEST_ASSERT(ws->getSolutionStatus() == Tempus::Status::FAILED);
  // Not resetting test as it is the last test.
}

// ************************************************************
// ************************************************************
// Test Integrator creation from ParameterList and ModelEvaluator.
TEUCHOS_UNIT_TEST(IntegratorBasic, PL_ME_Creation)
{
  // 1) Setup default Integrator
  RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::createIntegratorBasic<double>();

  // 2) Setup the ParameterList
  //    - Start with the default Tempus PL
  //    - Add Stepper PL
  RCP<ParameterList> tempusPL =
      Teuchos::rcp_const_cast<ParameterList>(integrator->getValidParameters());

  tempusPL->sublist("Default Integrator").set("Stepper Name", "Demo Stepper");
  RCP<ParameterList> stepperPL = Teuchos::parameterList();
  stepperPL->set("Stepper Type", "Forward Euler");
  tempusPL->set("Demo Stepper", *stepperPL);

  // 3) Create integrator from non-member function
  auto model = Teuchos::rcp(new Tempus_Test::SinCosModel<double>());
  integrator = Tempus::createIntegratorBasic<double>(tempusPL, model);

  // Test the ParameterList
  auto testPL = integrator->getValidParameters();
  // Write out ParameterList to rebaseline test.
  // writeParameterListToXmlFile(*testPL,"Tempus_IntegratorBasic_ref2-test.xml");

  // Read params from reference .xml file
  RCP<ParameterList> referencePL =
      Teuchos::getParametersFromXmlFile("Tempus_IntegratorBasic_ref.xml");

  bool pass = haveSameValuesSorted(*testPL, *referencePL, true);
  if (!pass) {
    out << std::endl;
    out << "testPL      -------------- \n"
        << *testPL << std::endl;
    out << "referencePL -------------- \n"
        << *referencePL << std::endl;
  }
  TEST_ASSERT(pass)
}

}  // namespace Tempus_Unit_Test
