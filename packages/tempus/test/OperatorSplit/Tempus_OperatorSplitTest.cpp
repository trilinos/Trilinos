//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_StepperOperatorSplit.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"

#include "../TestModels/VanDerPol_IMEX_ExplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEX_ImplicitModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <vector>

namespace Tempus_Test {

using Teuchos::getParametersFromXmlFile;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::sublist;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(OperatorSplit, ConstructingFromDefaults)
{
  double dt = 0.05;

  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_OperatorSplit_VanDerPol.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the explicit VanDerPol ModelEvaluator
  RCP<ParameterList> vdpmPL = sublist(pList, "VanDerPolModel", true);
  RCP<const Thyra::ModelEvaluator<double>> explicitModel =
      rcp(new VanDerPol_IMEX_ExplicitModel<double>(vdpmPL));

  // Setup the implicit VanDerPol ModelEvaluator (reuse vdpmPL)
  RCP<const Thyra::ModelEvaluator<double>> implicitModel =
      rcp(new VanDerPol_IMEX_ImplicitModel<double>(vdpmPL));

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperOperatorSplit<double>());

  auto subStepper1 =
      Tempus::createStepperForwardEuler(explicitModel, Teuchos::null);
  auto subStepper2 =
      Tempus::createStepperBackwardEuler(implicitModel, Teuchos::null);

  stepper->addStepper(subStepper1);
  stepper->addStepper(subStepper2);
  stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL =
      pl->sublist("Demo Integrator").sublist("Time Step Control");
  timeStepControl->setInitIndex(tscPL.get<int>("Initial Time Index"));
  timeStepControl->setInitTime(tscPL.get<double>("Initial Time"));
  timeStepControl->setFinalTime(tscPL.get<double>("Final Time"));
  timeStepControl->setInitTimeStep(dt);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC = stepper->getModel()->getNominalValues();
  auto icX      = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
  auto icXDot   = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot());
  auto icState  = Tempus::createSolutionStateX(icX, icXDot);
  icState->setTime(timeStepControl->getInitTime());
  icState->setIndex(timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);
  icState->setOrder(stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Setup Integrator -----------------------------------------
  RCP<Tempus::IntegratorBasic<double>> integrator =
      Tempus::createIntegratorBasic<double>();
  integrator->setStepper(stepper);
  integrator->setTimeStepControl(timeStepControl);
  integrator->setSolutionHistory(solutionHistory);
  // integrator->setObserver(...);
  integrator->initialize();

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)

  // Test if at 'Final Time'
  double time      = integrator->getTime();
  double timeFinal = pl->sublist("Demo Integrator")
                         .sublist("Time Step Control")
                         .get<double>("Final Time");
  TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

  // Time-integrated solution and the exact solution
  RCP<Thyra::VectorBase<double>> x = integrator->getX();

  // Check the order and intercept
  out << "  Stepper = " << stepper->description() << std::endl;
  out << "  =========================" << std::endl;
  out << "  Computed solution: " << get_ele(*(x), 0) << "   "
      << get_ele(*(x), 1) << std::endl;
  out << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), -2.223910, 1.0e-4);
  TEST_FLOATING_EQUALITY(get_ele(*(x), 1), 0.565441, 1.0e-4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(OperatorSplit, VanDerPol)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 4;  // 8 for Error plot
  double dt                = 0.1;
  double time              = 0.0;
  for (int n = 0; n < nTimeStepSizes; n++) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_OperatorSplit_VanDerPol.xml");

    // Setup the explicit VanDerPol ModelEvaluator
    RCP<ParameterList> vdpmPL = sublist(pList, "VanDerPolModel", true);
    auto explicitModel        = rcp(new VanDerPol_IMEX_ExplicitModel<double>(vdpmPL));

    // Setup the implicit VanDerPol ModelEvaluator (reuse vdpmPL)
    auto implicitModel = rcp(new VanDerPol_IMEX_ImplicitModel<double>(vdpmPL));

    // Setup vector of models
    std::vector<RCP<const Thyra::ModelEvaluator<double>>> models;
    models.push_back(explicitModel);
    models.push_back(implicitModel);

    // Set the step size
    dt /= 2;
    if (n == nTimeStepSizes - 1) dt /= 10.0;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);
    integrator = Tempus::createIntegratorBasic<double>(pl, models);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time             = integrator->getTime();
    double timeFinal = pl->sublist("Demo Integrator")
                           .sublist("Time Step Control")
                           .get<double>("Final Time");
    double tol = 100.0 * std::numeric_limits<double>::epsilon();
    TEST_FLOATING_EQUALITY(time, timeFinal, tol);

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(implicitModel->get_x_space());
    Thyra::copy(*(integrator->getX()), solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(implicitModel->get_x_space());
    Thyra::copy(*(integrator->getXDot()), solutionDot.ptr());
    solutionsDot.push_back(solutionDot);

    // Output finest temporal solution for plotting
    // This only works for ONE MPI process
    if ((n == 0) || (n == nTimeStepSizes - 1)) {
      std::string fname = "Tempus_OperatorSplit_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_OperatorSplit_VanDerPol.dat";
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      writeSolution(fname, solutionHistory);
      // solutionHistory->printHistory("medium");
    }
  }

  // Check the order and intercept
  double xSlope                        = 0.0;
  double xDotSlope                     = 0.0;
  RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
  double order                         = stepper->getOrder();
  writeOrderError("Tempus_OperatorSplit_VanDerPol-Error.dat", stepper, StepSize,
                  solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                  xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, order, 0.05);
  TEST_FLOATING_EQUALITY(xDotSlope, order, 0.05);  //=order at small dt
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 1.27294, 1.0e-4);
  TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 12.7102, 1.0e-4);

  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
