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
#include "Tempus_WrapperModelEvaluatorPairPartIMEX_Basic.hpp"
#include "Tempus_StepperIMEX_RK_Partition.hpp"

#include "../TestModels/VanDerPol_IMEX_ExplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEXPart_ImplicitModel.hpp"
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
TEUCHOS_UNIT_TEST(IMEX_RK_Partitioned, ConstructingFromDefaults)
{
  double dt = 0.025;

  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_IMEX_RK_VanDerPol.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the explicit VanDerPol ModelEvaluator
  RCP<ParameterList> vdpmPL   = sublist(pList, "VanDerPolModel", true);
  const bool useProductVector = true;
  auto explicitModel =
      rcp(new VanDerPol_IMEX_ExplicitModel<double>(vdpmPL, useProductVector));

  // Setup the implicit VanDerPol ModelEvaluator (reuse vdpmPL)
  auto implicitModel =
      rcp(new VanDerPol_IMEXPart_ImplicitModel<double>(vdpmPL));

  // Setup the IMEX Pair ModelEvaluator
  const int numExplicitBlocks = 1;
  const int parameterIndex    = 4;
  auto model                  = rcp(new Tempus::WrapperModelEvaluatorPairPartIMEX_Basic<double>(
      explicitModel, implicitModel, numExplicitBlocks, parameterIndex));

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperIMEX_RK_Partition<double>());
  stepper->setModel(model);
  stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL =
      pl->sublist("Default Integrator").sublist("Time Step Control");
  timeStepControl->setInitIndex(tscPL.get<int>("Initial Time Index"));
  timeStepControl->setInitTime(tscPL.get<double>("Initial Time"));
  timeStepControl->setFinalTime(tscPL.get<double>("Final Time"));
  timeStepControl->setInitTimeStep(dt);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC   = model->getNominalValues();
  auto icSolution = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
  auto icState    = Tempus::createSolutionStateX(icSolution);
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
  double timeFinal = pl->sublist("Default Integrator")
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
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), 1.810210, 1.0e-4);
  TEST_FLOATING_EQUALITY(get_ele(*(x), 1), -0.754602, 1.0e-4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK_Partitioned, VanDerPol)
{
  std::vector<std::string> stepperTypes;
  stepperTypes.push_back("Partitioned IMEX RK 1st order");
  stepperTypes.push_back("Partitioned IMEX RK SSP2");
  stepperTypes.push_back("Partitioned IMEX RK ARS 233");
  stepperTypes.push_back("General Partitioned IMEX RK");

  std::vector<double> stepperOrders;
  stepperOrders.push_back(1.07964);
  stepperOrders.push_back(2.00408);
  stepperOrders.push_back(2.70655);
  stepperOrders.push_back(2.00211);

  std::vector<double> stepperErrors;
  stepperErrors.push_back(0.0046423);
  stepperErrors.push_back(0.0154534);
  stepperErrors.push_back(0.000298908);
  stepperErrors.push_back(0.0071546);

  std::vector<double> stepperInitDt;
  stepperInitDt.push_back(0.0125);
  stepperInitDt.push_back(0.05);
  stepperInitDt.push_back(0.05);
  stepperInitDt.push_back(0.05);

  std::vector<std::string>::size_type m;
  for (m = 0; m != stepperTypes.size(); m++) {
    std::string stepperType = stepperTypes[m];
    std::string stepperName = stepperTypes[m];
    std::replace(stepperName.begin(), stepperName.end(), ' ', '_');
    std::replace(stepperName.begin(), stepperName.end(), '/', '.');

    RCP<Tempus::IntegratorBasic<double>> integrator;
    std::vector<RCP<Thyra::VectorBase<double>>> solutions;
    std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
    std::vector<double> StepSize;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;

    const int nTimeStepSizes = 3;  // 6 for error plot
    double dt                = stepperInitDt[m];
    double time              = 0.0;
    for (int n = 0; n < nTimeStepSizes; n++) {
      // Read params from .xml file
      RCP<ParameterList> pList =
          getParametersFromXmlFile("Tempus_IMEX_RK_VanDerPol.xml");

      // Setup the explicit VanDerPol ModelEvaluator
      RCP<ParameterList> vdpmPL   = sublist(pList, "VanDerPolModel", true);
      const bool useProductVector = true;
      auto explicitModel          = rcp(
                   new VanDerPol_IMEX_ExplicitModel<double>(vdpmPL, useProductVector));

      // Setup the implicit VanDerPol ModelEvaluator (reuse vdpmPL)
      auto implicitModel =
          rcp(new VanDerPol_IMEXPart_ImplicitModel<double>(vdpmPL));

      // Setup the IMEX Pair ModelEvaluator
      const int numExplicitBlocks = 1;
      const int parameterIndex    = 4;
      auto model =
          rcp(new Tempus::WrapperModelEvaluatorPairPartIMEX_Basic<double>(
              explicitModel, implicitModel, numExplicitBlocks, parameterIndex));

      // Set the Stepper
      RCP<ParameterList> pl = sublist(pList, "Tempus", true);

      if (stepperType == "General Partitioned IMEX RK") {
        // use the appropriate stepper sublist
        pl->sublist("Default Integrator")
            .set("Stepper Name", "General IMEX RK");
      }
      else {
        pl->sublist("Default Stepper").set("Stepper Type", stepperType);
      }

      // Set the step size
      if (n == nTimeStepSizes - 1)
        dt /= 10.0;
      else
        dt /= 2;

      // Setup the Integrator and reset initial time step
      pl->sublist("Default Integrator")
          .sublist("Time Step Control")
          .set("Initial Time Step", dt);
      integrator = Tempus::createIntegratorBasic<double>(pl, model);

      // Integrate to timeMax
      bool integratorStatus = integrator->advanceTime();
      TEST_ASSERT(integratorStatus)

      // Test if at 'Final Time'
      time             = integrator->getTime();
      double timeFinal = pl->sublist("Default Integrator")
                             .sublist("Time Step Control")
                             .get<double>("Final Time");
      double tol = 100.0 * std::numeric_limits<double>::epsilon();
      TEST_FLOATING_EQUALITY(time, timeFinal, tol);

      // Store off the final solution and step size
      StepSize.push_back(dt);
      auto solution = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(integrator->getX()), solution.ptr());
      solutions.push_back(solution);
      auto solutionDot = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(integrator->getXDot()), solutionDot.ptr());
      solutionsDot.push_back(solutionDot);

      // Output finest temporal solution for plotting
      // This only works for ONE MPI process
      if ((n == 0) || (n == nTimeStepSizes - 1)) {
        std::string fname = "Tempus_" + stepperName + "_VanDerPol-Ref.dat";
        if (n == 0) fname = "Tempus_" + stepperName + "_VanDerPol.dat";
        RCP<const SolutionHistory<double>> solutionHistory =
            integrator->getSolutionHistory();
        writeSolution(fname, solutionHistory);
      }
    }

    // Check the order and intercept
    double xSlope                        = 0.0;
    double xDotSlope                     = 0.0;
    RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
    // double order = stepper->getOrder();

    // xDot not yet available for DIRK methods, e.g., are not calc. and zero.
    solutionsDot.clear();

    writeOrderError("Tempus_" + stepperName + "_VanDerPol-Error.dat", stepper,
                    StepSize, solutions, xErrorNorm, xSlope, solutionsDot,
                    xDotErrorNorm, xDotSlope, out);

    TEST_FLOATING_EQUALITY(xSlope, stepperOrders[m], 0.02);
    TEST_FLOATING_EQUALITY(xErrorNorm[0], stepperErrors[m], 1.0e-4);
    // TEST_FLOATING_EQUALITY( xDotSlope,              1.74898, 0.02 );
    // TEST_FLOATING_EQUALITY( xDotErrorNorm[0],        1.0038, 1.0e-4 );
  }
  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
