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
#include "Tempus_IntegratorObserverSubcycling.hpp"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"
#include "Tempus_StepperSubcycling.hpp"
#include "Tempus_StepperOperatorSplit.hpp"
#include "Tempus_TimeStepControlStrategyConstant.hpp"
#include "Tempus_TimeStepControlStrategyBasicVS.hpp"

#include "../TestModels/SinCosModel.hpp"
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
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Subcycling, ParameterList)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_Subcycling_SinCos.xml");

  // std::ofstream ftmp("PL.txt");
  // pList->print(ftmp);
  // ftmp.close();

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  auto model                = rcp(new SinCosModel<double>(scm_pl));

  RCP<ParameterList> tempusPL = sublist(pList, "Tempus", true);

  // Test constructor IntegratorBasic(tempusPL, model)
  {
    RCP<Tempus::IntegratorBasic<double>> integrator =
        Tempus::createIntegratorBasic<double>(tempusPL, model);

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Demo Stepper", true);
    RCP<const ParameterList> defaultPL =
        integrator->getStepper()->getValidParameters();

    bool pass = haveSameValuesSorted(*stepperPL, *defaultPL, true);
    if (!pass) {
      out << std::endl;
      out << "stepperPL -------------- \n"
          << *stepperPL << std::endl;
      out << "defaultPL -------------- \n"
          << *defaultPL << std::endl;
    }
    TEST_ASSERT(pass)
  }

  // Test constructor IntegratorBasic(model, stepperType)
  {
    RCP<Tempus::IntegratorBasic<double>> integrator =
        Tempus::createIntegratorBasic<double>(model,
                                              std::string("Forward Euler"));

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Demo Stepper", true);
    RCP<const ParameterList> defaultPL =
        integrator->getStepper()->getValidParameters();

    bool pass = haveSameValuesSorted(*stepperPL, *defaultPL, true);
    if (!pass) {
      out << std::endl;
      out << "stepperPL -------------- \n"
          << *stepperPL << std::endl;
      out << "defaultPL -------------- \n"
          << *defaultPL << std::endl;
    }
    TEST_ASSERT(pass)
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Subcycling, ConstructingFromDefaults)
{
  double dt = 0.4;

  // Setup the SinCosModel ------------------------------------
  auto model   = rcp(new SinCosModel<double>());
  auto modelME = rcp_dynamic_cast<const Thyra::ModelEvaluator<double>>(model);

  // Setup Stepper for field solve ----------------------------
  auto stepper   = rcp(new Tempus::StepperSubcycling<double>());
  auto stepperFE = Tempus::createStepperForwardEuler(modelME, Teuchos::null);
  stepper->setSubcyclingStepper(stepperFE);

  stepper->setSubcyclingMinTimeStep(0.1);
  stepper->setSubcyclingInitTimeStep(0.1);
  stepper->setSubcyclingMaxTimeStep(0.1);
  stepper->setSubcyclingMaxFailures(10);
  stepper->setSubcyclingMaxConsecFailures(5);
  stepper->setSubcyclingScreenOutputIndexInterval(1);
  stepper->setSubcyclingIntegratorObserver(
      Teuchos::rcp(new Tempus::IntegratorObserverSubcycling<double>()));
  stepper->setSubcyclingPrintDtChanges(true);

  // Set subcycling strategy.
  auto subStrategy =
      rcp(new Tempus::TimeStepControlStrategyConstant<double>(dt));
  stepper->setSubcyclingTimeStepControlStrategy(subStrategy);

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
  timeStepControl->setInitIndex(0);
  timeStepControl->setFinalIndex(10);
  timeStepControl->setInitTime(0.0);
  timeStepControl->setFinalTime(1.0);
  timeStepControl->setInitTimeStep(dt);

  // Set TimeStepControl strategy.
  auto strategy = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
  strategy->initialize();
  timeStepControl->setTimeStepControlStrategy(strategy);

  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  auto inArgsIC   = stepper->getModel()->getNominalValues();
  auto icSolution = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
  auto icState    = Tempus::createSolutionStateX(icSolution);
  icState->setTime(timeStepControl->getInitTime());
  icState->setIndex(timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);                           // dt for ICs are indicated by zero.
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Ensure ICs are consistent and stepper memory is set (e.g., xDot is set).
  stepper->setInitialConditions(solutionHistory);
  stepper->initialize();

  // Setup Integrator -----------------------------------------
  RCP<Tempus::IntegratorBasic<double>> integrator =
      Tempus::createIntegratorBasic<double>();
  integrator->setStepper(stepper);
  integrator->setTimeStepControl(timeStepControl);
  integrator->setSolutionHistory(solutionHistory);
  integrator->setScreenOutputIndexInterval(1);
  // integrator->setObserver(...);
  integrator->initialize();

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)

  // Test if at 'Final Time'
  double time      = integrator->getTime();
  double timeFinal = 1.0;
  TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

  // Time-integrated solution and the exact solution
  RCP<Thyra::VectorBase<double>> x = integrator->getX();
  RCP<const Thyra::VectorBase<double>> x_exact =
      model->getExactSolution(time).get_x();

  // Calculate the error
  RCP<Thyra::VectorBase<double>> xdiff = x->clone_v();
  Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));

  // Check the order and intercept
  out << "  Stepper = " << stepper->description() << std::endl;
  out << "  =========================" << std::endl;
  out << "  Exact solution   : " << get_ele(*(x_exact), 0) << "   "
      << get_ele(*(x_exact), 1) << std::endl;
  out << "  Computed solution: " << get_ele(*(x), 0) << "   "
      << get_ele(*(x), 1) << std::endl;
  out << "  Difference       : " << get_ele(*(xdiff), 0) << "   "
      << get_ele(*(xdiff), 1) << std::endl;
  out << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), 0.882508, 1.0e-4);
  TEST_FLOATING_EQUALITY(get_ele(*(x), 1), 0.570790, 1.0e-4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Subcycling, SinCosAdapt)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;

  double dt = 0.05;

  // Setup the SinCosModel
  const int nTimeStepSizes       = 2;
  std::string output_file_string = "Tempus_Subcycling_SinCos";
  std::string output_file_name   = output_file_string + ".dat";
  std::string err_out_file_name  = output_file_string + "-Error.dat";
  double time                    = 0.0;
  for (int n = 0; n < nTimeStepSizes; n++) {
    dt /= 2;

    // Setup the SinCosModel ------------------------------------
    auto model   = rcp(new SinCosModel<double>());
    auto modelME = rcp_dynamic_cast<const Thyra::ModelEvaluator<double>>(model);

    // Setup Stepper for field solve ----------------------------
    auto stepper   = rcp(new Tempus::StepperSubcycling<double>());
    auto stepperFE = Tempus::createStepperForwardEuler(modelME, Teuchos::null);
    stepper->setSubcyclingStepper(stepperFE);

    stepper->setSubcyclingMinTimeStep(dt / 10.0);
    stepper->setSubcyclingInitTimeStep(dt / 10.0);
    stepper->setSubcyclingMaxTimeStep(dt);
    stepper->setSubcyclingMaxFailures(10);
    stepper->setSubcyclingMaxConsecFailures(5);
    stepper->setSubcyclingScreenOutputIndexInterval(1);
    // stepper->setSubcyclingIntegratorObserver(
    //   Teuchos::rcp(new Tempus::IntegratorObserverSubcycling<double>()));
    // stepper->setSubcyclingPrintDtChanges   (true);

    // Set variable strategy.
    auto strategy = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
    strategy->setMinEta(0.02);
    strategy->setMaxEta(0.04);
    strategy->initialize();
    stepper->setSubcyclingTimeStepControlStrategy(strategy);

    // Setup TimeStepControl ------------------------------------
    auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
    timeStepControl->setInitIndex(0);
    timeStepControl->setInitTime(0.0);
    timeStepControl->setFinalTime(1.0);
    timeStepControl->setMinTimeStep(dt);
    timeStepControl->setInitTimeStep(dt);
    timeStepControl->setMaxTimeStep(dt);
    timeStepControl->initialize();

    // Setup initial condition SolutionState --------------------
    auto inArgsIC = stepper->getModel()->getNominalValues();
    auto icSolution =
        rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
    auto icState = Tempus::createSolutionStateX(icSolution);
    icState->setTime(timeStepControl->getInitTime());
    icState->setIndex(timeStepControl->getInitIndex());
    icState->setTimeStep(0.0);                           // dt for ICs are zero.
    icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

    // Setup SolutionHistory ------------------------------------
    auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
    solutionHistory->setName("Forward States");
    solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
    solutionHistory->setStorageLimit(2);
    solutionHistory->addState(icState);

    // Ensure ICs are consistent and stepper memory is set (e.g., xDot is set).
    stepper->setInitialConditions(solutionHistory);
    stepper->initialize();

    // Setup Integrator -----------------------------------------
    integrator = Tempus::createIntegratorBasic<double>();
    integrator->setStepper(stepper);
    integrator->setTimeStepControl(timeStepControl);
    integrator->setSolutionHistory(solutionHistory);
    integrator->setScreenOutputIndexInterval(10);
    // integrator->setObserver(...);
    integrator->initialize();

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time             = integrator->getTime();
    double timeFinal = 1.0;
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double>> x = integrator->getX();
    RCP<const Thyra::VectorBase<double>> x_exact =
        model->getExactSolution(time).get_x();

    //// Plot sample solution and exact solution
    // if (n == 0) {
    //   std::ofstream ftmp(output_file_name);
    //   //Warning: the following assumes serial run
    //   FILE *gold_file = fopen("Tempus_Subcycling_SinCos_AdaptDt_gold.dat",
    //   "r"); RCP<const SolutionHistory<double> > solutionHistory =
    //     integrator->getSolutionHistory();
    //   RCP<const Thyra::VectorBase<double> > x_exact_plot;
    //   for (int i=0; i<solutionHistory->getNumStates(); i++) {
    //     char time_gold_char[100];
    //     fgets(time_gold_char, 100, gold_file);
    //     double time_gold;
    //     sscanf(time_gold_char, "%lf", &time_gold);
    //     RCP<const SolutionState<double> > solutionState =
    //     (*solutionHistory)[i]; double time_i = solutionState->getTime();
    //     //Throw error if time does not match time in gold file to specified
    //     tolerance TEST_FLOATING_EQUALITY( time_i, time_gold, 1.0e-5 );
    //     RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
    //     x_exact_plot = model->getExactSolution(time_i).get_x();
    //     ftmp << time_i << "   "
    //          << get_ele(*(x_plot), 0) << "   "
    //          << get_ele(*(x_plot), 1) << "   "
    //          << get_ele(*(x_exact_plot), 0) << "   "
    //          << get_ele(*(x_exact_plot), 1) << std::endl;
    //   }
    //   ftmp.close();
    // }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()), solution.ptr());
    solutions.push_back(solution);
    if (n == nTimeStepSizes - 1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solutionExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),
                  solutionExact.ptr());
      solutions.push_back(solutionExact);
    }
  }

  // Check the order and intercept
  if (nTimeStepSizes > 1) {
    double xSlope    = 0.0;
    double xDotSlope = 0.0;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;
    RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
    // double order = stepper->getOrder();
    writeOrderError("Tempus_BDF2_SinCos-Error.dat", stepper, StepSize,
                    solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                    xDotSlope, out);

    TEST_FLOATING_EQUALITY(xSlope, 1.00137, 0.01);
    // TEST_FLOATING_EQUALITY( xDotSlope,            1.95089, 0.01 );
    TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.00387948, 1.0e-4);
    // TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.000197325, 1.0e-4 );
  }

  Teuchos::TimeMonitor::summarize();
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Subcycling, VanDerPolOperatorSplit)
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
    // Set the step size
    dt /= 2;
    if (n == nTimeStepSizes - 1) dt /= 10.0;

    // Setup the explicit and implicit VanDerPol ModelEvaluators
    auto tmpModel = rcp(new VanDerPol_IMEX_ExplicitModel<double>());
    auto pl       = Teuchos::rcp_const_cast<Teuchos::ParameterList>(
        tmpModel->getValidParameters());
    pl->set("Coeff epsilon", 0.1);
    RCP<const Thyra::ModelEvaluator<double>> explicitModel =
        rcp(new VanDerPol_IMEX_ExplicitModel<double>(pl));
    RCP<const Thyra::ModelEvaluator<double>> implicitModel =
        rcp(new VanDerPol_IMEX_ImplicitModel<double>(pl));

    // Setup Steppers for field solve ---------------------------

    // Explicit Subcycling Stepper
    auto stepperSC = rcp(new Tempus::StepperSubcycling<double>());
    auto stepperFE =
        Tempus::createStepperForwardEuler(explicitModel, Teuchos::null);
    stepperFE->setUseFSAL(false);
    stepperFE->initialize();
    stepperSC->setSubcyclingStepper(stepperFE);

    stepperSC->setSubcyclingMinTimeStep(0.00001);
    stepperSC->setSubcyclingInitTimeStep(dt / 10.0);
    stepperSC->setSubcyclingMaxTimeStep(dt / 10.0);
    stepperSC->setSubcyclingMaxFailures(10);
    stepperSC->setSubcyclingMaxConsecFailures(5);
    stepperSC->setSubcyclingScreenOutputIndexInterval(1);
    // stepper->setSubcyclingIntegratorObserver(
    //   Teuchos::rcp(new Tempus::IntegratorObserverSubcycling<double>()));
    // stepperSC->setSubcyclingPrintDtChanges   (true);

    auto strategySC = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
    strategySC->setMinEta(0.000001);
    strategySC->setMaxEta(0.01);
    strategySC->initialize();
    stepperSC->setSubcyclingTimeStepControlStrategy(strategySC);

    // Implicit Stepper
    auto stepperBE =
        Tempus::createStepperBackwardEuler(implicitModel, Teuchos::null);

    // Operator-Split Stepper
    auto stepper = rcp(new Tempus::StepperOperatorSplit<double>());
    stepper->addStepper(stepperSC);
    stepper->addStepper(stepperBE);

    // Setup TimeStepControl ------------------------------------
    auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
    timeStepControl->setInitIndex(0);
    timeStepControl->setInitTime(0.0);
    // timeStepControl->setFinalIndex(2);
    timeStepControl->setFinalTime(2.0);
    timeStepControl->setMinTimeStep(0.000001);
    timeStepControl->setInitTimeStep(dt);
    timeStepControl->setMaxTimeStep(dt);

    // timeStepControl->setInitTimeStep(dt/2.0);
    // timeStepControl->setMaxTimeStep (dt);
    // auto strategy = rcp(new
    // Tempus::TimeStepControlStrategyBasicVS<double>());
    // strategy->setMinEta(1.0e-6);
    // strategy->setMaxEta(5.0);
    // strategy->initialize();
    // timeStepControl->getTimeStepControlStrategy()->clearObservers();
    // timeStepControl->getTimeStepControlStrategy()->addStrategy(strategy);

    timeStepControl->initialize();

    // Setup initial condition SolutionState --------------------
    auto inArgsIC = stepper->getModel()->getNominalValues();
    auto icX      = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
    auto icXDot =
        rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot());
    auto icState = Tempus::createSolutionStateX(icX, icXDot);
    icState->setTime(timeStepControl->getInitTime());
    icState->setIndex(timeStepControl->getInitIndex());
    icState->setTimeStep(0.0);  // dt for ICs are zero.
    icState->setOrder(stepper->getOrder());
    icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

    // Setup SolutionHistory ------------------------------------
    auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
    solutionHistory->setName("Forward States");
    solutionHistory->setStorageType(Tempus::STORAGE_TYPE_UNLIMITED);
    // solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
    solutionHistory->setStorageLimit(3);
    solutionHistory->addState(icState);

    // Ensure ICs are consistent and stepper memory is set (e.g., xDot is set).
    stepperSC->setInitialConditions(solutionHistory);
    stepper->initialize();

    // Setup Integrator -----------------------------------------
    integrator = Tempus::createIntegratorBasic<double>();
    integrator->setStepper(stepper);
    integrator->setTimeStepControl(timeStepControl);
    integrator->setSolutionHistory(solutionHistory);
    integrator->setScreenOutputIndexInterval(10);
    // integrator->setObserver(...);
    integrator->initialize();

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time             = integrator->getTime();
    double timeFinal = 2.0;
    double tol       = 100.0 * std::numeric_limits<double>::epsilon();
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
      std::string fname = "Tempus_Subcycling_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_Subcycling_VanDerPol.dat";
      writeSolution(fname, integrator->getSolutionHistory());
      // solutionHistory->printHistory("medium");
    }
  }

  // Check the order and intercept
  double xSlope                        = 0.0;
  double xDotSlope                     = 0.0;
  RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
  // double order = stepper->getOrder();
  writeOrderError("Tempus_Subcycling_VanDerPol-Error.dat", stepper, StepSize,
                  solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                  xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, 1.25708, 0.05);
  TEST_FLOATING_EQUALITY(xDotSlope, 1.20230, 0.05);
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.37156, 1.0e-4);
  TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 3.11651, 1.0e-4);

  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
