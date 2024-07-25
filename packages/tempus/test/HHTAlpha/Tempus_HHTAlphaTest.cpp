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

#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_StepperHHTAlpha.hpp"

#include "../TestModels/HarmonicOscillatorModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <limits>
#include <sstream>
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
TEUCHOS_UNIT_TEST(HHTAlpha, BallParabolic)
{
  // Tolerance to check if test passed
  double tolerance = 1.0e-14;
  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_HHTAlpha_BallParabolic.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  auto model                = rcp(new HarmonicOscillatorModel<double>(hom_pl));

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  RCP<Tempus::IntegratorBasic<double>> integrator =
      Tempus::createIntegratorBasic<double>(pl, model);

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
  RCP<const Thyra::VectorBase<double>> x_exact =
      model->getExactSolution(time).get_x();

  // Plot sample solution and exact solution
  std::ofstream ftmp("Tempus_HHTAlpha_BallParabolic.dat");
  ftmp.precision(16);
  RCP<const SolutionHistory<double>> solutionHistory =
      integrator->getSolutionHistory();
  bool passed = true;
  double err  = 0.0;
  RCP<const Thyra::VectorBase<double>> x_exact_plot;
  for (int i = 0; i < solutionHistory->getNumStates(); i++) {
    RCP<const SolutionState<double>> solutionState = (*solutionHistory)[i];
    double time_i                                  = solutionState->getTime();
    RCP<const Thyra::VectorBase<double>> x_plot    = solutionState->getX();
    x_exact_plot                                   = model->getExactSolution(time_i).get_x();
    ftmp << time_i << "   " << get_ele(*(x_plot), 0) << "   "
         << get_ele(*(x_exact_plot), 0) << std::endl;
    if (abs(get_ele(*(x_plot), 0) - get_ele(*(x_exact_plot), 0)) > err)
      err = abs(get_ele(*(x_plot), 0) - get_ele(*(x_exact_plot), 0));
  }
  ftmp.close();
  out << "Max error = " << err << "\n \n";
  if (err > tolerance) passed = false;

  TEUCHOS_TEST_FOR_EXCEPTION(
      !passed, std::logic_error,
      "\n Test failed!  Max error = " << err << " > tolerance = " << tolerance
                                      << "\n!");
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(HHTAlpha, ConstructingFromDefaults)
{
  double dt = 0.05;

  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_HHTAlpha_SinCos_SecondOrder.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  auto model                = rcp(new HarmonicOscillatorModel<double>(hom_pl));

  // Setup Stepper for field solve ----------------------------
  auto stepper = rcp(new Tempus::StepperHHTAlpha<double>());
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
  auto inArgsIC = model->getNominalValues();
  auto icX      = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
  auto icXDot   = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot());
  auto icXDotDot =
      rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot_dot());
  auto icState = Tempus::createSolutionStateX(icX, icXDot, icXDotDot);
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
  RCP<const Thyra::VectorBase<double>> x_exact =
      model->getExactSolution(time).get_x();

  // Calculate the error
  RCP<Thyra::VectorBase<double>> xdiff = x->clone_v();
  Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));

  // Check the order and intercept
  out << "  Stepper = " << stepper->description() << std::endl;
  out << "  =========================" << std::endl;
  out << "  Exact solution   : " << get_ele(*(x_exact), 0) << std::endl;
  out << "  Computed solution: " << get_ele(*(x), 0) << std::endl;
  out << "  Difference       : " << get_ele(*(xdiff), 0) << std::endl;
  out << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), 0.144918, 1.0e-4);
}

// ************************************************************
TEUCHOS_UNIT_TEST(HHTAlpha, SinCos_SecondOrder)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 7;
  double time              = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_HHTAlpha_SinCos_SecondOrder.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  auto model                = rcp(new HarmonicOscillatorModel<double>(hom_pl));

  // Get k and m coefficients from model, needed for computing energy
  double k = hom_pl->get<double>("x coeff k");
  double m = hom_pl->get<double>("Mass coeff m");

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Set initial time step = 2*dt specified in input file (for convergence
  // study)
  //
  double dt = pl->sublist("Default Integrator")
                  .sublist("Time Step Control")
                  .get<double>("Initial Time Step");
  dt *= 2.0;

  for (int n = 0; n < nTimeStepSizes; n++) {
    // Perform time-step refinement
    dt /= 2;
    out << "\n \n time step #" << n << " (out of " << nTimeStepSizes - 1
        << "), dt = " << dt << "\n";
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
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double>> x = integrator->getX();
    RCP<const Thyra::VectorBase<double>> x_exact =
        model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution at most-refined resolution
    if (n == nTimeStepSizes - 1) {
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      writeSolution("Tempus_HHTAlpha_SinCos_SecondOrder.dat", solutionHistory);

      auto solnHistExact = rcp(new Tempus::SolutionHistory<double>());
      for (int i = 0; i < solutionHistory->getNumStates(); i++) {
        double time_i = (*solutionHistory)[i]->getTime();
        auto state    = Tempus::createSolutionStateX(
               rcp_const_cast<Thyra::VectorBase<double>>(
                model->getExactSolution(time_i).get_x()),
               rcp_const_cast<Thyra::VectorBase<double>>(
                model->getExactSolution(time_i).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution("Tempus_HHTAlpha_SinCos_SecondOrder-Ref.dat",
                    solnHistExact);

      // Problem specific output
      {
        std::ofstream ftmp("Tempus_HHTAlpha_SinCos_SecondOrder-Energy.dat");
        ftmp.precision(16);
        RCP<const Thyra::VectorBase<double>> x_exact_plot;
        for (int i = 0; i < solutionHistory->getNumStates(); i++) {
          RCP<const SolutionState<double>> solutionState =
              (*solutionHistory)[i];
          double time_i                               = solutionState->getTime();
          RCP<const Thyra::VectorBase<double>> x_plot = solutionState->getX();
          RCP<const Thyra::VectorBase<double>> x_dot_plot =
              solutionState->getXDot();
          x_exact_plot = model->getExactSolution(time_i).get_x();
          // kinetic energy = 0.5*m*xdot*xdot
          double ke = Thyra::dot(*x_dot_plot, *x_dot_plot);
          ke *= 0.5 * m;
          // potential energy = 0.5*k*x*x
          double pe = Thyra::dot(*x_plot, *x_plot);
          pe *= 0.5 * k;
          double te = ke + pe;
          // Output to file the following:
          //[time, x computed, x exact, xdot computed, ke, pe, te]
          ftmp << time_i << "   " << get_ele(*(x_plot), 0) << "   "
               << get_ele(*(x_exact_plot), 0) << "   "
               << get_ele(*(x_dot_plot), 0) << "   " << ke << "   " << pe
               << "   " << te << std::endl;
        }
        ftmp.close();
      }
    }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()), solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXDot()), solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes - 1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solutionExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),
                  solutionExact.ptr());
      solutions.push_back(solutionExact);
      auto solutionDotExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDotExact.ptr());
      solutionsDot.push_back(solutionDotExact);
    }
  }

  // Check the order and intercept
  double xSlope                        = 0.0;
  double xDotSlope                     = 0.0;
  RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
  double order                         = stepper->getOrder();
  writeOrderError("Tempus_HHTAlpha_SinCos_SecondOrder-Error.dat", stepper,
                  StepSize, solutions, xErrorNorm, xSlope, solutionsDot,
                  xDotErrorNorm, xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, order, 0.02);
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.00644755, 1.0e-4);
  TEST_FLOATING_EQUALITY(xDotSlope, order, 0.01);
  TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 0.104392, 1.0e-4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(HHTAlpha, SinCos_FirstOrder)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 7;
  double time              = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_HHTAlpha_SinCos_FirstOrder.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  auto model                = rcp(new HarmonicOscillatorModel<double>(hom_pl));

  // Get k and m coefficients from model, needed for computing energy
  double k = hom_pl->get<double>("x coeff k");
  double m = hom_pl->get<double>("Mass coeff m");

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Set initial time step = 2*dt specified in input file (for convergence
  // study)
  //
  double dt = pl->sublist("Default Integrator")
                  .sublist("Time Step Control")
                  .get<double>("Initial Time Step");
  dt *= 2.0;

  for (int n = 0; n < nTimeStepSizes; n++) {
    // Perform time-step refinement
    dt /= 2;
    out << "\n \n time step #" << n << " (out of " << nTimeStepSizes - 1
        << "), dt = " << dt << "\n";
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
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double>> x = integrator->getX();
    RCP<const Thyra::VectorBase<double>> x_exact =
        model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution at most-refined resolution
    if (n == nTimeStepSizes - 1) {
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      writeSolution("Tempus_HHTAlpha_SinCos_FirstOrder.dat", solutionHistory);

      auto solnHistExact = rcp(new Tempus::SolutionHistory<double>());
      for (int i = 0; i < solutionHistory->getNumStates(); i++) {
        double time_i = (*solutionHistory)[i]->getTime();
        auto state    = Tempus::createSolutionStateX(
               rcp_const_cast<Thyra::VectorBase<double>>(
                model->getExactSolution(time_i).get_x()),
               rcp_const_cast<Thyra::VectorBase<double>>(
                model->getExactSolution(time_i).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution("Tempus_HHTAlpha_SinCos_FirstOrder-Ref.dat", solnHistExact);

      // Problem specific output
      {
        std::ofstream ftmp("Tempus_HHTAlpha_SinCos_FirstOrder-Energy.dat");
        ftmp.precision(16);
        RCP<const Thyra::VectorBase<double>> x_exact_plot;
        for (int i = 0; i < solutionHistory->getNumStates(); i++) {
          RCP<const SolutionState<double>> solutionState =
              (*solutionHistory)[i];
          double time_i                               = solutionState->getTime();
          RCP<const Thyra::VectorBase<double>> x_plot = solutionState->getX();
          RCP<const Thyra::VectorBase<double>> x_dot_plot =
              solutionState->getXDot();
          x_exact_plot = model->getExactSolution(time_i).get_x();
          // kinetic energy = 0.5*m*xdot*xdot
          double ke = Thyra::dot(*x_dot_plot, *x_dot_plot);
          ke *= 0.5 * m;
          // potential energy = 0.5*k*x*x
          double pe = Thyra::dot(*x_plot, *x_plot);
          pe *= 0.5 * k;
          double te = ke + pe;
          // Output to file the following:
          //[time, x computed, x exact, xdot computed, ke, pe, te]
          ftmp << time_i << "   " << get_ele(*(x_plot), 0) << "   "
               << get_ele(*(x_exact_plot), 0) << "   "
               << get_ele(*(x_dot_plot), 0) << "   " << ke << "   " << pe
               << "   " << te << std::endl;
        }
        ftmp.close();
      }
    }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()), solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXDot()), solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes - 1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solutionExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),
                  solutionExact.ptr());
      solutions.push_back(solutionExact);
      auto solutionDotExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDotExact.ptr());
      solutionsDot.push_back(solutionDotExact);
    }
  }

  // Check the order and intercept
  double xSlope                        = 0.0;
  double xDotSlope                     = 0.0;
  RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
  // double order = stepper->getOrder();
  writeOrderError("Tempus_HHTAlpha_SinCos_FirstOrder-Error.dat", stepper,
                  StepSize, solutions, xErrorNorm, xSlope, solutionsDot,
                  xDotErrorNorm, xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, 0.977568, 0.02);
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.048932, 1.0e-4);
  TEST_FLOATING_EQUALITY(xDotSlope, 1.2263, 0.01);
  TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 0.393504, 1.0e-4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(HHTAlpha, SinCos_CD)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 7;
  double time              = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_HHTAlpha_SinCos_ExplicitCD.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  auto model                = rcp(new HarmonicOscillatorModel<double>(hom_pl));

  // Get k and m coefficients from model, needed for computing energy
  double k = hom_pl->get<double>("x coeff k");
  double m = hom_pl->get<double>("Mass coeff m");

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Set initial time step = 2*dt specified in input file (for convergence
  // study)
  //
  double dt = pl->sublist("Default Integrator")
                  .sublist("Time Step Control")
                  .get<double>("Initial Time Step");
  dt *= 2.0;

  for (int n = 0; n < nTimeStepSizes; n++) {
    // Perform time-step refinement
    dt /= 2;
    out << "\n \n time step #" << n << " (out of " << nTimeStepSizes - 1
        << "), dt = " << dt << "\n";
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
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double>> x = integrator->getX();
    RCP<const Thyra::VectorBase<double>> x_exact =
        model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution at most-refined resolution
    if (n == nTimeStepSizes - 1) {
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      writeSolution("Tempus_HHTAlpha_SinCos_ExplicitCD.dat", solutionHistory);

      auto solnHistExact = rcp(new Tempus::SolutionHistory<double>());
      for (int i = 0; i < solutionHistory->getNumStates(); i++) {
        double time_i = (*solutionHistory)[i]->getTime();
        auto state    = Tempus::createSolutionStateX(
               rcp_const_cast<Thyra::VectorBase<double>>(
                model->getExactSolution(time_i).get_x()),
               rcp_const_cast<Thyra::VectorBase<double>>(
                model->getExactSolution(time_i).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution("Tempus_HHTAlpha_SinCos_ExplicitCD-Ref.dat", solnHistExact);

      // Problem specific output
      {
        std::ofstream ftmp("Tempus_HHTAlpha_SinCos_ExplicitCD-Energy.dat");
        ftmp.precision(16);
        RCP<const Thyra::VectorBase<double>> x_exact_plot;
        for (int i = 0; i < solutionHistory->getNumStates(); i++) {
          RCP<const SolutionState<double>> solutionState =
              (*solutionHistory)[i];
          double time_i                               = solutionState->getTime();
          RCP<const Thyra::VectorBase<double>> x_plot = solutionState->getX();
          RCP<const Thyra::VectorBase<double>> x_dot_plot =
              solutionState->getXDot();
          x_exact_plot = model->getExactSolution(time_i).get_x();
          // kinetic energy = 0.5*m*xdot*xdot
          double ke = Thyra::dot(*x_dot_plot, *x_dot_plot);
          ke *= 0.5 * m;
          // potential energy = 0.5*k*x*x
          double pe = Thyra::dot(*x_plot, *x_plot);
          pe *= 0.5 * k;
          double te = ke + pe;
          // Output to file the following:
          //[time, x computed, x exact, xdot computed, ke, pe, te]
          ftmp << time_i << "   " << get_ele(*(x_plot), 0) << "   "
               << get_ele(*(x_exact_plot), 0) << "   "
               << get_ele(*(x_dot_plot), 0) << "   " << ke << "   " << pe
               << "   " << te << std::endl;
        }
        ftmp.close();
      }
    }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()), solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXDot()), solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes - 1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solutionExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),
                  solutionExact.ptr());
      solutions.push_back(solutionExact);
      auto solutionDotExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDotExact.ptr());
      solutionsDot.push_back(solutionDotExact);
    }
  }

  // Check the order and intercept
  double xSlope                        = 0.0;
  double xDotSlope                     = 0.0;
  RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
  double order                         = stepper->getOrder();
  writeOrderError("Tempus_HHTAlpha_SinCos_ExplicitCD-Error.dat", stepper,
                  StepSize, solutions, xErrorNorm, xSlope, solutionsDot,
                  xDotErrorNorm, xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, order, 0.02);
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.00451069, 1.0e-4);
  TEST_FLOATING_EQUALITY(xDotSlope, order, 0.01);
  TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 0.0551522, 1.0e-4);
}

}  // namespace Tempus_Test
