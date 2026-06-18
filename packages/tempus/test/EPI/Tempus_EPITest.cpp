//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Teuchos_LocalTestingHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DetachedVectorView.hpp"

#include "Tempus_IntegratorBasic.hpp"

#include "Tempus_StepperEPI.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/ReactionModel.hpp"
#include "../TestModels/NonAutoSrcModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestModels/LotkaVolterraModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"
#include "Thyra_VectorStdOps_decl.hpp"

#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

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
TEUCHOS_UNIT_TEST(EPI, SinCos)
{
  // Run EPI integrator logic with different PhiEvaluator configurations
  std::vector<std::string> xml_cases = {
    "Taylor",
    "Leja"
  };

  for (const auto& xml_case : xml_cases ){
    RCP<Tempus::IntegratorBasic<double>> integrator;
    std::vector<RCP<Thyra::VectorBase<double>>> solutions;
    std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
    std::vector<double> StepSize;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;
    const int nTimeStepSizes = 4;
    double dt                = 2.0;
    double time              = 0.0;
    int expected_order;
    for (int n = 0; n < nTimeStepSizes; n++) {
      // Read params from .xml file
      RCP<ParameterList> pList =
          getParametersFromXmlFile("Tempus_EPI_SinCos.xml");

      // Setup the SinCosModel
      RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
      // RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
      auto model = rcp(new SinCosModel<double>(scm_pl));

      dt /= 2;

      // Setup the Integrator and reset initial time step
      RCP<ParameterList> pl = sublist(pList, "Tempus", true);
      pl->sublist("Demo Integrator")
          .sublist("Time Step Control")
          .set("Initial Time Step", dt);

      // Update the PhiEvaluator parameters based on the PhiEvaluator type
      auto& phiList = pl->sublist("Demo Stepper").sublist("PhiEvaluator");
      if (xml_case == "Leja") {
        phiList.set("PhiEvaluator Type", "Leja")
               .set("Expansion Order", 30)
               .set("Leja DD Method", 2)
               .set("leja_a", -1.0)
               .set("leja_c", 0.001);
      }
      else if (xml_case == "Taylor") {
        phiList.remove("Leja DD Method", false);
        phiList.remove("leja_a", false);
        phiList.remove("leja_c", false);

        phiList.set("PhiEvaluator Type", "Taylor")
               .set("Expansion Order", 20);
      }

      integrator = Tempus::createIntegratorBasic<double>(pl, model);
      // Initial Conditions
      // During the Integrator construction, the initial SolutionState
      // is set by default to model->getNominalVales().get_x().  However,
      // the application can set it also by integrator->initializeSolutionHistory.
      RCP<Thyra::VectorBase<double>> x0 =
          model->getNominalValues().get_x()->clone_v();

      integrator->initializeSolutionHistory(0.0, x0);
      integrator->initialize();

      // Integrate to timeMax
      bool integratorStatus = integrator->advanceTime();
      TEST_ASSERT(integratorStatus)

      // Test PhysicsState
      RCP<Tempus::PhysicsState<double>> physicsState =
          integrator->getSolutionHistory()->getCurrentState()->getPhysicsState();
      TEST_EQUALITY(physicsState->getName(), "Tempus::PhysicsState");

      // Test if at 'Final Time'
      time             = integrator->getTime();
      double timeFinal = pl->sublist("Demo Integrator")
                            .sublist("Time Step Control")
                            .get<double>("Final Time");
      TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

      // Time-integrated solution and the exact solution
      RCP<Thyra::VectorBase<double>> x = integrator->getX();
      RCP<const Thyra::VectorBase<double>> x_exact =
          model->getExactSolution(time).get_x();

      // Plot sample solution and exact solution
      if (n == 0) {
        RCP<const SolutionHistory<double>> solutionHistory =
            integrator->getSolutionHistory();
        writeSolution("Tempus_EPI_SinCos_" + xml_case + ".dat", solutionHistory);

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
        writeSolution("Tempus_EPI_SinCos-Ref_" + xml_case + ".dat", solnHistExact);
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

    // compute difference between expected and computed
    auto gold = solutions[solutions.size()-1];
    for (int i=0; i<solutions.size()-1; ++i) {
      auto calc = solutions[i];
      Thyra::Vp_StV(calc.ptr(), -1.0, *gold);
      TEST_COMPARE(calc->norm_2(), <=, 1e-8);
    }

    Teuchos::TimeMonitor::summarize();
  }
}

TEUCHOS_UNIT_TEST(EPI, Reaction)
{
  // Run EPI integrator logic with different PhiEvaluator configurations
  std::vector<std::string> xml_cases = {
    "Taylor",
    "Leja"
  };

  for (const auto& xml_case : xml_cases ){
    RCP<Tempus::IntegratorBasic<double>> integrator;
    std::vector<RCP<Thyra::VectorBase<double>>> solutions;
    std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
    std::vector<double> StepSize;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;
    const int nTimeStepSizes = 7;
    double dt                = 0.2;
    double time              = 0.0;
    // int expected_order;
    for (int n = 0; n < nTimeStepSizes; n++) {
      // Read params from .xml file
      RCP<ParameterList> pList = getParametersFromXmlFile("Tempus_EPI_Reaction.xml");

      // Setup the ReactionModel
      RCP<ParameterList> scm_pl = sublist(pList, "ReactionModel", true);
      // RCP<ReactionModel<double> > model = sineCosineModel(scm_pl);
      auto model = rcp(new ReactionModel<double>(scm_pl));

      dt /= 2;

      // Setup the Integrator and reset initial time step
      RCP<ParameterList> pl = sublist(pList, "Tempus", true);
      pl->sublist("Demo Integrator")
          .sublist("Time Step Control")
          .set("Initial Time Step", dt);
      
      // Update the PhiEvaluator parameters based on the PhiEvaluator type
      if (xml_case == "Leja") {
        pl->sublist("Demo Stepper").sublist("PhiEvaluator")
            .set("PhiEvaluator Type", "Leja")
            .set("Expansion Order", 50)
            .set("leja_tol", 1e-12)
            .set("leja_a", -1.0)
            .set("leja_c", 0.001);
      }
      else if (xml_case == "Taylor") {
        pl->sublist("Demo Stepper").sublist("PhiEvaluator")
            .set("PhiEvaluator Type", "Taylor")
            .set("Expansion Order", 16);
      }

      integrator = Tempus::createIntegratorBasic<double>(pl, model);
      // Initial Conditions
      // During the Integrator construction, the initial SolutionState
      // is set by default to model->getNominalVales().get_x().  However,
      // the application can set it also by integrator->initializeSolutionHistory.
      RCP<Thyra::VectorBase<double>> x0 =
          model->getNominalValues().get_x()->clone_v();

      integrator->initializeSolutionHistory(0.0, x0);
      integrator->initialize();

      // Integrate to timeMax
      bool integratorStatus = integrator->advanceTime();
      TEST_ASSERT(integratorStatus)

      // Test PhysicsState
      RCP<Tempus::PhysicsState<double>> physicsState =
          integrator->getSolutionHistory()->getCurrentState()->getPhysicsState();
      TEST_EQUALITY(physicsState->getName(), "Tempus::PhysicsState");

      // Test if at 'Final Time'
      time             = integrator->getTime();
      double timeFinal = pl->sublist("Demo Integrator")
                            .sublist("Time Step Control")
                            .get<double>("Final Time");
      TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

      // Time-integrated solution and the exact solution
      RCP<Thyra::VectorBase<double>> x = integrator->getX();
      RCP<const Thyra::VectorBase<double>> x_exact =
          model->getExactSolution(time).get_x();

      // Plot sample solution and exact solution
      if (n == nTimeStepSizes - 1) {
        RCP<const SolutionHistory<double>> solutionHistory =
            integrator->getSolutionHistory();
        writeSolution("Tempus_EPI_Reaction_" + xml_case + ".dat", solutionHistory);
      // solutionHistory->printHistory("high");

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
        writeSolution("Tempus_EPI_Reaction_" + xml_case + "-Ref.dat", solnHistExact);
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

    // Check against exact solution
    double xSlope                        = 0.0;
    double xDotSlope                     = 0.0;
    RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
    writeOrderError("Tempus_EPI_Reaction_" + xml_case + "-Error.dat", stepper, StepSize,
                    solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                    xDotSlope, out);

    for (int i=0; i < nTimeStepSizes - 1; ++i) {
      // Linear problem, expect near exact solution from exp. integrator
      // for all time step sizes.
      TEST_COMPARE(xErrorNorm[i], <=, 1.0e-10);
    }

    Teuchos::TimeMonitor::summarize();
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, VanDerPol)
{
  // Run EPI integrator logic with different PhiEvaluator configurations
  std::vector<std::string> xml_cases = {
    "Taylor",
    "Leja"
  };

  for (const auto& xml_case : xml_cases) {
    RCP<Tempus::IntegratorBasic<double>> integrator;
    std::vector<RCP<Thyra::VectorBase<double>>> solutions;
    std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
    std::vector<double> StepSize;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;
    const int nTimeStepSizes = 4;  // 8 for Error plot
    double dt                = 1.0;
    for (int n = 0; n < nTimeStepSizes; n++) {
      // Read params from .xml file
      RCP<ParameterList> pList =
          getParametersFromXmlFile("Tempus_EPI_VanDerPol.xml");

      // Setup the VanDerPolModel
      RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
      auto model                 = rcp(new VanDerPolModel<double>(vdpm_pl));

      // Set the step size
      dt /= 2;
      if (n == nTimeStepSizes - 1) dt /= 5.0;

      // Setup the Integrator and reset initial time step
      RCP<ParameterList> pl = sublist(pList, "Tempus", true);
      pl->sublist("Demo Integrator")
          .sublist("Time Step Control")
          .set("Initial Time Step", dt);

      // Update the PhiEvaluator parameters based on the PhiEvaluator type
      auto& phiList = pl->sublist("Demo Stepper").sublist("PhiEvaluator");
      if (xml_case == "Leja") {
        phiList.set("PhiEvaluator Type", "Leja")
               .set("Expansion Order", 30)
               .set("Leja DD Method", 2)
               .set("leja_tol", 1e-8)
               .set("leja_a", -1.0)
               .set("leja_c", 0.5);
      }
      else if (xml_case == "Taylor") {
        phiList.remove("Leja DD Method", false);
        phiList.remove("leja_tol", false);
        phiList.remove("leja_a", false);
        phiList.remove("leja_c", false);

        phiList.set("PhiEvaluator Type", "Taylor")
               .set("Expansion Order", 30);
      }

      integrator = Tempus::createIntegratorBasic<double>(pl, model);

      // Integrate to timeMax
      bool integratorStatus = integrator->advanceTime();
      TEST_ASSERT(integratorStatus)

      // Test if at 'Final Time'
      double time      = integrator->getTime();
      double timeFinal = pl->sublist("Demo Integrator")
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
        std::string fname = "Tempus_EPI_VanDerPol_" + xml_case + "-Ref.dat";
        if (n == 0) fname = "Tempus_EPI_VanDerPol_" + xml_case + ".dat";
        RCP<const SolutionHistory<double>> solutionHistory =
            integrator->getSolutionHistory();
        writeSolution(fname, solutionHistory);
      }
    }

    // Check the order and intercept
    double xSlope                        = 0.0;
    double xDotSlope                     = 0.0;
    RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
    double order                         = stepper->getOrder();

    // xDot not yet available for EPI methods, e.g., are not calc. and zero.
    solutionsDot.clear();

    writeOrderError("Tempus_EPI_VanDerPol_" + xml_case + "-Error.dat", stepper, StepSize,
                    solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                    xDotSlope, out);

    TEST_FLOATING_EQUALITY(xSlope, order, 0.2);
    TEST_COMPARE(xErrorNorm[0], <=, .1);

    Teuchos::TimeMonitor::summarize();
  }
}

TEUCHOS_UNIT_TEST(EPI, NonAutoSrc)
{
  // Run EPI integrator:
  // ii = 0 -> without nonautonomous correction (Epsilon for RHS finite difference<=0)
  // ii = 1 -> with nonautonomous correction (Epsilon for RHS finite difference>0)
  for (int ii = 0; ii < 2; ++ii) {
    std::string xml_case = "Tempus_EPI_NonAutoSrc";
    RCP<Tempus::IntegratorBasic<double>> integrator;
    std::vector<RCP<Thyra::VectorBase<double>>> solutions;
    std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
    std::vector<double> StepSize;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;
    const int nTimeStepSizes = 7;
    double dt                = 0.2;
    double time              = 0.0;
    for (int n = 0; n < nTimeStepSizes; n++) {
      // Read params from .xml file
      RCP<ParameterList> pList = getParametersFromXmlFile(xml_case + ".xml");

      // Setup the NonAutoSrcModel
      RCP<ParameterList> scm_pl = sublist(pList, "NonAutoSrcModel", true);
      auto model = rcp(new NonAutoSrcModel<double>(scm_pl));

      dt /= 2;

      // Setup the Integrator and reset initial time step
      RCP<ParameterList> pl = sublist(pList, "Tempus", true);
      pl->sublist("Demo Integrator")
          .sublist("Time Step Control")
          .set("Initial Time Step", dt);
      if (ii == 0)
        pl->sublist("Demo Stepper")
          .set("Epsilon for RHS finite difference", -1.0);
      else
        pl->sublist("Demo Stepper")
            .set("Epsilon for RHS finite difference", 1e-4);

      integrator = Tempus::createIntegratorBasic<double>(pl, model);
      // Initial Conditions
      // During the Integrator construction, the initial SolutionState
      // is set by default to model->getNominalVales().get_x().  However,
      // the application can set it also by integrator->initializeSolutionHistory.
      RCP<Thyra::VectorBase<double>> x0 =
          model->getNominalValues().get_x()->clone_v();

      integrator->initializeSolutionHistory(0.0, x0);
      integrator->initialize();

      // Integrate to timeMax
      bool integratorStatus = integrator->advanceTime();
      TEST_ASSERT(integratorStatus)

      // Test PhysicsState
      RCP<Tempus::PhysicsState<double>> physicsState =
          integrator->getSolutionHistory()->getCurrentState()->getPhysicsState();
      TEST_EQUALITY(physicsState->getName(), "Tempus::PhysicsState");

      // Test if at 'Final Time'
      time             = integrator->getTime();
      double timeFinal = pl->sublist("Demo Integrator")
                            .sublist("Time Step Control")
                            .get<double>("Final Time");
      TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

      // Time-integrated solution and the exact solution
      RCP<Thyra::VectorBase<double>> x = integrator->getX();
      RCP<const Thyra::VectorBase<double>> x_exact =
          model->getExactSolution(time).get_x();

      // Plot sample solution and exact solution
      if (n == nTimeStepSizes - 1) {
        RCP<const SolutionHistory<double>> solutionHistory =
            integrator->getSolutionHistory();
        writeSolution(xml_case + ".dat", solutionHistory);
      // solutionHistory->printHistory("high");

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
        writeSolution(xml_case + "-Ref.dat", solnHistExact);
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

    // Check against exact solution
    double xSlope                        = 0.0;
    double xDotSlope                     = 0.0;
    RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
    writeOrderError(xml_case + "-Error.dat", stepper, StepSize,
                    solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                    xDotSlope, out);
    if (ii == 0)
    {
      // Without nonautonomous correction, expect order ~1.0
      TEST_FLOATING_EQUALITY(xSlope, 1.0, 0.2);
    }
    else
    {
      // With nonautonomous correction, expect full order
      TEST_FLOATING_EQUALITY(xSlope, stepper->getOrder(), 0.2);
    }
    for (int i=0; i < nTimeStepSizes - 1; ++i) {
      TEST_COMPARE(xErrorNorm[i], <=, 1.0e-10);
    }

    Teuchos::TimeMonitor::summarize();
  }
}

TEUCHOS_UNIT_TEST(EPI, LotkaVolterra)
{
  // Convergence study for EPI3 (third-order exponential propagation iterative).
  // A RK4 solution at dt_ref = 1e-3 is pre-computed and used
  // as the reference.

  // Compute the fine RK4 reference solution
  const double dt_ref = 1.0e-3;

  RCP<Thyra::VectorBase<double>> xRef;
  double timeFinal = 0.0;
  {
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_RK4_LotkaVolterra.xml");

    RCP<ParameterList> lvm_pl = sublist(pList, "LotkaVolterraModel", true);
    auto modelRef             = rcp(new LotkaVolterraModel<double>(lvm_pl));

    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt_ref);
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Maximum Time Step", dt_ref);

    auto integratorRef = Tempus::createIntegratorBasic<double>(pl, modelRef);

    RCP<Thyra::VectorBase<double>> x0 =
        modelRef->getNominalValues().get_x()->clone_v();
    integratorRef->initializeSolutionHistory(0.0, x0);
    integratorRef->initialize();

    bool status = integratorRef->advanceTime();
    TEST_ASSERT(status)

    timeFinal = pl->sublist("Demo Integrator")
                    .sublist("Time Step Control")
                    .get<double>("Final Time");
    TEST_FLOATING_EQUALITY(integratorRef->getTime(), timeFinal, 1.0e-14);

    xRef = Thyra::createMember(modelRef->get_x_space());
    Thyra::copy(*(integratorRef->getX()), xRef.ptr());

    out << "  RK4 reference computed at dt_ref = " << dt_ref << std::endl;
  }

  // EPI3 convergence study
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;

  const int nTimeStepSizes = 5;
  double dt                = 0.2;

  for (int n = 0; n < nTimeStepSizes; n++) {
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_EPI_LotkaVolterra.xml");

    RCP<ParameterList> lvm_pl = sublist(pList, "LotkaVolterraModel", true);
    auto model                = rcp(new LotkaVolterraModel<double>(lvm_pl));

    dt /= 2.0;

    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Maximum Time Step", dt);

    integrator = Tempus::createIntegratorBasic<double>(pl, model);

    RCP<Thyra::VectorBase<double>> x0 =
        model->getNominalValues().get_x()->clone_v();
    integrator->initializeSolutionHistory(0.0, x0);
    integrator->initialize();

    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    double time = integrator->getTime();
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Write trajectory for the coarsest run
    if (n == 0) {
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      writeSolution("Tempus_EPI_LotkaVolterra.dat", solutionHistory);
    }

    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()), solution.ptr());
    solutions.push_back(solution);
  }

  // Append the RK4 reference as the last entry.
  StepSize.push_back(0.0);
  solutions.push_back(xRef);

  // ----------------------------------------------------------
  // NOTE: EPI methods do not populate xDot
  // at the final state, so only the state convergence slope is checked.
  // EPI methods do not compute xDot at the final state; only check x slope.
  double xSlope                        = 0.0;
  RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
  double order                         = stepper->getOrder();
  writeOrderError("Tempus_EPI_LotkaVolterra-Error.dat", stepper, StepSize,
                  solutions, xErrorNorm, xSlope, out);

  std::cout << "Estimated Order: " << xSlope << std::endl;
  std::cout << "Expected  Order: " << order << std::endl;

  // ~order is expected, super-convergence OK
  TEST_COMPARE(xSlope, >=, order * 0.9);

  Teuchos::TimeMonitor::summarize();
}



}  // namespace Tempus_Test
