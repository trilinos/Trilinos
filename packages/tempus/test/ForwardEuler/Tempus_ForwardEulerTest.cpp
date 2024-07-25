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
#include "Thyra_DetachedVectorView.hpp"

#include "Tempus_IntegratorBasic.hpp"

#include "Tempus_StepperForwardEuler.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
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
TEUCHOS_UNIT_TEST(ForwardEuler, ParameterList)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_ForwardEuler_SinCos.xml");

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
TEUCHOS_UNIT_TEST(ForwardEuler, ConstructingFromDefaults)
{
  double dt = 0.1;
  std::vector<std::string> options;
  options.push_back("useFSAL=true");
  options.push_back("useFSAL=false");
  options.push_back("ICConsistency and Check");

  for (const auto& option : options) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_ForwardEuler_SinCos.xml");
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    // RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
    auto model = rcp(new SinCosModel<double>(scm_pl));

    // Setup Stepper for field solve ----------------------------
    auto stepper = rcp(new Tempus::StepperForwardEuler<double>());
    stepper->setModel(model);
    if (option == "useFSAL=true")
      stepper->setUseFSAL(true);
    else if (option == "useFSAL=false")
      stepper->setUseFSAL(false);
    else if (option == "ICConsistency and Check") {
      stepper->setICConsistency("Consistent");
      stepper->setICConsistencyCheck(true);
    }
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
    auto inArgsIC = model()->getNominalValues();
    auto icSolution =
        rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
    auto icState = Tempus::createSolutionStateX(icSolution);
    icState->setTime(timeStepControl->getInitTime());
    icState->setIndex(timeStepControl->getInitIndex());
    icState->setTimeStep(0.0);
    icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

    // Setup SolutionHistory ------------------------------------
    auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
    solutionHistory->setName("Forward States");
    solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
    solutionHistory->setStorageLimit(2);
    solutionHistory->addState(icState);

    // Ensure ICs are consistent and stepper memory is set (e.g., xDot is set).
    stepper->setInitialConditions(solutionHistory);

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
    RCP<const Thyra::VectorBase<double>> x_exact =
        model->getExactSolution(time).get_x();

    // Calculate the error
    RCP<Thyra::VectorBase<double>> xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));

    // Check the order and intercept
    out << "  Stepper = " << stepper->description() << "\n            with "
        << option << std::endl;
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
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, SinCos)
{
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
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_ForwardEuler_SinCos.xml");

    // std::ofstream ftmp("PL.txt");
    // pList->print(ftmp);
    // ftmp.close();

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
      writeSolution("Tempus_ForwardEuler_SinCos.dat", solutionHistory);

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
      writeSolution("Tempus_ForwardEuler_SinCos-Ref.dat", solnHistExact);
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
  writeOrderError("Tempus_ForwardEuler_SinCos-Error.dat", stepper, StepSize,
                  solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                  xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, order, 0.01);
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.051123, 1.0e-4);
  // xDot not yet available for Forward Euler.
  // TEST_FLOATING_EQUALITY( xDotSlope,            order, 0.01   );
  // TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.0486418, 1.0e-4 );

  Teuchos::TimeMonitor::summarize();
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, VanDerPol)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 7;  // 8 for Error plot
  double dt                = 0.2;
  for (int n = 0; n < nTimeStepSizes; n++) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_ForwardEuler_VanDerPol.xml");

    // Setup the VanDerPolModel
    RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
    auto model                 = rcp(new VanDerPolModel<double>(vdpm_pl));

    // Set the step size
    dt /= 2;
    if (n == nTimeStepSizes - 1) dt /= 10.0;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);
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
      std::string fname = "Tempus_ForwardEuler_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_ForwardEuler_VanDerPol.dat";
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
  writeOrderError("Tempus_ForwardEuler_VanDerPol-Error.dat", stepper, StepSize,
                  solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                  xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, order, 0.10);
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.387476, 1.0e-4);
  // xDot not yet available for Forward Euler.
  // TEST_FLOATING_EQUALITY( xDotSlope,       1.74898, 0.10   );
  // TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 1.0038, 1.0e-4 );

  Teuchos::TimeMonitor::summarize();
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, NumberTimeSteps)
{
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  // const int nTimeStepSizes = 7;
  // double dt = 0.2;
  // double order = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_ForwardEuler_NumberOfTimeSteps.xml");

  // Setup the VanDerPolModel
  RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
  auto model                 = rcp(new VanDerPolModel<double>(vdpm_pl));

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // dt = pl->sublist("Demo Integrator")
  //         .sublist("Time Step Control")
  //         .get<double>("Initial Time Step");
  const int numTimeSteps = pl->sublist("Demo Integrator")
                               .sublist("Time Step Control")
                               .get<int>("Number of Time Steps");

  RCP<Tempus::IntegratorBasic<double>> integrator =
      Tempus::createIntegratorBasic<double>(pl, model);

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)

  // check that the number of time steps taken is whats is set
  // in the parameter list
  TEST_EQUALITY(numTimeSteps, integrator->getIndex());
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, Variable_TimeSteps)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_ForwardEuler_VanDerPol.xml");

  // Setup the VanDerPolModel
  RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
  auto model                 = rcp(new VanDerPolModel<double>(vdpm_pl));

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Set parameters for this test.
  pl->sublist("Demo Integrator")
      .sublist("Time Step Control")
      .set("Initial Time Step", 0.01);

  pl->sublist("Demo Integrator")
      .sublist("Time Step Control")
      .sublist("Time Step Control Strategy")
      .set("Reduction Factor", 0.9);
  pl->sublist("Demo Integrator")
      .sublist("Time Step Control")
      .sublist("Time Step Control Strategy")
      .set("Amplification Factor", 1.15);
  pl->sublist("Demo Integrator")
      .sublist("Time Step Control")
      .sublist("Time Step Control Strategy")
      .set("Minimum Value Monitoring Function", 0.05);
  pl->sublist("Demo Integrator")
      .sublist("Time Step Control")
      .sublist("Time Step Control Strategy")
      .set("Maximum Value Monitoring Function", 0.1);

  pl->sublist("Demo Integrator")
      .sublist("Solution History")
      .set("Storage Type", "Static");
  pl->sublist("Demo Integrator")
      .sublist("Solution History")
      .set("Storage Limit", 3);

  RCP<Tempus::IntegratorBasic<double>> integrator =
      Tempus::createIntegratorBasic<double>(pl, model);

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)

  // Check 'Final Time'
  double time      = integrator->getTime();
  double timeFinal = pl->sublist("Demo Integrator")
                         .sublist("Time Step Control")
                         .get<double>("Final Time");
  TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

  // Check TimeStep size
  auto state = integrator->getCurrentState();
  double dt  = state->getTimeStep();
  TEST_FLOATING_EQUALITY(dt, 0.008310677297208358, 1.0e-12);

  // Check number of time steps taken
  const int numTimeSteps = 60;
  TEST_EQUALITY(numTimeSteps, integrator->getIndex());

  // Time-integrated solution and the reference solution
  RCP<Thyra::VectorBase<double>> x     = integrator->getX();
  RCP<Thyra::VectorBase<double>> x_ref = x->clone_v();
  {
    Thyra::DetachedVectorView<double> x_ref_view(*x_ref);
    x_ref_view[0] = -1.931946840284863;
    x_ref_view[1] = 0.645346748303107;
  }

  // Calculate the error
  RCP<Thyra::VectorBase<double>> xdiff = x->clone_v();
  Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_ref, -1.0, *(x));

  // Check the solution
  out << "  Stepper = ForwardEuler" << std::endl;
  out << "  =========================" << std::endl;
  out << "  Reference solution: " << get_ele(*(x_ref), 0) << "   "
      << get_ele(*(x_ref), 1) << std::endl;
  out << "  Computed solution : " << get_ele(*(x), 0) << "   "
      << get_ele(*(x), 1) << std::endl;
  out << "  Difference        : " << get_ele(*(xdiff), 0) << "   "
      << get_ele(*(xdiff), 1) << std::endl;
  out << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), get_ele(*(x_ref), 0), 1.0e-12);
  TEST_FLOATING_EQUALITY(get_ele(*(x), 1), get_ele(*(x_ref), 1), 1.0e-12);
}

}  // namespace Tempus_Test
