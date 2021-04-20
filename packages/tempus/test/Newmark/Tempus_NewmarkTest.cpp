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

#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_StepperNewmarkImplicitAForm.hpp"
#include "Tempus_StepperNewmarkImplicitDForm.hpp"

#include "../TestModels/HarmonicOscillatorModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"


#ifdef Tempus_ENABLE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

namespace Tempus_Test {

using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;


// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkExplicitAForm, BallParabolic)
{
  //Tolerance to check if test passed
  double tolerance = 1.0e-14;
  std::vector<std::string> options;
  options.push_back("useFSAL=true");
  options.push_back("useFSAL=false");
  options.push_back("ICConsistency and Check");

  for(const auto& option: options) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_NewmarkExplicitAForm_BallParabolic.xml");

    // Setup the HarmonicOscillatorModel
    RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
    RCP<HarmonicOscillatorModel<double> > model =
      Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl));

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    RCP<ParameterList> stepperPL = sublist(pl, "Default Stepper", true);
    stepperPL->remove("Zero Initial Guess");
    if (option == "useFSAL=true") stepperPL->set("Use FSAL", true);
    else if (option == "useFSAL=false") stepperPL->set("Use FSAL", false);
    else if (option == "ICConsistency and Check") {
      stepperPL->set("Initial Condition Consistency", "Consistent");
      stepperPL->set("Initial Condition Consistency Check", true);
    }

    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::createIntegratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

//   Test if at 'Final Time'
    double time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::VectorBase<double> > x_exact =
      model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution
    std::ofstream ftmp("Tempus_NewmarkExplicitAForm_BallParabolic.dat");
    ftmp.precision(16);
    RCP<const SolutionHistory<double> > solutionHistory =
      integrator->getSolutionHistory();
    bool passed = true;
    double err = 0.0;
    RCP<const Thyra::VectorBase<double> > x_exact_plot;
    for (int i=0; i<solutionHistory->getNumStates(); i++) {
      RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
      double time_i = solutionState->getTime();
      RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
      x_exact_plot = model->getExactSolution(time_i).get_x();
      ftmp << time_i << "   "
           << get_ele(*(x_plot), 0) << "   "
           << get_ele(*(x_exact_plot), 0) << std::endl;
      if (abs(get_ele(*(x_plot),0) - get_ele(*(x_exact_plot), 0)) > err)
        err = abs(get_ele(*(x_plot),0) - get_ele(*(x_exact_plot), 0));
    }
    ftmp.close();
    std::cout << "Max error = " << err << "\n \n";
    if (err > tolerance)
      passed = false;

    TEUCHOS_TEST_FOR_EXCEPTION(!passed, std::logic_error,
      "\n Test failed!  Max error = " << err << " > tolerance = " << tolerance << "\n!");

    // Check the order and intercept
    RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
    std::cout << "  Stepper = " << stepper->description()
              << "\n            with " << option << std::endl;
    std::cout << "  =========================" << std::endl;
    std::cout << "  Exact solution   : " << get_ele(*(x_exact), 0) << std::endl;
    std::cout << "  Computed solution: " << get_ele(*(x      ), 0) << std::endl;
    std::cout << "  =========================" << std::endl;
    TEST_ASSERT(std::abs(get_ele(*(x), 0)) < 1.0e-14);
  }
}


// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkExplicitAForm, SinCos)
{
  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 9;
  double time = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_NewmarkExplicitAForm_SinCos.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<HarmonicOscillatorModel<double> > model =
    Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl));


  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  RCP<ParameterList> stepperPL = sublist(pl, "Default Stepper", true);
  stepperPL->remove("Zero Initial Guess");

  //Set initial time step = 2*dt specified in input file (for convergence study)
  //
  double dt =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Initial Time Step");
  dt *= 2.0;

  for (int n=0; n<nTimeStepSizes; n++) {

    //Perform time-step refinement
    dt /= 2;
    std::cout << "\n \n time step #" << n << " (out of "
              << nTimeStepSizes-1 << "), dt = " << dt << "\n";
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    integrator = Tempus::createIntegratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Plot sample solution and exact solution
    if (n == 0) {
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      writeSolution("Tempus_NewmarkExplicitAForm_SinCos.dat", solutionHistory);

      RCP<Tempus::SolutionHistory<double> > solnHistExact =
        Teuchos::rcp(new Tempus::SolutionHistory<double>());
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        double time_i = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double> > state =
          Tempus::createSolutionStateX(
          rcp_const_cast<Thyra::VectorBase<double> > (
            model->getExactSolution(time_i).get_x()),
          rcp_const_cast<Thyra::VectorBase<double> > (
            model->getExactSolution(time_i).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution("Tempus_NewmarkExplicitAForm_SinCos-Ref.dat",solnHistExact);
    }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXDot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes-1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solutionExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),solutionExact.ptr());
      solutions.push_back(solutionExact);
      auto solutionDotExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDotExact.ptr());
      solutionsDot.push_back(solutionDotExact);
    }
  }

  // Check the order and intercept
  double xSlope = 0.0;
  double xDotSlope = 0.0;
  RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
  double order = stepper->getOrder();
  writeOrderError("Tempus_NewmarkExplicitAForm_SinCos-Error.dat",
                  stepper, StepSize,
                  solutions,    xErrorNorm,    xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope);

  TEST_FLOATING_EQUALITY( xSlope,              order, 0.02   );
  TEST_FLOATING_EQUALITY( xErrorNorm[0],   0.0157928, 1.0e-4 );
  TEST_FLOATING_EQUALITY( xDotSlope,           order, 0.02   );
  TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.233045, 1.0e-4 );

  Teuchos::TimeMonitor::summarize();
}


// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkExplicitAForm, HarmonicOscillatorDamped)
{
  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 9;
  double time = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<HarmonicOscillatorModel<double> > model =
    Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl));


  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  RCP<ParameterList> stepperPL = sublist(pl, "Default Stepper", true);
  stepperPL->remove("Zero Initial Guess");

  //Set initial time step = 2*dt specified in input file (for convergence study)
  //
  double dt =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Initial Time Step");
  dt *= 2.0;

  for (int n=0; n<nTimeStepSizes; n++) {

    //Perform time-step refinement
    dt /= 2;
    std::cout << "\n \n time step #" << n << " (out of "
              << nTimeStepSizes-1 << "), dt = " << dt << "\n";
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    integrator = Tempus::createIntegratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Plot sample solution and exact solution
    if (n == 0) {
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      writeSolution("Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped.dat", solutionHistory);

      RCP<Tempus::SolutionHistory<double> > solnHistExact =
        Teuchos::rcp(new Tempus::SolutionHistory<double>());
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        double time_i = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double> > state =
          Tempus::createSolutionStateX(
          rcp_const_cast<Thyra::VectorBase<double> > (
            model->getExactSolution(time_i).get_x()),
          rcp_const_cast<Thyra::VectorBase<double> > (
            model->getExactSolution(time_i).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution("Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped-Ref.dat", solnHistExact);
    }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXDot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes-1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solutionExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),solutionExact.ptr());
      solutions.push_back(solutionExact);
      auto solutionDotExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDotExact.ptr());
      solutionsDot.push_back(solutionDotExact);
    }
  }

  // Check the order and intercept
  double xSlope = 0.0;
  double xDotSlope = 0.0;
  RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
  //double order = stepper->getOrder();
  writeOrderError("Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped-Error.dat",
                  stepper, StepSize,
                  solutions,    xErrorNorm,    xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope);

  TEST_FLOATING_EQUALITY( xSlope,           1.060930, 0.01   );
  TEST_FLOATING_EQUALITY( xErrorNorm[0],    0.508229, 1.0e-4 );
  TEST_FLOATING_EQUALITY( xDotSlope,        1.019300, 0.01   );
  TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.172900, 1.0e-4 );

  Teuchos::TimeMonitor::summarize();
}


// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, ConstructingFromDefaults)
{
  double dt = 1.0;
  std::vector<std::string> options;
  options.push_back("Default Parameters");
  options.push_back("ICConsistency and Check");

  for(const auto& option: options) {

    // Read params from .xml file
    RCP<ParameterList> pList = getParametersFromXmlFile(
      "Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.xml");
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);

    // Setup the HarmonicOscillatorModel
    RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
    auto model = Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl));
    auto modelME = rcp_dynamic_cast<const Thyra::ModelEvaluator<double>>(model);

    // Setup Stepper for field solve ----------------------------
    RCP<Tempus::StepperNewmarkImplicitAForm<double> > stepper =
      Tempus::createStepperNewmarkImplicitAForm(modelME, Teuchos::null);
    if (option == "ICConsistency and Check") {
      stepper->setICConsistency("Consistent");
      stepper->setICConsistencyCheck(true);
    }
    stepper->initialize();

    // Setup TimeStepControl ------------------------------------
    RCP<Tempus::TimeStepControl<double> > timeStepControl =
      Teuchos::rcp(new Tempus::TimeStepControl<double>());
    ParameterList tscPL = pl->sublist("Default Integrator")
                             .sublist("Time Step Control");
    timeStepControl->setInitIndex(tscPL.get<int>   ("Initial Time Index"));
    timeStepControl->setInitTime (tscPL.get<double>("Initial Time"));
    timeStepControl->setFinalTime(tscPL.get<double>("Final Time"));
    timeStepControl->setInitTimeStep(dt);
    timeStepControl->initialize();

    // Setup initial condition SolutionState --------------------
    using Teuchos::rcp_const_cast;
    auto inArgsIC = model->getNominalValues();
    RCP<Thyra::VectorBase<double> > icX =
      rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
    RCP<Thyra::VectorBase<double> > icXDot =
      rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot());
    RCP<Thyra::VectorBase<double> > icXDotDot =
      rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot_dot());
    RCP<Tempus::SolutionState<double> > icState =
      Tempus::createSolutionStateX(icX, icXDot, icXDotDot);
    icState->setTime    (timeStepControl->getInitTime());
    icState->setIndex   (timeStepControl->getInitIndex());
    icState->setTimeStep(0.0);
    icState->setOrder   (stepper->getOrder());
    icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

    // Setup SolutionHistory ------------------------------------
    RCP<Tempus::SolutionHistory<double> > solutionHistory =
      Teuchos::rcp(new Tempus::SolutionHistory<double>());
    solutionHistory->setName("Forward States");
    solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
    solutionHistory->setStorageLimit(2);
    solutionHistory->addState(icState);

    // Ensure ICs are consistent.
    stepper->setInitialConditions(solutionHistory);

    // Setup Integrator -----------------------------------------
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::createIntegratorBasic<double>();
    integrator->setStepper(stepper);
    integrator->setTimeStepControl(timeStepControl);
    integrator->setSolutionHistory(solutionHistory);
    //integrator->setObserver(...);
    integrator->initialize();


    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)


    // Test if at 'Final Time'
    double time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::VectorBase<double> > x_exact =
      model->getExactSolution(time).get_x();

    // Calculate the error
    RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));

    // Check the order and intercept
    std::cout << "  Stepper = " << stepper->description()
              << "\n            with " << option << std::endl;
    std::cout << "  =========================" << std::endl;
    std::cout << "  Exact solution   : " << get_ele(*(x_exact), 0) << std::endl;
    std::cout << "  Computed solution: " << get_ele(*(x      ), 0) << std::endl;
    std::cout << "  Difference       : " << get_ele(*(xdiff  ), 0) << std::endl;
    std::cout << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY(get_ele(*(x), 0), -0.222222, 1.0e-4 );
  }
}


// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitDForm, Constructing_From_Defaults)
{
  double dt = 1.0;

  // Read params from .xml file
  RCP<ParameterList> pList = getParametersFromXmlFile(
    "Tempus_NewmarkImplicitDForm_HarmonicOscillator_Damped_SecondOrder.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  auto model = Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl, true));
  auto modelME = rcp_dynamic_cast<const Thyra::ModelEvaluator<double>>(model);

  // Setup Stepper for field solve ----------------------------
  auto stepper = Tempus::createStepperNewmarkImplicitDForm(modelME, Teuchos::null);

  // Setup TimeStepControl ------------------------------------
  RCP<Tempus::TimeStepControl<double> > timeStepControl =
    Teuchos::rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL = pl->sublist("Default Integrator")
                           .sublist("Time Step Control");
  timeStepControl->setInitIndex(tscPL.get<int>   ("Initial Time Index"));
  timeStepControl->setInitTime (tscPL.get<double>("Initial Time"));
  timeStepControl->setFinalTime(tscPL.get<double>("Final Time"));
  timeStepControl->setInitTimeStep(dt);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  using Teuchos::rcp_const_cast;
  auto inArgsIC = model->getNominalValues();
  RCP<Thyra::VectorBase<double> > icX =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  RCP<Thyra::VectorBase<double> > icXDot =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot());
  RCP<Thyra::VectorBase<double> > icXDotDot =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot_dot());
  RCP<Tempus::SolutionState<double> > icState =
    Tempus::createSolutionStateX(icX, icXDot, icXDotDot);
  icState->setTime    (timeStepControl->getInitTime());
  icState->setIndex   (timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);
  icState->setOrder   (stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  RCP<Tempus::SolutionHistory<double> > solutionHistory =
    Teuchos::rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Setup Integrator -----------------------------------------
  RCP<Tempus::IntegratorBasic<double> > integrator =
    Tempus::createIntegratorBasic<double>();
  integrator->setStepper(stepper);
  integrator->setTimeStepControl(timeStepControl);
  integrator->setSolutionHistory(solutionHistory);
  //integrator->setObserver(...);
  integrator->initialize();


  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)


  // Test if at 'Final Time'
  double time = integrator->getTime();
  double timeFinal =pl->sublist("Default Integrator")
     .sublist("Time Step Control").get<double>("Final Time");
  TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

  // Time-integrated solution and the exact solution
  RCP<Thyra::VectorBase<double> > x = integrator->getX();
  RCP<const Thyra::VectorBase<double> > x_exact =
    model->getExactSolution(time).get_x();

  // Calculate the error
  RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
  Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));

  // Check the order and intercept
  std::cout << "  Stepper = " << stepper->description() << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Exact solution   : " << get_ele(*(x_exact), 0) << std::endl;
  std::cout << "  Computed solution: " << get_ele(*(x      ), 0) << std::endl;
  std::cout << "  Difference       : " << get_ele(*(xdiff  ), 0) << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), -0.222222, 1.0e-4 );
}


// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, HarmonicOscillatorDamped_SecondOrder)
{
  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 10;
  double time = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<HarmonicOscillatorModel<double> > model =
    Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl));


  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  //Set initial time step = 2*dt specified in input file (for convergence study)
  //
  double dt =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Initial Time Step");
  dt *= 2.0;

  for (int n=0; n<nTimeStepSizes; n++) {

    //Perform time-step refinement
    dt /= 2;
    std::cout << "\n \n time step #" << n << " (out of "
              << nTimeStepSizes-1 << "), dt = " << dt << "\n";
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    integrator = Tempus::createIntegratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Plot sample solution and exact solution
    if (n == 0) {
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      writeSolution("Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.dat", solutionHistory);

      RCP<Tempus::SolutionHistory<double> > solnHistExact =
        Teuchos::rcp(new Tempus::SolutionHistory<double>());
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        double time_i = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double> > state =
          Tempus::createSolutionStateX(
          rcp_const_cast<Thyra::VectorBase<double> > (
            model->getExactSolution(time_i).get_x()),
          rcp_const_cast<Thyra::VectorBase<double> > (
            model->getExactSolution(time_i).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution("Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder-Ref.dat", solnHistExact);
    }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXDot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes-1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solutionExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),solutionExact.ptr());
      solutions.push_back(solutionExact);
      auto solutionDotExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDotExact.ptr());
      solutionsDot.push_back(solutionDotExact);
    }
  }

  // Check the order and intercept
  double xSlope = 0.0;
  double xDotSlope = 0.0;
  RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
  double order = stepper->getOrder();
  writeOrderError("Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder_SinCos-Error.dat",
                  stepper, StepSize,
                  solutions,    xErrorNorm,    xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope);

  TEST_FLOATING_EQUALITY( xSlope,               order, 0.01   );
  TEST_FLOATING_EQUALITY( xErrorNorm[0],    0.0484483, 1.0e-4 );
  TEST_FLOATING_EQUALITY( xDotSlope,            order, 0.01   );
  TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.0484483, 1.0e-4 );

  Teuchos::TimeMonitor::summarize();
}


// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitDForm, HarmonicOscillatorDamped_SecondOrder)
{
  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 10;
  double time = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_NewmarkImplicitDForm_HarmonicOscillator_Damped_SecondOrder.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<HarmonicOscillatorModel<double> > model =
    Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl, true));


  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  //Set initial time step = 2*dt specified in input file (for convergence study)
  //
  double dt =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Initial Time Step");
  dt *= 2.0;

  for (int n=0; n<nTimeStepSizes; n++) {

    //Perform time-step refinement
    dt /= 2;
    std::cout << "\n \n time step #" << n << " (out of "
              << nTimeStepSizes-1 << "), dt = " << dt << "\n";
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    integrator = Tempus::createIntegratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Plot sample solution and exact solution
    if (n == 0) {
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      writeSolution("Tempus_NewmarkImplicitDForm_HarmonicOscillator_Damped_SecondOrder.dat", solutionHistory);

      RCP<Tempus::SolutionHistory<double> > solnHistExact =
        Teuchos::rcp(new Tempus::SolutionHistory<double>());
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        double time_i = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double> > state =
          Tempus::createSolutionStateX(
          rcp_const_cast<Thyra::VectorBase<double> > (
            model->getExactSolution(time_i).get_x()),
          rcp_const_cast<Thyra::VectorBase<double> > (
            model->getExactSolution(time_i).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution("Tempus_NewmarkImplicitDForm_HarmonicOscillator_Damped_SecondOrder-Ref.dat", solnHistExact);
    }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXDot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes-1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solutionExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),solutionExact.ptr());
      solutions.push_back(solutionExact);
      auto solutionDotExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDotExact.ptr());
      solutionsDot.push_back(solutionDotExact);
    }
  }

  // Check the order and intercept
  double xSlope = 0.0;
  double xDotSlope = 0.0;
  RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
  double order = stepper->getOrder();
  writeOrderError("Tempus_NewmarkImplicitDForm_HarmonicOscillator_Damped_SecondOrder_SinCos-Error.dat",
                  stepper, StepSize,
                  solutions,    xErrorNorm,    xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope);

  TEST_FLOATING_EQUALITY( xSlope,               order, 0.01   );
  TEST_FLOATING_EQUALITY( xErrorNorm[0],    0.0484483, 1.0e-4 );
  TEST_FLOATING_EQUALITY( xDotSlope,            order, 0.01   );
  TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.0484483, 1.0e-4 );

  Teuchos::TimeMonitor::summarize();
}


// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, HarmonicOscillatorDamped_FirstOrder)
{
  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 10;
  double time = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_FirstOrder.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<HarmonicOscillatorModel<double> > model =
    Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl));


  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  //Set initial time step = 2*dt specified in input file (for convergence study)
  //
  double dt =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Initial Time Step");
  dt *= 2.0;

  for (int n=0; n<nTimeStepSizes; n++) {

    //Perform time-step refinement
    dt /= 2;
    std::cout << "\n \n time step #" << n << " (out of "
              << nTimeStepSizes-1 << "), dt = " << dt << "\n";
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    integrator = Tempus::createIntegratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Plot sample solution and exact solution
    if (n == 0) {
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      writeSolution("Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_FirstOrder.dat", solutionHistory);

      RCP<Tempus::SolutionHistory<double> > solnHistExact =
        Teuchos::rcp(new Tempus::SolutionHistory<double>());
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        double time_i = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double> > state =
          Tempus::createSolutionStateX(
          rcp_const_cast<Thyra::VectorBase<double> > (
            model->getExactSolution(time_i).get_x()),
          rcp_const_cast<Thyra::VectorBase<double> > (
            model->getExactSolution(time_i).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution("Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_FirstOrder-Ref.dat", solnHistExact);
    }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXDot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes-1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solutionExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),solutionExact.ptr());
      solutions.push_back(solutionExact);
      auto solutionDotExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDotExact.ptr());
      solutionsDot.push_back(solutionDotExact);
    }
  }

  // Check the order and intercept
  double xSlope = 0.0;
  double xDotSlope = 0.0;
  RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
  double order = stepper->getOrder();
  writeOrderError("Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_FirstOrder-Error.dat",
                  stepper, StepSize,
                  solutions,    xErrorNorm,    xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope);

  TEST_FLOATING_EQUALITY( xSlope,               order, 0.02   );
  TEST_FLOATING_EQUALITY( xErrorNorm[0],    0.0224726, 1.0e-4 );
  TEST_FLOATING_EQUALITY( xDotSlope,            order, 0.02   );
  TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.0122223, 1.0e-4 );

  Teuchos::TimeMonitor::summarize();
}


}
