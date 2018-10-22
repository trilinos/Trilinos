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
#include "Tempus_StepperNewmarkImplicitAForm.hpp"

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
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

//IKT, 3/22/17: comment out any of the following
//if you wish not to build/run all the test cases.
#define TEST_BALL_PARABOLIC
#define TEST_SINCOS_EXPLICIT
#define TEST_HARMONIC_OSCILLATOR_DAMPED_EXPLICIT
#define TEST_HARMONIC_OSCILLATOR_DAMPED_CTOR
#define TEST_HARMONIC_OSCILLATOR_DAMPED


#ifdef TEST_BALL_PARABOLIC
// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkExplicitAForm, BallParabolic)
{
  //Tolerance to check if test passed
  double tolerance = 1.0e-14;
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

  RCP<Tempus::IntegratorBasic<double> > integrator =
    Tempus::integratorBasic<double>(pl, model);

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
    double time = solutionState->getTime();
    RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
    x_exact_plot = model->getExactSolution(time).get_x();
    ftmp << time << "   "
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
}
#endif


#ifdef TEST_SINCOS_EXPLICIT
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
    integrator = Tempus::integratorBasic<double>(pl, model);

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
        double time = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double> > state =
          Teuchos::rcp(new Tempus::SolutionState<double>(
            model->getExactSolution(time).get_x(),
            model->getExactSolution(time).get_x_dot()));
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
    Thyra::copy(*(integrator->getXdot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes-1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solution = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),solution.ptr());
      solutions.push_back(solution);
      auto solutionDot = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDot.ptr());
      solutionsDot.push_back(solutionDot);
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
#endif

#ifdef TEST_HARMONIC_OSCILLATOR_DAMPED_EXPLICIT
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
    integrator = Tempus::integratorBasic<double>(pl, model);

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
        double time = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double> > state =
          Teuchos::rcp(new Tempus::SolutionState<double>(
            model->getExactSolution(time).get_x(),
            model->getExactSolution(time).get_x_dot()));
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
    Thyra::copy(*(integrator->getXdot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes-1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solution = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),solution.ptr());
      solutions.push_back(solution);
      auto solutionDot = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDot.ptr());
      solutionsDot.push_back(solutionDot);
    }
  }

  // Check the order and intercept
  double xSlope = 0.0;
  double xDotSlope = 0.0;
  RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
  double order = stepper->getOrder();
  writeOrderError("Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped-Error.dat",
                  stepper, StepSize,
                  solutions,    xErrorNorm,    xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope);

  TEST_FLOATING_EQUALITY( xSlope,              order, 0.01   );
  TEST_FLOATING_EQUALITY( xErrorNorm[0],    0.617129, 1.0e-4 );
  TEST_FLOATING_EQUALITY( xDotSlope,           order, 0.01   );
  TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.270671, 1.0e-4 );

  Teuchos::TimeMonitor::summarize();
}
#endif


#ifdef TEST_HARMONIC_OSCILLATOR_DAMPED_CTOR
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, ConstructingFromDefaults)
{
  double dt = 1.0;

  // Read params from .xml file
  RCP<ParameterList> pList = getParametersFromXmlFile(
    "Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<HarmonicOscillatorModel<double> > model =
    Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl));

  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperNewmarkImplicitAForm<double> > stepper =
    Teuchos::rcp(new Tempus::StepperNewmarkImplicitAForm<double>(model));

  // Setup TimeStepControl ------------------------------------
  RCP<Tempus::TimeStepControl<double> > timeStepControl =
    Teuchos::rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL = pl->sublist("Default Integrator")
                           .sublist("Time Step Control");
  timeStepControl->setStepType (tscPL.get<std::string>("Integrator Step Type"));
  timeStepControl->setInitIndex(tscPL.get<int>   ("Initial Time Index"));
  timeStepControl->setInitTime (tscPL.get<double>("Initial Time"));
  timeStepControl->setFinalTime(tscPL.get<double>("Final Time"));
  timeStepControl->setInitTimeStep(dt);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  using Teuchos::rcp_const_cast;
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
    stepper->getModel()->getNominalValues();
  RCP<Thyra::VectorBase<double> > icX =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  RCP<Thyra::VectorBase<double> > icXDot =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot());
  RCP<Thyra::VectorBase<double> > icXDotDot =
    rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x_dot_dot());
  RCP<Tempus::SolutionState<double> > icState =
      Teuchos::rcp(new Tempus::SolutionState<double>(icX, icXDot, icXDotDot));
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
    Tempus::integratorBasic<double>();
  integrator->setStepperWStepper(stepper);
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
#endif


#ifdef TEST_HARMONIC_OSCILLATOR_DAMPED
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
    integrator = Tempus::integratorBasic<double>(pl, model);

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
        double time = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double> > state =
          Teuchos::rcp(new Tempus::SolutionState<double>(
            model->getExactSolution(time).get_x(),
            model->getExactSolution(time).get_x_dot()));
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
    Thyra::copy(*(integrator->getXdot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes-1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solution = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),solution.ptr());
      solutions.push_back(solution);
      auto solutionDot = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDot.ptr());
      solutionsDot.push_back(solutionDot);
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
    integrator = Tempus::integratorBasic<double>(pl, model);

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
        double time = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double> > state =
          Teuchos::rcp(new Tempus::SolutionState<double>(
            model->getExactSolution(time).get_x(),
            model->getExactSolution(time).get_x_dot()));
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
    Thyra::copy(*(integrator->getXdot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes-1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solution = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),solution.ptr());
      solutions.push_back(solution);
      auto solutionDot = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDot.ptr());
      solutionsDot.push_back(solutionDot);
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
#endif
}
