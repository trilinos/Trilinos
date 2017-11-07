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

#include <vector>
#include <sstream>
#include <limits>

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
#define TEST_HARMONIC_OSCILLATOR_DAMPED
#define TEST_HARMONIC_OSCILLATOR_DAMPED_EXPLICIT
#define TEST_SINCOS_EXPLICIT
#define TEST_BALL_PARABOLIC


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
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 9;
  double order = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_NewmarkExplicitAForm_SinCos.xml");

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
    std::cout << "\n \n time step #" << n << " (out of " << nTimeStepSizes-1 << "), dt = " << dt << "\n"; 
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(pl, model);
    order = integrator->getStepper()->getOrder();

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

    // Plot sample solution and exact solution at most-refined resolution 
    if (n == nTimeStepSizes-1) {
      std::ofstream ftmp("Tempus_NewmarkExplicitAForm_SinCos.dat");
      ftmp.precision(16); 
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      RCP<const Thyra::VectorBase<double> > x_exact_plot;
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        double time = solutionState->getTime();
        RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
        x_exact_plot = model->getExactSolution(time).get_x();
        ftmp << time << "   "
             << get_ele(*(x_plot), 0) << "   "
             << get_ele(*(x_exact_plot), 0) << std::endl;
      }
      ftmp.close();
    }

    // Calculate the error
    RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));
    StepSize.push_back(dt);
    const double L2norm = Thyra::norm_2(*xdiff);
    ErrorNorm.push_back(L2norm);
  }

  // Check the order and intercept
  double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
  std::cout << "  Stepper = Newmark Explicit a-Form" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Expected order: " << order << std::endl;
  std::cout << "  Observed order: " << slope << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope, order, 0.01 );
  //TEST_FLOATING_EQUALITY( ErrorNorm[0], 0.0486418, 1.0e-4 );

  std::ofstream ftmp("Tempus_NewmarkExplicitAForm-Error_SinCos.dat");
  ftmp.precision(16); 
  double error0 = 0.8*ErrorNorm[0];
  for (int n=0; n<nTimeStepSizes; n++) {
    ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
         << error0*(StepSize[n]/StepSize[0]) << std::endl;
  }
  ftmp.close();
}
#endif

#ifdef TEST_HARMONIC_OSCILLATOR_DAMPED_EXPLICIT
// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkExplicitAForm, HarmonicOscillatorDamped)
{
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 9;
  double order = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped.xml");

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
    std::cout << "\n \n time step #" << n << " (out of " << nTimeStepSizes-1 << "), dt = " << dt << "\n";
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(pl, model);
    order = integrator->getStepper()->getOrder();

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

    // Plot sample solution and exact solution at most-refined resolution 
    if (n == nTimeStepSizes-1) {
      std::ofstream ftmp("Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped.dat");
      ftmp.precision(16);
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      RCP<const Thyra::VectorBase<double> > x_exact_plot;
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        double time = solutionState->getTime();
        RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
        x_exact_plot = model->getExactSolution(time).get_x();
        ftmp << time << "   "
             << get_ele(*(x_plot), 0) << "   "
             << get_ele(*(x_exact_plot), 0) << std::endl;
      }
      ftmp.close();
    }

    // Calculate the error
    RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));
    StepSize.push_back(dt);
    const double L2norm = Thyra::norm_2(*xdiff);
    ErrorNorm.push_back(L2norm);
  }

  // Check the order and intercept
  double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
  std::cout << "  Stepper = Newmark Implicit a-Form" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Expected order: " << order << std::endl;
  std::cout << "  Observed order: " << slope << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope, order, 0.01 );
  //TEST_FLOATING_EQUALITY( ErrorNorm[0], 0.0486418, 1.0e-4 );

  std::ofstream ftmp("Tempus_NewmarkExplicitAForm-Error_HarmonicOscillator_Damped.dat");
  ftmp.precision(16);
  double error0 = 0.8*ErrorNorm[0];
  for (int n=0; n<nTimeStepSizes; n++) {
    ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
         << error0*(StepSize[n]/StepSize[0]) << std::endl;
  }
  ftmp.close();
}
#endif


#ifdef TEST_HARMONIC_OSCILLATOR_DAMPED
// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, HarmonicOscillatorDamped_SecondOrder)
{
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 10;
  double order = 0.0;

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
    std::cout << "\n \n time step #" << n << " (out of " << nTimeStepSizes-1 << "), dt = " << dt << "\n"; 
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(pl, model);
    order = integrator->getStepper()->getOrder();

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

    // Plot sample solution and exact solution at most-refined resolution 
    if (n == nTimeStepSizes-1) {
      std::ofstream ftmp("Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.dat");
      ftmp.precision(16); 
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      RCP<const Thyra::VectorBase<double> > x_exact_plot;
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        double time = solutionState->getTime();
        RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
        x_exact_plot = model->getExactSolution(time).get_x();
        ftmp << time << "   "
             << get_ele(*(x_plot), 0) << "   "
             << get_ele(*(x_exact_plot), 0) << std::endl;
      }
      ftmp.close();
    }

    // Calculate the error
    RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));
    StepSize.push_back(dt);
    const double L2norm = Thyra::norm_2(*xdiff);
    ErrorNorm.push_back(L2norm);
  }

  // Check the order and intercept
  double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
  std::cout << "  Stepper = Newmark Implicit a-Form" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Expected order: " << order << std::endl;
  std::cout << "  Observed order: " << slope << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope, order, 0.01 );
  //TEST_FLOATING_EQUALITY( ErrorNorm[0], 0.0486418, 1.0e-4 );

  std::ofstream ftmp("Tempus_Newmark-Error_HarmonicOscillator_Damped_SecondOrder.dat");
  ftmp.precision(16); 
  double error0 = 0.8*ErrorNorm[0];
  for (int n=0; n<nTimeStepSizes; n++) {
    ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
         << error0*(StepSize[n]/StepSize[0]) << std::endl;
  }
  ftmp.close();
}

// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicitAForm, HarmonicOscillatorDamped_FirstOrder)
{
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 10;
  double order = 0.0;

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
    std::cout << "\n \n time step #" << n << " (out of " << nTimeStepSizes-1 << "), dt = " << dt << "\n"; 
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(pl, model);
    order = integrator->getStepper()->getOrder();

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

    // Plot sample solution and exact solution at most-refined resolution 
    if (n == nTimeStepSizes-1) {
      std::ofstream ftmp("Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_FirstOrder.dat");
      ftmp.precision(16); 
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      RCP<const Thyra::VectorBase<double> > x_exact_plot;
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        double time = solutionState->getTime();
        RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
        x_exact_plot = model->getExactSolution(time).get_x();
        ftmp << time << "   "
             << get_ele(*(x_plot), 0) << "   "
             << get_ele(*(x_exact_plot), 0) << std::endl;
      }
      ftmp.close();
    }

    // Calculate the error
    RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));
    StepSize.push_back(dt);
    const double L2norm = Thyra::norm_2(*xdiff);
    ErrorNorm.push_back(L2norm);
  }

  // Check the order and intercept
  double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
  std::cout << "  Stepper = Newmark Implicit a-Form" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Expected order: " << order << std::endl;
  std::cout << "  Observed order: " << slope << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope, order, 0.02 );
  //TEST_FLOATING_EQUALITY( ErrorNorm[0], 0.0486418, 1.0e-4 );

  std::ofstream ftmp("Tempus_Newmark-Error_HarmonicOscillator_Damped_FirstOrder.dat");
  ftmp.precision(16); 
  double error0 = 0.8*ErrorNorm[0];
  for (int n=0; n<nTimeStepSizes; n++) {
    ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
         << error0*(StepSize[n]/StepSize[0]) << std::endl;
  }
  ftmp.close();
}
#endif
}
