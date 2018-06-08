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
#include "Teuchos_DefaultComm.hpp"

#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_StepperBDF2.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/CDR_Model.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"

#ifdef Tempus_ENABLE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

//IKT, 11/20/17: comment out any of the following
//if you wish not to build/run all the test cases.
#define TEST_PARAMETERLIST
#define TEST_CONSTRUCTING_FROM_DEFAULTS
#define TEST_SINCOS
#define TEST_SINCOS_ADAPT
#define TEST_CDR
#define TEST_VANDERPOL

namespace Tempus_Test {

using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;


#ifdef TEST_PARAMETERLIST
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, ParameterList)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_BDF2_SinCos.xml");

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  RCP<SinCosModel<double> > model =
    Teuchos::rcp(new SinCosModel<double> (scm_pl));

  RCP<ParameterList> tempusPL  = sublist(pList, "Tempus", true);

  // Test constructor IntegratorBasic(tempusPL, model)
  {
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(tempusPL, model);

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
    // Remove Start Up Stepper for comparison
    stepperPL->remove("Start Up Stepper Name");
    stepperPL->remove("Default Start Up Stepper");
    RCP<ParameterList> defaultPL =
      integrator->getStepper()->getDefaultParameters();
    TEST_ASSERT(haveSameValues(*stepperPL, *defaultPL, true))
  }

  // Test constructor IntegratorBasic(model, stepperType)
  {
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(model, "BDF2");

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
    RCP<ParameterList> defaultPL =
      integrator->getStepper()->getDefaultParameters();

    //std::cout << std::endl;
    //std::cout << "stepperPL ----------------- \n" << *stepperPL << std::endl;
    //std::cout << "defaultPL ----------------- \n" << *defaultPL << std::endl;
    TEST_ASSERT(haveSameValues(*stepperPL, *defaultPL, true))
  }
}
#endif // TEST_PARAMETERLIST


#ifdef TEST_CONSTRUCTING_FROM_DEFAULTS
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, ConstructingFromDefaults)
{
  double dt = 0.1;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_BDF2_SinCos.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  //RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
  RCP<SinCosModel<double> > model =
    Teuchos::rcp(new SinCosModel<double>(scm_pl));

  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperBDF2<double> > stepper =
    Teuchos::rcp(new Tempus::StepperBDF2<double>(model));
  //{
  //  // Setup a linear NOX solve
  //  RCP<ParameterList> sPL = stepper->getNonconstParameterList();
  //  std::string solverName = sPL->get<std::string>("Solver Name");
  //  RCP<ParameterList> solverPL = Teuchos::sublist(sPL, solverName, true);
  //  stepper->setSolver(solverPL);
  //  stepper->initialize();
  //}

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
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
    stepper->getModel()->getNominalValues();
  RCP<Thyra::VectorBase<double> > icSolution =
    Teuchos::rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  RCP<Tempus::SolutionState<double> > icState =
      Teuchos::rcp(new Tempus::SolutionState<double>(icSolution));
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
  solutionHistory->setStorageLimit(3);
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
  std::cout << "  Stepper = BDF2" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Exact solution   : " << get_ele(*(x_exact), 0) << "   "
                                       << get_ele(*(x_exact), 1) << std::endl;
  std::cout << "  Computed solution: " << get_ele(*(x      ), 0) << "   "
                                       << get_ele(*(x      ), 1) << std::endl;
  std::cout << "  Difference       : " << get_ele(*(xdiff  ), 0) << "   "
                                       << get_ele(*(xdiff  ), 1) << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), 0.839732, 1.0e-4 );
  TEST_FLOATING_EQUALITY(get_ele(*(x), 1), 0.542663, 1.0e-4 );
}
#endif // TEST_CONSTRUCTING_FROM_DEFAULTS


#ifdef TEST_SINCOS
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, SinCos)
{
  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;

  // Read params from .xml file
  RCP<ParameterList> pList = getParametersFromXmlFile("Tempus_BDF2_SinCos.xml");
  //Set initial time step = 2*dt specified in input file (for convergence study)
  //
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  double dt = pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Initial Time Step");
  dt *= 2.0;

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  const int nTimeStepSizes = scm_pl->get<int>("Number of Time Step Sizes", 7);
  std::string output_file_string =
                    scm_pl->get<std::string>("Output File Name", "Tempus_BDF2_SinCos");
  std::string output_file_name = output_file_string + ".dat";
  std::string ref_out_file_name = output_file_string + "-Ref.dat";
  std::string err_out_file_name = output_file_string + "-Error.dat";
  double time = 0.0;
  for (int n=0; n<nTimeStepSizes; n++) {

    //std::ofstream ftmp("PL.txt");
    //pList->print(ftmp);
    //ftmp.close();

    RCP<SinCosModel<double> > model =
      Teuchos::rcp(new SinCosModel<double>(scm_pl));

    dt /= 2;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    integrator = Tempus::integratorBasic<double>(pl, model);

    // Initial Conditions
    // During the Integrator construction, the initial SolutionState
    // is set by default to model->getNominalVales().get_x().  However,
    // the application can set it also by integrator->setInitialState.
    RCP<Thyra::VectorBase<double> > x0 =
      model->getNominalValues().get_x()->clone_v();
    integrator->setInitialState(0.0, x0);

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
      writeSolution(output_file_name, solutionHistory);

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
      writeSolution(ref_out_file_name, solnHistExact);
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
  if (nTimeStepSizes > 1) {
    double xSlope = 0.0;
    double xDotSlope = 0.0;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;
    RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
    double order = stepper->getOrder();
    writeOrderError(err_out_file_name,
                    stepper, StepSize,
                    solutions,    xErrorNorm,    xSlope,
                    solutionsDot, xDotErrorNorm, xDotSlope);

    TEST_FLOATING_EQUALITY( xSlope,                 order, 0.01   );
    TEST_FLOATING_EQUALITY( xDotSlope,              order, 0.01   );
    TEST_FLOATING_EQUALITY( xErrorNorm[0],    5.13425e-05, 1.0e-4 );
    TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 5.13425e-05, 1.0e-4 );
  }

  Teuchos::TimeMonitor::summarize();
}
#endif // TEST_SINCOS


#ifdef TEST_SINCOS_ADAPT
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, SinCosAdapt)
{
  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_BDF2_SinCos_AdaptDt.xml");
  //Set initial time step = 2*dt specified in input file (for convergence study)
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  double dt = pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Initial Time Step");
  dt *= 2.0;

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  const int nTimeStepSizes = scm_pl->get<int>("Number of Time Step Sizes", 7);
  std::string output_file_string =
    scm_pl->get<std::string>("Output File Name", "Tempus_BDF2_SinCos");
  std::string output_file_name = output_file_string + ".dat";
  std::string err_out_file_name = output_file_string + "-Error.dat";
  double time = 0.0;
  for (int n=0; n<nTimeStepSizes; n++) {

    //std::ofstream ftmp("PL.txt");
    //pList->print(ftmp);
    //ftmp.close();

    RCP<SinCosModel<double> > model =
      Teuchos::rcp(new SinCosModel<double>(scm_pl));

    dt /= 2;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt/4.0);
    // Ensure time step does not get larger than the initial time step size,
    // as that would mess up the convergence rates.
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Maximum Time Step", dt);
    // Ensure time step does not get too small and therefore too many steps.
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Minimum Time Step", dt/4.0);
    // For the SinCos problem eta is directly related to dt
    pl->sublist("Default Integrator")
       .sublist("Time Step Control")
       .sublist("Time Step Control Strategy")
       .sublist("basic_vs")
       .set("Minimum Value Monitoring Function", dt*0.99);
    integrator = Tempus::integratorBasic<double>(pl, model);

    // Initial Conditions
    // During the Integrator construction, the initial SolutionState
    // is set by default to model->getNominalVales().get_x().  However,
    // the application can set it also by integrator->setInitialState.
    RCP<Thyra::VectorBase<double> > x0 =
      model->getNominalValues().get_x()->clone_v();
    integrator->setInitialState(0.0, x0);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::VectorBase<double> > x_exact =
      model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution
    if (n == 0) {
      std::ofstream ftmp(output_file_name);
      //Warning: the following assumes serial run
      FILE *gold_file = fopen("Tempus_BDF2_SinCos_AdaptDt_gold.dat", "r");
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      RCP<const Thyra::VectorBase<double> > x_exact_plot;
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        char time_gold_char[100];
        fgets(time_gold_char, 100, gold_file);
        double time_gold;
        sscanf(time_gold_char, "%lf", &time_gold);
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        double time = solutionState->getTime();
        //Throw error if time does not match time in gold file to specified tolerance
        TEST_FLOATING_EQUALITY( time, time_gold, 1.0e-5 );
        RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
        x_exact_plot = model->getExactSolution(time).get_x();
        ftmp << time << "   "
             << get_ele(*(x_plot), 0) << "   "
             << get_ele(*(x_plot), 1) << "   "
             << get_ele(*(x_exact_plot), 0) << "   "
             << get_ele(*(x_exact_plot), 1) << std::endl;
      }
      ftmp.close();
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
  if (nTimeStepSizes > 1) {
    double xSlope = 0.0;
    double xDotSlope = 0.0;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;
    RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
    //double order = stepper->getOrder();
    writeOrderError("Tempus_BDF2_SinCos-Error.dat",
                    stepper, StepSize,
                    solutions,    xErrorNorm,    xSlope,
                    solutionsDot, xDotErrorNorm, xDotSlope);

    TEST_FLOATING_EQUALITY( xSlope,               1.95089, 0.01 );
    TEST_FLOATING_EQUALITY( xDotSlope,            1.95089, 0.01 );
    TEST_FLOATING_EQUALITY( xErrorNorm[0],    0.000197325, 1.0e-4 );
    TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.000197325, 1.0e-4 );
  }

  Teuchos::TimeMonitor::summarize();
}
#endif // TEST_SINCOS_ADAPT


#ifdef TEST_CDR
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, CDR)
{
  // Create a communicator for Epetra objects
  RCP<Epetra_Comm> comm;
#ifdef Tempus_ENABLE_MPI
  comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_BDF2_CDR.xml");
  //Set initial time step = 2*dt specified in input file (for convergence study)
  //
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  double dt = pl->sublist("Demo Integrator")
       .sublist("Time Step Control").get<double>("Initial Time Step");
  dt *= 2.0;
  RCP<ParameterList> model_pl = sublist(pList, "CDR Model", true);

  const int nTimeStepSizes = model_pl->get<int>("Number of Time Step Sizes", 5);

  for (int n=0; n<nTimeStepSizes; n++) {

    // Create CDR Model
    const int num_elements = model_pl->get<int>("num elements");
    const double left_end = model_pl->get<double>("left end");
    const double right_end = model_pl->get<double>("right end");
    const double a_convection = model_pl->get<double>("a (convection)");
    const double k_source = model_pl->get<double>("k (source)");

    RCP<Tempus_Test::CDR_Model<double>> model =
      Teuchos::rcp(new Tempus_Test::CDR_Model<double>(comm,
                                                      num_elements,
                                                      left_end,
                                                      right_end,
                                                      a_convection,
                                                      k_source));

    // Set the factory
    ::Stratimikos::DefaultLinearSolverBuilder builder;

    Teuchos::RCP<Teuchos::ParameterList> p =
      Teuchos::rcp(new Teuchos::ParameterList);
    p->set("Linear Solver Type", "Belos");
    p->set("Preconditioner Type", "None");
    builder.setParameterList(p);

    Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = builder.createLinearSolveStrategy("");

    model->set_W_factory(lowsFactory);

    // Set the step size
    dt /= 2;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    integrator = Tempus::integratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    double time = integrator->getTime();
    double timeFinal =pl->sublist("Demo Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    double tol = 100.0 * std::numeric_limits<double>::epsilon();
    TEST_FLOATING_EQUALITY(time, timeFinal, tol);

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXdot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);

    // Output finest temporal solution for plotting
    // This only works for ONE MPI process
    if ((n == nTimeStepSizes-1) && (comm->NumProc() == 1)) {
      std::ofstream ftmp("Tempus_BDF2_CDR.dat");
      ftmp << "TITLE=\"BDF2 Solution to CDR\"\n"
           << "VARIABLES=\"z\",\"T\"\n";
      const double dx = std::fabs(left_end-right_end) /
                        static_cast<double>(num_elements);
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      int nStates = solutionHistory->getNumStates();
      for (int i=0; i<nStates; i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        RCP<const Thyra::VectorBase<double> > x = solutionState->getX();
        double ttime = solutionState->getTime();
        ftmp << "ZONE T=\"Time="<<ttime<<"\", I="
             <<num_elements+1<<", F=BLOCK\n";
        for (int j = 0; j < num_elements+1; j++) {
          const double x_coord = left_end + static_cast<double>(j) * dx;
          ftmp << x_coord << "   ";
        }
        ftmp << std::endl;
        for (int j=0; j<num_elements+1; j++) ftmp << get_ele(*x, j) << "   ";
        ftmp << std::endl;
      }
      ftmp.close();
    }
  }

  // Check the order and intercept
  if (nTimeStepSizes > 2) {
    double xSlope = 0.0;
    double xDotSlope = 0.0;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;
    RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
    double order = stepper->getOrder();
    writeOrderError("Tempus_BDF2_CDR-Error.dat",
                    stepper, StepSize,
                    solutions,    xErrorNorm,    xSlope,
                    solutionsDot, xDotErrorNorm, xDotSlope);
    TEST_FLOATING_EQUALITY( xSlope, order, 0.35 );
    TEST_COMPARE(xSlope, >, 0.95);
    TEST_FLOATING_EQUALITY( xDotSlope, order, 0.35 );
    TEST_COMPARE(xDotSlope, >, 0.95);

    TEST_FLOATING_EQUALITY( xErrorNorm[0],    0.0145747, 1.0e-4 );
    TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.0563621, 1.0e-4 );
  }

  // Write fine mesh solution at final time
  // This only works for ONE MPI process
  if (comm->NumProc() == 1) {
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_BDF2_CDR.xml");
    RCP<ParameterList> model_pl = sublist(pList, "CDR Model", true);
    const int num_elements = model_pl->get<int>("num elements");
    const double left_end = model_pl->get<double>("left end");
    const double right_end = model_pl->get<double>("right end");

    const Thyra::VectorBase<double>& x = *(solutions[solutions.size()-1]);

    std::ofstream ftmp("Tempus_BDF2_CDR-Solution.dat");
    for (int n = 0; n < num_elements+1; n++) {
      const double dx = std::fabs(left_end-right_end) /
                        static_cast<double>(num_elements);
      const double x_coord = left_end + static_cast<double>(n) * dx;
      ftmp << x_coord << "   " <<  Thyra::get_ele(x,n) << std::endl;
    }
    ftmp.close();
  }

  Teuchos::TimeMonitor::summarize();
}
#endif // TEST_CDR


#ifdef TEST_VANDERPOL
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, VanDerPol)
{
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_BDF2_VanDerPol.xml");
  //Set initial time step = 2*dt specified in input file (for convergence study)
  //
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  double dt = pl->sublist("Demo Integrator")
       .sublist("Time Step Control").get<double>("Initial Time Step");
  dt *= 2.0;

  RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
  const int nTimeStepSizes = vdpm_pl->get<int>("Number of Time Step Sizes", 3);
  //const int nTimeStepSizes = 5;
  double order = 0.0;

  for (int n=0; n<nTimeStepSizes; n++) {

    // Setup the VanDerPolModel
    RCP<VanDerPolModel<double> > model =
      Teuchos::rcp(new VanDerPolModel<double>(vdpm_pl));

    // Set the step size
    dt /= 2;
    if (n == nTimeStepSizes-1) dt /= 10.0;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(pl, model);
    order = integrator->getStepper()->getOrder();

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    double time = integrator->getTime();
    double timeFinal =pl->sublist("Demo Integrator")
      .sublist("Time Step Control").get<double>("Final Time");
    double tol = 100.0 * std::numeric_limits<double>::epsilon();
    TEST_FLOATING_EQUALITY(time, timeFinal, tol);

    // Store off the final solution and step size
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    StepSize.push_back(dt);

    // Output finest temporal solution for plotting
    // This only works for ONE MPI process
    if ((n == 0) or (n == nTimeStepSizes-1)) {
      std::string fname = "Tempus_BDF2_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_BDF2_VanDerPol.dat";
      std::ofstream ftmp(fname);
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      int nStates = solutionHistory->getNumStates();
      for (int i=0; i<nStates; i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        RCP<const Thyra::VectorBase<double> > x = solutionState->getX();
        double ttime = solutionState->getTime();
        ftmp << ttime << "   " << get_ele(*x, 0) << "   " << get_ele(*x, 1)
             << std::endl;
      }
      ftmp.close();
    }
  }

  // Calculate the error - use the most temporally refined mesh for
  // the reference solution.
  auto ref_solution = solutions[solutions.size()-1];
  std::vector<double> StepSizeCheck;
  for (std::size_t i=0; i < (solutions.size()-1); ++i) {
    auto tmp = solutions[i];
    Thyra::Vp_StV(tmp.ptr(), -1.0, *ref_solution);
    const double L2norm = Thyra::norm_2(*tmp);
    StepSizeCheck.push_back(StepSize[i]);
    ErrorNorm.push_back(L2norm);
  }

  if (nTimeStepSizes > 2) {
    // Check the order and intercept
    double slope = computeLinearRegressionLogLog<double>(StepSizeCheck,ErrorNorm);
    std::cout << "  Stepper = BDF2" << std::endl;
    std::cout << "  =========================" << std::endl;
    std::cout << "  Expected order: " << order << std::endl;
    std::cout << "  Observed order: " << slope << std::endl;
    std::cout << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY( slope, order, 0.10 );
    out << "\n\n ** Slope on BDF2 Method = " << slope
        << "\n" << std::endl;
  }

  // Write error data
  {
    std::ofstream ftmp("Tempus_BDF2_VanDerPol-Error.dat");
    double error0 = 0.8*ErrorNorm[0];
    for (std::size_t n = 0; n < StepSizeCheck.size(); n++) {
      ftmp << StepSizeCheck[n]  << "   " << ErrorNorm[n] << "   "
           << error0*(pow(StepSize[n]/StepSize[0],order)) << std::endl;
    }
    ftmp.close();
  }

  Teuchos::TimeMonitor::summarize();
}
#endif // TEST_VANDERPOL

} // namespace Tempus_Test
