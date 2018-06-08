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

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_StepperOperatorSplit.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"

#include "../TestModels/VanDerPol_IMEX_ExplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEX_ImplicitModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <vector>

namespace Tempus_Test {

using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

// Comment out any of the following tests to exclude from build/run.
#define TEST_CONSTRUCTING_FROM_DEFAULTS
#define TEST_VANDERPOL


#ifdef TEST_CONSTRUCTING_FROM_DEFAULTS
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
  RCP<VanDerPol_IMEX_ExplicitModel<double> > explicitModel =
    Teuchos::rcp(new VanDerPol_IMEX_ExplicitModel<double>(vdpmPL));

  // Setup the implicit VanDerPol ModelEvaluator (reuse vdpmPL)
  RCP<VanDerPol_IMEX_ImplicitModel<double> > implicitModel =
    Teuchos::rcp(new VanDerPol_IMEX_ImplicitModel<double>(vdpmPL));


  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperOperatorSplit<double> > stepper =
    Teuchos::rcp(new Tempus::StepperOperatorSplit<double>());

  RCP<Tempus::StepperForwardEuler<double> > subStepper1 =
    Teuchos::rcp(new Tempus::StepperForwardEuler<double>(explicitModel));
  RCP<Tempus::StepperBackwardEuler<double> > subStepper2 =
    Teuchos::rcp(new Tempus::StepperBackwardEuler<double>(implicitModel));

  stepper->addStepper(subStepper1);
  stepper->addStepper(subStepper2);
  stepper->initialize();


  // Setup TimeStepControl ------------------------------------
  RCP<Tempus::TimeStepControl<double> > timeStepControl =
    Teuchos::rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL = pl->sublist("Demo Integrator")
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
  double timeFinal =pl->sublist("Demo Integrator")
     .sublist("Time Step Control").get<double>("Final Time");
  TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

  // Time-integrated solution and the exact solution
  RCP<Thyra::VectorBase<double> > x = integrator->getX();

  // Check the order and intercept
  std::cout << "  Stepper = " << stepper->description() << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Computed solution: " << get_ele(*(x      ), 0) << "   "
                                       << get_ele(*(x      ), 1) << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), -2.223910, 1.0e-4);
  TEST_FLOATING_EQUALITY(get_ele(*(x), 1),  0.565441, 1.0e-4);
}
#endif // TEST_CONSTRUCTING_FROM_DEFAULTS


#ifdef TEST_VANDERPOL
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(OperatorSplit, VanDerPol)
{
  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 4;  // 8 for Error plot
  double dt = 0.1;
  double time = 0.0;
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_OperatorSplit_VanDerPol.xml");

    // Setup the explicit VanDerPol ModelEvaluator
    RCP<ParameterList> vdpmPL = sublist(pList, "VanDerPolModel", true);
    RCP<VanDerPol_IMEX_ExplicitModel<double> > explicitModel =
      Teuchos::rcp(new VanDerPol_IMEX_ExplicitModel<double>(vdpmPL));

    // Setup the implicit VanDerPol ModelEvaluator (reuse vdpmPL)
    RCP<VanDerPol_IMEX_ImplicitModel<double> > implicitModel =
      Teuchos::rcp(new VanDerPol_IMEX_ImplicitModel<double>(vdpmPL));

    // Setup vector of models
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<double> > > models;
    models.push_back(explicitModel);
    models.push_back(implicitModel);

    // Set the step size
    dt /= 2;
    if (n == nTimeStepSizes-1) dt /= 10.0;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    integrator = Tempus::integratorBasic<double>(pl, models);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time = integrator->getTime();
    double timeFinal =pl->sublist("Demo Integrator")
      .sublist("Time Step Control").get<double>("Final Time");
    double tol = 100.0 * std::numeric_limits<double>::epsilon();
    TEST_FLOATING_EQUALITY(time, timeFinal, tol);

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(implicitModel->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(implicitModel->get_x_space());
    Thyra::copy(*(integrator->getXdot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);

    // Output finest temporal solution for plotting
    // This only works for ONE MPI process
    if ((n == 0) or (n == nTimeStepSizes-1)) {
      std::string fname = "Tempus_OperatorSplit_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_OperatorSplit_VanDerPol.dat";
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      writeSolution(fname, solutionHistory);
    }
  }

  // Check the order and intercept
  double xSlope = 0.0;
  double xDotSlope = 0.0;
  RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
  double order = stepper->getOrder();
  writeOrderError("Tempus_OperatorSplit_VanDerPol-Error.dat",
                  stepper, StepSize,
                  solutions,    xErrorNorm,    xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope);

  TEST_FLOATING_EQUALITY( xSlope,            order, 0.05 );
  TEST_FLOATING_EQUALITY( xDotSlope,         order, 0.05 );//=order at small dt
  TEST_FLOATING_EQUALITY( xErrorNorm[0],    1.27294, 1.0e-4 );
  TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 12.7102, 1.0e-4 );

  Teuchos::TimeMonitor::summarize();
}
#endif // TEST_VANDERPOL

} // namespace Tempus_Test
