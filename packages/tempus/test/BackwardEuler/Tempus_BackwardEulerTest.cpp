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
#include "Tempus_IntegratorForwardSensitivity.hpp"
#include "Tempus_IntegratorPseudoTransientForwardSensitivity.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/CDR_Model.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestModels/SteadyQuadraticModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

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


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, ParameterList)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_BackwardEuler_SinCos.xml");

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
    // Remove Predictor for comparison
    stepperPL->remove("Predictor Name");
    stepperPL->remove("Default Predictor");
    RCP<ParameterList> defaultPL =
      integrator->getStepper()->getDefaultParameters();
    TEST_ASSERT(haveSameValues(*stepperPL,*defaultPL))
  }

  // Test constructor IntegratorBasic(model, stepperType)
  {
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(model, "Backward Euler");

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
    RCP<ParameterList> defaultPL =
      integrator->getStepper()->getDefaultParameters();

    TEST_ASSERT(haveSameValues(*stepperPL,*defaultPL))
  }
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, SinCos)
{
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 7;
  double dt = 0.2;
  double order = 0.0;
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_BackwardEuler_SinCos.xml");

    //std::ofstream ftmp("PL.txt");
    //pList->print(ftmp);
    //ftmp.close();

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    //RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
    RCP<SinCosModel<double> > model =
      Teuchos::rcp(new SinCosModel<double>(scm_pl));

    dt /= 2;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(pl, model);
    order = integrator->getStepper()->getOrder();

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
    double time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::VectorBase<double> > x_exact =
      model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution
    if (n == 0) {
      std::ofstream ftmp("Tempus_BackwardEuler_SinCos.dat");
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
             << get_ele(*(x_plot), 1) << "   "
             << get_ele(*(x_exact_plot), 0) << "   "
             << get_ele(*(x_exact_plot), 1) << std::endl;
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
  std::cout << "  Stepper = BackwardEuler" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Expected order: " << order << std::endl;
  std::cout << "  Observed order: " << slope << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope, order, 0.01 );
  TEST_FLOATING_EQUALITY( ErrorNorm[0], 0.0486418, 1.0e-4 );

  std::ofstream ftmp("Tempus_BackwardEuler_SinCos-Error.dat");
  double error0 = 0.8*ErrorNorm[0];
  for (int n=0; n<nTimeStepSizes; n++) {
    ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
         << error0*(StepSize[n]/StepSize[0]) << std::endl;
  }
  ftmp.close();
}

// ************************************************************
// ************************************************************
void test_sincos_fsa(const bool use_combined_method,
                     const bool use_dfdp_as_tangent,
                     Teuchos::FancyOStream &out, bool &success)
{
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 7;
  double dt = 0.2;
  double order = 0.0;
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();
  Teuchos::RCP<Teuchos::FancyOStream> my_out =
    Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  my_out->setProcRankAndSize(comm->getRank(), comm->getSize());
  my_out->setOutputToRootOnly(0);
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_BackwardEuler_SinCos.xml");

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    scm_pl->set("Use DfDp as Tangent", use_dfdp_as_tangent);
    RCP<SinCosModel<double> > model =
      Teuchos::rcp(new SinCosModel<double>(scm_pl));

    dt /= 2;

    // Setup sensitivities
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    ParameterList& sens_pl = pl->sublist("Sensitivities");
    if (use_combined_method)
      sens_pl.set("Sensitivity Method", "Combined");
    else {
      sens_pl.set("Sensitivity Method", "Staggered");
      sens_pl.set("Reuse State Linear Solver", true);
    }
    sens_pl.set("Use DfDp as Tangent", use_dfdp_as_tangent);

    // Setup the Integrator and reset initial time step
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    RCP<Tempus::IntegratorForwardSensitivity<double> > integrator =
      Tempus::integratorForwardSensitivity<double>(pl, model);
    order = integrator->getStepper()->getOrder();

    // Initial Conditions
    double t0 = pl->sublist("Default Integrator")
      .sublist("Time Step Control").get<double>("Initial Time");
    RCP<const Thyra::VectorBase<double> > x0 =
      model->getExactSolution(t0).get_x();
    const int num_param = model->get_p_space(0)->dim();
    RCP<Thyra::MultiVectorBase<double> > DxDp0 =
      Thyra::createMembers(model->get_x_space(), num_param);
    for (int i=0; i<num_param; ++i)
      Thyra::assign(DxDp0->col(i).ptr(),
                    *(model->getExactSensSolution(i, t0).get_x()));
    integrator->setInitialState(t0, x0, Teuchos::null, Teuchos::null,
                                DxDp0, Teuchos::null, Teuchos::null);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    double time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<const Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::MultiVectorBase<double> > DxDp = integrator->getDxDp();
    RCP<const Thyra::VectorBase<double> > x_exact =
      model->getExactSolution(time).get_x();
    RCP<Thyra::MultiVectorBase<double> > DxDp_exact =
      Thyra::createMembers(model->get_x_space(), num_param);
    for (int i=0; i<num_param; ++i)
      Thyra::assign(DxDp_exact->col(i).ptr(),
                    *(model->getExactSensSolution(i, time).get_x()));

    // Plot sample solution and exact solution
    if (comm->getRank() == 0 && n == nTimeStepSizes-1) {
      typedef Thyra::DefaultMultiVectorProductVector<double> DMVPV;

      std::ofstream ftmp("Tempus_BackwardEuler_SinCos_Sens.dat");
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      RCP< Thyra::MultiVectorBase<double> > DxDp_exact_plot =
        Thyra::createMembers(model->get_x_space(), num_param);
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        double time = solutionState->getTime();
        RCP<const DMVPV> x_prod_plot =
          Teuchos::rcp_dynamic_cast<const DMVPV>(solutionState->getX());
        RCP<const Thyra::VectorBase<double> > x_plot =
          x_prod_plot->getMultiVector()->col(0);
        RCP<const Thyra::MultiVectorBase<double> > DxDp_plot =
          x_prod_plot->getMultiVector()->subView(Teuchos::Range1D(1,num_param));
        RCP<const Thyra::VectorBase<double> > x_exact_plot =
          model->getExactSolution(time).get_x();
        for (int j=0; j<num_param; ++j)
          Thyra::assign(DxDp_exact_plot->col(j).ptr(),
                        *(model->getExactSensSolution(j, time).get_x()));
        ftmp << std::fixed << std::setprecision(7)
             << time
             << std::setw(11) << get_ele(*(x_plot), 0)
             << std::setw(11) << get_ele(*(x_plot), 1);
        for (int j=0; j<num_param; ++j)
          ftmp << std::setw(11) << get_ele(*(DxDp_plot->col(j)), 0)
               << std::setw(11) << get_ele(*(DxDp_plot->col(j)), 1);
        ftmp << std::setw(11) << get_ele(*(x_exact_plot), 0)
             << std::setw(11) << get_ele(*(x_exact_plot), 1);
        for (int j=0; j<num_param; ++j)
          ftmp << std::setw(11) << get_ele(*(DxDp_exact_plot->col(j)), 0)
               << std::setw(11) << get_ele(*(DxDp_exact_plot->col(j)), 1);
        ftmp << std::endl;
      }
      ftmp.close();
    }

    // Calculate the error
    RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
    RCP<Thyra::MultiVectorBase<double> > DxDpdiff = DxDp->clone_mv();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));
    Thyra::V_VmV(DxDpdiff.ptr(), *DxDp_exact, *DxDp);
    StepSize.push_back(dt);
    double L2norm = Thyra::norm_2(*xdiff);
    L2norm *= L2norm;
    Teuchos::Array<double> L2norm_DxDp(num_param);
    Thyra::norms_2(*DxDpdiff, L2norm_DxDp());
    for (int i=0; i<num_param; ++i)
      L2norm += L2norm_DxDp[i]*L2norm_DxDp[i];
    L2norm = std::sqrt(L2norm);
    ErrorNorm.push_back(L2norm);

    *my_out << " n = " << n << " dt = " << dt << " error = " << L2norm
            << std::endl;

  }

  // Check the order and intercept
  double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
  *my_out << "  Stepper = BackwardEuler" << std::endl;
  *my_out << "  =========================" << std::endl;
  *my_out << "  Expected order: " << order << std::endl;
  *my_out << "  Observed order: " << slope << std::endl;
  *my_out << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope, order, 0.015 );
  TEST_FLOATING_EQUALITY( ErrorNorm[0], 0.163653, 1.0e-4 );

  if (comm->getRank() == 0) {
    std::ofstream ftmp("Tempus_BackwardEuler_SinCos_Sens-Error.dat");
    double error0 = 0.8*ErrorNorm[0];
    for (int n=0; n<nTimeStepSizes; n++) {
      ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
           << error0*(StepSize[n]/StepSize[0]) << std::endl;
    }
    ftmp.close();
  }

}

TEUCHOS_UNIT_TEST(BackwardEuler, SinCos_Combined_FSA)
{
  test_sincos_fsa(true, false, out, success);
}

TEUCHOS_UNIT_TEST(BackwardEuler, SinCos_Combined_FSA_Tangent)
{
  test_sincos_fsa(true, true, out, success);
}

TEUCHOS_UNIT_TEST(BackwardEuler, SinCos_Staggered_FSA)
{
  test_sincos_fsa(false, false, out, success);
}

TEUCHOS_UNIT_TEST(BackwardEuler, SinCos_Staggered_FSA_Tangent)
{
  test_sincos_fsa(false, true, out, success);
}

// ************************************************************
// ************************************************************
void test_pseudotransient_fsa(const bool use_dfdp_as_tangent,
                              Teuchos::FancyOStream &out, bool &success)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_BackwardEuler_SteadyQuadratic.xml");

  // Setup the SteadyQuadraticModel
  RCP<ParameterList> scm_pl = sublist(pList, "SteadyQuadraticModel", true);
  scm_pl->set("Use DfDp as Tangent", use_dfdp_as_tangent);
  RCP<SteadyQuadraticModel<double> > model =
    Teuchos::rcp(new SteadyQuadraticModel<double>(scm_pl));

  // Setup sensitivities
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  ParameterList& sens_pl = pl->sublist("Sensitivities");
  sens_pl.set("Use DfDp as Tangent", use_dfdp_as_tangent);
  sens_pl.set("Reuse State Linear Solver", true);
  sens_pl.set("Force W Update", true); // Have to do this because for this
  // model the solver seems to be overwriting the matrix

  // Setup the Integrator
  RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<double> > integrator =
    Tempus::integratorPseudoTransientForwardSensitivity<double>(pl, model);

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus);

  // Test if at 'Final Time'
  double time = integrator->getTime();
  double timeFinal =pl->sublist("Default Integrator")
    .sublist("Time Step Control").get<double>("Final Time");
  TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

  // Time-integrated solution and the exact solution
  RCP<const Thyra::VectorBase<double> > x_vec = integrator->getX();
  RCP<const Thyra::MultiVectorBase<double> > DxDp_vec = integrator->getDxDp();
  const double x = Thyra::get_ele(*x_vec, 0);
  const double dxdb = Thyra::get_ele(*(DxDp_vec->col(0)), 0);
  const double x_exact = model->getSteadyStateSolution();
  const double dxdb_exact = model->getSteadyStateSolutionSensitivity();

  TEST_FLOATING_EQUALITY( x,    x_exact,    1.0e-6 );
  TEST_FLOATING_EQUALITY( dxdb, dxdb_exact, 1.0e-6 );
}

TEUCHOS_UNIT_TEST(BackwardEuler, SteadyQuadratic_PseudoTransient_FSA)
{
  test_pseudotransient_fsa(false, out, success);
}

TEUCHOS_UNIT_TEST(BackwardEuler, SteadyQuadratic_PseudoTransient_FSA_Tangent)
{
  test_pseudotransient_fsa(true, out, success);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, CDR)
{
  // Create a communicator for Epetra objects
  RCP<Epetra_Comm> comm;
#ifdef Tempus_ENABLE_MPI
  comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 5;
  double dt = 0.2;
  double order = 0.0;
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_BackwardEuler_CDR.xml");

    // Create CDR Model
    RCP<ParameterList> model_pl = sublist(pList, "CDR Model", true);
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
    if ((n == nTimeStepSizes-1) && (comm->NumProc() == 1)) {
      std::ofstream ftmp("Tempus_BackwardEuler_CDR.dat");
      ftmp << "TITLE=\"Backward Euler Solution to CDR\"\n"
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

  // Calculate the error - use the most temporally refined mesh for
  // the base solution.
  auto base_solution = solutions[solutions.size()-1];
  std::vector<double> StepSizeCheck;
  for (std::size_t i=0; i < (solutions.size()-1); ++i) {
    auto tmp = solutions[i];
    Thyra::Vp_StV(tmp.ptr(), -1.0, *base_solution);
    const double L2norm = Thyra::norm_2(*tmp);
    StepSizeCheck.push_back(StepSize[i]);
    ErrorNorm.push_back(L2norm);
  }

  // Check the order and intercept
  double slope = computeLinearRegressionLogLog<double>(StepSizeCheck,ErrorNorm);
  std::cout << "  Stepper = BackwardEuler" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Expected order: " << order << std::endl;
  std::cout << "  Observed order: " << slope << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope, order, 0.35 );
  TEST_COMPARE(slope, >, 0.95);
  out << "\n\n ** Slope on Backward Euler Method = " << slope
      << "\n" << std::endl;

  // Write error data
  {
    std::ofstream ftmp("Tempus_BackwardEuler_CDR-Error.dat");
    double error0 = 0.8*ErrorNorm[0];
    for (std::size_t n = 0; n < StepSizeCheck.size(); n++) {
      ftmp << StepSizeCheck[n]  << "   " << ErrorNorm[n] << "   "
           << error0*(StepSize[n]/StepSize[0]) << std::endl;
    }
    ftmp.close();
  }

  // Write fine mesh solution at final time
  // This only works for ONE MPI process
  if (comm->NumProc() == 1) {
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_BackwardEuler_CDR.xml");
    RCP<ParameterList> model_pl = sublist(pList, "CDR Model", true);
    const int num_elements = model_pl->get<int>("num elements");
    const double left_end = model_pl->get<double>("left end");
    const double right_end = model_pl->get<double>("right end");

    const Thyra::VectorBase<double>& x = *(solutions[solutions.size()-1]);

    std::ofstream ftmp("Tempus_BackwardEuler_CDR-Solution.dat");
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

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, VanDerPol)
{
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 4;
  double dt = 0.05;
  double order = 0.0;
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_BackwardEuler_VanDerPol.xml");

    // Setup the VanDerPolModel
    RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
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
      std::string fname = "Tempus_BackwardEuler_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_BackwardEuler_VanDerPol.dat";
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

  // Check the order and intercept
  double slope = computeLinearRegressionLogLog<double>(StepSizeCheck,ErrorNorm);
  std::cout << "  Stepper = BackwardEuler" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Expected order: " << order << std::endl;
  std::cout << "  Observed order: " << slope << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope, order, 0.10 );
  out << "\n\n ** Slope on Backward Euler Method = " << slope
      << "\n" << std::endl;

  // Write error data
  {
    std::ofstream ftmp("Tempus_BackwardEuler_VanDerPol-Error.dat");
    double error0 = 0.8*ErrorNorm[0];
    for (std::size_t n = 0; n < StepSizeCheck.size(); n++) {
      ftmp << StepSizeCheck[n]  << "   " << ErrorNorm[n] << "   "
           << error0*(StepSize[n]/StepSize[0]) << std::endl;
    }
    ftmp.close();
  }

  Teuchos::TimeMonitor::summarize();
}


} // namespace Tempus_Test
