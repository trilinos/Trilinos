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

#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"
#include "Thyra_VectorStdOps_decl.hpp"

#ifdef TEMPUS_ENABLE_EPETRA_STACK
#include "../TestModels/CDR_Model.hpp"
#ifdef Tempus_ENABLE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#endif
#ifdef TEMPUS_ENABLE_TPETRA_STACK
#include "../TestModels/CDR_Model_Tpetra.hpp"
#include "Tpetra_Core.hpp"
#endif
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
template <typename SC, typename Model, typename Comm>
void CDR_Test(const Comm& comm, const int commSize, Teuchos::FancyOStream& out,
              bool& success, const std::string& caseName, const bool lumped)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 4;
  double dt                = 0.0005;
  if (caseName == "Taylor") dt /= 2;
  for (int n = 0; n < nTimeStepSizes; n++) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_EPI_CDR.xml");

    // Create CDR Model
    RCP<ParameterList> model_pl = sublist(pList, "CDR Model", true);
    const auto num_elements     = model_pl->get<int>("num elements");
    const auto left_end         = model_pl->get<SC>("left end");
    const auto right_end        = model_pl->get<SC>("right end");
    const auto a_convection     = model_pl->get<SC>("a (convection)");
    const auto k_source         = model_pl->get<SC>("k (source)");

    auto model = rcp(new Model(comm, num_elements, left_end, right_end,
                               a_convection, k_source));

    // Set the factory
    ::Stratimikos::DefaultLinearSolverBuilder builder;

    auto p = rcp(new ParameterList);
    p->set("Linear Solver Type", "Belos");
    p->set("Preconditioner Type", "None");
    p->sublist("Linear Solver Types").sublist("Belos")
        .set("Solver Type", "Pseudo Block GMRES");
    p->sublist("Linear Solver Types").sublist("Belos")
        .sublist("Solver Types").sublist("Pseudo Block GMRES")
        .set("Convergence Tolerance", 1e-13);
    builder.setParameterList(p);

    auto lowsFactory = builder.createLinearSolveStrategy("");

    model->set_W_factory(lowsFactory);

    dt /= 2;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);

    auto& phiList = pl->sublist("Demo Stepper").sublist("PhiEvaluator");
    if (caseName == "Leja") {
      phiList.set("PhiEvaluator Type", "Leja")
             .set("Expansion Order", 100)
             .set("Leja DD Method", 2)
             .set("leja_tol", 1.e-12)
             .set("leja_a", -50000.0)
             .set("leja_c", 1000.0);
    }
    else if (caseName == "Taylor") {
      phiList.remove("Leja DD Method", false);
      phiList.remove("leja_tol", false);
      phiList.remove("leja_a", false);
      phiList.remove("leja_c", false);

      phiList.set("PhiEvaluator Type", "Taylor")
             .set("Expansion Order", 100);
    }

    phiList.set("Lump Mass Matrix", lumped);

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
    if ((n == nTimeStepSizes - 1) && (commSize == 1)) {
      std::ofstream ftmp("Tempus_EPI_CDR.dat");
      ftmp << "TITLE=\"Exponential Euler Solution to CDR\"\n"
           << "VARIABLES=\"z\",\"T\"\n";
      const double dx =
          std::fabs(left_end - right_end) / static_cast<double>(num_elements);
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      int nStates = solutionHistory->getNumStates();
      for (int i = 0; i < nStates; i++) {
        RCP<const SolutionState<double>> solutionState = (*solutionHistory)[i];
        RCP<const Thyra::VectorBase<double>> x         = solutionState->getX();
        double ttime                                   = solutionState->getTime();
        ftmp << "ZONE T=\"Time=" << ttime << "\", I=" << num_elements + 1
             << ", F=BLOCK\n";
        for (int j = 0; j < num_elements + 1; j++) {
          const double x_coord = left_end + static_cast<double>(j) * dx;
          ftmp << x_coord << "   ";
        }
        ftmp << std::endl;
        for (int j = 0; j < num_elements + 1; j++)
          ftmp << get_ele(*x, j) << "   ";
        ftmp << std::endl;
      }
      ftmp.close();
    }
  }

  // Check the order and intercept
  double xSlope                        = 0.0;
  double xDotSlope                     = 0.0;
  RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
  writeOrderError("Tempus_EPI_CDR_" + caseName + "-Error.dat", stepper, StepSize,
                  solutions, xErrorNorm, xSlope, out);

  TEST_COMPARE(std::abs(xErrorNorm[0]), <=, 1.0e-10);
  TEST_COMPARE(std::abs(xErrorNorm[nTimeStepSizes - 2]), <=, 1.e-12);

  // ---------------------------------------------------------------
  // Run SDIRK 3 Stage 4th order at the finest EPI step as a baseline
  // and verify that the fine EPI solution is close to it.
  // ---------------------------------------------------------------
  RCP<Thyra::VectorBase<double>> xSDIRK;
  {
    const double dt_sdirk = StepSize[nTimeStepSizes - 1] / 4.;  // finest EPI dt / 4

    RCP<ParameterList> pListS =
        getParametersFromXmlFile("Tempus_EPI_CDR.xml");

    RCP<ParameterList> model_plS = sublist(pListS, "CDR Model", true);
    const auto num_elemS  = model_plS->get<int>("num elements");
    const auto left_endS  = model_plS->get<SC>("left end");
    const auto right_endS = model_plS->get<SC>("right end");
    const auto a_convS    = model_plS->get<SC>("a (convection)");
    const auto k_srcS     = model_plS->get<SC>("k (source)");

    auto modelSDIRK = rcp(new Model(comm, num_elemS, left_endS, right_endS,
                                    a_convS, k_srcS));

    ::Stratimikos::DefaultLinearSolverBuilder builderS;
    auto pS = rcp(new ParameterList);
    pS->set("Linear Solver Type", "Belos");
    pS->set("Preconditioner Type", "None");
    builderS.setParameterList(pS);
    modelSDIRK->set_W_factory(builderS.createLinearSolveStrategy(""));

    RCP<ParameterList> plS = sublist(pListS, "Tempus", true);
    plS->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt_sdirk);
    plS->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Maximum Time Step", dt_sdirk);

    plS->sublist("Demo Stepper").set("Stepper Type", "SDIRK 3 Stage 4th order");

    plS->sublist("Demo Stepper").remove("EPI Order", false);
    plS->sublist("Demo Stepper").remove("PhiEvaluator", false);

    auto integratorSDIRK =
        Tempus::createIntegratorBasic<double>(plS, modelSDIRK);
    bool sdirkStatus = integratorSDIRK->advanceTime();
    TEST_ASSERT(sdirkStatus);

    double timeS      = integratorSDIRK->getTime();
    double timeFinalS = plS->sublist("Demo Integrator")
                            .sublist("Time Step Control")
                            .get<double>("Final Time");

    double tol = 100.0 * std::numeric_limits<double>::epsilon();
    TEST_FLOATING_EQUALITY(timeS, timeFinalS, tol);

    xSDIRK = Thyra::createMember(modelSDIRK->get_x_space());
    Thyra::copy(*(integratorSDIRK->getX()), xSDIRK.ptr());

    out << "  SDIRK 3 Stage 4th order reference computed at dt = "
        << dt_sdirk << std::endl;
  }

  // Compare finest EPI solution to SDIRK 4th-order reference
  {
    auto xDiff = Thyra::createMember(xSDIRK->space());
    Thyra::V_StVpStV(xDiff.ptr(), 1.0, *xSDIRK,
                     -1.0, *(solutions[nTimeStepSizes - 1]));
    double sdirkDiffNorm = Thyra::norm_2(*xDiff);
    double epiNorm = Thyra::norm_2(*(solutions[nTimeStepSizes - 1]));
    double sdirkNorm = Thyra::norm_2(*xSDIRK);
    out << "  ||EPI_fine||_2 = " << epiNorm << std::endl;
    out << "  ||SDIRK_fine||_2 = " << sdirkNorm << std::endl;
    out << "  ||EPI_fine - SDIRK_fine||_2 = " << sdirkDiffNorm << std::endl;
    TEST_COMPARE(sdirkDiffNorm, <=, 1.0e-2);
  }

  // Write fine mesh solution at final time
  // This only works for ONE MPI process
  if (commSize == 1) {
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_EPI_CDR.xml");
    RCP<ParameterList> model_pl = sublist(pList, "CDR Model", true);
    const int num_elements      = model_pl->get<int>("num elements");
    const double left_end       = model_pl->get<double>("left end");
    const double right_end      = model_pl->get<double>("right end");

    const Thyra::VectorBase<double>& x = *(solutions[solutions.size() - 1]);

    std::ofstream ftmp("Tempus_EPI_CDR_" + caseName + "-Solution.dat");
    for (int n = 0; n < num_elements + 1; n++) {
      const double dx =
          std::fabs(left_end - right_end) / static_cast<double>(num_elements);
      const double x_coord = left_end + static_cast<double>(n) * dx;
      ftmp << x_coord << "   " << Thyra::get_ele(x, n) << "   " << Thyra::get_ele(*xSDIRK, n) << std::endl;
    }
    ftmp.close();
  }

  Teuchos::TimeMonitor::summarize();
}

#ifdef TEMPUS_ENABLE_EPETRA_STACK
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, CDR_Leja_lumped)
{
  // Create a communicator for Epetra objects
  RCP<Epetra_Comm> comm;
#ifdef Tempus_ENABLE_MPI
  comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = rcp(new Epetra_SerialComm);
#endif

  CDR_Test<double, Tempus_Test::CDR_Model<double>>(comm, comm->NumProc(), out,
                                                   success, "Leja", true);

  std::cout << "Running EPETRA" << std::endl;
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, CDR_Leja)
{
  // Create a communicator for Epetra objects
  RCP<Epetra_Comm> comm;
#ifdef Tempus_ENABLE_MPI
  comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = rcp(new Epetra_SerialComm);
#endif

  CDR_Test<double, Tempus_Test::CDR_Model<double>>(comm, comm->NumProc(), out,
                                                   success, "Leja", false);

  std::cout << "Running EPETRA" << std::endl;
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, CDR_Taylor_lumped)
{
  // Create a communicator for Epetra objects
  RCP<Epetra_Comm> comm;
#ifdef Tempus_ENABLE_MPI
  comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = rcp(new Epetra_SerialComm);
#endif

  CDR_Test<double, Tempus_Test::CDR_Model<double>>(comm, comm->NumProc(), out,
                                                   success, "Taylor", true);

  std::cout << "Running EPETRA" << std::endl;
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, CDR_Taylor)
{
  // Create a communicator for Epetra objects
  RCP<Epetra_Comm> comm;
#ifdef Tempus_ENABLE_MPI
  comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = rcp(new Epetra_SerialComm);
#endif

  CDR_Test<double, Tempus_Test::CDR_Model<double>>(comm, comm->NumProc(), out,
                                                   success, "Taylor", false);

  std::cout << "Running EPETRA" << std::endl;
}
#endif

#ifdef TEMPUS_ENABLE_TPETRA_STACK
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, CDR_Tpetra_Leja_lumped)
{
  // Get default Tpetra template types
  using SC   = Tpetra::Vector<>::scalar_type;
  using LO   = Tpetra::Vector<>::local_ordinal_type;
  using GO   = Tpetra::Vector<>::global_ordinal_type;
  using Node = Tpetra::Vector<>::node_type;

  auto comm = Tpetra::getDefaultComm();

  CDR_Test<SC, Tempus_Test::CDR_Model_Tpetra<SC, LO, GO, Node>>(
     comm, comm->getSize(), out, success, "Leja", true);

  std::cout << "Running TPETRA" << std::endl;
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, CDR_Tpetra_Leja)
{
  // Get default Tpetra template types
  using SC   = Tpetra::Vector<>::scalar_type;
  using LO   = Tpetra::Vector<>::local_ordinal_type;
  using GO   = Tpetra::Vector<>::global_ordinal_type;
  using Node = Tpetra::Vector<>::node_type;

  auto comm = Tpetra::getDefaultComm();

  CDR_Test<SC, Tempus_Test::CDR_Model_Tpetra<SC, LO, GO, Node>>(
     comm, comm->getSize(), out, success, "Leja", false);

  std::cout << "Running TPETRA" << std::endl;
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, CDR_Tpetra_Taylor_lumped)
{
  // Get default Tpetra template types
  using SC   = Tpetra::Vector<>::scalar_type;
  using LO   = Tpetra::Vector<>::local_ordinal_type;
  using GO   = Tpetra::Vector<>::global_ordinal_type;
  using Node = Tpetra::Vector<>::node_type;

  auto comm = Tpetra::getDefaultComm();

  CDR_Test<SC, Tempus_Test::CDR_Model_Tpetra<SC, LO, GO, Node>>(
     comm, comm->getSize(), out, success, "Taylor", true);

  std::cout << "Running TPETRA" << std::endl;
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, CDR_Tpetra_Taylor)
{
  // Get default Tpetra template types
  using SC   = Tpetra::Vector<>::scalar_type;
  using LO   = Tpetra::Vector<>::local_ordinal_type;
  using GO   = Tpetra::Vector<>::global_ordinal_type;
  using Node = Tpetra::Vector<>::node_type;

  auto comm = Tpetra::getDefaultComm();

  CDR_Test<SC, Tempus_Test::CDR_Model_Tpetra<SC, LO, GO, Node>>(
     comm, comm->getSize(), out, success, "Taylor", false);

  std::cout << "Running TPETRA" << std::endl;
}
#endif



}  // namespace Tempus_Test
