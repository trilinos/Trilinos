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
#include "Teuchos_DefaultComm.hpp"

#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_StepperBDF2.hpp"

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
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

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
template <typename SC, typename Model, typename Comm>
void CDR_Test(const Comm& comm, const int commSize, Teuchos::FancyOStream& out,
              bool& success)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;

  // Read params from .xml file
  auto pList = getParametersFromXmlFile("Tempus_BDF2_CDR.xml");
  // Set initial time step = 2*dt specified in input file (for convergence
  // study)
  //
  auto pl = sublist(pList, "Tempus", true);
  auto dt = pl->sublist("Demo Integrator")
                .sublist("Time Step Control")
                .get<double>("Initial Time Step");
  dt *= 2.0;
  auto model_pl = sublist(pList, "CDR Model", true);

  const auto nTimeStepSizes =
      model_pl->get<int>("Number of Time Step Sizes", 5);

  for (int n = 0; n < nTimeStepSizes; n++) {
    // Create CDR Model
    const int num_elements  = model_pl->get<int>("num elements");
    const auto left_end     = model_pl->get<SC>("left end");
    const auto right_end    = model_pl->get<SC>("right end");
    const auto a_convection = model_pl->get<SC>("a (convection)");
    const auto k_source     = model_pl->get<SC>("k (source)");

    auto model = rcp(new Model(comm, num_elements, left_end, right_end,
                               a_convection, k_source));

    // Set the factory
    ::Stratimikos::DefaultLinearSolverBuilder builder;

    auto p = rcp(new ParameterList);
    p->set("Linear Solver Type", "Belos");
    p->set("Preconditioner Type", "None");
    builder.setParameterList(p);

    auto lowsFactory = builder.createLinearSolveStrategy("");

    model->set_W_factory(lowsFactory);

    // Set the step size
    dt /= 2;

    // Setup the Integrator and reset initial time step
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
    if ((n == nTimeStepSizes - 1) && (commSize == 1)) {
      std::ofstream ftmp("Tempus_BDF2_CDR.dat");
      ftmp << "TITLE=\"BDF2 Solution to CDR\"\n"
           << "VARIABLES=\"z\",\"T\"\n";
      const auto dx =
          std::fabs(left_end - right_end) / static_cast<double>(num_elements);
      auto solutionHistory = integrator->getSolutionHistory();
      int nStates          = solutionHistory->getNumStates();
      for (int i = 0; i < nStates; i++) {
        auto solutionState = (*solutionHistory)[i];
        auto x             = solutionState->getX();
        auto ttime         = solutionState->getTime();
        ftmp << "ZONE T=\"Time=" << ttime << "\", I=" << num_elements + 1
             << ", F=BLOCK\n";
        for (int j = 0; j < num_elements + 1; j++) {
          const auto x_coord = left_end + static_cast<double>(j) * dx;
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
  if (nTimeStepSizes > 2) {
    double xSlope    = 0.0;
    double xDotSlope = 0.0;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;
    auto stepper = integrator->getStepper();
    auto order   = stepper->getOrder();
    writeOrderError("Tempus_BDF2_CDR-Error.dat", stepper, StepSize, solutions,
                    xErrorNorm, xSlope, solutionsDot, xDotErrorNorm, xDotSlope,
                    out);
    TEST_FLOATING_EQUALITY(xSlope, order, 0.35);
    TEST_COMPARE(xSlope, >, 0.95);
    TEST_FLOATING_EQUALITY(xDotSlope, order, 0.35);
    TEST_COMPARE(xDotSlope, >, 0.95);

    TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.0145747, 1.0e-4);
    TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 0.0563621, 1.0e-4);
  }

  // Write fine mesh solution at final time
  // This only works for ONE MPI process
  if (commSize == 1) {
    auto pListCDR           = getParametersFromXmlFile("Tempus_BDF2_CDR.xml");
    auto model_pl_CDR       = sublist(pListCDR, "CDR Model", true);
    const auto num_elements = model_pl_CDR->get<int>("num elements");
    const auto left_end     = model_pl_CDR->get<double>("left end");
    const auto right_end    = model_pl_CDR->get<double>("right end");

    const auto& x = *(solutions[solutions.size() - 1]);

    std::ofstream ftmp("Tempus_BDF2_CDR-Solution.dat");
    for (int n = 0; n < num_elements + 1; n++) {
      const auto dx =
          std::fabs(left_end - right_end) / static_cast<double>(num_elements);
      const auto x_coord = left_end + static_cast<double>(n) * dx;
      ftmp << x_coord << "   " << Thyra::get_ele(x, n) << std::endl;
    }
    ftmp.close();
  }

  Teuchos::TimeMonitor::summarize();
}

#ifdef TEMPUS_ENABLE_EPETRA_STACK
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, CDR)
{
  // Create a communicator for Epetra objects
  RCP<Epetra_Comm> comm;
#ifdef Tempus_ENABLE_MPI
  comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = rcp(new Epetra_SerialComm);
#endif

  CDR_Test<double, Tempus_Test::CDR_Model<double>>(comm, comm->NumProc(), out,
                                                   success);
}
#endif

#ifdef TEMPUS_ENABLE_TPETRA_STACK
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, CDR_Tpetra)
{
  // Get default Tpetra template types
  using SC   = Tpetra::Vector<>::scalar_type;
  using LO   = Tpetra::Vector<>::local_ordinal_type;
  using GO   = Tpetra::Vector<>::global_ordinal_type;
  using Node = Tpetra::Vector<>::node_type;

  auto comm = Tpetra::getDefaultComm();

  CDR_Test<SC, Tempus_Test::CDR_Model_Tpetra<SC, LO, GO, Node>>(
      comm, comm->getSize(), out, success);
}
#endif

}  // namespace Tempus_Test
