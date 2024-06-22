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

#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

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
TEUCHOS_UNIT_TEST(BDF2, VanDerPol)
{
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;

  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_BDF2_VanDerPol.xml");
  // Set initial time step = 2*dt specified in input file (for convergence
  // study)
  //
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  double dt             = pl->sublist("Demo Integrator")
                  .sublist("Time Step Control")
                  .get<double>("Initial Time Step");
  dt *= 2.0;

  RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
  const int nTimeStepSizes   = vdpm_pl->get<int>("Number of Time Step Sizes", 3);
  // const int nTimeStepSizes = 5;
  double order = 0.0;

  for (int n = 0; n < nTimeStepSizes; n++) {
    // Setup the VanDerPolModel
    auto model = rcp(new VanDerPolModel<double>(vdpm_pl));

    // Set the step size
    dt /= 2;
    if (n == nTimeStepSizes - 1) dt /= 10.0;

    // Setup the Integrator and reset initial time step
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);
    RCP<Tempus::IntegratorBasic<double>> integrator =
        Tempus::createIntegratorBasic<double>(pl, model);
    order = integrator->getStepper()->getOrder();

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
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()), solution.ptr());
    solutions.push_back(solution);
    StepSize.push_back(dt);

    // Output finest temporal solution for plotting
    // This only works for ONE MPI process
    if ((n == 0) || (n == nTimeStepSizes - 1)) {
      std::string fname = "Tempus_BDF2_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_BDF2_VanDerPol.dat";
      std::ofstream ftmp(fname);
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      int nStates = solutionHistory->getNumStates();
      for (int i = 0; i < nStates; i++) {
        RCP<const SolutionState<double>> solutionState = (*solutionHistory)[i];
        RCP<const Thyra::VectorBase<double>> x         = solutionState->getX();
        double ttime                                   = solutionState->getTime();
        ftmp << ttime << "   " << get_ele(*x, 0) << "   " << get_ele(*x, 1)
             << std::endl;
      }
      ftmp.close();
    }
  }

  // Calculate the error - use the most temporally refined mesh for
  // the reference solution.
  auto ref_solution = solutions[solutions.size() - 1];
  std::vector<double> StepSizeCheck;
  for (std::size_t i = 0; i < (solutions.size() - 1); ++i) {
    auto tmp = solutions[i];
    Thyra::Vp_StV(tmp.ptr(), -1.0, *ref_solution);
    const double L2norm = Thyra::norm_2(*tmp);
    StepSizeCheck.push_back(StepSize[i]);
    ErrorNorm.push_back(L2norm);
  }

  if (nTimeStepSizes > 2) {
    // Check the order and intercept
    double slope =
        computeLinearRegressionLogLog<double>(StepSizeCheck, ErrorNorm);
    out << "  Stepper = BDF2" << std::endl;
    out << "  =========================" << std::endl;
    out << "  Expected order: " << order << std::endl;
    out << "  Observed order: " << slope << std::endl;
    out << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY(slope, order, 0.10);
    out << "\n\n ** Slope on BDF2 Method = " << slope << "\n"
        << std::endl;
  }

  // Write error data
  {
    std::ofstream ftmp("Tempus_BDF2_VanDerPol-Error.dat");
    double error0 = 0.8 * ErrorNorm[0];
    for (std::size_t n = 0; n < StepSizeCheck.size(); n++) {
      ftmp << StepSizeCheck[n] << "   " << ErrorNorm[n] << "   "
           << error0 * (pow(StepSize[n] / StepSize[0], order)) << std::endl;
    }
    ftmp.close();
  }

  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
