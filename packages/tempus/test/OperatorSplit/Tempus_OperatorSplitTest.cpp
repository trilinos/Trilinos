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


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(OperatorSplit, VanDerPol)
{
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 4;  // 8 for Error plot
  double dt = 0.1;
  double order = 0.0;
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
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(pl, models);
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
    auto solution = Thyra::createMember(implicitModel->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    StepSize.push_back(dt);

    // Output finest temporal solution for plotting
    // This only works for ONE MPI process
    if ((n == 0) or (n == nTimeStepSizes-1)) {
      std::string fname = "Tempus_OperatorSplit_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_OperatorSplit_VanDerPol.dat";
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
  std::cout << "  Stepper = OperatorSplit" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Expected order: " << order << std::endl;
  std::cout << "  Observed order: " << slope << std::endl;
  std::cout << "  =========================" << std::endl;
  //TEST_FLOATING_EQUALITY( slope, order, 0.05 );
  TEST_FLOATING_EQUALITY( slope, 0.956, 0.02 );

  // Write error data
  {
    std::ofstream ftmp("Tempus_OperatorSplit_VanDerPol-Error.dat");
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
