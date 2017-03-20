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

#include "../TestModels/BallParabolicModel.hpp"
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



// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkImplicit, BallParabolic)
{
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 7;
  double order = 0.0;
  double dt = 0.4;
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_Newmark_BallParabolic.xml");

    //std::ofstream ftmp("PL.txt");
    //pList->print(ftmp);
    //ftmp.close();

    // Setup the BallParabolicModel
    RCP<ParameterList> scm_pl = sublist(pList, "BallParabolicModel", true);
    //RCP<BallParabolicModel<double> > model = sineCosineModel(scm_pl);
    RCP<BallParabolicModel<double> > model =
      Teuchos::rcp(new BallParabolicModel<double>(scm_pl));


    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    dt /= 2;
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

    // Plot sample solution and exact solution
    if (n == 0) {
      std::ofstream ftmp("Tempus_Newmark_BallParabolic.dat");
      RCP<SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      RCP<const Thyra::VectorBase<double> > x_exact_plot;
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        RCP<SolutionState<double> > solutionState = (*solutionHistory)[i];
        double time = solutionState->getTime();
        RCP<Thyra::VectorBase<double> > x_plot = solutionState->getX();
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
  std::cout << "  Stepper = Newmark" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Expected order: " << order << std::endl;
  std::cout << "  Observed order: " << slope << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope, order, 0.01 );
  TEST_FLOATING_EQUALITY( ErrorNorm[0], 0.0486418, 1.0e-4 );

  std::ofstream ftmp("Tempus_Newmark-Error_BallParabolic.dat");
  double error0 = 0.8*ErrorNorm[0];
  for (int n=0; n<nTimeStepSizes; n++) {
    ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
         << error0*(StepSize[n]/StepSize[0]) << std::endl;
  }
  ftmp.close();
}
}
