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

#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_StepperNewmarkImplicitAForm.hpp"
#include "Tempus_StepperNewmarkImplicitDForm.hpp"

#include "../TestModels/HarmonicOscillatorModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

namespace Tempus_Test {

using Teuchos::getParametersFromXmlFile;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkExplicitAForm, BallParabolic)
{
  // Tolerance to check if test passed
  double tolerance = 1.0e-14;
  std::vector<std::string> options;
  options.push_back("useFSAL=true");
  options.push_back("useFSAL=false");
  options.push_back("ICConsistency and Check");

  for (const auto& option : options) {
    // Read params from .xml file
    RCP<ParameterList> pList = getParametersFromXmlFile(
        "Tempus_Test_NewmarkExplicitAForm_BallParabolic.xml");

    // Setup the HarmonicOscillatorModel
    RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
    RCP<HarmonicOscillatorModel<double> > model =
        Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl));

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl        = sublist(pList, "Tempus", true);
    RCP<ParameterList> stepperPL = sublist(pl, "Default Stepper", true);
    stepperPL->remove("Zero Initial Guess");
    if (option == "useFSAL=true")
      stepperPL->set("Use FSAL", true);
    else if (option == "useFSAL=false")
      stepperPL->set("Use FSAL", false);
    else if (option == "ICConsistency and Check") {
      stepperPL->set("Initial Condition Consistency", "Consistent");
      stepperPL->set("Initial Condition Consistency Check", true);
    }

    RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::createIntegratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    //   Test if at 'Final Time'
    double time      = integrator->getTime();
    double timeFinal = pl->sublist("Default Integrator")
                           .sublist("Time Step Control")
                           .get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::VectorBase<double> > x_exact =
        model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution
    std::ofstream ftmp("Tempus_Test_NewmarkExplicitAForm_BallParabolic.dat");
    ftmp.precision(16);
    RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
    bool passed = true;
    double err  = 0.0;
    RCP<const Thyra::VectorBase<double> > x_exact_plot;
    for (int i = 0; i < solutionHistory->getNumStates(); i++) {
      RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
      double time_i                                   = solutionState->getTime();
      RCP<const Thyra::VectorBase<double> > x_plot    = solutionState->getX();
      x_exact_plot                                    = model->getExactSolution(time_i).get_x();
      ftmp << time_i << "   " << get_ele(*(x_plot), 0) << "   "
           << get_ele(*(x_exact_plot), 0) << std::endl;
      if (abs(get_ele(*(x_plot), 0) - get_ele(*(x_exact_plot), 0)) > err)
        err = abs(get_ele(*(x_plot), 0) - get_ele(*(x_exact_plot), 0));
    }
    ftmp.close();
    out << "Max error = " << err << "\n \n";
    if (err > tolerance) passed = false;

    TEUCHOS_TEST_FOR_EXCEPTION(
        !passed, std::logic_error,
        "\n Test failed!  Max error = " << err << " > tolerance = " << tolerance
                                        << "\n!");

    // Check the order and intercept
    RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
    out << "  Stepper = " << stepper->description() << "\n            with "
        << option << std::endl;
    out << "  =========================" << std::endl;
    out << "  Exact solution   : " << get_ele(*(x_exact), 0) << std::endl;
    out << "  Computed solution: " << get_ele(*(x), 0) << std::endl;
    out << "  =========================" << std::endl;
    TEST_ASSERT(std::abs(get_ele(*(x), 0)) < 1.0e-14);
  }
}

}  // namespace Tempus_Test
