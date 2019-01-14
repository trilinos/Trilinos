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

#include "Tempus_IntegratorBasic.hpp"

#include "../TestModels/SinCosModel.hpp"

#include <vector>

namespace Tempus_Test {

using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

// Test Integrator construction from ParameterList and ModelEvaluator.
TEUCHOS_UNIT_TEST(IntegratorBasic, PL_ME_Construction)
{
  // 1) Setup the ParameterList (here we start with params from .xml file)
  RCP<ParameterList> pl = getParametersFromXmlFile("Tempus_default.xml");

  // 2) Setup the ModelEvaluator
  RCP<SinCosModel<double> > model = Teuchos::rcp(new SinCosModel<double> ());

  // 3) Setup the Integrator
  RCP<ParameterList> tempusPL = sublist(pl, "Tempus", true);
  RCP<Tempus::IntegratorBasic<double> > integrator =
    Tempus::integratorBasic<double>(tempusPL, model);

  // Test the ParameterList
  RCP<ParameterList> testPL = integrator->getTempusParameterList();
  // Write out ParameterList to rebaseline test.
  //writeParameterListToXmlFile(*testPL, "Tempus_IntegratorBasic_ref-test.xml");

  // Read params from reference .xml file
  RCP<ParameterList> referencePL =
    getParametersFromXmlFile("Tempus_IntegratorBasic_ref.xml");

  bool pass = haveSameValues(*testPL, *referencePL, true);
  if (!pass) {
    std::cout << std::endl;
    std::cout << "testPL      -------------- \n" << *testPL << std::endl;
    std::cout << "referencePL -------------- \n" << *referencePL << std::endl;
  }
  TEST_ASSERT(pass)
}


// Test integator construction, and then setParameterList, setStepper, and
// initialization.
TEUCHOS_UNIT_TEST(IntegratorBasic, Construction)
{
  // 1) Setup the Integrator
  RCP<Tempus::IntegratorBasic<double> > integrator =
    Tempus::integratorBasic<double>();

  // 2) Setup the ParameterList
  //    - Start with the default Tempus PL
  //    - Add Stepper PL
  RCP<ParameterList> tempusPL = integrator->getTempusParameterList();

  tempusPL->sublist("Default Integrator").set("Stepper Name", "Demo Stepper");
  RCP<ParameterList> stepperPL = Teuchos::parameterList();
  stepperPL->set("Stepper Type", "Forward Euler");
  tempusPL->set("Demo Stepper", *stepperPL);

  integrator->setTempusParameterList(tempusPL);

  // 3) Setup the Stepper
  RCP<SinCosModel<double> > model = Teuchos::rcp(new SinCosModel<double> ());
  integrator->setStepper(model);

  // 4) Initialize integrator
  integrator->initialize();

  // Test the ParameterList
  RCP<ParameterList> testPL = integrator->getTempusParameterList();
  // Write out ParameterList to rebaseline test.
  //writeParameterListToXmlFile(*testPL,"Tempus_IntegratorBasic_ref2-test.xml");

  // Read params from reference .xml file
  RCP<ParameterList> referencePL =
    getParametersFromXmlFile("Tempus_IntegratorBasic_ref2.xml");

  bool pass = haveSameValues(*testPL, *referencePL, true);
  if (!pass) {
    std::cout << std::endl;
    std::cout << "testPL      -------------- \n" << *testPL << std::endl;
    std::cout << "referencePL -------------- \n" << *referencePL << std::endl;
  }
  TEST_ASSERT(pass)
}

} // namespace Tempus_Test
