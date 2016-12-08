#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Tempus_IntegratorBasic.hpp"

#include "../TestModels/SinCosModel.hpp"

#include <vector>

using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(IntegratorBasic, DefaultConstruction)
{
  // Read params from .xml file
  RCP<ParameterList> tempusPL = getParametersFromXmlFile("Tempus_default.xml");

  // Setup the SinCosModel
  RCP<SinCosModel<double> > model = Teuchos::rcp(new SinCosModel<double> ());

  // Setup the Integrator
  RCP<ParameterList> integratorPL = sublist(tempusPL, "Tempus", true);
  RCP<Tempus::IntegratorBasic<double> > integrator =
    Tempus::integratorBasic<double>(integratorPL, model);

  RCP<const ParameterList> defaultIntegratorPL = integrator->getParameterList();
  // Write out ParameterList to rebaseline test.
  //writeParameterListToXmlFile(*defaultIntegratorPL,
  //                            "Tempus_IntegratorBasic_ref.xml");

  // Read params from reference .xml file
  RCP<ParameterList> refIntegratorPL =
    getParametersFromXmlFile("Tempus_IntegratorBasic_ref.xml");

  TEST_ASSERT(*defaultIntegratorPL == *refIntegratorPL)

  Teuchos::TimeMonitor::summarize();
}

} // namespace Tempus_Test
