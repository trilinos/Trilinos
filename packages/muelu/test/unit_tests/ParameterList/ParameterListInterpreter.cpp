#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "MueLu_TestHelpers.hpp"

#include "MueLu_ParameterListInterpreter.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap; // TODO: remove

  TEUCHOS_UNIT_TEST(ParameterListInterpreter, SetParameterList)
  {
    std::string  hierarchyConfigurationFiles[] = {"Smoothed-Aggregation.xml"};
    int         nHierarchyConfigurationFiles = 1;

    for(int i=0; i< nHierarchyConfigurationFiles; i++) {
      out << "Processing file: " << hierarchyConfigurationFiles[i] << std::endl;
      ParameterListInterpreter f("ParameterList/ParameterListInterpreter/" + hierarchyConfigurationFiles[i]);

      RCP<Hierarchy> h = f.CreateHierarchy();

      //       h->GetLevel(0).Set<Operator>
      //       f.SetupHierarchy(*h);

      //TODO: check results of Iterate()
    }
  }
  
} // namespace MueLuTests
