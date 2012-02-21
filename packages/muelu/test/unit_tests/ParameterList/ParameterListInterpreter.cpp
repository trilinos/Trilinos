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

    //TODO: this test can be done at compilation time
#if !defined(HAVE_MUELU_IFPACK)  or !defined(HAVE_MUELU_AMESOS)
    if (TestHelpers::Parameters::getLib() == Xpetra::UseEpetra) {
      out << "Test skipped (dependencies not available)" << std::endl;
      return;
    }
#endif

#if !defined(HAVE_MUELU_IFPACK2) or !defined(HAVE_MUELU_AMESOS2)
    if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
      out << "Test skipped (dependencies not available)" << std::endl;
      return;
    }
#endif

    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99);

    std::string  hierarchyConfigurationFiles[] = {"Smoothed-Aggregation.xml", "Smoothed-Aggregation2.xml", "Smoothed-Aggregation3.xml"};
    int         nHierarchyConfigurationFiles = 3;

    for(int i=0; i< nHierarchyConfigurationFiles; i++) {
      out << "Processing file: " << hierarchyConfigurationFiles[i] << std::endl;
      ParameterListInterpreter mueluFactory("ParameterList/ParameterListInterpreter/" + hierarchyConfigurationFiles[i]);

      RCP<Hierarchy> H = mueluFactory.CreateHierarchy();
      H->GetLevel(0)->Set<RCP<Operator> >("A", A);
     
      mueluFactory.SetupHierarchy(*H);

      //TODO: check results of Iterate()
    }
  }
  
} // namespace MueLuTests
