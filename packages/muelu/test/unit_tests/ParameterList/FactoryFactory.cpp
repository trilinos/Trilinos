#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "MueLu_TestHelpers.hpp"

#include "MueLu_FactoryFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap; // TODO: remove

  // This is not a real unit test, because output of BuildFactory is not verified. But anyway, it still useful.
  TEUCHOS_UNIT_TEST(FactoryFactory, BuildFactory)
  {
    std::string  factoryConfigurationFiles[] = {"RAPFactory.xml", "SaPFactory.xml", "TentativePFactory.xml", "TrilinosSmoother.xml", "UCAggregationFactory.xml", "DirectSolver.xml"};
    int         nFactoryConfigurationFiles = 6;

    for(int i=0; i< nFactoryConfigurationFiles; i++) {
      out << "Processing file: " << factoryConfigurationFiles[i] << std::endl;
      Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile("ParameterList/FactoryFactory/" + factoryConfigurationFiles[i]);

      for (Teuchos::ParameterList::ConstIterator param = paramList->begin(); param != paramList->end(); ++param) {
        Teuchos::OSTab tab(out);
        const std::string             & paramName  = paramList->name(param);
        const Teuchos::ParameterEntry & paramValue = paramList->entry(param);

        out << "Building object '" << paramName << "'" << std::endl;

        const FactoryMap factoryMapIn;
        FactoryFactory().BuildFactory(paramValue, factoryMapIn);
      }
    }
  }
  
} // namespace MueLuTests


