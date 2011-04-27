#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_TransPFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"  
#include "MueLu_UseShortNames.hpp"  

namespace {

  TEUCHOS_UNIT_TEST(TransPFactory, Test0)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<TransPFactory> transPFact = rcp(new TransPFactory);
    TEUCHOS_TEST_EQUALITY(transPFact != Teuchos::null, true, out, success);

    out << *transPFact << std::endl;

  }

}//namespace <anonymous>

