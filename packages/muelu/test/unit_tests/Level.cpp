#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_TestHelpers.hpp"

#include "MueLu_Level.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(Level, SetCoreData)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(2); //can be an empty operator

    RCP<DefaultFactoryHandlerBase> defaultFactHandler = rcp(new DefaultFactoryHandler());
    Level firstLevel(defaultFactHandler);
    firstLevel.Set("A",A);
    RCP<Operator> newA = firstLevel.Get< RCP<Operator> >("A");
    TEUCHOS_TEST_EQUALITY(newA, A, out, success);
    firstLevel.Set("R", A);
    TEUCHOS_TEST_EQUALITY(firstLevel.Get< RCP<Operator> >("R"), A, out, success); //TODO from JG: must be tested using another matrix !
    firstLevel.Set("P", A);
    TEUCHOS_TEST_EQUALITY(firstLevel.Get< RCP<Operator> >("P"), A, out, success);
    firstLevel.SetLevelID(42);
    TEUCHOS_TEST_EQUALITY(firstLevel.GetLevelID(), 42, out, success); //TODO: test default value of LevelID

    /*
      RCP<Smoother> preSmoo = Smoother<Scalar, LO, GO, Node, LMO>();
      TEUCHOS_TEST_EQUALITY(firstLevel.Get< RCP<SmootherPrototype> >("PreSmoother"), preSmoo, out, success);
      //RCP<Smoother> postSmoo = Smoother<Scalar, LO, GO, Map, Operator>();
      */


    //out << firstLevel << std::endl;
    /*
      out << "Testing copy ctor" << std::endl;
      Level secondLevel(firstLevel);
      //out << secondLevel << std::endl;
      */
  }

}//namespace MueLuTests

