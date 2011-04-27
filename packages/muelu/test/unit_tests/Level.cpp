#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;

  TEUCHOS_UNIT_TEST(Level, SetCoreData)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(2); //can be an empty operator

    Level firstLevel;
    firstLevel.SetA(A);
    TEUCHOS_TEST_EQUALITY(firstLevel.GetA(), A, out, success);
    firstLevel.SetR(A);
    TEUCHOS_TEST_EQUALITY(firstLevel.GetR(), A, out, success); //TODO from JG: must be tested using another matrix !
    firstLevel.SetP(A);
    TEUCHOS_TEST_EQUALITY(firstLevel.GetP(), A, out, success);
    firstLevel.SetLevelID(42);
    TEUCHOS_TEST_EQUALITY(firstLevel.GetLevelID(), 42, out, success); //TODO: test default value of LevelID

    /*
      RCP<Smoother> preSmoo = Smoother<Scalar, LO, GO, Node, LMO>();
      TEUCHOS_TEST_EQUALITY(firstLevel.GetPreSmoother(), preSmoo, out, success);
      //RCP<Smoother> postSmoo = Smoother<Scalar, LO, GO, Map, Operator>();
      */


    //out << firstLevel << std::endl;
    /*
      out << "Testing copy ctor" << std::endl;
      Level secondLevel(firstLevel);
      //out << secondLevel << std::endl;
      */
  }

}//namespace <anonymous>

