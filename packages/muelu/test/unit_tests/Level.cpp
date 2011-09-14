#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_config.hpp"
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Version.hpp"



#include "MueLu_Utilities.hpp"


#include "MueLu_TestHelpers.hpp"

#include "MueLu_Level.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"


namespace MueLuTests {

  TEUCHOS_UNIT_TEST(Level, SetCoreData)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level aLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSingleLevelHierarchy(aLevel);

    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(2); //can be an empty operator

    aLevel.Request("Hitchhiker's Guide");
    TEST_EQUALITY(aLevel.IsRequested("Hitchhiker's Guide"), true);
    aLevel.Set("Hitchhiker's Guide",42);
    int fff = aLevel.Get<int>("Hitchhiker's Guide");
    TEST_EQUALITY(fff, 42);
    TEST_EQUALITY(aLevel.IsRequested("PI"), false);
    aLevel.Request("PI");
    TEST_EQUALITY(aLevel.IsRequested("PI"), true);
    aLevel.Set("PI",3.14159265);
    double ggg = aLevel.Get<double>("PI");
    TEST_EQUALITY(ggg, 3.14159265);
    TEST_EQUALITY(aLevel.IsAvailable("PI"), true);
    aLevel.Release("PI");
    TEST_EQUALITY(aLevel.IsAvailable("PI"), false);
    aLevel.Request("Hello MueLu");
    aLevel.Set("Hello MueLu", std::string("Greetings to MueMat"));
    std::string hhh = aLevel.Get<std::string>("Hello MueLu");
    TEST_EQUALITY(hhh, "Greetings to MueMat");
    aLevel.Request("A");
    aLevel.Set("A",A);
    RCP<Operator> newA = aLevel.Get< RCP<Operator> >("A");
    TEST_EQUALITY(newA, A);
    aLevel.Request("R");
    aLevel.Set("R", A);
    RCP<Operator> newR = aLevel.Get< RCP<Operator> >("R");
    TEST_EQUALITY(newR, A); //TODO from JG: must be tested using another matrix !
    aLevel.Request("P");
    aLevel.Set("P", A);
    RCP<Operator> newP = aLevel.Get< RCP<Operator> >("P");
    TEST_EQUALITY(newP, A);
    aLevel.SetLevelID(42);
    TEST_EQUALITY(aLevel.GetLevelID(), 42); //TODO: test default value of LevelID

    //aLevel.print(std::cout);

  }

} // namespace MueLuTests

