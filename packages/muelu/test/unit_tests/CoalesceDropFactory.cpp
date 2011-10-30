#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<CoalesceDropFactory> coalesceDropFact = rcp(new CoalesceDropFactory());
    TEST_EQUALITY(coalesceDropFact != Teuchos::null, true);

  } //Constructor

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, Build)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level fineLevel;
    TestHelpers::Factory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(fineLevel);

    RCP<Operator> A = TestHelpers::Factory<SC,LO,GO,NO,LMO>::Build1DPoisson(36);
    fineLevel.Set("A", A);

    CoalesceDropFactory coalesceDropFact;
    coalesceDropFact.Build(fineLevel);
    //FIXME how do we verify that this is correct?
  } //Build

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, PreDrop)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level fineLevel;
    TestHelpers::Factory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(fineLevel);

    RCP<Operator> A = TestHelpers::Factory<SC,LO,GO,NO,LMO>::Build1DPoisson(3);
    fineLevel.Set("A", A);
    A->describe(out,Teuchos::VERB_EXTREME);

    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetVerbLevel(MueLu::Extreme);
    dropFact.SetPreDropFunction(rcp(new PreDropFunctionConstVal(0.00001)));

    fineLevel.Request("Graph", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<Graph> graph = fineLevel.Get<RCP<Graph> >("Graph", &dropFact);

    std::cout << graph->GetDomainMap()->getGlobalNumElements() << std::endl;
    graph->print(out, MueLu::Debug);

//    TEST_EQUALITY(1 == 0, true);

  } //PreDrop

} // namespace MueLuTests

