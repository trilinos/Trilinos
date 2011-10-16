/*
 * ThresholdAFilterFactory.cpp
 *
 *  Created on: 16.10.2011
 *      Author: tobias
 */


#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_config.hpp"
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Utilities.hpp"

#include "MueLu_NoFactory.hpp"

#include "MueLu_TestHelpers.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_ThresholdAFilterFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"


namespace MueLuTests {

  TEUCHOS_UNIT_TEST(ThresholdAFilterFactory, Basic)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level aLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSingleLevelHierarchy(aLevel);

    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(20); //can be an empty operator

    RCP<ThresholdAFilterFactory> AfilterFactory0 = rcp(new ThresholdAFilterFactory("A",NULL,0.1)); // keep all
    RCP<ThresholdAFilterFactory> AfilterFactory1 = rcp(new ThresholdAFilterFactory("A",NULL,1.1)); // keep only diagonal
    RCP<ThresholdAFilterFactory> AfilterFactory2 = rcp(new ThresholdAFilterFactory("A",NULL,3)); // keep only diagonal

    aLevel.Set("A",A);

    aLevel.Request("A",AfilterFactory0.get());
    AfilterFactory0->Build(aLevel);
    TEST_EQUALITY(aLevel.IsAvailable("A",AfilterFactory0.get()), true);
    RCP<Operator> A0 = aLevel.Get< RCP<Operator> >("A",AfilterFactory0.get());
    aLevel.Release("A",AfilterFactory0.get());
    TEST_EQUALITY(aLevel.IsAvailable("A",AfilterFactory0.get()), false);
    TEST_EQUALITY(A0->getNodeNumEntries(), A->getNodeNumEntries());
    TEST_EQUALITY(A0->getGlobalNumEntries(), A->getGlobalNumEntries());

    aLevel.Request("A",AfilterFactory1.get());
    AfilterFactory1->Build(aLevel);
    TEST_EQUALITY(aLevel.IsAvailable("A",AfilterFactory1.get()), true);
    RCP<Operator> A1 = aLevel.Get< RCP<Operator> >("A",AfilterFactory1.get());
    aLevel.Release("A",AfilterFactory1.get());
    TEST_EQUALITY(aLevel.IsAvailable("A",AfilterFactory1.get()), false);
    TEST_EQUALITY(A1->getGlobalNumEntries(), A1->getGlobalNumRows());

    aLevel.Request("A",AfilterFactory2.get());
    AfilterFactory2->Build(aLevel);
    TEST_EQUALITY(aLevel.IsAvailable("A",AfilterFactory2.get()), true);
    RCP<Operator> A2 = aLevel.Get< RCP<Operator> >("A",AfilterFactory2.get());
    aLevel.Release("A",AfilterFactory2.get());
    TEST_EQUALITY(aLevel.IsAvailable("A",AfilterFactory2.get()), false);
    TEST_EQUALITY(A2->getGlobalNumEntries(), A2->getGlobalNumRows());


  }
} // end MueLuTests namespace

