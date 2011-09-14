#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_RAPFactory.hpp"
#include "MueLu_DefaultFactoryHandler.hpp"
#include "MueLu_UCAggregationFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  // Little utility to generate aggregates.
  RCP<Aggregates> gimmeAggregates(RCP<Operator> const &A)
  {
    UCAggregationFactory aggFact;
    aggFact.SetMinNodesPerAggregate(3);
    aggFact.SetMaxNeighAlreadySelected(0);
    aggFact.SetOrdering(MueLu::AggOptions::NATURAL);
    aggFact.SetPhase3AggCreation(0.5);

    Level level;
    TestHelpers::Factory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(level);
    level.Request("A",NULL);	// use default factory handler
    level.Set("A",A,NULL);
    level.Request("Aggregates", &aggFact);
    aggFact.DeclareInput(level);
    aggFact.Build(level);
    RCP<Aggregates> aggregates = level.Get<RCP<Aggregates> >("Aggregates",Teuchos::rcp(&aggFact,false)); // fix me
    level.Release("Aggregates", &aggFact);
    return aggregates;
  }  //gimmeAggregates

  TEUCHOS_UNIT_TEST(Aggregates, JustAggregation)
  {
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(15);
    RCP<Aggregates> aggregates = gimmeAggregates(A);
    TEST_EQUALITY(aggregates != Teuchos::null, true);
  }

///////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST(Aggregates, GetNumAggregates)
  {
      out << "version: " << MueLu::Version() << std::endl;

      RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(36);
      RCP<const Map> rowmap = A->getRowMap();
      RCP<Aggregates> aggregates = gimmeAggregates(A);
      GO numAggs = aggregates->GetNumAggregates();
      RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

      ArrayRCP<GO> aggSizes  = aggregates->ComputeAggregateSizes();
      bool foundAggNotSize3=false;
      for (int i=0; i<aggSizes.size(); ++i)
        if (aggSizes[i] != 3) {
          foundAggNotSize3=true;
          break;
        }

      switch (comm->getSize()) {

        case 1 :
           TEST_EQUALITY(numAggs, 12);
           TEST_EQUALITY(foundAggNotSize3, false);
           break;

        case 2:
           TEST_EQUALITY(numAggs, 6);
           TEST_EQUALITY(foundAggNotSize3, false);
           break;

        case 3:
           TEST_EQUALITY(numAggs, 4);
           TEST_EQUALITY(foundAggNotSize3, false);
           break;

        case 4:
           TEST_EQUALITY(numAggs, 3);
           TEST_EQUALITY(foundAggNotSize3, false);
           break;

        default:
           std::string msg = "Only 1-4 MPI processes are supported.";
           //throw(MueLu::Exceptions::NotImplemented(msg));
           out << msg << std::endl;
           break;
      }

      ArrayRCP< ArrayRCP<GO> > aggToRowMap(numAggs);
      aggregates->ComputeAggregateToRowMap(aggToRowMap);
      int root = out.getOutputToRootOnly();
      out.setOutputToRootOnly(-1);
      for (int j=0; j<comm->getSize(); ++j) {
        if (comm->getRank() == j) {
            out << "++ pid " << j << " ++" << std::endl;
            out << "   num local DOFs = " << rowmap->getNodeNumElements() << std::endl;
          for (int i=0; i< aggToRowMap.size(); ++i) {
            out << "   aggregate " << i << ": ";
            for (int k=0; k< aggToRowMap[i].size(); ++k)
              out << aggToRowMap[i][k] << " ";
            out << std::endl;
          }
        }
        comm->barrier();
      }
      out.setOutputToRootOnly(root);

  } //GetNumAggregates


} // namespace MueLuTests
