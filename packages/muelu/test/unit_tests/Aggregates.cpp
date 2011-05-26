#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UCAggregationFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;

  // Little utility to generate aggregates.
  RCP<Aggregates> gimmeAggregates(RCP<Operator> const &A)
  {
    RCP<Graph> graph = rcp(new Graph(A->getCrsGraph(), "someGraphLabel"));

    MueLu::AggregationOptions aggOptions;

    aggOptions.SetPrintFlag(0);
    aggOptions.SetMinNodesPerAggregate(3);
    aggOptions.SetMaxNeighAlreadySelected(0);
    aggOptions.SetOrdering(MueLu::AggOptions::NATURAL);
    aggOptions.SetPhase3AggCreation(0.5);

    UCAggregationFactory aggFact(aggOptions);
    RCP<Aggregates> aggregates = aggFact.Build(*graph);
    return aggregates;
  }  //gimmeAggregates

  TEUCHOS_UNIT_TEST(Aggregates, JustAggregation)
  {
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(15);
    RCP<Aggregates> aggregates = gimmeAggregates(A);
  }

///////////////////////////////////////////////////////////////////////////

  TEUCHOS_UNIT_TEST(Aggregates, GetNumAggregates)
  {
      out << "version: " << MueLu::Version() << std::endl;

      //RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(36);
      RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(39);
      RCP<const Map> rowmap = A->getRowMap();
      RCP<Aggregates> aggregates = gimmeAggregates(A);
      GO numAggs = aggregates->GetNumAggregates();
      RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();

      switch (comm->getSize()) {

        case 1 :
           TEUCHOS_TEST_EQUALITY(numAggs, 12, out, success);
           break;

        case 2:
           TEUCHOS_TEST_EQUALITY(numAggs, 6, out, success);
           break;

        case 3:
           TEUCHOS_TEST_EQUALITY(numAggs, 4, out, success);
           break;

        case 4:
           TEUCHOS_TEST_EQUALITY(numAggs, 3, out, success);
           break;

        default:
           std::string msg = "Only 1-4 MPI processes are supported.";
           //throw(MueLu::Exceptions::NotImplemented(msg));
           out << msg << endl;
           break;
      }

      ArrayRCP<GO> aggSizes  = aggregates->ComputeAggregateSizes();
      bool foundAggNotSize3=false;
      for (int i=0; i<aggSizes.size(); ++i)
        if (aggSizes[i] != 3) {
          foundAggNotSize3=true;
          break;
        }

      TEUCHOS_TEST_EQUALITY(foundAggNotSize3, false, out, success);

      ArrayRCP< ArrayRCP<GO> > aggToRowMap(numAggs);
      aggregates->ComputeAggregateToRowMap(aggToRowMap);
      for (int j=0; j<comm->getSize(); ++j) {
        if (comm->getRank() == j) {
            std::cout << "++ pid " << j << " ++" << std::endl;
            std::cout << "   num local DOFs = " << rowmap->getNodeNumElements() << std::endl;
          for (int i=0; i< aggToRowMap.size(); ++i) {
            std::cout << "   aggregate " << i << ": ";
            for (int k=0; k< aggToRowMap[i].size(); ++k)
              std::cout << aggToRowMap[i][k] << " ";
            std::cout << std::endl;
          }
        }
        comm->barrier();
      }

  } //GetNumAggregates


} // namespace <anonymous>
