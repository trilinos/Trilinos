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
      // RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(15);
//       RCP<Aggregates> aggregates = gimmeAggregates(A);
  }

} // namespace <anonymous>
