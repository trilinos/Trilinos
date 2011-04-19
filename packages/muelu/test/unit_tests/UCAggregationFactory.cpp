#include "Teuchos_UnitTestHarness.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UCAggregationFactory.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UCAggregationFactory, Constructor, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps)
  {
#include "MueLu_UseShortNamesOrdinal.hpp"

    out << "version: " << MueLu::Version() << std::endl;
    
    RCP<UCAggregationFactory> aggFact = rcp(new UCAggregationFactory());
    TEUCHOS_TEST_EQUALITY(aggFact != Teuchos::null, true, out, success);
  } // Constructor

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UCAggregationFactory, Build, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps)
  {
#include "MueLu_UseShortNamesOrdinal.hpp"
    
    out << "version: " << MueLu::Version() << std::endl;

    UCAggregationFactory aggFact;

    MueLu::AggregationOptions aggOptions;

    aggOptions.SetPrintFlag(6);
    aggOptions.SetMinNodesPerAggregate(2);
    aggOptions.SetMaxNeighAlreadySelected(5);
    aggOptions.SetOrdering(2);
    aggOptions.SetPhase3AggCreation(0.5);

    RCP<CrsOperator> Op = MueLu::UnitTest::create_1d_poisson_matrix<double,LO,GO>(16);

    RCP<Graph> graph = rcp(new Graph(Op->getCrsGraph(), "someGraphLabel"));
    RCP<Aggregates> aggregates = aggFact.Build(*graph,aggOptions);

    RCP<LOVector> Final_;
    Final_ = LOVectorFactory::Build( aggregates->GetVertex2AggId()->getMap() );

    ArrayRCP<LO> Final = Final_->getDataNonConst(0);
    ArrayRCP<const LO> vertex2AggId = aggregates->GetVertex2AggId()->getData(0);
    ArrayRCP<const LO> procWinner   = aggregates->GetProcWinner()->getData(0);

    for (size_t i=0; i<aggregates->GetVertex2AggId()->getMap()->getNodeNumElements(); i++)
      Final[i] = vertex2AggId[i] + procWinner[i]*1000;

    cout << *Final_ << endl;

  } // Build

  // 
  // INSTANTIATIONS
  //

typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<ScalarType,LocalOrdinal,Node>::SparseOps LocalMatOps;
  
#define UNIT_TEST_GROUP_4(LOCALORDINAL, GLOBALORDINAL, NODE, LOCALMATOPS) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UCAggregationFactory, Constructor, LOCALORDINAL, GLOBALORDINAL, NODE, LOCALMATOPS) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UCAggregationFactory, Build,       LOCALORDINAL, GLOBALORDINAL, NODE, LOCALMATOPS)
    
#define UNIT_TEST_GROUP_2(LOCALORDINAL, GLOBALORDINAL) \
  UNIT_TEST_GROUP_4(LOCALORDINAL, GLOBALORDINAL,  Node, LocalMatOps)

  UNIT_TEST_GROUP_2(int, int)
  //UNIT_TEST_GROUP_2(long, int)
  //UNIT_TEST_GROUP_2(long, long)

} // namespace <anonymous>
