#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UCAggregationFactory.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;

  MUELU_UNIT_TEST_TEMPLATE_5_DECL(UCAggregationFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps)
  {
#include "MueLu_UseShortNames.hpp"

    out << "version: " << MueLu::Version() << std::endl;
    
    RCP<UCAggregationFactory> aggFact = rcp(new UCAggregationFactory());
    TEUCHOS_TEST_EQUALITY(aggFact != Teuchos::null, true, out, success);
  } // Constructor
  MUELU_UNIT_TEST_DECL_END()

  MUELU_UNIT_TEST_TEMPLATE_5_DECL(UCAggregationFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps)
  {
    //    typedef double Scalar;
#include "MueLu_UseShortNames.hpp"
    
    out << "version: " << MueLu::Version() << std::endl;

    UCAggregationFactory aggFact;

    MueLu::AggregationOptions aggOptions;

    aggOptions.SetPrintFlag(6);
    aggOptions.SetMinNodesPerAggregate(2);
    aggOptions.SetMaxNeighAlreadySelected(5);
    aggOptions.SetOrdering(2);
    aggOptions.SetPhase3AggCreation(0.5);

    RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(16);

    RCP<Graph> graph = rcp(new Graph(Op->getCrsGraph(), "someGraphLabel"));

    // This test must be fixed after Jonathan modifications
    success=false;

    //    RCP<Aggregates> aggregates = aggFact.Build(*graph,aggOptions);

//     RCP<LOVector> Final_;
//     Final_ = LOVectorFactory::Build( aggregates->GetVertex2AggId()->getMap() );

//     ArrayRCP<LO> Final = Final_->getDataNonConst(0);
//     ArrayRCP<const LO> vertex2AggId = aggregates->GetVertex2AggId()->getData(0);
//     ArrayRCP<const LO> procWinner   = aggregates->GetProcWinner()->getData(0);

//     for (size_t i=0; i<aggregates->GetVertex2AggId()->getMap()->getNodeNumElements(); i++)
//       Final[i] = vertex2AggId[i] + procWinner[i]*1000;

//     cout << *Final_ << endl;

  } // Build
  MUELU_UNIT_TEST_DECL_END()

  // 
  // INSTANTIATIONS
  //

  typedef double Scalar;                             // Scalar is not relevant for this test
  typedef Kokkos::DefaultNode::DefaultNodeType Node; // Kokkos Node is not relevant for this test   

  typedef long int LongInt;                          // macros dislike parameters with space...
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  typedef long long int LongLongInt;  
#endif
  
#define UNIT_TEST_GROUP_5(SC, LO, GO, NO, LMO)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(UCAggregationFactory, Constructor, SC, LO, GO, NO, LMO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(UCAggregationFactory, Build,       SC, LO, GO, NO, LMO)

#define UNIT_TEST_GROUP_2(LO, GO)                                      \
  typedef Kokkos::DefaultKernels<Scalar,LO,Node>::SparseOps LMO ## LO; \
  UNIT_TEST_GROUP_5(Scalar, LO, GO, Node, LMO ## LO)

  UNIT_TEST_GROUP_2(int, int)
  UNIT_TEST_GROUP_2(int, LongInt)
  UNIT_TEST_GROUP_2(LongInt, LongInt)

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  UNIT_TEST_GROUP_2(LongInt, LongLongInt)
#endif

} // namespace <anonymous>
