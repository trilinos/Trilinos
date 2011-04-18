#include "Teuchos_UnitTestHarness.hpp"
#include "test_helpers.hpp"
#include "Cthulhu_Example.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_UCAggregationFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames_Graph.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(AggregationFactory, Constructor)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  //RCP<UCAggregationFactory> aggFact= rcp(new UCAggregationFactory()); //FIXME
  RCP<MueLu::UCAggregationFactory<LO,GO,NO,LMO> > aggFact;
  aggFact= rcp(new MueLu::UCAggregationFactory<LO,GO,NO,LMO>());
  TEUCHOS_TEST_EQUALITY(aggFact != Teuchos::null, true, out, success);
} //Constructor

TEUCHOS_UNIT_TEST(UCAggregationFactory, Build)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;

  out << "version: " << MueLu::Version() << std::endl;

  //UCAggregationFactory aggFact; //FIXME
  MueLu::UCAggregationFactory<LO,GO,NO,LMO> aggFact;


  MueLu::AggregationOptions aggOptions;

  //these are pulled right from the AggregationExample file.
  int printFlag=6;
  aggOptions.SetPrintFlag(printFlag);
  aggOptions.SetMinNodesPerAggregate(2);
  aggOptions.SetMaxNeighAlreadySelected(5);
  aggOptions.SetOrdering(2);
  aggOptions.SetPhase3AggCreation(0.5);

  RCP<CrsOperator> Op = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(16);

  RCP<Graph> graph = rcp(new Graph(Op->getCrsGraph(), "someGraphLabel"));
  RCP<Aggregates> aggregates = aggFact.Build(*graph,aggOptions);

  RCP<Cthulhu::Vector<int> > Final_;
  Final_ = Cthulhu::VectorFactory<int>::Build( aggregates->GetVertex2AggId()->getMap() );

  ArrayRCP<int> Final = Final_->getDataNonConst(0);
  ArrayRCP<const int> vertex2AggId = aggregates->GetVertex2AggId()->getData(0);
  ArrayRCP<const int> procWinner   = aggregates->GetProcWinner()->getData(0);

  for (size_t i=0; i<aggregates->GetVertex2AggId()->getMap()->getNodeNumElements(); i++)
    Final[i] = vertex2AggId[i] + procWinner[i]*1000;

  cout << *Final_ << endl;

} //Build


}//namespace <anonymous>

