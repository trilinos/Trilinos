#include "Teuchos_UnitTestHarness.hpp"
#include "test_helpers.hpp"
#include "Cthulhu.hpp"
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

TEUCHOS_UNIT_TEST(UCAggregationFactory, GetSetMethods)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  //UCAggregationFactory aggFact; //FIXME
  MueLu::UCAggregationFactory<LO,GO,NO,LMO> aggFact;
} //GetSetMethods


}//namespace <anonymous>

