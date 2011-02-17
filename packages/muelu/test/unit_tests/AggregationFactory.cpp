#include "Teuchos_UnitTestHarness.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_AggregationFactory.hpp"

#include "Cthulhu.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

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

  RCP<AggregationFactory> aggFact= rcp(new AggregationFactory());
  //RCP<MueLu::AggregationFactory<LO,GO,NO,LMO> > aggFact= rcp(new AggregationFactory());
  TEUCHOS_TEST_EQUALITY(aggFact != Teuchos::null, true, out, success);
} //Constructor

TEUCHOS_UNIT_TEST(AggregationFactory, GetSetMethods)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  AggregationFactory aggFact;
  std::string alg="Uncoupled";
  aggFact.SetAlgorithm("Uncoupled");
  TEUCHOS_TEST_EQUALITY( aggFAct.GetAlgorithm(), alg, out, success);
} //GetSetMethods


}//namespace <anonymous>

