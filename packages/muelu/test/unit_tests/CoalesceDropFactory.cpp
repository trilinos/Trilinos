#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_CoalesceDropFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(CoalesceDropFactory, Constructor)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<CoalesceDropFactory> cdFact = rcp(new CoalesceDropFactory);
  TEUCHOS_TEST_EQUALITY(cdFact != Teuchos::null, true, out, success);

} //Constructor

TEUCHOS_UNIT_TEST(CoalesceDropFactory, Build)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  CoalesceDropFactory cdFact;
  RCP<Operator> A = MueLu::TestHelpers::Factory<SC,LO,GO,NO,LMO>::Build1DPoisson(36);
  Level fineLevel;
  fineLevel.SetA(A);

  cdFact.Build(fineLevel);
  //FIXME how do we verify that this is correct?
} //Build

}//namespace <anonymous>

