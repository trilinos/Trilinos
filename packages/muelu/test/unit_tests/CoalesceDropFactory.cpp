#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_CoalesceDropFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<CoalesceDropFactory> cdFact = rcp(new CoalesceDropFactory);
    TEUCHOS_TEST_EQUALITY(cdFact != Teuchos::null, true, out, success);

  } //Constructor

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, Build)
  {
    out << "version: " << MueLu::Version() << std::endl;

    CoalesceDropFactory cdFact;
    RCP<Operator> A = MueLu::TestHelpers::Factory<SC,LO,GO,NO,LMO>::Build1DPoisson(36);
    Level fineLevel;
    fineLevel.SetA(A);

    cdFact.Build(fineLevel);
    //FIXME how do we verify that this is correct?
  } //Build

}//namespace <anonymous>

