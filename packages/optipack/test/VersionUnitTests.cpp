
#include "OptiPack_Version.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


TEUCHOS_UNIT_TEST( UnitTests, default )
{
  TEST_EQUALITY_CONST(OptiPack::OptiPack_Version(), "OptiPack Dev in Trilinos Dev");
}


} // namespace


