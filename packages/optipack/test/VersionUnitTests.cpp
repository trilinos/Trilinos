
#include "OptiPack_Version.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Trilinos_version.h"


namespace {


TEUCHOS_UNIT_TEST( Version, default )
{
  TEST_EQUALITY_CONST(OptiPack::OptiPack_Version(),
    ("OptiPack in Trilinos " TRILINOS_VERSION_STRING));
}


} // namespace


