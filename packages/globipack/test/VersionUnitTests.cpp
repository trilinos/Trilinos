
#include "GlobiPack_Version.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Trilinos_version.h"


namespace {


TEUCHOS_UNIT_TEST( Version, default )
{
  TEST_EQUALITY_CONST(GlobiPack::GlobiPack_Version(),
    ("GlobiPack in Trilinos " TRILINOS_VERSION_STRING));
}


} // namespace


