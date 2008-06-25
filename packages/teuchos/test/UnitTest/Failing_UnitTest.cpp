#include "Teuchos_UnitTestHarness.hpp"


namespace {


TEUCHOS_UNIT_TEST( Int,  BadAssignment )
{
  int i1 = 4;
  int i2 = i1 + 1;
  TEST_EQUALITY( i2, i1 );
}


} // namespace
