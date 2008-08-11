#include "Teuchos_UnitTestHarness.hpp"


namespace {


TEUCHOS_UNIT_TEST( Int, Basic )
{
  int i1 = 5;
  TEST_EQUALITY_CONST( i1, 5 );
}


TEUCHOS_UNIT_TEST( Int, Assignment )
{
  int i1 = 4;
  int i2 = i1;
  TEST_EQUALITY( i2, i1 );
}


} // namespace
