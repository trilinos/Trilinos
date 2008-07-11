#include "Teuchos_UnitTestHarness.hpp"


namespace {


TEUCHOS_UNIT_TEST( std_map, basic )
{
  std::map<std::string, int> myMap;
  TEST_THROW( *myMap.begin(), std::exception );
}


} // namespace
