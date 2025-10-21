#include <Teuchos_UnitTestHarness.hpp>

#include <array>

TEUCHOS_UNIT_TEST( PR14546, AllocationIssue )
{
  std::array<int, 2800000> array;
}
