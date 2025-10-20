#include <Kokkos_Core.hpp>
#include <Teuchos_UnitTestHarness.hpp>

class A {
public:
  Kokkos::Array<int,7> activeDims_;
};

TEUCHOS_UNIT_TEST( PR14546, AllocationIssue )
{
  using BigArray = Kokkos::Array< A, 400000 >;
  BigArray vectorComponents;
}
