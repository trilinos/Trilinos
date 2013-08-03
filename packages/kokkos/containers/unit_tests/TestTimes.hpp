#ifndef KOKKOS_CONTAINERS_UNITTESTS_TEST_TIMES_HPP
#define KOKKOS_CONTAINERS_UNITTESTS_TEST_TIMES_HPP

namespace Test {

struct map_test_times
{
  map_test_times()
    : construct(0.0)
    , santity_check(0.0)
    , insert(0.0)
    , find(0.0)
  {}

  double construct;
  double santity_check;
  double insert;
  double find;
};

} // namespace Test

#endif //KOKKOS_CONTAINERS_UNITTESTS_TEST_TIMES_HPP
