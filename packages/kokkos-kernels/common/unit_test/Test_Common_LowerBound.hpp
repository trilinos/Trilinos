//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

/// \file Test_Common_LowerBound.hpp
/// \brief Tests lower bounds search routines

#include <KokkosKernels_LowerBound.hpp>
#include <KokkosKernels_ExecSpaceUtils.hpp>

template <typename Ordinal>
size_t std_lower_bound(const std::vector<Ordinal> &haystack, const Ordinal needle) {
  const auto it = std::lower_bound(haystack.begin(), haystack.end(), needle);
  return it - haystack.begin();
}

/*! \brief count the number of incorrect values */
template <typename HaystackView>
struct ThreadLowerBoundFunctor {
  using hv_value_type = typename HaystackView::non_const_value_type;
  using hv_size_type  = typename HaystackView::size_type;

  ThreadLowerBoundFunctor(const hv_size_type &expected, const HaystackView &haystack, const hv_value_type &needle)
      : expected_(expected), haystack_(haystack), needle_(needle) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t i, int &lerrCount) const {
    if (0 == i) {
      hv_size_type idx = KokkosKernels::lower_bound_thread(haystack_, needle_);
      if (idx != expected_) {
        Kokkos::printf("%s:%d thread %d expected %d got %d\n", __FILE__, __LINE__, int(i), int(expected_), int(idx));
        ++lerrCount;
      }
    }
  }

  hv_size_type expected_;
  HaystackView haystack_;
  hv_value_type needle_;
};

template <typename T, typename Device>
void test_lower_bound_thread(const std::vector<T> &_haystack, const T &_needle) {
  using execution_space = typename Device::execution_space;
  using Policy          = Kokkos::RangePolicy<execution_space>;
  using view_t          = Kokkos::View<T *, Device>;
  using u_const_view_t  = Kokkos::View<const T *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using size_type       = typename u_const_view_t::size_type;

  // get expected value
  const size_type expected = std_lower_bound(_haystack, _needle);

  // create device views of input data
  u_const_view_t uhaystack(_haystack.data(), _haystack.size());
  view_t haystack("haystack", uhaystack.size());
  Kokkos::deep_copy(haystack, uhaystack);

  // test lower_bound search
  int errCount;
  // run a single thread
  Kokkos::parallel_reduce(Policy(0, 1), ThreadLowerBoundFunctor(expected, haystack, _needle), errCount);

  EXPECT_EQ(0, errCount);
}

/*! \brief count the number of incorrect values */
template <typename Member, typename HaystackView>
struct TeamLowerBoundFunctor {
  using hv_value_type = typename HaystackView::non_const_value_type;
  using hv_size_type  = typename HaystackView::size_type;

  TeamLowerBoundFunctor(const hv_size_type &expected, const HaystackView &haystack, const hv_value_type &needle)
      : expected_(expected), haystack_(haystack), needle_(needle) {}

  KOKKOS_INLINE_FUNCTION void operator()(const Member &handle, int &lerrCount) const {
    hv_size_type idx = KokkosKernels::lower_bound_team(handle, haystack_, needle_);
    if (idx != expected_) {
      Kokkos::printf("%s:%d thread %d expected %d got %d\n", __FILE__, __LINE__, int(handle.team_rank()),
                     int(expected_), int(idx));
      ++lerrCount;
    }
  }

  hv_size_type expected_;
  HaystackView haystack_;
  hv_value_type needle_;
};

template <typename T, typename Device>
void test_lower_bound_team(const std::vector<T> &_haystack, const T _needle) {
  using execution_space = typename Device::execution_space;
  using Policy          = Kokkos::TeamPolicy<execution_space>;
  using Member          = typename Policy::member_type;
  using view_t          = Kokkos::View<T *, Device>;
  using u_const_view_t  = Kokkos::View<const T *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using size_type       = typename u_const_view_t::size_type;

  // get expected value
  const size_type expected = std_lower_bound(_haystack, _needle);

  // create device views of input data
  u_const_view_t uhaystack(_haystack.data(), _haystack.size());
  view_t haystack("haystack", uhaystack.size());
  Kokkos::deep_copy(haystack, uhaystack);

  // test lower_bound search
  const int leagueSize = 1;
  const int teamSize   = KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>() ? 64 : 1;
  int errCount;
  Kokkos::parallel_reduce(Policy(leagueSize, teamSize),
                          TeamLowerBoundFunctor<Member, view_t>(expected, haystack, _needle), errCount);

  EXPECT_EQ(0, errCount);
}

template <typename T, typename Device>
void test_lower_bound(const std::vector<T> &haystack, const T needle) {
  test_lower_bound_thread<T, Device>(haystack, needle);
  test_lower_bound_team<T, Device>(haystack, needle);
}

template <typename T>
T randn(T n) {
  if constexpr (std::is_floating_point_v<T>) {
    return T(rand()) / T(RAND_MAX) * n;
  } else {
    return T(rand()) % n;
  }
}

/* define specific and random lower-bound test cases
 */
template <typename T, typename Device>
void test_lower_bound() {
  test_lower_bound<T, Device>({}, T(0));
  test_lower_bound<T, Device>({}, T(1));
  test_lower_bound<T, Device>({}, T(-1));

  test_lower_bound<T, Device>({0}, T(0));
  test_lower_bound<T, Device>({0}, T(1));
  test_lower_bound<T, Device>({0}, T(-1));

  test_lower_bound<T, Device>({1}, T(0));
  test_lower_bound<T, Device>({1}, T(1));
  test_lower_bound<T, Device>({1}, T(-1));

  test_lower_bound<T, Device>({T(-1)}, T(0));
  test_lower_bound<T, Device>({T(-1)}, T(1));
  test_lower_bound<T, Device>({T(-1)}, T(-1));

  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(-1));
  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(0));
  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(1));
  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(2));
  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(2.4));
  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(2.5));
  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(2.6));
  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(3));
  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(4));
  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(5));
  test_lower_bound<T, Device>({0, 1, T(2.5), 3, 4, 5}, T(6));

  auto randn = [](T n) -> T {
    T ret;
    if constexpr (std::is_floating_point_v<T>) {
      ret = T(rand()) / T(RAND_MAX) * n;
    } else {
      ret = T(rand()) % n;
    }
    return ret;
  };

  T maxEntry         = 20;
  const int numTests = 100;
  for (int n = 0; n < numTests; ++n) {
    for (size_t sz : {10, 100, 1000}) {
      // generate a sorted random vector
      std::vector<T> haystack;
      for (size_t i = 0; i < sz; ++i) {
        haystack.push_back(randn(maxEntry));
      }
      std::sort(haystack.begin(), haystack.end());

      // generate a random value to search for
      const T needle = randn(maxEntry);

      // do the test
      test_lower_bound<T, Device>(haystack, needle);
    }
  }
}

#define EXECUTE_TEST(T, DEVICE) \
  TEST_F(TestCategory, common##_##lower_bound##_##T##_##DEVICE) { test_lower_bound<T, DEVICE>(); }

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(int, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(int64_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_ORDINAL_INT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(size_t, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_FLOAT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_DOUBLE)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, TestDevice)
#endif

#undef EXECUTE_TEST
