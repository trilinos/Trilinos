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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "KokkosGraph_Merge.hpp"

#include <vector>

/* test diagonal search on the host
 */
template <typename ordinal_type>
void test_diagonal_search() {
  using iota_type    = KokkosKernels::Impl::Iota<ordinal_type, size_t>;
  using view_type    = Kokkos::View<ordinal_type *, Kokkos::HostSpace>;
  using um_view_type = Kokkos::View<ordinal_type *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using DSR = KokkosGraph::Impl::DiagonalSearchResult<typename view_type::size_type, typename iota_type::size_type>;

  {
    /*   a->
        x 2 2 2 2 4
      b  *       :
      | 0|       :
      v  *       :
        1|       :
         *-*-*-*-*
        2        |
                 *
        3        |
                 *-*
    */

    iota_type iota(4);
    view_type v("", 5);
    v(0) = 2;
    v(1) = 2;
    v(2) = 2;
    v(3) = 2;
    v(4) = 4;
#define KK_TEST_DIAG(_d, _ai, _bi)                             \
  {                                                            \
    DSR dsr = KokkosGraph::Impl::diagonal_search(v, iota, _d); \
    EXPECT_EQ(static_cast<decltype(dsr.ai)>(_ai), dsr.ai);     \
    EXPECT_EQ(static_cast<decltype(dsr.bi)>(_bi), dsr.bi);     \
  }
    KK_TEST_DIAG(0, 0, 0);
    KK_TEST_DIAG(1, 0, 1);
    KK_TEST_DIAG(2, 0, 2);
    KK_TEST_DIAG(3, 1, 2);
    KK_TEST_DIAG(4, 2, 2);
    KK_TEST_DIAG(5, 3, 2);
    KK_TEST_DIAG(6, 4, 2);
    KK_TEST_DIAG(7, 4, 3);
    KK_TEST_DIAG(8, 4, 4);
    KK_TEST_DIAG(9, 5, 4);
#undef KK_TEST_DIAG
  }

  {
    /*  From Odeh, Green, Mwassi, Shmueli, Birk
        Merge Path - Parallel Merging Made Simple
        b ->
      a   3 5 12 22 45 64 69 82
      | 17
      v 29
        35
        73
        86
        90
        95
        99
    */

    std::vector<ordinal_type> _a{17, 29, 35, 75, 86, 90, 95, 99};
    std::vector<ordinal_type> _b{3, 5, 12, 22, 45, 64, 69, 82};
    um_view_type a(_a.data(), _a.size());
    um_view_type b(_b.data(), _b.size());
#define KK_TEST_DIAG(_d, _ai, _bi)                           \
  {                                                          \
    auto dsr = KokkosGraph::Impl::diagonal_search(a, b, _d); \
    EXPECT_EQ(static_cast<decltype(dsr.ai)>(_ai), dsr.ai);   \
    EXPECT_EQ(static_cast<decltype(dsr.bi)>(_bi), dsr.bi);   \
  }
    KK_TEST_DIAG(0, 0, 0);
    KK_TEST_DIAG(1, 0, 1);
    KK_TEST_DIAG(2, 0, 2);
    KK_TEST_DIAG(3, 0, 3);
    KK_TEST_DIAG(4, 1, 3);
    KK_TEST_DIAG(5, 1, 4);
    KK_TEST_DIAG(6, 2, 4);
    KK_TEST_DIAG(7, 3, 4);
    KK_TEST_DIAG(8, 3, 5);
    KK_TEST_DIAG(9, 3, 6);
    KK_TEST_DIAG(10, 3, 7);
    KK_TEST_DIAG(11, 4, 7);
    KK_TEST_DIAG(12, 4, 8);
    KK_TEST_DIAG(13, 5, 8);
    KK_TEST_DIAG(14, 6, 8);
    KK_TEST_DIAG(15, 7, 8);
    KK_TEST_DIAG(16, 8, 8);
#undef KK_TEST_DIAG
  }
}

template <typename Ordinal>
std::vector<Ordinal> vector_merge(const std::vector<Ordinal> &a, const std::vector<Ordinal> &b) {
  size_t ai = 0;
  size_t bi = 0;
  std::vector<Ordinal> c;
  while (ai < a.size() && bi < b.size()) {
    if (a[ai] < b[bi]) {
      c.push_back(a[ai++]);
    } else {
      c.push_back(b[bi++]);
    }
  }
  while (ai < a.size()) {
    c.push_back(a[ai++]);
  }
  while (bi < b.size()) {
    c.push_back(b[bi++]);
  }
  return c;
}

/* compare expected and actual
   if errors, print a, b, expected, and actual
*/
template <typename Ordinal>
bool vector_compare(const std::vector<Ordinal> &a, const std::vector<Ordinal> &b, const std::vector<Ordinal> &exp,
                    const std::vector<Ordinal> &act) {
  size_t err = 0;
  for (size_t i = 0; i < act.size() || i < exp.size(); ++i) {
    if (i < exp.size() && i < act.size()) {
      if (act[i] != exp[i]) {
        ++err;  // value mimatch
      }
    } else {
      ++err;  // size mismatch
    }
  }
  if (err != 0) {
    std::cerr << "i\ta\tb\texp.\tact." << std::endl;
    std::cerr << "-\t-\t-\t----\t----" << std::endl;
    for (size_t i = 0; i < exp.size() || i < act.size(); ++i) {
      std::cerr << i << "\t";
      if (i < a.size()) {
        std::cerr << a[i] << "\t";
      } else {
        std::cerr << " "
                  << "\t";
      }

      if (i < b.size()) {
        std::cerr << b[i] << "\t";
      } else {
        std::cerr << " "
                  << "\t";
      }

      if (i < exp.size()) {
        std::cerr << exp[i] << "\t";
      } else {
        std::cerr << " "
                  << "\t";
      }

      std::cerr << "\t";
      if (i < act.size()) {
        std::cerr << act[i];
      } else {
        std::cerr << " ";
      }

      std::cerr << "\n";
    }
  }
  return 0 == err;
}

template <typename ViewType>
struct ThreadMergeIntoFunctor {
  ThreadMergeIntoFunctor(const ViewType &c, const ViewType &a, const ViewType &b) : c_(c), a_(a), b_(b) {}

  // one thread does the whole merge
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t t) const {
    if (0 == t) {
      KokkosGraph::merge_into_thread(c_, a_, b_);
    }
  }

  ViewType c_;
  ViewType a_;
  ViewType b_;
};

template <typename Ordinal, typename Device>
std::vector<Ordinal> run_thread_merge_into(const std::vector<Ordinal> &_a, const std::vector<Ordinal> &_b) {
  using execution_space = typename Device::execution_space;
  using Policy          = Kokkos::RangePolicy<execution_space>;
  using view_t          = Kokkos::View<Ordinal *, Device>;
  using u_const_view_t  = Kokkos::View<const Ordinal *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using u_view_t        = Kokkos::View<Ordinal *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // create device views of input data
  u_const_view_t ua(_a.data(), _a.size());
  u_const_view_t ub(_b.data(), _b.size());
  view_t a("a", ua.size());
  Kokkos::deep_copy(a, ua);
  view_t b("b", ub.size());
  Kokkos::deep_copy(b, ub);

  // device view of output data
  view_t c("c", a.size() + b.size());

  // run merge
  Kokkos::parallel_for(Policy(0, 1) /*one thread will merge*/, ThreadMergeIntoFunctor(c, a, b));

  // return output to host
  std::vector<Ordinal> _c(c.size());
  u_view_t uc(_c.data(), _c.size());
  Kokkos::deep_copy(uc, c);
  Kokkos::fence();

  return _c;
}

template <typename ViewType>
struct TeamMergeIntoFunctor {
  TeamMergeIntoFunctor(const ViewType &c, const ViewType &a, const ViewType &b) : c_(c), a_(a), b_(b) {}

  template <typename Member>
  KOKKOS_INLINE_FUNCTION void operator()(const Member &handle) const {
    KokkosGraph::merge_into_team(handle, c_, a_, b_);
  }

  ViewType c_;
  ViewType a_;
  ViewType b_;
};

template <typename Ordinal, typename Device>
std::vector<Ordinal> run_merge_into_team(const std::vector<Ordinal> &_a, const std::vector<Ordinal> &_b) {
  using execution_space = typename Device::execution_space;
  using Policy          = Kokkos::TeamPolicy<execution_space>;
  using view_t          = Kokkos::View<Ordinal *, Device>;
  using u_const_view_t  = Kokkos::View<const Ordinal *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using u_view_t        = Kokkos::View<Ordinal *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // create device views of input data
  u_const_view_t ua(_a.data(), _a.size());
  u_const_view_t ub(_b.data(), _b.size());
  view_t a("a", ua.size());
  Kokkos::deep_copy(a, ua);
  view_t b("b", ub.size());
  Kokkos::deep_copy(b, ub);

  // device view of output data
  view_t c("c", a.size() + b.size());

  // one team, small for non-GPU spaces
  const int leagueSize = 1;
  int teamSize;
  Policy policy;
  if constexpr (KokkosKernels::Impl::is_gpu_exec_space_v<execution_space>) {
    teamSize = 64;
  } else {
    teamSize = 1;
  }

  // run merge
  Kokkos::parallel_for(Policy(leagueSize, teamSize), TeamMergeIntoFunctor(c, a, b));

  // return output to host
  std::vector<Ordinal> _c(c.size());
  u_view_t uc(_c.data(), _c.size());
  Kokkos::deep_copy(uc, c);
  Kokkos::fence();

  return _c;
}

template <typename Ordinal, typename Device>
std::vector<Ordinal> run_global_merge(const std::vector<Ordinal> &_a, const std::vector<Ordinal> &_b) {
  using view_t         = Kokkos::View<Ordinal *, Device>;
  using u_const_view_t = Kokkos::View<const Ordinal *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using u_view_t       = Kokkos::View<Ordinal *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // create device views of input data
  u_const_view_t ua(_a.data(), _a.size());
  u_const_view_t ub(_b.data(), _b.size());
  view_t a("a", ua.size());
  Kokkos::deep_copy(a, ua);
  view_t b("b", ub.size());
  Kokkos::deep_copy(b, ub);

  // run merge
  view_t c = KokkosGraph::merge<view_t>(a, b);

  // return output to host
  std::vector<Ordinal> _c(c.size());
  u_view_t uc(_c.data(), _c.size());
  Kokkos::deep_copy(uc, c);
  Kokkos::fence();

  return _c;
}

template <typename Ordinal, typename Device>
void test_merge_into_thread(const std::vector<Ordinal> &a, const std::vector<Ordinal> &b) {
  // expected values
  std::vector<Ordinal> c_exp = vector_merge(a, b);

  // result under test
  std::vector<Ordinal> c_act = run_thread_merge_into<Ordinal, Device>(a, b);

  EXPECT_TRUE(vector_compare(a, b, c_exp, c_act));
}

template <typename Ordinal, typename Device>
void test_merge_into_team(const std::vector<Ordinal> &a, const std::vector<Ordinal> &b) {
  // expected values
  std::vector<Ordinal> c_exp = vector_merge(a, b);

  // result under test
  std::vector<Ordinal> c_act = run_merge_into_team<Ordinal, Device>(a, b);

  EXPECT_TRUE(vector_compare(a, b, c_exp, c_act));
}

template <typename Ordinal, typename Device>
void test_global_merge(const std::vector<Ordinal> &a, const std::vector<Ordinal> &b) {
  // expected values
  std::vector<Ordinal> c_exp = vector_merge(a, b);

  // result under test
  std::vector<Ordinal> c_act = run_global_merge<Ordinal, Device>(a, b);

  EXPECT_TRUE(vector_compare(a, b, c_exp, c_act));
}

/* test merging `a` and `b` on Device
 */
template <typename T, typename Device>
void test_merge(const std::vector<T> &a, const std::vector<T> &b) {
  test_merge_into_thread<T, Device>(a, b);
  test_merge_into_team<T, Device>(a, b);
  test_global_merge<T, Device>(a, b);
}

/* test both merge a/b and merge b/a
 */
template <typename T, typename Device>
void test_merge_commutes(const std::vector<T> &a, const std::vector<T> &b) {
  test_merge<T, Device>(a, b);
  test_merge<T, Device>(b, a);
}

/* define specific and random merge test cases on the device
 */
template <typename T, typename Device>
void test_merge() {
  test_merge_commutes<T, Device>({}, {});
  test_merge_commutes<T, Device>({}, {1});
  test_merge_commutes<T, Device>({1}, {1});
  test_merge_commutes<T, Device>({1}, {2});
  test_merge_commutes<T, Device>({1, 1}, {1, 1});
  test_merge_commutes<T, Device>({1, 2}, {3, 4});
  test_merge_commutes<T, Device>({1, 3}, {2, 4});
  test_merge_commutes<T, Device>({1, 2, 3}, {4, 5, 6});
  test_merge_commutes<T, Device>({1, 2, 3, 5, 6}, {4});
  test_merge_commutes<T, Device>({1, 2, 3, 4, 6}, {5});
  test_merge_commutes<T, Device>({1, 1, 1, 2, 2}, {2});

  // some random cases
  for (auto aLen : {1ul, 10ul, 100ul}) {
    for (auto bLen : {1ul, 10ul, 100ul}) {
      std::vector<T> a, b;
      for (size_t i = 0; i < aLen; ++i) {
        a.push_back(rand() % std::max(aLen, bLen));
      }
      for (size_t i = 0; i < bLen; ++i) {
        b.push_back(rand() % std::max(aLen, bLen));
      }
      std::sort(a.begin(), a.end());
      std::sort(b.begin(), b.end());
      test_merge_commutes<T, Device>(a, b);
    }
  }
}

#define EXECUTE_TEST(ORDINAL, DEVICE)                            \
  TEST_F(TestCategory, graph##_##merge##_##ORDINAL##_##DEVICE) { \
    test_diagonal_search<ORDINAL>();                             \
    test_merge<ORDINAL, DEVICE>();                               \
  }

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

#if (defined(KOKKOSKERNELS_INST_SCALAR_FLOAT)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(float, TestDevice)
#endif

#if (defined(KOKKOSKERNELS_INST_SCALAR_DOUBLE)) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
EXECUTE_TEST(double, TestDevice)
#endif

#undef EXECUTE_TEST
