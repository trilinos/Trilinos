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

/*! \file

This test file was motivated by an observation in the SpGEMM on SYCL that
strange values were coming out of the pointer alignment functions, causing
Kokkos::atomic_add to be a no-op or write 0. The Kokkos Kernels alignPtrTo
function was updated with the one of four implementations that was observed to
work on SYCL (even though all four in here should be okay.)

TEST_FN 0-3 are various implemetations, and TEST_FN 4 is testing Kokkos Kernels
implementation. The tests are written to PASS for the observed SYCL behavor -
i.e., that TEST_FN 1,4 produce aligned pointers, and the others do not (even
though they should). If the other functions start working on SYCL, then this
test will "fail", and the Kokkos Kernels implementation should be updated with
one of the now-working (and faster?) implementations.
*/

#ifndef TEST_COMMON_ALIGNPTRTO_HPP
#define TEST_COMMON_ALIGNPTRTO_HPP

#include <type_traits>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_Utils.hpp>

namespace {

// the original Kokkos Kernels implementation
template <typename T, typename InPtr>
KOKKOS_INLINE_FUNCTION T *f0(InPtr p) {
  std::uintptr_t ptrVal = reinterpret_cast<std::uintptr_t>(p);
  return reinterpret_cast<T *>((ptrVal + alignof(T) - 1) & (~(alignof(T) - 1)));
}

// an implementation that works for SYCL
template <typename T, typename InPtr>
KOKKOS_INLINE_FUNCTION T *f1(InPtr p) {
  std::uintptr_t ptrVal = reinterpret_cast<std::uintptr_t>(p);
  while (ptrVal % alignof(T)) {
    ++ptrVal;
  }
  return reinterpret_cast<T *>(ptrVal);
}

// another valid implementation
template <typename T, typename InPtr>
KOKKOS_INLINE_FUNCTION T *f2(InPtr p) {
  std::uintptr_t ptrVal = reinterpret_cast<std::uintptr_t>(p);
  return reinterpret_cast<T *>((ptrVal + alignof(T) - 1) / alignof(T) * alignof(T));
}

// the way GCC does it (roughly)
template <typename T, typename InPtr>
KOKKOS_INLINE_FUNCTION T *f3(InPtr p) {
  std::uintptr_t ptrVal = reinterpret_cast<std::uintptr_t>(p);
  return reinterpret_cast<T *>((ptrVal - uint64_t(1) + alignof(T)) & -alignof(T));
}

// Function to be executed by each team
template <int TEST_FN, typename Results>
struct TeamFunction {
  TeamFunction() = default;
  TeamFunction(const Results &results) : results_(results) {}

  template <typename Team>
  KOKKOS_INLINE_FUNCTION void operator()(const Team &team) const {
    // get an "aligned" pointer to scratch memory
    char *shmem = (char *)(team.team_shmem().get_shmem(team.team_size() * sizeof(double)));
    double *vals;
    if constexpr (0 == TEST_FN) {
      vals = f0<double>(shmem);
    } else if constexpr (1 == TEST_FN) {
      vals = f1<double>(shmem);
    } else if constexpr (2 == TEST_FN) {
      vals = f2<double>(shmem);
    } else if constexpr (3 == TEST_FN) {
      vals = f3<double>(shmem);
    } else if constexpr (4 == TEST_FN) {
      vals = KokkosKernels::Impl::alignPtrTo<double>(shmem);
    } else {
      static_assert(std::is_void_v<Results>, "Unexpected test function");
    }

    const size_t i = team.team_rank();
    double val     = team.team_rank();
    vals[i]        = 0;  // zero shared memory
    Kokkos::atomic_add(&vals[i], val);
#if 0  // debugging
    Kokkos::printf("%s:%i result(%lu) += %f yielded %f\n", __FILE__, __LINE__, i, val, vals[i]);
#endif

    results_(i) = vals[i];
  }

  size_t team_shmem_size(int team_size) const { return team_size * sizeof(double); }

  Results results_;
};

// use atomic add to set result(i) = i
template <int TEST_FN, typename Device>
void test_alignPtrTo() {
  using MemorySpace = typename Device::memory_space;
  using ExecSpace   = typename Device::execution_space;
  using TestView    = Kokkos::View<double *, MemorySpace>;
  using TestPolicy  = Kokkos::TeamPolicy<ExecSpace>;
  const int teamSize =
      TestPolicy(1, Kokkos::AUTO).team_size_max(TeamFunction<TEST_FN, TestView>(), Kokkos::ParallelForTag{});

  ExecSpace space;

  TestView results("TestView", teamSize);
  TestPolicy policy(space, 1, teamSize);
  Kokkos::parallel_for("test alignment", policy, TeamFunction<TEST_FN, TestView>(results));

  int errs;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecSpace>(space, 0, teamSize),
      KOKKOS_LAMBDA(int i, int &lerr) { lerr += (results(i) != i); }, errs);

// if SYCL is enabled, only TEST_FN 1 and 4 should work
#if defined(KOKKOS_ENABLE_SYCL)
  if constexpr (std::is_same_v<ExecSpace, Kokkos::Experimental::SYCL>) {
    if constexpr ((1 == TEST_FN) || (4 == TEST_FN)) {
      EXPECT_EQ(0, errs);
    } else {
      EXPECT_NE(0, errs);
    }
  } else {
    EXPECT_EQ(0, errs);
  }
#else
  EXPECT_EQ(0, errs);
#endif
}

TEST_F(TestCategory, common_AlignPtrTo_0) { test_alignPtrTo<0, TestDevice>(); }
TEST_F(TestCategory, common_AlignPtrTo_1) { test_alignPtrTo<1, TestDevice>(); }
TEST_F(TestCategory, common_AlignPtrTo_2) { test_alignPtrTo<2, TestDevice>(); }
TEST_F(TestCategory, common_AlignPtrTo_3) { test_alignPtrTo<3, TestDevice>(); }
TEST_F(TestCategory, common_AlignPtrTo_kk) { test_alignPtrTo<4, TestDevice>(); }

}  // anonymous namespace

#endif  // TEST_COMMON_ALIGNPTRTO
