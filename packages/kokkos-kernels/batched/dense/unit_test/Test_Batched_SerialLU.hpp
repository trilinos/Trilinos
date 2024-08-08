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
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

// #include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {

template <typename DeviceType, typename ViewType, typename AlgoTagType>
struct Functor_TestBatchedSerialLU {
  using execution_space = typename DeviceType::execution_space;
  ViewType _a;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedSerialLU(const ViewType &a) : _a(a) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());

    for (int i = 0; i < static_cast<int>(aa.extent(0)); ++i) aa(i, i) += 10.0;

    SerialLU<AlgoTagType>::invoke(aa);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::SerialLU");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, _a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ViewType, typename AlgoTagType>
void impl_test_batched_lu(const int N, const int BlkSize) {
  typedef typename ViewType::value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;

  /// randomized input testing views
  ViewType a0("a0", N, BlkSize, BlkSize), a1("a1", N, BlkSize, BlkSize);

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  Kokkos::fill_random(a0, random, value_type(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(a1, a0);

  Functor_TestBatchedSerialLU<DeviceType, ViewType, Algo::LU::Unblocked>(a0).run();
  Functor_TestBatchedSerialLU<DeviceType, ViewType, AlgoTagType>(a1).run();

  Kokkos::fence();

  /// for comparison send it to host
  typename ViewType::HostMirror a0_host = Kokkos::create_mirror_view(a0);
  typename ViewType::HostMirror a1_host = Kokkos::create_mirror_view(a1);

  Kokkos::deep_copy(a0_host, a0);
  Kokkos::deep_copy(a1_host, a1);

  /// check b0 = b1 ; this eps is about 10^-14
  typedef typename ats::mag_type mag_type;
  mag_type sum(1), diff(0);
  const mag_type eps = 1.0e3 * ats::epsilon();

  for (int k = 0; k < N; ++k)
    for (int i = 0; i < BlkSize; ++i)
      for (int j = 0; j < BlkSize; ++j) {
        sum += ats::abs(a0_host(k, i, j));
        diff += ats::abs(a0_host(k, i, j) - a1_host(k, i, j));
      }
  EXPECT_NEAR_KK(diff / sum, 0, eps);
}
}  // namespace Test

template <typename DeviceType, typename ValueType, typename AlgoTagType>
int test_batched_lu() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType> ViewType;
    Test::impl_test_batched_lu<DeviceType, ViewType, AlgoTagType>(0, 10);
    for (int i = 0; i < 10; ++i) {
      // printf("Testing: LayoutLeft,  Blksize %d\n", i);
      Test::impl_test_batched_lu<DeviceType, ViewType, AlgoTagType>(1024, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType> ViewType;
    Test::impl_test_batched_lu<DeviceType, ViewType, AlgoTagType>(0, 10);
    for (int i = 0; i < 10; ++i) {
      // printf("Testing: LayoutLeft,  Blksize %d\n", i);
      Test::impl_test_batched_lu<DeviceType, ViewType, AlgoTagType>(1024, i);
    }
  }
#endif

  return 0;
}
