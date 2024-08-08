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
/// \author Vinh Dang (vqdang@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

// #include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"
#include "KokkosBatched_InverseLU_Decl.hpp"
// #include "KokkosBatched_InverseLU_Serial_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace SerialInverseLU {

template <typename TA, typename TB>
struct ParamTag {
  typedef TA transA;
  typedef TB transB;
};

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
struct Functor_BatchedSerialGemm {
  using execution_space = typename DeviceType::execution_space;
  ViewType _a, _b, _c;

  ScalarType _alpha, _beta;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGemm(const ScalarType alpha, const ViewType &a, const ViewType &b, const ScalarType beta,
                            const ViewType &c)
      : _a(a), _b(b), _c(c), _alpha(alpha), _beta(beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(_c, k, Kokkos::ALL(), Kokkos::ALL());

    for (int i = 0; i < static_cast<int>(aa.extent(0)); ++i) aa(i, i) += 10.0;

    SerialGemm<typename ParamTagType::transA, typename ParamTagType::transB, AlgoTagType>::invoke(_alpha, aa, bb, _beta,
                                                                                                  cc);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::SerialInverseLU");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, _c.extent(0));
    Kokkos::parallel_for((name + "::GemmFunctor").c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ViewType, typename AlgoTagType>
struct Functor_BatchedSerialLU {
  using execution_space = typename DeviceType::execution_space;
  ViewType _a;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialLU(const ViewType &a) : _a(a) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());

    for (int i = 0; i < static_cast<int>(aa.extent(0)); ++i) aa(i, i) += 10.0;

    SerialLU<AlgoTagType>::invoke(aa);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::SerialInverseLU");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, _a.extent(0));
    Kokkos::parallel_for((name + "::LUFunctor").c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename AViewType, typename WViewType, typename AlgoTagType>
struct Functor_TestBatchedSerialInverseLU {
  using execution_space = typename DeviceType::execution_space;
  AViewType _a;
  WViewType _w;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedSerialInverseLU(const AViewType &a, const WViewType &w) : _a(a), _w(w) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto ww = Kokkos::subview(_w, k, Kokkos::ALL());

    KokkosBatched::SerialInverseLU<AlgoTagType>::invoke(aa, ww);
  }

  inline void run() {
    typedef typename AViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::SerialInverseLU");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, _a.extent(0));
    Kokkos::parallel_for((name + "::InverseLUFunctor").c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename AViewType, typename WViewType, typename AlgoTagType>
void impl_test_batched_inverselu(const int N, const int BlkSize) {
  typedef typename AViewType::value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;

  /// randomized input testing views
  AViewType a0("a0", N, BlkSize, BlkSize);
  AViewType a1("a1", N, BlkSize, BlkSize);
  WViewType w("w", N, BlkSize * BlkSize);
  AViewType c0("c0", N, BlkSize, BlkSize);

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  Kokkos::fill_random(a0, random, value_type(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(a1, a0);
  Kokkos::deep_copy(w, value_type(0.0));

  Functor_BatchedSerialLU<DeviceType, AViewType, AlgoTagType>(a1).run();

  Functor_TestBatchedSerialInverseLU<DeviceType, AViewType, WViewType, AlgoTagType>(a1, w).run();

  value_type alpha = 1.0, beta = 0.0;
  typedef SerialInverseLU::ParamTag<Trans::NoTranspose, Trans::NoTranspose> param_tag_type;

  Functor_BatchedSerialGemm<DeviceType, AViewType, value_type, param_tag_type, AlgoTagType>(alpha, a0, a1, beta, c0)
      .run();

  Kokkos::fence();

  /// for comparison send it to host
  typename AViewType::HostMirror c0_host = Kokkos::create_mirror_view(c0);

  Kokkos::deep_copy(c0_host, c0);

  /// check identity matrix ; this eps is about 10^-14
  typedef typename ats::mag_type mag_type;
  mag_type sum_diag(0), sum_all(0), sum_diag_ref(N * BlkSize);
  const mag_type eps = 1.0e3 * ats::epsilon();

  for (int k = 0; k < N; ++k)
    for (int i = 0; i < BlkSize; ++i)
      for (int j = 0; j < BlkSize; ++j) {
        sum_all += ats::abs(c0_host(k, i, j));
        if (i == j) sum_diag += ats::abs(c0_host(k, i, j));
      }
  EXPECT_NEAR_KK(sum_all - sum_diag, 0, eps);
  EXPECT_NEAR_KK(sum_diag - sum_diag_ref, 0, eps);
}
}  // namespace SerialInverseLU
}  // namespace Test

template <typename DeviceType, typename ValueType, typename AlgoTagType>
int test_batched_inverselu() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType> AViewType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> WViewType;
    Test::SerialInverseLU::impl_test_batched_inverselu<DeviceType, AViewType, WViewType, AlgoTagType>(0, 10);
    for (int i = 0; i < 10; ++i) {
      Test::SerialInverseLU::impl_test_batched_inverselu<DeviceType, AViewType, WViewType, AlgoTagType>(1024, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType> AViewType;
    typedef Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType> WViewType;
    Test::SerialInverseLU::impl_test_batched_inverselu<DeviceType, AViewType, WViewType, AlgoTagType>(0, 10);
    for (int i = 0; i < 10; ++i) {
      Test::SerialInverseLU::impl_test_batched_inverselu<DeviceType, AViewType, WViewType, AlgoTagType>(1024, i);
    }
  }
#endif

  return 0;
}
