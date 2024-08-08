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
#include "KokkosBatched_SolveLU_Decl.hpp"
// #include "KokkosBatched_SolveLU_Serial_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace SerialSolveLU {

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
    std::string name_region("KokkosBatched::Test::SerialSolveLU");
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
    std::string name_region("KokkosBatched::Test::SerialSolveLU");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, _a.extent(0));
    Kokkos::parallel_for((name + "::LUFunctor").c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ViewType, typename TransType, typename AlgoTagType>
struct Functor_TestBatchedSerialSolveLU {
  using execution_space = typename DeviceType::execution_space;
  ViewType _a;
  ViewType _b;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedSerialSolveLU(const ViewType &a, const ViewType &b) : _a(a), _b(b) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialSolveLU<TransType, AlgoTagType>::invoke(aa, bb);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::SerialSolveLU");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, _a.extent(0));
    Kokkos::parallel_for((name + "::SolveLUFunctor").c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ViewType, typename AlgoTagType>
void impl_test_batched_solvelu(const int N, const int BlkSize) {
  typedef typename ViewType::value_type value_type;
  typedef Kokkos::ArithTraits<value_type> ats;

  /// randomized input testing views
  ViewType a0("a0", N, BlkSize, BlkSize);
  ViewType a1("a1", N, BlkSize, BlkSize);
  ViewType b("b", N, BlkSize, 5);
  ViewType x0("x0", N, BlkSize, 5);
  // ViewType a0_T("a0_T", N, BlkSize, BlkSize);
  // ViewType b_T ("b_T",  N, BlkSize, 5 );

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  Kokkos::fill_random(a0, random, value_type(1.0));
  Kokkos::fill_random(x0, random, value_type(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(a1, a0);
  // Kokkos::deep_copy(a0_T, a0);

  value_type alpha = 1.0, beta = 0.0;
  typedef ParamTag<Trans::NoTranspose, Trans::NoTranspose> param_tag_type;

  Functor_BatchedSerialGemm<DeviceType, ViewType, value_type, param_tag_type, AlgoTagType>(alpha, a0, x0, beta, b)
      .run();

  Functor_BatchedSerialLU<DeviceType, ViewType, AlgoTagType>(a1).run();

  Functor_TestBatchedSerialSolveLU<DeviceType, ViewType, Trans::NoTranspose, AlgoTagType>(a1, b).run();

  Kokkos::fence();

  // //Transpose
  // typedef ParamTag<Trans::Transpose,Trans::NoTranspose> param_tag_type_T;

  // Functor_BatchedSerialGemm<DeviceType,ViewType,value_type,
  //   param_tag_type_T,AlgoTagType>(alpha, a0_T, x0, beta, b_T).run();

  // Functor_TestBatchedSerialSolveLU<DeviceType,ViewType,Trans::Transpose,AlgoTagType>(a1,b_T).run();

  // Kokkos::fence();

  /// for comparison send it to host
  typename ViewType::HostMirror x0_host = Kokkos::create_mirror_view(x0);
  typename ViewType::HostMirror b_host  = Kokkos::create_mirror_view(b);
  // typename ViewType::HostMirror b_T_host = Kokkos::create_mirror_view(b_T);

  Kokkos::deep_copy(x0_host, x0);
  Kokkos::deep_copy(b_host, b);
  // Kokkos::deep_copy(b_T_host, b_T);

  /// check x0 = b ; this eps is about 10^-14
  typedef typename ats::mag_type mag_type;
  mag_type sum(1), diff(0);
  const mag_type eps = 1.0e3 * ats::epsilon();

  for (int k = 0; k < N; ++k)
    for (int i = 0; i < BlkSize; ++i)
      for (int j = 0; j < 5; ++j) {
        sum += ats::abs(x0_host(k, i, j));
        diff += ats::abs(x0_host(k, i, j) - b_host(k, i, j));
      }
  // printf("NoTranspose -- N=%d, BlkSize=%d, sum=%f, diff=%f\n", N, BlkSize,
  // sum, diff);
  EXPECT_NEAR_KK(diff / sum, 0.0, eps);

  /// check x0 = b_T ; this eps is about 10^-14
  // mag_type sum_T(1), diff_T(0);

  // for (int k=0;k<N;++k)
  //   for (int i=0;i<BlkSize;++i)
  //     for (int j=0;j<5;++j) {
  //       sum_T  += ats::abs(x0_host(k,i,j));
  //       diff_T += ats::abs(x0_host(k,i,j)-b_T_host(k,i,j));
  //     }
  // //printf("Transpose -- N=%d, BlkSize=%d, sum=%f, diff=%f\n", N, BlkSize,
  // sum_T, diff_T); EXPECT_NEAR_KK( diff_T/sum_T, 0.0, eps);
}
}  // namespace SerialSolveLU
}  // namespace Test

template <typename DeviceType, typename ValueType, typename AlgoTagType>
int test_batched_solvelu() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType> ViewType;
    Test::SerialSolveLU::impl_test_batched_solvelu<DeviceType, ViewType, AlgoTagType>(0, 10);
    for (int i = 0; i < 10; ++i) {
      Test::SerialSolveLU::impl_test_batched_solvelu<DeviceType, ViewType, AlgoTagType>(1024, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType> ViewType;
    Test::SerialSolveLU::impl_test_batched_solvelu<DeviceType, ViewType, AlgoTagType>(0, 10);
    for (int i = 0; i < 10; ++i) {
      Test::SerialSolveLU::impl_test_batched_solvelu<DeviceType, ViewType, AlgoTagType>(1024, i);
    }
  }
#endif

  return 0;
}
