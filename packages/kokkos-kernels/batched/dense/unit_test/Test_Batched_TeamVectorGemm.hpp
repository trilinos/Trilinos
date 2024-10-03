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
#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_TeamVector_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace TeamVectorGemm {

template <typename TA, typename TB>
struct ParamTag {
  typedef TA transA;
  typedef TB transB;
};

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
struct Functor_TestBatchedTeamVector {
  using execution_space = typename DeviceType::execution_space;
  ViewType _a, _b, _c;

  ScalarType _alpha, _beta;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamVector(const ScalarType alpha, const ViewType &a, const ViewType &b, const ScalarType beta,
                                const ViewType &c)
      : _a(a), _b(b), _c(c), _alpha(alpha), _beta(beta) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const ParamTagType &, const MemberType &member) const {
    const int k = member.league_rank();

    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(_c, k, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::TeamVectorGemm<MemberType, typename ParamTagType::transA, typename ParamTagType::transB,
                                  AlgoTagType>::invoke(member, _alpha, aa, bb, _beta, cc);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::TeamVector");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    const int league_size = _c.extent(0);
    Kokkos::TeamPolicy<execution_space, ParamTagType> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_teamvectorgemm(const int N, const int matAdim1, const int matAdim2, const int matBdim1,
                                      const int matBdim2, const int matCdim1, const int matCdim2) {
  using transA          = typename ParamTagType::transA;
  using transB          = typename ParamTagType::transB;
  using execution_space = typename DeviceType::execution_space;
  using value_type      = typename ViewType::value_type;
  using ats             = Kokkos::ArithTraits<value_type>;

  /// randomized input testing views
  ScalarType alpha = ScalarType(1.5), beta = ScalarType(3.0);

  ViewType a_expected("a_expected", N, matAdim1, matAdim2), a_actual("a_actual", N, matAdim1, matAdim2),
      b_expected("b_expected", N, matBdim1, matBdim2), b_actual("b_actual", N, matBdim1, matBdim2),
      c_expected("c_expected", N, matCdim1, matCdim2), c_actual("c_actual", N, matCdim1, matCdim2);

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);

  Kokkos::fill_random(a_expected, random, value_type(1.0));
  Kokkos::fill_random(b_expected, random, value_type(1.0));
  Kokkos::fill_random(c_expected, random, value_type(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(a_actual, a_expected);
  Kokkos::deep_copy(b_actual, b_expected);
  Kokkos::deep_copy(c_actual, c_expected);

  // Functor_TestBatchedTeamVector<DeviceType,ViewType,ScalarType,
  //   ParamTagType,Algo::Gemm::Unblocked>(alpha, a_expected, b_expected,
  //   beta, c_expected).run();
  Functor_BatchedVanillaGEMM<ViewType, ViewType, ViewType, execution_space> vgemm;
  vgemm.A_t = std::is_same<transA, Trans::Transpose>::value;
  vgemm.B_t = std::is_same<transB, Trans::Transpose>::value;
  vgemm.A_c = vgemm.B_c = false;
  vgemm.A               = a_expected;
  vgemm.B               = b_expected;
  vgemm.C               = c_expected;
  vgemm.alpha           = alpha;
  vgemm.beta            = beta;
  vgemm.run();  // Compute c_expected

  Functor_TestBatchedTeamVector<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(alpha, a_actual, b_actual,
                                                                                             beta, c_actual)
      .run();

  Kokkos::fence();

  typename ViewType::HostMirror c_expected_host = Kokkos::create_mirror_view(c_expected);
  typename ViewType::HostMirror c_actual_host   = Kokkos::create_mirror_view(c_actual);

  // Copy to host for comparison
  Kokkos::deep_copy(c_expected_host, c_expected);
  Kokkos::deep_copy(c_actual_host, c_actual);

  using mag_type = typename ats::mag_type;
  mag_type sum(1), diff(0);

  mag_type eps = ats::epsilon();

  eps *= std::is_same<value_type, Kokkos::Experimental::half_t>::value ||
                 std::is_same<value_type, Kokkos::Experimental::bhalf_t>::value
             ? 4
             : 1e3;

  for (int k = 0; k < N; ++k)
    for (int i = 0; i < matCdim1; ++i)
      for (int j = 0; j < matCdim2; ++j) {
        sum += ats::abs(c_expected_host(k, i, j));
        diff += ats::abs(c_expected_host(k, i, j) - c_actual_host(k, i, j));
      }
  EXPECT_NEAR_KK(diff / sum, 0, eps);
}
}  // namespace TeamVectorGemm
}  // namespace Test

// void (*impl_test)(const int, const int, const int, const int, const int,
// const int, const int)
template <typename DeviceType, typename ValueType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_teamvectorgemm() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutLeft, DeviceType> ViewType;
    Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(
        0, 10, 10, 10, 10, 10, 10);
    for (int i = 0; i < 10; ++i) {
      // printf("Testing: LayoutLeft,  Blksize %d\n", i);
      Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType,
                                                             AlgoTagType>(1024, i, i, i, i, i, i);
    }
    for (int i = 0; i < 10; ++i) {
      // printf("Testing: LayoutLeft,  Blksize %d\n", i);
      int dimM = i;
      int dimN = 2 * i;
      int dimK = 3 * i;
      if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>::value) &&
          (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>::value)) {
        Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType,
                                                               AlgoTagType>(1024, dimM, dimK, dimK, dimN, dimM, dimN);
      }
      if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>::value) &&
          (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::Transpose>::value)) {
        Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType,
                                                               AlgoTagType>(1024, dimM, dimK, dimN, dimK, dimM, dimN);
      }
      if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::Transpose>::value) &&
          (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>::value)) {
        Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType,
                                                               AlgoTagType>(1024, dimK, dimM, dimK, dimN, dimM, dimN);
      }
      if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::Transpose>::value) &&
          (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::Transpose>::value)) {
        Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType,
                                                               AlgoTagType>(1024, dimK, dimM, dimN, dimK, dimM, dimN);
      }
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    typedef Kokkos::View<ValueType ***, Kokkos::LayoutRight, DeviceType> ViewType;
    Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(
        0, 10, 10, 10, 10, 10, 10);
    for (int i = 0; i < 10; ++i) {
      // printf("Testing: LayoutRight, Blksize %d\n", i);
      Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType,
                                                             AlgoTagType>(1024, i, i, i, i, i, i);
    }
    for (int i = 0; i < 10; ++i) {
      // printf("Testing: LayoutLeft,  Blksize %d\n", i);
      int dimM = i;
      int dimN = 2 * i;
      int dimK = 3 * i;
      if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>::value) &&
          (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>::value)) {
        Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType,
                                                               AlgoTagType>(1024, dimM, dimK, dimK, dimN, dimM, dimN);
      }
      if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>::value) &&
          (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::Transpose>::value)) {
        Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType,
                                                               AlgoTagType>(1024, dimM, dimK, dimN, dimK, dimM, dimN);
      }
      if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::Transpose>::value) &&
          (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>::value)) {
        Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType,
                                                               AlgoTagType>(1024, dimK, dimM, dimK, dimN, dimM, dimN);
      }
      if ((std::is_same<typename ParamTagType::transA, KokkosBatched::Trans::Transpose>::value) &&
          (std::is_same<typename ParamTagType::transB, KokkosBatched::Trans::Transpose>::value)) {
        Test::TeamVectorGemm::impl_test_batched_teamvectorgemm<DeviceType, ViewType, ScalarType, ParamTagType,
                                                               AlgoTagType>(1024, dimK, dimM, dimN, dimK, dimM, dimN);
      }
    }
  }
#endif

  return 0;
}
