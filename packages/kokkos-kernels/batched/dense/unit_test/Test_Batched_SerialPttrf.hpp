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
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Pttrf.hpp"
#include "Test_Batched_DenseUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace Pttrf {

template <typename DeviceType, typename DViewType, typename EViewType, typename AlgoTagType>
struct Functor_BatchedSerialPttrf {
  using execution_space = typename DeviceType::execution_space;
  DViewType _d;
  EViewType _e;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialPttrf(const DViewType &d, const EViewType &e) : _d(d), _e(e) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    auto dd = Kokkos::subview(_d, k, Kokkos::ALL());
    auto ee = Kokkos::subview(_e, k, Kokkos::ALL());

    info += KokkosBatched::SerialPttrf<AlgoTagType>::invoke(dd, ee);
  }

  inline int run() {
    using value_type = typename DViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPttrf");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, _d.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename ScalarType, typename AViewType, typename BViewType, typename CViewType,
          typename ArgTransB>
struct Functor_BatchedSerialGemm {
  using execution_space = typename DeviceType::execution_space;
  AViewType _a;
  BViewType _b;
  CViewType _c;
  ScalarType _alpha, _beta;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGemm(const ScalarType alpha, const AViewType &a, const BViewType &b, const ScalarType beta,
                            const CViewType &c)
      : _a(a), _b(b), _c(c), _alpha(alpha), _beta(beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(_c, k, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialGemm<Trans::NoTranspose, ArgTransB, Algo::Gemm::Unblocked>::invoke(_alpha, aa, bb, _beta, cc);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPttrf");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, _a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename ScalarType, typename LayoutType, typename AlgoTagType>
/// \brief Implementation details of batched pttrf test for random matrix
///
/// \param N [in] Batch size of matrix A
/// \param BlkSize [in] Block size of matrix A
void impl_test_batched_pttrf(const int N, const int BlkSize) {
  using ats            = typename Kokkos::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType     = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  View3DType A("A", N, BlkSize, BlkSize), A_reconst("A_reconst", N, BlkSize, BlkSize);
  View3DType EL("EL", N, BlkSize, BlkSize), EU("EU", N, BlkSize, BlkSize), D("D", N, BlkSize, BlkSize),
      LD("LD", N, BlkSize, BlkSize), L("L", N, BlkSize, BlkSize), I("I", N, BlkSize, BlkSize);
  RealView2DType d("d", N, BlkSize),  // Diagonal components
      ones(Kokkos::view_alloc("ones", Kokkos::WithoutInitializing), N, BlkSize);
  View2DType e_upper("e_upper", N, BlkSize - 1), e_lower("e_lower", N,
                                                         BlkSize - 1);  // upper and lower diagonal components

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  RealType realRandStart, realRandEnd;
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, realRandStart, realRandEnd);
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);

  // Add BlkSize to ensure positive definiteness
  Kokkos::fill_random(d, rand_pool, realRandStart + BlkSize, realRandEnd + BlkSize);
  Kokkos::fill_random(e_upper, rand_pool, randStart, randEnd);

  auto h_e_upper = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), e_upper);
  auto h_e_lower = Kokkos::create_mirror_view(e_lower);

  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize - 1; i++) {
      // Fill the lower diagonal with conjugate of the upper diagonal
      h_e_lower(ib, i) = Kokkos::ArithTraits<ScalarType>::conj(h_e_upper(ib, i));
    }
  }

  Kokkos::deep_copy(e_lower, h_e_lower);
  Kokkos::deep_copy(ones, RealType(1.0));

  // Reconstruct Tridiagonal matrix A
  // A = D + EL + EU
  create_diagonal_matrix(e_lower, EL, -1);
  create_diagonal_matrix(e_upper, EU, 1);
  create_diagonal_matrix(d, D);
  create_diagonal_matrix(ones, I);

  // Matrix matrix addition by Gemm
  // D + EU by D * I + EU (result stored in EU)
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, D, I,
                                                                                                            1.0, EU)
      .run();

  // Copy EL to A
  Kokkos::deep_copy(A, EL);

  // EU + EL by EU * I + A (result stored in A)
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, EU, I,
                                                                                                            1.0, A)
      .run();

  // Factorize matrix A -> L * D * L**H
  // d and e are updated by pttrf
  auto info = Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e_lower).run();

  Kokkos::fence();

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  EXPECT_EQ(info, 0);
#endif

  // Reconstruct L and D from factorized matrix
  create_diagonal_matrix(e_lower, EL, -1);
  create_diagonal_matrix(d, D);

  // Copy I to L
  Kokkos::deep_copy(L, I);

  // EL + I by EL * I + L (result stored in L)
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, EL, I,
                                                                                                            1.0, L)
      .run();

  // Reconstruct A by L*D*L**H
  // Gemm to compute L*D -> LD
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, L, D,
                                                                                                            0.0, LD)
      .run();

  // FIXME: We should use SerialGemm Trans::ConjTranspose.
  // For the moment, we compute the complex conjugate of L and
  // then use Trans::Transpose.
  // Gemm to compute (L*D)*L**H -> A_reconst
  // Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType,
  //                          View3DType, Trans::ConjTranspose>(1.0, LD, L, 0.0,
  //                                                            A_reconst)
  //    .run();

  // Compute the complex conjugate of L
  // L -> conj(L)
  auto h_L = Kokkos::create_mirror_view(L);
  Kokkos::deep_copy(h_L, L);
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      for (int j = 0; j < BlkSize; j++) {
        h_L(ib, i, j) = Kokkos::ArithTraits<ScalarType>::conj(h_L(ib, i, j));
      }
    }
  }
  Kokkos::deep_copy(L, h_L);

  // Gemm to compute (L*D)*(conj(L))**T -> A_reconst
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::Transpose>(
      1.0, LD, L, 0.0, A_reconst)
      .run();

  Kokkos::fence();

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  auto h_A         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  auto h_A_reconst = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_reconst);

  // Check A = L*D*L**H
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      for (int j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_A_reconst(ib, i, j), h_A(ib, i, j), eps);
      }
    }
  }
}

template <typename DeviceType, typename ScalarType, typename LayoutType, typename AlgoTagType>
/// \brief Implementation details of batched pttrf test for early return
///        BlkSize must be 0 or 1
///
/// \param N [in] Batch size of matrix A
/// \param BlkSize [in] Block size of matrix A
void impl_test_batched_pttrf_quick_return(const int N, const int BlkSize) {
  using ats            = typename Kokkos::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;

  if (BlkSize > 1) return;

  const int BlkSize_minus_1 = BlkSize > 0 ? BlkSize - 1 : 0;

  RealView2DType d("d", N, BlkSize), d2("d2", N, BlkSize);  // Diagonal components
  View2DType e("e", N,
               BlkSize_minus_1);  // lower diagonal components

  const RealType reference_value = 4.0;

  Kokkos::deep_copy(d, reference_value);
  Kokkos::deep_copy(d2, -reference_value);
  Kokkos::deep_copy(e, ScalarType(1.0));

  // Factorize matrix A -> L * D * L**H
  // d and e are updated by pttrf
  // Early return if BlkSize is 0 or 1
  auto info = Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e).run();

  // For negative values, info should be 1 for BlkSize = 1
  auto info2 = Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d2, e).run();

  Kokkos::fence();

  int expected_info2 = BlkSize == 0 ? 0 : N;
  EXPECT_EQ(info, 0);
  EXPECT_EQ(info2, expected_info2);

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  auto h_d  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d);
  auto h_d2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d2);

  // Check if d is unchanged
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_d(ib, i), reference_value, eps);
      EXPECT_NEAR_KK(h_d2(ib, i), -reference_value, eps);
    }
  }
}

template <typename DeviceType, typename ScalarType, typename LayoutType, typename AlgoTagType>
/// \brief Implementation details of batched pttrf test
///
/// \param N [in] Batch size of matrix A
/// \param BlkSize [in] Block size of matrix A
void impl_test_batched_pttrf_analytical(const int N, const int BlkSize) {
  using ats            = typename Kokkos::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType     = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  View3DType A("A", N, BlkSize, BlkSize), A_reconst("A_reconst", N, BlkSize, BlkSize);
  View3DType EL("EL", N, BlkSize, BlkSize), EU("EU", N, BlkSize, BlkSize), D("D", N, BlkSize, BlkSize),
      LD("LD", N, BlkSize, BlkSize), L("L", N, BlkSize, BlkSize), I("I", N, BlkSize, BlkSize);
  RealView2DType d(Kokkos::view_alloc("d", Kokkos::WithoutInitializing), N,
                   BlkSize),  // Diagonal components
      ones(Kokkos::view_alloc("ones", Kokkos::WithoutInitializing), N, BlkSize);
  View2DType e(Kokkos::view_alloc("e", Kokkos::WithoutInitializing), N,
               BlkSize - 1);  // Upper and lower diagonal components (identical)

  Kokkos::deep_copy(d, RealType(4.0));
  Kokkos::deep_copy(e, ScalarType(1.0));
  Kokkos::deep_copy(ones, RealType(1.0));

  // Reconstruct Tridiaonal matrix A
  // A = D + EL + EU
  create_diagonal_matrix(e, EL, -1);
  create_diagonal_matrix(e, EU, 1);
  create_diagonal_matrix(d, D);
  create_diagonal_matrix(ones, I);

  // Matrix matrix addition by Gemm
  // D + EU by D * I + EU (result stored in EU)
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, D, I,
                                                                                                            1.0, EU)
      .run();

  // Copy EL to A
  Kokkos::deep_copy(A, EL);

  // EU + EL by EU * I + A (result stored in A)
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, EU, I,
                                                                                                            1.0, A)
      .run();

  // Factorize matrix A -> L * D * L**T
  // d and e are updated by pttrf
  auto info = Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e).run();

  Kokkos::fence();

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  EXPECT_EQ(info, 0);
#endif

  // Reconstruct L and D from factorized matrix
  create_diagonal_matrix(e, EL, -1);
  create_diagonal_matrix(d, D);

  // Copy I to L
  Kokkos::deep_copy(L, I);

  // EL + I by EL * I + L (result stored in L)
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, EL, I,
                                                                                                            1.0, L)
      .run();

  // Reconstruct A by L*D*L**T
  // Gemm to compute L*D -> LD
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, L, D,
                                                                                                            0.0, LD)
      .run();

  // Gemm to compute (L*D)*L**T -> A_reconst
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::Transpose>(
      1.0, LD, L, 0.0, A_reconst)
      .run();

  Kokkos::fence();

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  auto h_A         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  auto h_A_reconst = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_reconst);

  // Check A = L*D*L.T
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      for (int j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_A_reconst(ib, i, j), h_A(ib, i, j), eps);
      }
    }
  }
}

}  // namespace Pttrf
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename AlgoTagType>
int test_batched_pttrf() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    for (int i = 0; i < 2; i++) {
      Test::Pttrf::impl_test_batched_pttrf_quick_return<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Pttrf::impl_test_batched_pttrf_quick_return<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
    for (int i = 2; i < 10; i++) {
      Test::Pttrf::impl_test_batched_pttrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Pttrf::impl_test_batched_pttrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
      Test::Pttrf::impl_test_batched_pttrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Pttrf::impl_test_batched_pttrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    for (int i = 0; i < 2; i++) {
      Test::Pttrf::impl_test_batched_pttrf_quick_return<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Pttrf::impl_test_batched_pttrf_quick_return<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
    for (int i = 2; i < 10; i++) {
      Test::Pttrf::impl_test_batched_pttrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Pttrf::impl_test_batched_pttrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
      Test::Pttrf::impl_test_batched_pttrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Pttrf::impl_test_batched_pttrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
  }
#endif

  return 0;
}
