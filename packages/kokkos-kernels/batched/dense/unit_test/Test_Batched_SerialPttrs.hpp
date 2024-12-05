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
#include "KokkosBatched_Pttrs.hpp"
#include "Test_Batched_DenseUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace Pttrs {

template <typename U>
struct ParamTag {
  using uplo = U;
};

template <typename DeviceType, typename DViewType, typename EViewType, typename AlgoTagType>
struct Functor_BatchedSerialPttrf {
  using execution_space = typename DeviceType::execution_space;
  DViewType _d;
  EViewType _e;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialPttrf(const DViewType &d, const EViewType &e) : _d(d), _e(e) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto dd = Kokkos::subview(_d, k, Kokkos::ALL());
    auto ee = Kokkos::subview(_e, k, Kokkos::ALL());

    KokkosBatched::SerialPttrf<AlgoTagType>::invoke(dd, ee);
  }

  inline void run() {
    using value_type = typename EViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPttrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, _d.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename DViewType, typename EViewType, typename BViewType, typename ParamTagType,
          typename AlgoTagType>
struct Functor_BatchedSerialPttrs {
  using execution_space = typename DeviceType::execution_space;
  DViewType _d;
  EViewType _e;
  BViewType _b;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialPttrs(const DViewType &d, const EViewType &e, const BViewType &b) : _d(d), _e(e), _b(b) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k, int &info) const {
    auto dd = Kokkos::subview(_d, k, Kokkos::ALL());
    auto ee = Kokkos::subview(_e, k, Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL());

    info += KokkosBatched::SerialPttrs<typename ParamTagType::uplo, AlgoTagType>::invoke(dd, ee, bb);
  }

  inline int run() {
    using value_type = typename BViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPttrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, _d.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
struct Functor_BatchedSerialGemv {
  using execution_space = typename DeviceType::execution_space;
  AViewType _a;
  xViewType _x;
  yViewType _y;
  ScalarType _alpha, _beta;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGemv(const ScalarType alpha, const AViewType &a, const xViewType &x, const ScalarType beta,
                            const yViewType &y)
      : _a(a), _x(x), _y(y), _alpha(alpha), _beta(beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto xx = Kokkos::subview(_x, k, Kokkos::ALL());
    auto yy = Kokkos::subview(_y, k, Kokkos::ALL());

    KokkosBlas::SerialGemv<Trans::NoTranspose, Algo::Gemv::Unblocked>::invoke(_alpha, aa, xx, _beta, yy);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPttrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, _a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
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
    std::string name_region("KokkosBatched::Test::SerialPttrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, _a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
/// \brief Implementation details of batched pttrs test
///        Confirm A * x = b, where
///        A: [[4, 1],
///            [1, 4]]
///        b: [1, 1]
///        x: [1/5, 1/5]
///
///        This corresponds to the following system of equations:
///        4 x0 +   x1 = 1
///          x0 + 4 x1 = 1
///
/// \param N [in] Batch size of RHS (banded matrix can also be batched matrix)
/// \param k [in] Number of superdiagonals or subdiagonals of matrix A
/// \param BlkSize [in] Block size of matrix A
void impl_test_batched_pttrs_analytical(const int N) {
  using ats            = typename Kokkos::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;

  constexpr int BlkSize = 2;
  RealView2DType d(Kokkos::view_alloc("d", Kokkos::WithoutInitializing), N,
                   BlkSize);  // Diagonal components
  View2DType e(Kokkos::view_alloc("e", Kokkos::WithoutInitializing), N,
               BlkSize - 1);  // Upper and lower diagonal components (identical)
  View2DType x(Kokkos::view_alloc("x", Kokkos::WithoutInitializing), N,
               BlkSize);  // Solutions

  Kokkos::deep_copy(d, RealType(4.0));
  Kokkos::deep_copy(e, ScalarType(1.0));
  Kokkos::deep_copy(x, ScalarType(1.0));  // This initialy stores b

  // Factorize matrix A -> L * D * L**H
  // d and e are updated by pttrf
  Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e).run();

  // pttrs (Note, d and e must be factorized by pttrf)
  auto info =
      Functor_BatchedSerialPttrs<DeviceType, RealView2DType, View2DType, View2DType, ParamTagType, AlgoTagType>(d, e, x)
          .run();

  Kokkos::fence();

  EXPECT_EQ(info, 0);

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  auto h_x = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);

  // Check x = [1/5, 1/5]
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), ScalarType(1.0 / 5.0), eps);
    }
  }
}

template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
/// \brief Implementation details of batched pttrs test
///
/// \param N [in] Batch size of RHS (banded matrix can also be batched matrix)
/// \param k [in] Number of superdiagonals or subdiagonals of matrix A
/// \param BlkSize [in] Block size of matrix A
void impl_test_batched_pttrs(const int N, const int BlkSize) {
  using ats            = typename Kokkos::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType     = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  View3DType A("A", N, BlkSize, BlkSize), EL("EL", N, BlkSize, BlkSize), EU("EU", N, BlkSize, BlkSize),
      D("D", N, BlkSize, BlkSize), I("I", N, BlkSize, BlkSize);
  RealView2DType d(Kokkos::view_alloc("d", Kokkos::WithoutInitializing), N,
                   BlkSize),  // Diagonal components
      ones(Kokkos::view_alloc("ones", Kokkos::WithoutInitializing), N, BlkSize);
  View2DType e_upper("e_upper", N, BlkSize - 1), e_lower("e_lower", N,
                                                         BlkSize - 1);            // upper and lower diagonal components
  View2DType x("x", N, BlkSize), b_ref("x_ref", N, BlkSize), b("b", N, BlkSize);  // Solutions

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  RealType realRandStart, realRandEnd;
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, realRandStart, realRandEnd);
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);

  // Add BlkSize to ensure positive definiteness
  Kokkos::fill_random(d, rand_pool, realRandStart + BlkSize, realRandEnd + BlkSize);
  Kokkos::fill_random(e_upper, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);

  auto h_e_upper = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), e_upper);
  auto h_e_lower = Kokkos::create_mirror_view(e_lower);

  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize - 1; i++) {
      // Fill the lower diagonal with conjugate of the upper diagonal
      h_e_lower(ib, i) = Kokkos::ArithTraits<ScalarType>::conj(h_e_upper(ib, i));
    }
  }

  Kokkos::deep_copy(e_lower, h_e_lower);
  Kokkos::deep_copy(b_ref, x);
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

  Kokkos::fence();

  int info = 0;
  if (std::is_same_v<typename ParamTagType::uplo, KokkosBatched::Uplo::Upper>) {
    // Factorize matrix A -> U**H * D * U
    // d and e are updated by pttrf
    Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e_upper).run();

    // pttrs (Note, d and e must be factorized by pttrf)
    info = Functor_BatchedSerialPttrs<DeviceType, RealView2DType, View2DType, View2DType, ParamTagType, AlgoTagType>(
               d, e_upper, x)
               .run();
  } else {
    // Factorize matrix A -> L * D * L**H
    // d and e are updated by pttrf
    Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e_lower).run();

    // pttrs (Note, d and e must be factorized by pttrf)
    info = Functor_BatchedSerialPttrs<DeviceType, RealView2DType, View2DType, View2DType, ParamTagType, AlgoTagType>(
               d, e_lower, x)
               .run();
  }

  Kokkos::fence();

  EXPECT_EQ(info, 0);

  // Gemv to compute b = A * x, this should be identical to b_ref
  Functor_BatchedSerialGemv<DeviceType, ScalarType, View3DType, View2DType, View2DType>(1.0, A, x, 0.0, b).run();

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  auto h_b     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b);
  auto h_b_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b_ref);

  // Check A * x = b_ref
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_b(ib, i), h_b_ref(ib, i), eps);
    }
  }
}

template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
/// \brief Implementation details of batched pttrs test for early return
///        BlkSize must be 0 or 1
///
/// \param N [in] Batch size of RHS (banded matrix can also be batched matrix)
/// \param k [in] Number of superdiagonals or subdiagonals of matrix A
/// \param BlkSize [in] Block size of matrix A
void impl_test_batched_pttrs_quick_return(const int N, const int BlkSize) {
  using ats            = typename Kokkos::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;

  if (BlkSize > 1) return;

  const int BlkSize_minus_1 = BlkSize > 0 ? BlkSize - 1 : 0;

  RealView2DType d(Kokkos::view_alloc("d", Kokkos::WithoutInitializing), N,
                   BlkSize);  // Diagonal components
  View2DType e(Kokkos::view_alloc("e", Kokkos::WithoutInitializing), N,
               BlkSize_minus_1);  // lower diagonal components
  View2DType x("x", N, BlkSize);  // Solutions

  const RealType reference_value = 4.0;

  Kokkos::deep_copy(d, reference_value);
  Kokkos::deep_copy(e, ScalarType(1.0));
  Kokkos::deep_copy(x, ScalarType(1.0));

  // Factorize matrix A -> U**H * D * U or L * D * L**H
  // d and e are updated by pttrf
  Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e).run();

  // pttrs (Note, d and e must be factorized by pttrf)
  auto info =
      Functor_BatchedSerialPttrs<DeviceType, RealView2DType, View2DType, View2DType, ParamTagType, AlgoTagType>(d, e, x)
          .run();

  Kokkos::fence();

  EXPECT_EQ(info, 0);

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  auto h_x = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);

  // Check x = x_ref
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), ScalarType(1.0 / reference_value), eps);
    }
  }
}

}  // namespace Pttrs
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_pttrs() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Pttrs::impl_test_batched_pttrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Pttrs::impl_test_batched_pttrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 2; i++) {
      Test::Pttrs::impl_test_batched_pttrs_quick_return<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
          1, i);
      Test::Pttrs::impl_test_batched_pttrs_quick_return<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
          2, i);
    }
    for (int i = 2; i < 10; i++) {
      Test::Pttrs::impl_test_batched_pttrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, i);
      Test::Pttrs::impl_test_batched_pttrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Pttrs::impl_test_batched_pttrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Pttrs::impl_test_batched_pttrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 2; i++) {
      Test::Pttrs::impl_test_batched_pttrs_quick_return<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
          1, i);
      Test::Pttrs::impl_test_batched_pttrs_quick_return<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
          2, i);
    }
    for (int i = 2; i < 10; i++) {
      Test::Pttrs::impl_test_batched_pttrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, i);
      Test::Pttrs::impl_test_batched_pttrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, i);
    }
  }
#endif

  return 0;
}
