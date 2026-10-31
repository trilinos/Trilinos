// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <concepts>
#include <type_traits>
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Syr2k.hpp>

#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Syr2k {

template <typename U, typename T, typename M>
struct ParamTag {
  using uplo  = U;
  using trans = T;
  using mode  = M;
};

template <typename DeviceType, typename ParamTagType, typename AlphaType, typename AViewType, typename BViewType,
          typename BetaType, typename CViewType>
struct Functor_BatchedSyr2k {
  using execution_space = typename DeviceType::execution_space;
  using member_type     = typename Kokkos::TeamPolicy<execution_space>::member_type;
  using ArgMode         = typename ParamTagType::mode;
  using ArgUplo         = typename ParamTagType::uplo;
  using ArgTrans        = typename ParamTagType::trans;
  AViewType m_A;
  BViewType m_B;
  CViewType m_C;
  AlphaType m_alpha;
  BetaType m_beta;

  Functor_BatchedSyr2k(const AlphaType alpha, const AViewType &A, const BViewType &B, const BetaType beta,
                       const CViewType &C)
      : m_A(A), m_B(B), m_C(C), m_alpha(alpha), m_beta(beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type &member, int &info) const {
    const int k = member.league_rank();
    auto sub_A  = Kokkos::subview(m_A, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_B  = Kokkos::subview(m_B, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_C  = Kokkos::subview(m_C, k, Kokkos::ALL(), Kokkos::ALL());

    if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        info += KokkosBatched::SerialSyr2k<ArgUplo, ArgTrans>::invoke(m_alpha, sub_A, sub_B, m_beta, sub_C);
      });
    } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Team>) {
      info += KokkosBatched::TeamSyr2k<member_type, ArgUplo, ArgTrans>::invoke(member, m_alpha, sub_A, sub_B, m_beta,
                                                                               sub_C);
    } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::TeamVector>) {
      info += KokkosBatched::TeamVectorSyr2k<member_type, ArgUplo, ArgTrans>::invoke(member, m_alpha, sub_A, sub_B,
                                                                                     m_beta, sub_C);
    }
  }

  inline int run() {
    using value_type        = typename AViewType::non_const_value_type;
    std::string name_region = std::same_as<ArgMode, KokkosBatched::Mode::Serial> ? "KokkosBatched::Test::SerialSyr2k"
                              : std::same_as<ArgMode, KokkosBatched::Mode::Team>
                                  ? "KokkosBatched::Test::TeamSyr2k"
                                  : "KokkosBatched::Test::TeamVectorSyr2k";
    const std::string name_value_type = TestUtils::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    const int league_size = m_A.extent_int(0);

    Kokkos::TeamPolicy<execution_space> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename ArgMode, typename ArgTransA, typename ArgTransB, typename ScalarType,
          typename AViewType, typename BViewType, typename CViewType>
struct Functor_BatchedGemm {
  using execution_space = typename DeviceType::execution_space;

  AViewType m_a;
  BViewType m_b;
  CViewType m_c;
  ScalarType m_alpha;
  ScalarType m_beta;

  Functor_BatchedGemm(const ScalarType alpha, const AViewType &a, const BViewType &b, const ScalarType beta,
                      const CViewType &c)
      : m_a(a), m_b(b), m_c(c), m_alpha(alpha), m_beta(beta) {}

  KOKKOS_INLINE_FUNCTION void operator()(const int k) const {
    auto aa = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(m_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(m_c, k, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialGemm<ArgTransA, ArgTransB, KokkosBatched::Algo::Gemm::Unblocked>::invoke(m_alpha, aa, bb,
                                                                                                  m_beta, cc);
  }

  inline void run() {
    using value_type        = typename CViewType::non_const_value_type;
    std::string name_region = std::same_as<ArgMode, KokkosBatched::Mode::Serial> ? "KokkosBatched::Test::SerialSyr2k"
                              : std::same_as<ArgMode, KokkosBatched::Mode::Team>
                                  ? "KokkosBatched::Test::TeamSyr2k"
                                  : "KokkosBatched::Test::TeamVectorSyr2k";
    const std::string name_value_type = TestUtils::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());

    Kokkos::RangePolicy<execution_space> policy(0, m_b.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

/// \brief Implementation details of batched syr2k analytical test
///        to confirm C:= alpha*A*B^T + alpha*B*A^T + beta*C or C:= alpha*A^T*B + alpha*B^T*A + beta*C is computed
///        correctly A and C are deliberately chosen to have non symmetric values to confirm both cases. alpha = 1.5,
///        beta = 1.2
///        4x4 matrix (upper)
///        A: [[1,  -3, -2,  0],
///            [3,  3, -1, -2],
///            [2,  1,  9, 5],
///            [0,  2, -5, 27]]
///        B: [[2, 1, 0, -1],
///            [-1, 4, 2, 0],
///            [3, -2, 1, 5],
///            [1, 0, -3, 2]]
///        C = A
///
///        ArgTrans == NoTranspose/ConjNoTranspose
///        C: = alpha*A*B^T + alpha*B*A^T + beta * C
///        ArgUplo == Upper
///        Ref: [[ -1.8 -12.6   8.1 -27. ]
///              [  3.   24.6  16.8  -2.4]
///              [  2.    1.  124.8 172.5]
///              [  0.    2.   -5.  239.4]]
///
///        ArgUplo == Lower
///        Ref: [[ -1.8  -3.   -2.    0. ]
///              [ -5.4  24.6  -1.   -2. ]
///              [ 12.9  19.2 124.8   5. ]
///              [-27.    2.4 160.5 239.4]]
///
///        ArgTrans == Transpose/ConjTranspose
///        C: = alpha*A^T*B + alpha*B^T*A + beta * C
///        ArgUplo == Upper
///        Ref: [[ 16.2   3.9  38.1  79.5]
///              [  3.   24.6 -35.7 -11.4]
///              [  2.    1.   76.8 -58.5]
///              [  0.    2.   -5.  269.4]]
///
///        ArgUplo == Lower
///        Ref: [[ 16.2  -3.   -2.    0. ]
///              [ 11.1  24.6  -1.   -2. ]
///              [ 42.9 -33.3  76.8   5. ]
///              [ 79.5  -6.6 -70.5 269.4]]
/// \param[in] Nb Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_syr2k_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using StridedView3DType = Kokkos::View<ScalarType ***, Kokkos::LayoutStride, DeviceType>;
  using ArgUplo           = typename ParamTagType::uplo;
  using ArgTrans          = typename ParamTagType::trans;

  const std::size_t BlkSize = 4;
  View3DType A("A", Nb, BlkSize, BlkSize), B("B", Nb, BlkSize, BlkSize);
  View3DType C("C", Nb, BlkSize, BlkSize), C_ref("C_ref", Nb, BlkSize, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, BlkSize, Nb * incx, BlkSize, Nb * incx * BlkSize};
  StridedView3DType A_s("A_s", layout), B_s("B_s", layout), C_s("C_s", layout);

  // Only filling a and b
  auto h_A     = Kokkos::create_mirror_view(A);
  auto h_B     = Kokkos::create_mirror_view(B);
  auto h_C_ref = Kokkos::create_mirror_view(C_ref);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    h_A(ib, 0, 0) = 1;
    h_A(ib, 0, 1) = -3;
    h_A(ib, 0, 2) = -2;
    h_A(ib, 0, 3) = 0;
    h_A(ib, 1, 0) = 3;
    h_A(ib, 1, 1) = 3;
    h_A(ib, 1, 2) = -1;
    h_A(ib, 1, 3) = -2;
    h_A(ib, 2, 0) = 2;
    h_A(ib, 2, 1) = 1;
    h_A(ib, 2, 2) = 9;
    h_A(ib, 2, 3) = 5;
    h_A(ib, 3, 0) = 0;
    h_A(ib, 3, 1) = 2;
    h_A(ib, 3, 2) = -5;
    h_A(ib, 3, 3) = 27;

    h_B(ib, 0, 0) = 2;
    h_B(ib, 0, 1) = 1;
    h_B(ib, 0, 2) = 0;
    h_B(ib, 0, 3) = -1;
    h_B(ib, 1, 0) = -1;
    h_B(ib, 1, 1) = 4;
    h_B(ib, 1, 2) = 2;
    h_B(ib, 1, 3) = 0;
    h_B(ib, 2, 0) = 3;
    h_B(ib, 2, 1) = -2;
    h_B(ib, 2, 2) = 1;
    h_B(ib, 2, 3) = 5;
    h_B(ib, 3, 0) = 1;
    h_B(ib, 3, 1) = 0;
    h_B(ib, 3, 2) = -3;
    h_B(ib, 3, 3) = 2;

    if (std::same_as<ArgTrans, KokkosBatched::Trans::NoTranspose> ||
        std::same_as<ArgTrans, KokkosBatched::Trans::ConjNoTranspose>) {
      if (std::same_as<ArgUplo, KokkosBatched::Uplo::Upper>) {
        h_C_ref(ib, 0, 0) = -1.8;
        h_C_ref(ib, 0, 1) = -12.6;
        h_C_ref(ib, 0, 2) = 8.1;
        h_C_ref(ib, 0, 3) = -27.0;
        h_C_ref(ib, 1, 0) = 3;
        h_C_ref(ib, 1, 1) = 24.6;
        h_C_ref(ib, 1, 2) = 16.8;
        h_C_ref(ib, 1, 3) = -2.4;
        h_C_ref(ib, 2, 0) = 2;
        h_C_ref(ib, 2, 1) = 1;
        h_C_ref(ib, 2, 2) = 124.8;
        h_C_ref(ib, 2, 3) = 172.5;
        h_C_ref(ib, 3, 0) = 0;
        h_C_ref(ib, 3, 1) = 2;
        h_C_ref(ib, 3, 2) = -5;
        h_C_ref(ib, 3, 3) = 239.4;
      } else {
        h_C_ref(ib, 0, 0) = -1.8;
        h_C_ref(ib, 0, 1) = -3;
        h_C_ref(ib, 0, 2) = -2;
        h_C_ref(ib, 0, 3) = 0;
        h_C_ref(ib, 1, 0) = -5.4;
        h_C_ref(ib, 1, 1) = 24.6;
        h_C_ref(ib, 1, 2) = -1;
        h_C_ref(ib, 1, 3) = -2;
        h_C_ref(ib, 2, 0) = 12.9;
        h_C_ref(ib, 2, 1) = 19.2;
        h_C_ref(ib, 2, 2) = 124.8;
        h_C_ref(ib, 2, 3) = 5;
        h_C_ref(ib, 3, 0) = -27;
        h_C_ref(ib, 3, 1) = 2.4;
        h_C_ref(ib, 3, 2) = 160.5;
        h_C_ref(ib, 3, 3) = 239.4;
      }
    } else {
      if (std::same_as<ArgUplo, KokkosBatched::Uplo::Upper>) {
        h_C_ref(ib, 0, 0) = 16.2;
        h_C_ref(ib, 0, 1) = 3.9;
        h_C_ref(ib, 0, 2) = 38.1;
        h_C_ref(ib, 0, 3) = 79.5;
        h_C_ref(ib, 1, 0) = 3;
        h_C_ref(ib, 1, 1) = 24.6;
        h_C_ref(ib, 1, 2) = -35.7;
        h_C_ref(ib, 1, 3) = -11.4;
        h_C_ref(ib, 2, 0) = 2;
        h_C_ref(ib, 2, 1) = 1;
        h_C_ref(ib, 2, 2) = 76.8;
        h_C_ref(ib, 2, 3) = -58.5;
        h_C_ref(ib, 3, 0) = 0;
        h_C_ref(ib, 3, 1) = 2;
        h_C_ref(ib, 3, 2) = -5;
        h_C_ref(ib, 3, 3) = 269.4;
      } else {
        h_C_ref(ib, 0, 0) = 16.2;
        h_C_ref(ib, 0, 1) = -3;
        h_C_ref(ib, 0, 2) = -2;
        h_C_ref(ib, 0, 3) = 0;
        h_C_ref(ib, 1, 0) = 11.1;
        h_C_ref(ib, 1, 1) = 24.6;
        h_C_ref(ib, 1, 2) = -1;
        h_C_ref(ib, 1, 3) = -2;
        h_C_ref(ib, 2, 0) = 42.9;
        h_C_ref(ib, 2, 1) = -33.3;
        h_C_ref(ib, 2, 2) = 76.8;
        h_C_ref(ib, 2, 3) = 5;
        h_C_ref(ib, 3, 0) = 79.5;
        h_C_ref(ib, 3, 1) = -6.6;
        h_C_ref(ib, 3, 2) = -70.5;
        h_C_ref(ib, 3, 3) = 269.4;
      }
    }
  }

  Kokkos::deep_copy(A, h_A);
  Kokkos::deep_copy(B, h_B);
  Kokkos::deep_copy(C, A);
  Kokkos::deep_copy(A_s, A);
  Kokkos::deep_copy(B_s, B);
  Kokkos::deep_copy(C_s, A);

  using AlphaType       = ScalarType;
  using BetaType        = std::conditional_t<std::same_as<ArgTrans, KokkosBatched::Trans::ConjNoTranspose> ||
                                          std::same_as<ArgTrans, KokkosBatched::Trans::ConjTranspose>,
                                      RealType, ScalarType>;
  const AlphaType alpha = 1.5;
  const BetaType beta   = 1.2;

  auto info = Functor_BatchedSyr2k<DeviceType, ParamTagType, AlphaType, View3DType, View3DType, BetaType, View3DType>(
                  alpha, A, B, beta, C)
                  .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // With strided views
  info = Functor_BatchedSyr2k<DeviceType, ParamTagType, AlphaType, StridedView3DType, StridedView3DType, BetaType,
                              StridedView3DType>(alpha, A_s, B_s, beta, C_s)
             .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  RealType eps = 1.0e2 * ats::epsilon();
  auto h_C     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C);

  // Check if C:= alpha*A*B^T + alpha*B*A^T + beta*C or C:= alpha*A^T*B + alpha*B^T*A + beta*C is computed correctly
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      for (std::size_t j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_C(ib, i, j), h_C_ref(ib, i, j), eps);
      }
    }
  }

  // Testing for strided views, reusing C
  Kokkos::deep_copy(C, C_s);
  Kokkos::deep_copy(h_C, C);

  // Check if C:= alpha*A*B^T + alpha*B*A^T + beta*C or C:= alpha*A^T*B + alpha*B^T*A + beta*C is computed correctly
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      for (std::size_t j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_C(ib, i, j), h_C_ref(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched syr2k test
///
/// \param[in] Nb Batch size of matrices
/// \param[in] M Number of rows of matrices
/// \param[in] N Number of columns of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_syr2k(const std::size_t Nb, const std::size_t M, const std::size_t N) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo    = typename ParamTagType::uplo;
  using ArgTrans   = typename ParamTagType::trans;

  constexpr bool is_trans = std::same_as<ArgTrans, KokkosBatched::Trans::Transpose> ||
                            std::same_as<ArgTrans, KokkosBatched::Trans::ConjTranspose>;

  // For NoTrans, A is N*M, C is NxN
  // For Trans, A is M*N, C is MxM
  const std::size_t n = is_trans ? M : N;
  const std::size_t k = is_trans ? N : M;

  View3DType A("A", Nb, n, k), B("B", Nb, n, k), C("C", Nb, N, N), C_alpha0("C_alpha0", Nb, N, N),
      C_beta0("C_beta0", Nb, N, N), C_ref("C_ref", Nb, N, N), C_alpha0_ref("C_alpha0_ref", Nb, N, N),
      C_beta0_ref("C_beta0_ref", Nb, N, N);

  // Create a random matrix A and x
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(B, rand_pool, randStart, randEnd);
  Kokkos::fill_random(C_ref, rand_pool, randStart, randEnd);

  // Make C Symmetric or Hermitian
  dense_to_symmetric<ArgUplo, KokkosBatched::Trans::ConjTranspose>(C_ref, C);

  Kokkos::deep_copy(C_ref, C);
  Kokkos::deep_copy(C_alpha0, C);
  Kokkos::deep_copy(C_beta0, C);
  Kokkos::deep_copy(C_alpha0_ref, C_alpha0);
  Kokkos::deep_copy(C_beta0_ref, C_beta0);

  using AlphaType  = ScalarType;
  using BetaType   = std::conditional_t<std::same_as<ArgTrans, KokkosBatched::Trans::ConjNoTranspose> ||
                                          std::same_as<ArgTrans, KokkosBatched::Trans::ConjTranspose>,
                                      RealType, ScalarType>;
  AlphaType alpha0 = 1.5;
  if constexpr (KokkosKernels::ArithTraits<ScalarType>::is_complex) {
    alpha0 = AlphaType(1.5, 1.8);
  }

  const AlphaType alpha = alpha0;
  const BetaType beta = 1.2, zero = KokkosKernels::ArithTraits<BetaType>::zero();

  const ScalarType cbeta = 1.2, czero = KokkosKernels::ArithTraits<ScalarType>::zero(),
                   cone = KokkosKernels::ArithTraits<ScalarType>::one();
  auto info0 = Functor_BatchedSyr2k<DeviceType, ParamTagType, AlphaType, View3DType, View3DType, BetaType, View3DType>(
                   alpha, A, B, beta, C)
                   .run();
  auto info1 = Functor_BatchedSyr2k<DeviceType, ParamTagType, AlphaType, View3DType, View3DType, BetaType, View3DType>(
                   zero, A, B, beta, C_alpha0)
                   .run();
  auto info2 = Functor_BatchedSyr2k<DeviceType, ParamTagType, AlphaType, View3DType, View3DType, BetaType, View3DType>(
                   alpha, A, B, zero, C_beta0)
                   .run();

  Kokkos::fence();

  if (is_trans) {
    if (M == 0 && N > 0) {
// Quick return case: M=0, N>0, this is not allowed in blas
#ifndef NDEBUG
      EXPECT_GT(info0, 0);
      EXPECT_GT(info1, 0);
      EXPECT_GT(info2, 0);
#else
      EXPECT_EQ(info0, 0);
      EXPECT_EQ(info1, 0);
      EXPECT_EQ(info2, 0);
#endif
      return;
    }
  }

  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);
  EXPECT_EQ(info2, 0);

  using ArgTransA = std::conditional_t<is_trans, typename ParamTagType::trans, KokkosBatched::Trans::NoTranspose>;
  using ArgTransB = std::conditional_t<
      is_trans, KokkosBatched::Trans::NoTranspose,
      std::conditional_t<std::same_as<typename ParamTagType::trans, KokkosBatched::Trans::NoTranspose>,
                         KokkosBatched::Trans::Transpose, KokkosBatched::Trans::ConjTranspose>>;

  // Make a reference with GEMM
  // C_ref = alpha * A * B^T + beta * C_ref if (ArgTrans == Trans::NoTranspose)
  // C_ref = alpha * A^T * B + beta * C_ref if (ArgTrans == Trans::Transpose)
  using ArgMode = typename ParamTagType::mode;
  Functor_BatchedGemm<DeviceType, ArgMode, ArgTransA, ArgTransB, ScalarType, View3DType, View3DType, View3DType>(
      alpha, A, B, cbeta, C_ref)
      .run();
  Functor_BatchedGemm<DeviceType, ArgMode, ArgTransA, ArgTransB, ScalarType, View3DType, View3DType, View3DType>(
      czero, A, B, cbeta, C_alpha0_ref)
      .run();
  Functor_BatchedGemm<DeviceType, ArgMode, ArgTransA, ArgTransB, ScalarType, View3DType, View3DType, View3DType>(
      alpha, A, B, czero, C_beta0_ref)
      .run();

  // C_ref = alpha * B * A^T + C_ref if (ArgTrans == Trans::NoTranspose)
  // C_ref = alpha * B^T * A + C_ref if (ArgTrans == Trans::Transpose)
  // Note: mathematically, alpha * A * B^H + conj(alpha) * B * A^H + beta * C would have real diagnoals if C is
  // Hermitian matrix.
  const ScalarType maybe_conj_alpha =
      (std::same_as<typename ParamTagType::trans, KokkosBatched::Trans::ConjTranspose> ||
       std::same_as<typename ParamTagType::trans, KokkosBatched::Trans::ConjNoTranspose>)
          ? KokkosKernels::ArithTraits<ScalarType>::conj(alpha)
          : alpha;
  Functor_BatchedGemm<DeviceType, ArgMode, ArgTransA, ArgTransB, ScalarType, View3DType, View3DType, View3DType>(
      maybe_conj_alpha, B, A, cone, C_ref)
      .run();
  Functor_BatchedGemm<DeviceType, ArgMode, ArgTransA, ArgTransB, ScalarType, View3DType, View3DType, View3DType>(
      czero, B, A, cone, C_alpha0_ref)
      .run();
  Functor_BatchedGemm<DeviceType, ArgMode, ArgTransA, ArgTransB, ScalarType, View3DType, View3DType, View3DType>(
      maybe_conj_alpha, B, A, cone, C_beta0_ref)
      .run();

  RealType eps = 1.0e1 * ats::epsilon();

  auto h_C            = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C);
  auto h_C_alpha0     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_alpha0);
  auto h_C_beta0      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_beta0);
  auto h_C_ref        = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_ref);
  auto h_C_alpha0_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_alpha0_ref);
  auto h_C_beta0_ref  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_beta0_ref);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    if (std::same_as<ArgUplo, KokkosBatched::Uplo::Lower>) {
      for (std::size_t i = 0; i < N; i++) {
        for (std::size_t j = 0; j < i + 1; j++) {
          EXPECT_NEAR_KK(h_C(ib, i, j), h_C_ref(ib, i, j), eps);
          EXPECT_NEAR_KK(h_C_alpha0(ib, i, j), h_C_alpha0_ref(ib, i, j), eps);
          EXPECT_NEAR_KK(h_C_beta0(ib, i, j), h_C_beta0_ref(ib, i, j), eps);
        }
      }
    } else {
      for (std::size_t i = 0; i < N; i++) {
        for (std::size_t j = i; j < N; j++) {
          EXPECT_NEAR_KK(h_C(ib, i, j), h_C_ref(ib, i, j), eps);
          EXPECT_NEAR_KK(h_C_alpha0(ib, i, j), h_C_alpha0_ref(ib, i, j), eps);
          EXPECT_NEAR_KK(h_C_beta0(ib, i, j), h_C_beta0_ref(ib, i, j), eps);
        }
      }
    }
  }
}

}  // namespace Syr2k
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType, typename LayoutType>
void test_batched_syr2k_layout() {
  Test::Syr2k::impl_test_batched_syr2k_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
  Test::Syr2k::impl_test_batched_syr2k_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
  for (int i = 1; i < 2; i++) {
    for (int j = 1; j < 2; j++) {
      Test::Syr2k::impl_test_batched_syr2k<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i, j);
      Test::Syr2k::impl_test_batched_syr2k<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i, j);
    }
  }
  Test::Syr2k::impl_test_batched_syr2k<DeviceType, ScalarType, LayoutType, ParamTagType>(1, 0, 0);
  Test::Syr2k::impl_test_batched_syr2k<DeviceType, ScalarType, LayoutType, ParamTagType>(1, 5, 0);
  Test::Syr2k::impl_test_batched_syr2k<DeviceType, ScalarType, LayoutType, ParamTagType>(1, 0, 5);
}

template <typename DeviceType, typename ScalarType, typename ParamTagType>
int test_batched_syr2k() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  test_batched_syr2k_layout<DeviceType, ScalarType, ParamTagType, Kokkos::LayoutLeft>();
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  test_batched_syr2k_layout<DeviceType, ScalarType, ParamTagType, Kokkos::LayoutRight>();
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_syr2k_l_nt_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_t_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_nc_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_c_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_nt_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_t_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_nc_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_c_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
// Team
TEST_F(TestCategory, test_batched_team_syr2k_l_nt_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_t_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_nc_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_c_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_nt_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_t_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_nc_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_c_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
// TeamVector
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_nt_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_t_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_nc_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_c_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_nt_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_t_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_nc_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_c_float) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_syr2k_l_nt_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_t_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_nc_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_c_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_nt_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_t_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_nc_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_c_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
// Team
TEST_F(TestCategory, test_batched_team_syr2k_l_nt_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_t_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_nc_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_c_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_nt_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_t_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_nc_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_c_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
// TeamVector
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_nt_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_t_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_nc_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_c_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_nt_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_t_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_nc_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_c_double) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, double, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_syr2k_l_nt_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_t_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_nc_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_c_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_nt_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_t_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_nc_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_c_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
// Team
TEST_F(TestCategory, test_batched_team_syr2k_l_nt_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_t_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_nc_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_c_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_nt_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_t_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_nc_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_c_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
// TeamVector
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_nt_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_t_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_nc_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_c_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_nt_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_t_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_nc_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_c_fcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_syr2k_l_nt_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_t_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_nc_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_l_c_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_nt_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_t_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_nc_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syr2k_u_c_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
// Team
TEST_F(TestCategory, test_batched_team_syr2k_l_nt_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_t_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_nc_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_l_c_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_nt_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_t_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_nc_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syr2k_u_c_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
// TeamVector
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_nt_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_t_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_nc_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_l_c_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_nt_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_t_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_nc_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syr2k_u_c_dcomplex) {
  using param_tag_type = ::Test::Syr2k::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syr2k<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
#endif
