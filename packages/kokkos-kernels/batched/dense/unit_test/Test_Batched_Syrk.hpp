// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Syrk.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Syrk {

template <typename U, typename T, typename M>
struct ParamTag {
  using uplo  = U;
  using trans = T;
  using mode  = M;
};

template <typename DeviceType, typename ParamTagType, typename ScalarType, typename AViewType, typename CViewType>
struct Functor_BatchedSyrk {
  using execution_space = typename DeviceType::execution_space;
  using member_type     = typename Kokkos::TeamPolicy<execution_space>::member_type;
  using ArgMode         = typename ParamTagType::mode;
  using ArgUplo         = typename ParamTagType::uplo;
  using ArgTrans        = typename ParamTagType::trans;
  AViewType m_A;
  CViewType m_C;
  ScalarType m_alpha, m_beta;

  Functor_BatchedSyrk(const ScalarType alpha, const AViewType &A, const ScalarType beta, const CViewType &C)
      : m_A(A), m_C(C), m_alpha(alpha), m_beta(beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type &member, int &info) const {
    const int k = member.league_rank();
    auto sub_C  = Kokkos::subview(m_C, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_A  = Kokkos::subview(m_A, k, Kokkos::ALL(), Kokkos::ALL());
    if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        info += KokkosBatched::SerialSyrk<ArgUplo, ArgTrans>::invoke(m_alpha, sub_A, m_beta, sub_C);
      });
    } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Team>) {
      info += KokkosBatched::TeamSyrk<member_type, ArgUplo, ArgTrans>::invoke(member, m_alpha, sub_A, m_beta, sub_C);
    } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::TeamVector>) {
      info +=
          KokkosBatched::TeamVectorSyrk<member_type, ArgUplo, ArgTrans>::invoke(member, m_alpha, sub_A, m_beta, sub_C);
    }
  }

  inline int run() {
    using value_type        = typename AViewType::non_const_value_type;
    std::string name_region = std::same_as<ArgMode, KokkosBatched::Mode::Serial> ? "KokkosBatched::Test::SerialSyrk"
                              : std::same_as<ArgMode, KokkosBatched::Mode::Team>
                                  ? "KokkosBatched::Test::TeamSyrk"
                                  : "KokkosBatched::Test::TeamVectorSyrk";
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

/// \brief Implementation details of batched syrk analytical test
///        to confirm C:= A*A**T + C or C:= A**T*A + C is computed correctly
///        A and C are deliberately chosen to have non symmetric values to confirm both cases.
///        alpha = 1.5, beta = 1.2
///        4x4 matrix (upper)
///        A: [[1,  -3, -2,  0],
///            [3,  3, -1, -2],
///            [2,  1,  9, 5],
///            [0,  2, -5, 27]]
///        C = A
///        ArgTrans == NoTranspose/ConjNoTranspose
///        C: = alpha * U * U**T + beta * C
///        Ref: [[22.2, -9.6, -30.9,   6.0],
///              [3,    38.1, -16.2, -66.9],
///              [2,       1, 177.3, 144.0],
///              [0,       2,   -5, 1169.4]]
///
///        C: = alpha * L**T * L + beta * C
///        Ref: [[ 22.2,    -3,    -2,      0],
///              [ -2.4,  38.1,    -1,     -2],
///              [-26.1, -13.8, 177.3,      5],
///              [  6.0, -62.1, 132.0, 1169.4]]
///
///        ArgTrans == Transpose/ConjTranspose
///        C: = alpha * U**T * U + beta * C
///        Ref: [[22.2,  8.4,  17.1,    6.0],
///              [3,    38.1,   1.8,   77.1],
///              [2,       1, 177.3, -126.0],
///              [0,       2,    -5, 1169.4]]
///
///        C: = alpha * L**T * L + beta * C
///        Ref: [[ 22.2,    -3,     -2,      0],
///              [ 15.6,  38.1,     -1,     -2],
///              [ 21.9,   4.2,  177.3,    5.0],
///              [  6.0,  81.9, -138.0, 1169.4]]
/// \param[in] Nb Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_syrk_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using StridedView3DType = Kokkos::View<ScalarType ***, Kokkos::LayoutStride, DeviceType>;
  using ArgUplo           = typename ParamTagType::uplo;
  using ArgTrans          = typename ParamTagType::trans;

  const std::size_t BlkSize = 4;
  View3DType A("A", Nb, BlkSize, BlkSize);
  View3DType C("C", Nb, BlkSize, BlkSize), C_ref("C_ref", Nb, BlkSize, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, BlkSize, Nb * incx, BlkSize, Nb * incx * BlkSize};
  StridedView3DType A_s("A_s", layout), C_s("C_s", layout);

  // Only filling a
  auto h_A     = Kokkos::create_mirror_view(A);
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

    if (std::same_as<ArgTrans, KokkosBatched::Trans::NoTranspose> ||
        std::same_as<ArgTrans, KokkosBatched::Trans::ConjNoTranspose>) {
      if (std::same_as<ArgUplo, KokkosBatched::Uplo::Upper>) {
        h_C_ref(ib, 0, 0) = 22.2;
        h_C_ref(ib, 0, 1) = -9.6;
        h_C_ref(ib, 0, 2) = -30.9;
        h_C_ref(ib, 0, 3) = 6;
        h_C_ref(ib, 1, 0) = 3;
        h_C_ref(ib, 1, 1) = 38.1;
        h_C_ref(ib, 1, 2) = -16.2;
        h_C_ref(ib, 1, 3) = -66.9;
        h_C_ref(ib, 2, 0) = 2;
        h_C_ref(ib, 2, 1) = 1;
        h_C_ref(ib, 2, 2) = 177.3;
        h_C_ref(ib, 2, 3) = 144;
        h_C_ref(ib, 3, 0) = 0;
        h_C_ref(ib, 3, 1) = 2;
        h_C_ref(ib, 3, 2) = -5;
        h_C_ref(ib, 3, 3) = 1169.4;
      } else {
        h_C_ref(ib, 0, 0) = 22.2;
        h_C_ref(ib, 0, 1) = -3;
        h_C_ref(ib, 0, 2) = -2;
        h_C_ref(ib, 0, 3) = 0;
        h_C_ref(ib, 1, 0) = -2.4;
        h_C_ref(ib, 1, 1) = 38.1;
        h_C_ref(ib, 1, 2) = -1;
        h_C_ref(ib, 1, 3) = -2;
        h_C_ref(ib, 2, 0) = -26.1;
        h_C_ref(ib, 2, 1) = -13.8;
        h_C_ref(ib, 2, 2) = 177.3;
        h_C_ref(ib, 2, 3) = 5;
        h_C_ref(ib, 3, 0) = 6;
        h_C_ref(ib, 3, 1) = -62.1;
        h_C_ref(ib, 3, 2) = 132;
        h_C_ref(ib, 3, 3) = 1169.4;
      }
    } else {
      if (std::same_as<ArgUplo, KokkosBatched::Uplo::Upper>) {
        h_C_ref(ib, 0, 0) = 22.2;
        h_C_ref(ib, 0, 1) = 8.4;
        h_C_ref(ib, 0, 2) = 17.1;
        h_C_ref(ib, 0, 3) = 6;
        h_C_ref(ib, 1, 0) = 3;
        h_C_ref(ib, 1, 1) = 38.1;
        h_C_ref(ib, 1, 2) = 1.8;
        h_C_ref(ib, 1, 3) = 77.1;
        h_C_ref(ib, 2, 0) = 2;
        h_C_ref(ib, 2, 1) = 1;
        h_C_ref(ib, 2, 2) = 177.3;
        h_C_ref(ib, 2, 3) = -126;
        h_C_ref(ib, 3, 0) = 0;
        h_C_ref(ib, 3, 1) = 2;
        h_C_ref(ib, 3, 2) = -5;
        h_C_ref(ib, 3, 3) = 1169.4;
      } else {
        h_C_ref(ib, 0, 0) = 22.2;
        h_C_ref(ib, 0, 1) = -3;
        h_C_ref(ib, 0, 2) = -2;
        h_C_ref(ib, 0, 3) = 0;
        h_C_ref(ib, 1, 0) = 15.6;
        h_C_ref(ib, 1, 1) = 38.1;
        h_C_ref(ib, 1, 2) = -1;
        h_C_ref(ib, 1, 3) = -2;
        h_C_ref(ib, 2, 0) = 21.9;
        h_C_ref(ib, 2, 1) = 4.2;
        h_C_ref(ib, 2, 2) = 177.3;
        h_C_ref(ib, 2, 3) = 5;
        h_C_ref(ib, 3, 0) = 6;
        h_C_ref(ib, 3, 1) = 81.9;
        h_C_ref(ib, 3, 2) = -138;
        h_C_ref(ib, 3, 3) = 1169.4;
      }
    }
  }

  Kokkos::deep_copy(A, h_A);
  Kokkos::deep_copy(C, A);
  Kokkos::deep_copy(A_s, A);
  Kokkos::deep_copy(C_s, A);

  const ScalarType alpha = 1.5, beta = 1.2;

  auto info =
      Functor_BatchedSyrk<DeviceType, ParamTagType, ScalarType, View3DType, View3DType>(alpha, A, beta, C).run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // With strided views
  info = Functor_BatchedSyrk<DeviceType, ParamTagType, ScalarType, StridedView3DType, StridedView3DType>(alpha, A_s,
                                                                                                         beta, C_s)
             .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  RealType eps = 1.0e1 * ats::epsilon();
  auto h_C     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C);

  // Check if C:= alpha * A * A**T + C or C:= alpha * A**T * A + C is computed correctly
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

  // Check if C:= alpha * A * A**T + C or C:= alpha * A**T * A + C is computed correctly
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      for (std::size_t j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_C(ib, i, j), h_C_ref(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched syrk test
///
/// \param[in] Nb Batch size of matrices
/// \param[in] M Number of rows of matrices
/// \param[in] N Number of columns of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_syrk(const std::size_t Nb, const std::size_t M, const std::size_t N) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo    = typename ParamTagType::uplo;

  const bool is_trans = std::same_as<typename ParamTagType::trans, KokkosBatched::Trans::Transpose> ||
                        std::same_as<typename ParamTagType::trans, KokkosBatched::Trans::ConjTranspose>;

  // For NoTrans, A is N*M, C is NxN
  // For Trans, A is M*N, C is MxM
  const std::size_t n = is_trans ? M : N;
  const std::size_t k = is_trans ? N : M;

  View3DType A("A", Nb, n, k), C("C", Nb, N, N), C_alpha0("C_alpha0", Nb, N, N), C_beta0("C_beta0", Nb, N, N),
      C_ref("C_ref", Nb, N, N), C_alpha0_ref("C_alpha0_ref", Nb, N, N), C_beta0_ref("C_beta0_ref", Nb, N, N);

  // Create a random matrix A and x
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(C, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(C_ref, C);
  Kokkos::deep_copy(C_alpha0, C);
  Kokkos::deep_copy(C_beta0, C);
  Kokkos::deep_copy(C_alpha0_ref, C_alpha0);
  Kokkos::deep_copy(C_beta0_ref, C_beta0);

  // When beta is zero
  const RealType alpha = 1.5, beta = 1.2, zero = KokkosKernels::ArithTraits<RealType>::zero();
  const ScalarType czero = KokkosKernels::ArithTraits<ScalarType>::zero();
  auto info0 =
      Functor_BatchedSyrk<DeviceType, ParamTagType, ScalarType, View3DType, View3DType>(alpha, A, beta, C).run();
  auto info1 =
      Functor_BatchedSyrk<DeviceType, ParamTagType, ScalarType, View3DType, View3DType>(zero, A, beta, C_alpha0).run();
  auto info2 =
      Functor_BatchedSyrk<DeviceType, ParamTagType, ScalarType, View3DType, View3DType>(alpha, A, zero, C_beta0).run();

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

  // Make a reference at host
  auto h_A            = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, A);
  auto h_C_ref        = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_ref);
  auto h_C_alpha0_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_alpha0_ref);
  auto h_C_beta0_ref  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_beta0_ref);

  // Note: ConjTranspose/ConjNoTranspose corresponds to {c,z}herk for Hermitian matrix
  constexpr bool is_conj = std::same_as<typename ParamTagType::trans, Trans::ConjTranspose> ||
                           std::same_as<typename ParamTagType::trans, Trans::ConjNoTranspose>;
  using Op     = std::conditional_t<is_conj, KokkosBlas::Impl::OpConj, KokkosBlas::Impl::OpID>;
  using Sym_Op = std::conditional_t<is_conj, KokkosBlas::Impl::OpReal, KokkosBlas::Impl::OpID>;

  Op op;
  Sym_Op sym_op;

  for (std::size_t ib = 0; ib < Nb; ib++) {
    if (std::same_as<ArgUplo, Uplo::Upper>) {
      for (std::size_t j = 0; j < N; j++) {
        for (std::size_t i = 0; i < j + 1; i++) {
          h_C_alpha0_ref(ib, i, j) = beta * h_C_alpha0_ref(ib, i, j);
        }
        h_C_alpha0_ref(ib, j, j) = sym_op(h_C_alpha0_ref(ib, j, j));
      }
    } else {
      for (std::size_t j = 0; j < N; j++) {
        for (std::size_t i = j; i < N; i++) {
          h_C_alpha0_ref(ib, i, j) = beta * h_C_alpha0_ref(ib, i, j);
        }
        h_C_alpha0_ref(ib, j, j) = sym_op(h_C_alpha0_ref(ib, j, j));
      }
    }

    if (!is_trans) {
      if (std::same_as<ArgUplo, Uplo::Upper>) {
        for (std::size_t j = 0; j < N; j++) {
          for (std::size_t i = 0; i < j + 1; i++) {
            h_C_beta0_ref(ib, i, j) = zero;
            h_C_ref(ib, i, j)       = beta * h_C_ref(ib, i, j);
          }
          for (std::size_t l = 0; l < M; l++) {
            if (h_A(ib, j, l) != czero) {
              auto temp = alpha * op(h_A(ib, j, l));
              for (std::size_t i = 0; i < j + 1; i++) {
                h_C_ref(ib, i, j) += temp * h_A(ib, i, l);
                h_C_beta0_ref(ib, i, j) += temp * h_A(ib, i, l);
              }
            }
          }
          h_C_ref(ib, j, j)       = sym_op(h_C_ref(ib, j, j));
          h_C_beta0_ref(ib, j, j) = sym_op(h_C_beta0_ref(ib, j, j));
        }
      } else {
        for (std::size_t j = 0; j < N; j++) {
          for (std::size_t i = j; i < N; i++) {
            h_C_beta0_ref(ib, i, j) = zero;
            h_C_ref(ib, i, j)       = beta * h_C_ref(ib, i, j);
          }
          for (std::size_t l = 0; l < M; l++) {
            if (h_A(ib, j, l) != czero) {
              auto temp = alpha * op(h_A(ib, j, l));
              for (std::size_t i = j; i < N; i++) {
                h_C_ref(ib, i, j) += temp * h_A(ib, i, l);
                h_C_beta0_ref(ib, i, j) += temp * h_A(ib, i, l);
              }
            }
          }
          h_C_ref(ib, j, j)       = sym_op(h_C_ref(ib, j, j));
          h_C_beta0_ref(ib, j, j) = sym_op(h_C_beta0_ref(ib, j, j));
        }
      }
    } else {
      if (std::same_as<ArgUplo, Uplo::Upper>) {
        for (std::size_t j = 0; j < N; j++) {
          for (std::size_t i = 0; i < j + 1; i++) {
            auto temp = czero;
            for (std::size_t l = 0; l < M; l++) {
              temp += op(h_A(ib, l, i)) * h_A(ib, l, j);
            }
            h_C_beta0_ref(ib, i, j) = alpha * temp;
            h_C_ref(ib, i, j)       = alpha * temp + beta * h_C_ref(ib, i, j);
          }
          h_C_ref(ib, j, j)       = sym_op(h_C_ref(ib, j, j));
          h_C_beta0_ref(ib, j, j) = sym_op(h_C_beta0_ref(ib, j, j));
        }
      } else {
        for (std::size_t j = 0; j < N; j++) {
          for (std::size_t i = j; i < N; i++) {
            auto temp = czero;
            for (std::size_t l = 0; l < M; l++) {
              temp += op(h_A(ib, l, i)) * h_A(ib, l, j);
            }
            h_C_beta0_ref(ib, i, j) = alpha * temp;
            h_C_ref(ib, i, j)       = alpha * temp + beta * h_C_ref(ib, i, j);
          }
          h_C_ref(ib, j, j)       = sym_op(h_C_ref(ib, j, j));
          h_C_beta0_ref(ib, j, j) = sym_op(h_C_beta0_ref(ib, j, j));
        }
      }
    }
  }

  RealType eps = 1.0e1 * ats::epsilon();

  auto h_C        = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C);
  auto h_C_alpha0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_alpha0);
  auto h_C_beta0  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_beta0);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < N; i++) {
      for (std::size_t j = 0; j < N; j++) {
        EXPECT_NEAR_KK(h_C(ib, i, j), h_C_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_C_alpha0(ib, i, j), h_C_alpha0_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_C_beta0(ib, i, j), h_C_beta0_ref(ib, i, j), eps);
      }
    }
  }
}

}  // namespace Syrk
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType>
int test_batched_syrk() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Syrk::impl_test_batched_syrk_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Syrk::impl_test_batched_syrk_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        Test::Syrk::impl_test_batched_syrk<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i, j);
        Test::Syrk::impl_test_batched_syrk<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i, j);
      }
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Syrk::impl_test_batched_syrk_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Syrk::impl_test_batched_syrk_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        Test::Syrk::impl_test_batched_syrk<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i, j);
        Test::Syrk::impl_test_batched_syrk<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i, j);
      }
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_syrk_l_nt_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_t_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_nc_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_c_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_nt_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_t_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_nc_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_c_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
// Team
TEST_F(TestCategory, test_batched_team_syrk_l_nt_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_t_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_nc_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_c_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_nt_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_t_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_nc_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_c_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
// TeamVector
TEST_F(TestCategory, test_batched_team_vector_syrk_l_nt_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_t_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_nc_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_c_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_nt_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_t_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_nc_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_c_float) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_syrk_l_nt_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_t_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_nc_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_c_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_nt_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_t_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_nc_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_c_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
// Team
TEST_F(TestCategory, test_batched_team_syrk_l_nt_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_t_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_nc_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_c_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_nt_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_t_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_nc_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_c_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
// TeamVector
TEST_F(TestCategory, test_batched_team_vector_syrk_l_nt_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_t_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_nc_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_c_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_nt_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_t_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_nc_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_c_double) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, double, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_syrk_l_nt_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_t_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_nc_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_c_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_nt_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_t_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_nc_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_c_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
// Team
TEST_F(TestCategory, test_batched_team_syrk_l_nt_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_t_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_nc_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_c_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_nt_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_t_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_nc_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_c_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
// TeamVector
TEST_F(TestCategory, test_batched_team_vector_syrk_l_nt_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_t_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_nc_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_c_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_nt_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_t_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_nc_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_c_fcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_syrk_l_nt_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_t_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_nc_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_l_c_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_nt_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_t_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_nc_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_syrk_u_c_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Serial>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
// Team
TEST_F(TestCategory, test_batched_team_syrk_l_nt_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_t_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_nc_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_l_c_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_nt_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_t_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_nc_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_syrk_u_c_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::Team>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
// TeamVector
TEST_F(TestCategory, test_batched_team_vector_syrk_l_nt_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_t_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_nc_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_l_c_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_nt_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_t_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_nc_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjNoTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_syrk_u_c_dcomplex) {
  using param_tag_type = ::Test::Syrk::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>;
  test_batched_syrk<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
#endif
