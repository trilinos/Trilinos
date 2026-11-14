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
#include <KokkosBatched_Symm.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Symm {

template <bool Conj, typename S, typename U, typename M>
struct ParamTag {
  static constexpr bool conj = Conj;
  using side                 = S;
  using uplo                 = U;
  using mode                 = M;
};

template <typename DeviceType, typename ScalarType, typename AViewType, typename BViewType, typename CViewType,
          typename ParamTagType>
struct Functor_BatchedSymm {
  using execution_space      = typename DeviceType::execution_space;
  using member_type          = typename Kokkos::TeamPolicy<execution_space>::member_type;
  using ArgSide              = typename ParamTagType::side;
  using ArgUplo              = typename ParamTagType::uplo;
  using ArgMode              = typename ParamTagType::mode;
  static constexpr bool conj = ParamTagType::conj;

  AViewType m_A;
  BViewType m_B;
  CViewType m_C;
  ScalarType m_alpha, m_beta;

  Functor_BatchedSymm(const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta,
                      const CViewType &C)
      : m_A(A), m_B(B), m_C(C), m_alpha(alpha), m_beta(beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type &member, int &info) const {
    const int k = member.league_rank();

    auto sub_A = Kokkos::subview(m_A, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_B = Kokkos::subview(m_B, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_C = Kokkos::subview(m_C, k, Kokkos::ALL(), Kokkos::ALL());

    if constexpr (std::same_as<ArgMode, KokkosBatched::Mode::Serial>) {
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        info += KokkosBatched::SerialSymm<ArgSide, ArgUplo, conj>::invoke(m_alpha, sub_A, sub_B, m_beta, sub_C);
      });
    } else if constexpr (std::same_as<ArgMode, KokkosBatched::Mode::Team>) {
      info += KokkosBatched::TeamSymm<member_type, ArgSide, ArgUplo, conj>::invoke(member, m_alpha, sub_A, sub_B,
                                                                                   m_beta, sub_C);
    } else if constexpr (std::same_as<ArgMode, KokkosBatched::Mode::TeamVector>) {
      info += KokkosBatched::TeamVectorSymm<member_type, ArgSide, ArgUplo, conj>::invoke(member, m_alpha, sub_A, sub_B,
                                                                                         m_beta, sub_C);
    }
  }

  inline int run() {
    using value_type        = typename AViewType::non_const_value_type;
    std::string name_region = std::same_as<ArgMode, KokkosBatched::Mode::Serial> ? "KokkosBatched::Test::SerialSymm"
                              : std::same_as<ArgMode, KokkosBatched::Mode::Team>
                                  ? "KokkosBatched::Test::TeamSymm"
                                  : "KokkosBatched::Test::TeamVectorSymm";
    const std::string name_value_type = KokkosKernels::ArithTraits<value_type>::name();
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

template <typename DeviceType, typename ArgMode, typename ScalarType, typename AViewType, typename BViewType,
          typename CViewType>
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

    using NoTrans = KokkosBatched::Trans::NoTranspose;
    KokkosBatched::SerialGemm<NoTrans, NoTrans, KokkosBatched::Algo::Gemm::Unblocked>::invoke(m_alpha, aa, bb, m_beta,
                                                                                              cc);
  }

  inline void run() {
    using value_type        = typename CViewType::non_const_value_type;
    std::string name_region = std::same_as<ArgMode, KokkosBatched::Mode::Serial> ? "KokkosBatched::Test::SerialSymm"
                              : std::same_as<ArgMode, KokkosBatched::Mode::Team>
                                  ? "KokkosBatched::Test::TeamSymm"
                                  : "KokkosBatched::Test::TeamVectorSymm";
    const std::string name_value_type = KokkosKernels::ArithTraits<value_type>::name();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());

    Kokkos::RangePolicy<execution_space> policy(0, m_b.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

/// \brief Implementation details of batched symm analytical test
///        to confirm
///        C := alpha*A*B + beta*C (Side::Left)
///        C := alpha*B*A + beta*C (Side::Right)
///        is computed correctly
///        alpha = 1.5, beta = 1.2
///        3x3 matrix
///        A: [[4, 1, 2],
///            [-1, 5, 3],
///            [7, 9, 6]]
///
///        Side::Left
///        B: [[1, 2],      C: [[1, 0],
///            [3, 4],          [0, 1],
///            [5, 6]]          [1, 1]]
///
///        Lower                   Upper
///        Ref: [[55.2, 69],     [[26.7, 36],
///              [88.5, 109.2],   [46.5, 61.2],
///              [97.2, 130.2]]   [62.7, 79.2]]
///
///        Side::Right
///        B: [[7, 8, 9],      C: [[1, 0, 0],
///            [10, 11, 12]]       [0, 1, 0]]
///
///        Lower                        Upper
///        Ref: [[125.7, 171, 262.5],     [[82.2, 111, 138],
///              [169.5, 230.7, 361.5]]    [112.5, 152.7, 187.5]]
///
/// \param Nb [in] Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_symm_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using StridedView3DType = Kokkos::View<ScalarType ***, Kokkos::LayoutStride, DeviceType>;
  using ArgSide           = typename ParamTagType::side;
  using ArgUplo           = typename ParamTagType::uplo;

  const std::size_t M = 2, N = 3;
  const std::size_t m = std::same_as<ArgSide, KokkosBatched::Side::Left> ? N : M;
  const std::size_t n = std::same_as<ArgSide, KokkosBatched::Side::Left> ? M : N;
  View3DType A("A", Nb, N, N), B("B", Nb, m, n), C("C", Nb, m, n), C_ref("C_ref", Nb, m, n);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, N, Nb * incx, N, Nb * incx * N},
      layout2{Nb, incx, m, Nb * incx, n, Nb * incx * m};
  StridedView3DType A_s("A_s", layout), B_s("B_s", layout2), C_s("C_s", layout2);

  // Initialization on host
  auto h_A     = Kokkos::create_mirror_view(A);
  auto h_B     = Kokkos::create_mirror_view(B);
  auto h_C     = Kokkos::create_mirror_view(C);
  auto h_C_ref = Kokkos::create_mirror_view(C_ref);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    h_A(ib, 0, 0) = 4;
    h_A(ib, 0, 1) = 1;
    h_A(ib, 0, 2) = 2;
    h_A(ib, 1, 0) = -1;
    h_A(ib, 1, 1) = 5;
    h_A(ib, 1, 2) = 3;
    h_A(ib, 2, 0) = 7;
    h_A(ib, 2, 1) = 9;
    h_A(ib, 2, 2) = 6;

    if (std::same_as<ArgSide, KokkosBatched::Side::Left>) {
      h_B(ib, 0, 0) = 1;
      h_B(ib, 0, 1) = 2;
      h_B(ib, 1, 0) = 3;
      h_B(ib, 1, 1) = 4;
      h_B(ib, 2, 0) = 5;
      h_B(ib, 2, 1) = 6;

      h_C(ib, 0, 0) = 1;
      h_C(ib, 0, 1) = 0;
      h_C(ib, 1, 0) = 0;
      h_C(ib, 1, 1) = 1;
      h_C(ib, 2, 0) = 1;
      h_C(ib, 2, 1) = 1;

      if (std::same_as<ArgUplo, KokkosBatched::Uplo::Lower>) {
        h_C_ref(ib, 0, 0) = 55.2;
        h_C_ref(ib, 0, 1) = 69;
        h_C_ref(ib, 1, 0) = 88.5;
        h_C_ref(ib, 1, 1) = 109.2;
        h_C_ref(ib, 2, 0) = 97.2;
        h_C_ref(ib, 2, 1) = 130.2;
      } else {
        h_C_ref(ib, 0, 0) = 26.7;
        h_C_ref(ib, 0, 1) = 36;
        h_C_ref(ib, 1, 0) = 46.5;
        h_C_ref(ib, 1, 1) = 61.2;
        h_C_ref(ib, 2, 0) = 62.7;
        h_C_ref(ib, 2, 1) = 79.2;
      }
    } else {
      h_B(ib, 0, 0) = 7;
      h_B(ib, 0, 1) = 8;
      h_B(ib, 0, 2) = 9;
      h_B(ib, 1, 0) = 10;
      h_B(ib, 1, 1) = 11;
      h_B(ib, 1, 2) = 12;

      h_C(ib, 0, 0) = 1;
      h_C(ib, 0, 1) = 0;
      h_C(ib, 0, 2) = 0;
      h_C(ib, 1, 0) = 0;
      h_C(ib, 1, 1) = 1;
      h_C(ib, 1, 2) = 0;

      if (std::same_as<ArgUplo, KokkosBatched::Uplo::Lower>) {
        h_C_ref(ib, 0, 0) = 125.7;
        h_C_ref(ib, 0, 1) = 171;
        h_C_ref(ib, 0, 2) = 262.5;
        h_C_ref(ib, 1, 0) = 169.5;
        h_C_ref(ib, 1, 1) = 230.7;
        h_C_ref(ib, 1, 2) = 361.5;
      } else {
        h_C_ref(ib, 0, 0) = 82.2;
        h_C_ref(ib, 0, 1) = 111;
        h_C_ref(ib, 0, 2) = 138;
        h_C_ref(ib, 1, 0) = 112.5;
        h_C_ref(ib, 1, 1) = 152.7;
        h_C_ref(ib, 1, 2) = 187.5;
      }
    }
  }

  Kokkos::deep_copy(A, h_A);
  Kokkos::deep_copy(B, h_B);
  Kokkos::deep_copy(C, h_C);

  // Deep copy to strided views
  Kokkos::deep_copy(A_s, A);
  Kokkos::deep_copy(B_s, B);
  Kokkos::deep_copy(C_s, C);
  const ScalarType alpha = 1.5, beta = 1.2;

  auto info = Functor_BatchedSymm<DeviceType, ScalarType, View3DType, View3DType, View3DType, ParamTagType>(alpha, A, B,
                                                                                                            beta, C)
                  .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // With strided views
  info = Functor_BatchedSymm<DeviceType, ScalarType, StridedView3DType, StridedView3DType, StridedView3DType,
                             ParamTagType>(alpha, A_s, B_s, beta, C_s)
             .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  RealType eps = 1.0e2 * ats::epsilon();
  Kokkos::deep_copy(h_C, C);

  // Check if C := alpha*A*B + beta*C or alpha*B*A + beta*C is correct
  for (std::size_t ib = 0; ib < h_C.extent(0); ib++) {
    for (std::size_t i = 0; i < h_C.extent(1); i++) {
      for (std::size_t j = 0; j < h_C.extent(2); j++) {
        EXPECT_NEAR_KK(h_C(ib, i, j), h_C_ref(ib, i, j), eps);
      }
    }
  }

  // Testing for strided views, reusing C
  Kokkos::deep_copy(C, C_s);
  Kokkos::deep_copy(h_C, C);

  // Check if C := alpha*A*B + beta*C or alpha*B*A + beta*C is correct
  for (std::size_t ib = 0; ib < h_C.extent(0); ib++) {
    for (std::size_t i = 0; i < h_C.extent(1); i++) {
      for (std::size_t j = 0; j < h_C.extent(2); j++) {
        EXPECT_NEAR_KK(h_C(ib, i, j), h_C_ref(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched symm test
///
/// \param Nb [in] Batch size of matrices
/// \param[in] M Number of rows of matrices
/// \param[in] N Number of columns of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_symm(const std::size_t Nb, const std::size_t M, const std::size_t N) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgSide    = typename ParamTagType::side;
  using ArgUplo    = typename ParamTagType::uplo;
  using ArgMode    = typename ParamTagType::mode;

  // A is N*N matrix
  // B is N*M matrix (Side::Left) or M*N matrix (Side::Right)
  // C is N*M matrix (Side::Left) or M*N matrix (Side::Right)

  const std::size_t m = std::same_as<ArgSide, KokkosBatched::Side::Left> ? N : M;
  const std::size_t n = std::same_as<ArgSide, KokkosBatched::Side::Left> ? M : N;

  View3DType A("A", Nb, N, N), A0("A0", Nb, N, N), B("B", Nb, m, n), C("C", Nb, m, n), C_ref("C_ref", Nb, m, n);

  // Create a random matrix A and x
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A0, rand_pool, randStart, randEnd);
  Kokkos::fill_random(B, rand_pool, randStart, randEnd);
  Kokkos::fill_random(C, rand_pool, randStart, randEnd);

  // Make A Symmetric or Hermitian
  using Trans =
      std::conditional_t<ParamTagType::conj, KokkosBatched::Trans::ConjTranspose, KokkosBatched::Trans::Transpose>;
  dense_to_symmetric<ArgUplo, Trans>(A0, A);
  Kokkos::deep_copy(C_ref, C);

  const RealType alpha = 1.5, beta = 1.2;
  auto info0 = Functor_BatchedSymm<DeviceType, ScalarType, View3DType, View3DType, View3DType, ParamTagType>(alpha, A,
                                                                                                             B, beta, C)
                   .run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);

  // Make a reference with GEMM
  if constexpr (std::same_as<ArgSide, KokkosBatched::Side::Left>) {
    // C_ref = alpha * A * B + beta * C
    Functor_BatchedGemm<DeviceType, ArgMode, ScalarType, View3DType, View3DType, View3DType>(alpha, A, B, beta, C_ref)
        .run();
  } else {
    // C_ref = alpha * B * A + beta * C
    Functor_BatchedGemm<DeviceType, ArgMode, ScalarType, View3DType, View3DType, View3DType>(alpha, B, A, beta, C_ref)
        .run();
  }

  RealType eps = 1.0e2 * ats::epsilon();

  auto h_C     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C);
  auto h_C_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, C_ref);

  // Check if C := alpha*A*B + beta*C or alpha*B*A + beta*C is correct
  for (std::size_t ib = 0; ib < h_C.extent(0); ib++) {
    for (std::size_t i = 0; i < h_C.extent(1); i++) {
      for (std::size_t j = 0; j < h_C.extent(2); j++) {
        EXPECT_NEAR_KK(h_C(ib, i, j), h_C_ref(ib, i, j), eps);
      }
    }
  }
}

}  // namespace Symm
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType, typename LayoutType>
void test_batched_symm_layout() {
  Test::Symm::impl_test_batched_symm_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
  Test::Symm::impl_test_batched_symm_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      Test::Symm::impl_test_batched_symm<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i, j);
      Test::Symm::impl_test_batched_symm<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i, j);
    }
  }
}

template <typename DeviceType, typename ScalarType, typename ParamTagType>
int test_batched_symm() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  test_batched_symm_layout<DeviceType, ScalarType, ParamTagType, Kokkos::LayoutLeft>();
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  test_batched_symm_layout<DeviceType, ScalarType, ParamTagType, Kokkos::LayoutRight>();
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_serial_symm_f_l_l_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_l_u_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_r_l_float) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_r_u_float) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_l_l_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_l_u_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_r_l_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_r_u_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_l_l_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_l_u_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_r_l_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_r_u_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_l_l_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_l_u_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_r_l_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_r_u_float) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_l_l_float) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_l_u_float) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_r_l_float) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_r_u_float) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_l_l_float) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_l_u_float) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_r_l_float) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_r_u_float) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_serial_symm_f_l_l_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_l_u_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_r_l_double) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_r_u_double) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_l_l_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_l_u_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_r_l_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_r_u_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_l_l_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_l_u_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_r_l_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_r_u_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_l_l_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_l_u_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_r_l_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_r_u_double) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_l_l_double) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_l_u_double) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_r_l_double) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_r_u_double) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_l_l_double) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_l_u_double) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_r_l_double) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_r_u_double) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, double, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_serial_symm_f_l_l_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_l_u_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_r_l_fcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_r_u_fcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_l_l_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_l_u_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_r_l_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_r_u_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_l_l_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_l_u_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_r_l_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_r_u_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_l_l_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_l_u_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_r_l_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_r_u_fcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_l_l_fcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_l_u_fcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_r_l_fcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_r_u_fcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_l_l_fcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_l_u_fcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_r_l_fcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_r_u_fcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_serial_symm_f_l_l_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_l_u_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_r_l_dcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_f_r_u_dcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_l_l_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_l_u_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_r_l_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symm_t_r_u_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Serial>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_l_l_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_l_u_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_r_l_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_f_r_u_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_l_l_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_l_u_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_r_l_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symm_t_r_u_dcomplex) {
  using param_tag_type =
      ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper, KokkosBatched::Mode::Team>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_l_l_dcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_l_u_dcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_r_l_dcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_f_r_u_dcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<false, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_l_l_dcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_l_u_dcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_r_l_dcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symm_t_r_u_dcomplex) {
  using param_tag_type = ::Test::Symm::ParamTag<true, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symm<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
#endif
