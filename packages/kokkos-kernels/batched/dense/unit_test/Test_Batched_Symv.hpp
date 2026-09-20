// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <concepts>
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Symv.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Symv {

template <typename U, typename T, typename M>
struct ParamTag {
  using uplo  = U;
  using trans = T;
  using mode  = M;
};

template <typename DeviceType, typename ScalarType, typename AViewType, typename XViewType, typename YViewType,
          typename ParamTagType>
struct Functor_BatchedSymv {
  using execution_space = typename DeviceType::execution_space;
  using member_type     = typename Kokkos::TeamPolicy<execution_space>::member_type;
  using ArgMode         = typename ParamTagType::mode;
  using ArgUplo         = typename ParamTagType::uplo;
  using ArgTrans        = typename ParamTagType::trans;

  AViewType m_A;
  XViewType m_x;
  YViewType m_y;
  ScalarType m_alpha, m_beta;

  Functor_BatchedSymv(const ScalarType alpha, const AViewType &A, const XViewType &x, const ScalarType beta,
                      const YViewType &y)
      : m_A(A), m_x(x), m_y(y), m_alpha(alpha), m_beta(beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type &member, int &info) const {
    const int k = member.league_rank();

    auto sub_A = Kokkos::subview(m_A, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_x = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto sub_y = Kokkos::subview(m_y, k, Kokkos::ALL());

    if constexpr (std::same_as<ArgMode, KokkosBatched::Mode::Serial>) {
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        info += KokkosBatched::SerialSymv<typename ParamTagType::uplo, typename ParamTagType::trans>::invoke(
            m_alpha, sub_A, sub_x, m_beta, sub_y);
      });
    } else if constexpr (std::same_as<ArgMode, KokkosBatched::Mode::Team>) {
      info += KokkosBatched::TeamSymv<member_type, typename ParamTagType::uplo, typename ParamTagType::trans>::invoke(
          member, m_alpha, sub_A, sub_x, m_beta, sub_y);
    } else if constexpr (std::same_as<ArgMode, KokkosBatched::Mode::TeamVector>) {
      info +=
          KokkosBatched::TeamVectorSymv<member_type, typename ParamTagType::uplo, typename ParamTagType::trans>::invoke(
              member, m_alpha, sub_A, sub_x, m_beta, sub_y);
    }
  }

  inline int run() {
    using value_type        = typename AViewType::non_const_value_type;
    std::string name_region = std::same_as<ArgMode, KokkosBatched::Mode::Serial> ? "KokkosBatched::Test::SerialSymv"
                              : std::same_as<ArgMode, KokkosBatched::Mode::Team>
                                  ? "KokkosBatched::Test::TeamSymv"
                                  : "KokkosBatched::Test::TeamVectorSymv";
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

/// \brief Implementation details of batched symv analytical test
///        to confirm y := alpha*op( A )*x + beta*y is computed correctly
///        alpha = 1.5, beta = 1.2
///        4x4 matrix
///        A: [[1,  -3, -2,  0],
///            [3,  3, -1, -2],
///            [2,  1,  9, 5],
///            [0,  2, -5, 27]]
///
///        x: [1, 2, 3, 4], y: [5, 6, 7, 8]
///
///        Upper
///        Ref: [-10.5  -4.8  72.9 188.1]
///
///        Lower
///        Ref: [ 25.5  37.2  24.9 155.1]
///
/// \param Nb [in] Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_symv_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo           = typename ParamTagType::uplo;

  const std::size_t BlkSize = 4;
  View3DType A("A", Nb, BlkSize, BlkSize);
  View2DType x("x", Nb, BlkSize), y("y", Nb, BlkSize), y_ref("y_ref", Nb, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, BlkSize, Nb * incx};
  StridedView2DType x_s("x_s", layout), y_s("y_s", layout);

  // Only filling x2, A2 and deep_copy from its subview
  auto h_A     = Kokkos::create_mirror_view(A);
  auto h_x     = Kokkos::create_mirror_view(x);
  auto h_y     = Kokkos::create_mirror_view(y);
  auto h_y_ref = Kokkos::create_mirror_view(y_ref);

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

    if (std::is_same_v<ArgUplo, KokkosBatched::Uplo::Upper>) {
      h_y_ref(ib, 0) = -10.5;
      h_y_ref(ib, 1) = -4.8;
      h_y_ref(ib, 2) = 72.9;
      h_y_ref(ib, 3) = 188.1;
    } else {
      h_y_ref(ib, 0) = 25.5;
      h_y_ref(ib, 1) = 37.2;
      h_y_ref(ib, 2) = 24.9;
      h_y_ref(ib, 3) = 155.1;
    }

    for (std::size_t j = 0; j < BlkSize; j++) {
      h_x(ib, j) = static_cast<ScalarType>(j + 1);
      h_y(ib, j) = static_cast<ScalarType>(j + 5);
    }
  }

  Kokkos::deep_copy(A, h_A);
  Kokkos::deep_copy(x, h_x);
  Kokkos::deep_copy(y, h_y);

  // Deep copy to strided views
  Kokkos::deep_copy(x_s, x);
  Kokkos::deep_copy(y_s, y);
  const ScalarType alpha = 1.5, beta = 1.2;

  auto info = Functor_BatchedSymv<DeviceType, ScalarType, View3DType, View2DType, View2DType, ParamTagType>(alpha, A, x,
                                                                                                            beta, y)
                  .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // With strided views
  info = Functor_BatchedSymv<DeviceType, ScalarType, View3DType, StridedView2DType, StridedView2DType, ParamTagType>(
             alpha, A, x_s, beta, y_s)
             .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  RealType eps = 1.0e2 * ats::epsilon();
  Kokkos::deep_copy(h_y, y);

  // Check if y:= alpha*op( A )*x + beta*y
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_y(ib, i), h_y_ref(ib, i), eps);
    }
  }

  // Testing for strided views, reusing y
  Kokkos::deep_copy(y, y_s);
  Kokkos::deep_copy(h_y, y);

  // Check if y:= alpha*op( A )*x + beta*y
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_y(ib, i), h_y_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched symv test
///
/// \param Nb [in] Batch size of matrices
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_symv(const std::size_t Nb, const std::size_t BlkSize) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View2DType = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo    = typename ParamTagType::uplo;
  using ArgTrans   = typename ParamTagType::trans;

  View3DType A("A", Nb, BlkSize, BlkSize);
  View2DType x("x", Nb, BlkSize), y("y", Nb, BlkSize), y_alpha0("y_alpha0", Nb, BlkSize),
      y_beta1("y_beta1", Nb, BlkSize), y_ref("y_ref", Nb, BlkSize), y_alpha0_ref("y_alpha0_ref", Nb, BlkSize),
      y_beta1_ref("y_beta1_ref", Nb, BlkSize);

  // Create a random matrix A and x
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);
  Kokkos::fill_random(y, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(y_alpha0, y);
  Kokkos::deep_copy(y_beta1, y);
  Kokkos::deep_copy(y_ref, y);
  Kokkos::deep_copy(y_alpha0_ref, y);
  Kokkos::deep_copy(y_beta1_ref, y);

  const RealType alpha = 1.5, beta = 1.2, zero = KokkosKernels::ArithTraits<RealType>::zero(),
                 one = KokkosKernels::ArithTraits<RealType>::one();
  auto info0 = Functor_BatchedSymv<DeviceType, ScalarType, View3DType, View2DType, View2DType, ParamTagType>(alpha, A,
                                                                                                             x, beta, y)
                   .run();

  // When alpha = 0
  auto info1 = Functor_BatchedSymv<DeviceType, ScalarType, View3DType, View2DType, View2DType, ParamTagType>(
                   zero, A, x, beta, y_alpha0)
                   .run();

  // When beta = 1
  auto info2 = Functor_BatchedSymv<DeviceType, ScalarType, View3DType, View2DType, View2DType, ParamTagType>(
                   alpha, A, x, one, y_beta1)
                   .run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);
  EXPECT_EQ(info2, 0);

  // Make a reference at host
  auto h_x            = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x);
  auto h_A            = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, A);
  auto h_y_ref        = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y_ref);
  auto h_y_alpha0_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y_alpha0_ref);
  auto h_y_beta1_ref  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y_beta1_ref);

  // Note: ConjTranspose corresponds to {c,z}hemv for Hermitian matrix
  using Op     = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  using Sym_op = std::conditional_t<std::same_as<ArgTrans, Trans::ConjTranspose>, KokkosBlas::Impl::OpReal,
                                    KokkosBlas::Impl::OpID>;
  Op op;
  Sym_op sym_op;
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      h_y_ref(ib, i)        = beta * h_y_ref(ib, i);
      h_y_alpha0_ref(ib, i) = beta * h_y_alpha0_ref(ib, i);
    }
    for (std::size_t j = 0; j < BlkSize; j++) {
      auto temp1 = alpha * h_x(ib, j);
      ScalarType temp2(0);
      if (std::is_same_v<ArgUplo, KokkosBatched::Uplo::Upper>) {
        for (std::size_t i = 0; i < j; i++) {
          h_y_ref(ib, i) += temp1 * h_A(ib, i, j);
          h_y_beta1_ref(ib, i) += temp1 * h_A(ib, i, j);
          temp2 += op(h_A(ib, i, j)) * h_x(ib, i);
        }
        h_y_ref(ib, j) += temp1 * sym_op(h_A(ib, j, j)) + alpha * temp2;
        h_y_beta1_ref(ib, j) += temp1 * sym_op(h_A(ib, j, j)) + alpha * temp2;
      } else {
        h_y_ref(ib, j) += temp1 * sym_op(h_A(ib, j, j));
        h_y_beta1_ref(ib, j) += temp1 * sym_op(h_A(ib, j, j));
        for (std::size_t i = j + 1; i < BlkSize; i++) {
          h_y_ref(ib, i) += temp1 * h_A(ib, i, j);
          h_y_beta1_ref(ib, i) += temp1 * h_A(ib, i, j);
          temp2 += op(h_A(ib, i, j)) * h_x(ib, i);
        }
        h_y_ref(ib, j) += alpha * temp2;
        h_y_beta1_ref(ib, j) += alpha * temp2;
      }
    }
  }

  RealType eps = 1.0e2 * ats::epsilon();

  auto h_y        = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y);
  auto h_y_alpha0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y_alpha0);
  auto h_y_beta1  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y_beta1);

  // Check if y := alpha*op( A )*x + beta*y
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t j = 0; j < BlkSize; j++) {
      EXPECT_NEAR_KK(h_y(ib, j), h_y_ref(ib, j), eps);
      EXPECT_NEAR_KK(h_y_alpha0(ib, j), h_y_alpha0_ref(ib, j), eps);
      EXPECT_NEAR_KK(h_y_beta1(ib, j), h_y_beta1_ref(ib, j), eps);
    }
  }
}

}  // namespace Symv
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType>
int test_batched_symv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Symv::impl_test_batched_symv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Symv::impl_test_batched_symv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Symv::impl_test_batched_symv<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Symv::impl_test_batched_symv<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Symv::impl_test_batched_symv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Symv::impl_test_batched_symv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Symv::impl_test_batched_symv<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Symv::impl_test_batched_symv<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_serial_symv_l_t_float) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_l_c_float) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_u_t_float) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_u_c_float) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_l_t_float) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_l_c_float) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_u_t_float) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_u_c_float) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_l_t_float) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_l_c_float) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_u_t_float) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_u_c_float) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_serial_symv_l_t_double) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_l_c_double) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_u_t_double) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_u_c_double) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_l_t_double) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_l_c_double) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_u_t_double) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_u_c_double) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_l_t_double) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_l_c_double) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_u_t_double) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_u_c_double) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, double, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_serial_symv_l_t_fcomplex) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_l_c_fcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_u_t_fcomplex) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_u_c_fcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_l_t_fcomplex) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_l_c_fcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_u_t_fcomplex) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_u_c_fcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_l_t_fcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_l_c_fcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_u_t_fcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_u_c_fcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_serial_symv_l_t_dcomplex) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_l_c_dcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_u_t_dcomplex) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_symv_u_c_dcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Serial>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_l_t_dcomplex) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_l_c_dcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_u_t_dcomplex) {
  using param_tag_type =
      ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose, KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_symv_u_c_dcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::Team>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_l_t_dcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::Transpose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_l_c_dcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_u_t_dcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::Transpose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_team_vector_symv_u_c_dcomplex) {
  using param_tag_type = ::Test::Symv::ParamTag<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::ConjTranspose,
                                                KokkosBatched::Mode::TeamVector>;
  test_batched_symv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
#endif
