// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Trmv.hpp>
#include <KokkosBlas2_gemv.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Trmv {

template <typename U, typename T, typename D, typename M>
struct ParamTag {
  using uplo  = U;
  using trans = T;
  using diag  = D;
  using mode  = M;
};

template <typename DeviceType, typename AViewType, typename XViewType, typename ScalarType, typename ParamTagType>
struct Functor_BatchedTrmv {
  using execution_space = typename DeviceType::execution_space;
  using member_type     = typename Kokkos::TeamPolicy<execution_space>::member_type;
  using ArgMode         = typename ParamTagType::mode;
  using ArgUplo         = typename ParamTagType::uplo;
  using ArgTrans        = typename ParamTagType::trans;
  using ArgDiag         = typename ParamTagType::diag;
  AViewType m_A;
  XViewType m_x;
  ScalarType m_alpha;

  Functor_BatchedTrmv(const ScalarType alpha, const AViewType &A, const XViewType &x)
      : m_A(A), m_x(x), m_alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type &member, int &info) const {
    const int k = member.league_rank();
    auto sub_A  = Kokkos::subview(m_A, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_x  = Kokkos::subview(m_x, k, Kokkos::ALL());

    if constexpr (std::same_as<ArgMode, KokkosBatched::Mode::Serial>) {
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        info += KokkosBatched::SerialTrmv<ArgUplo, ArgTrans, ArgDiag>::invoke(m_alpha, sub_A, sub_x);
      });
    }
  }

  inline int run() {
    using value_type                  = typename AViewType::non_const_value_type;
    std::string name_region           = "KokkosBatched::Test::SerialTrmv";
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

template <typename DeviceType, typename AViewType, typename xViewType, typename ScalarType, typename yViewType,
          typename ParamTagType>
struct Functor_BatchedGemv {
  using execution_space = typename DeviceType::execution_space;
  using ArgTrans        = typename ParamTagType::trans;
  AViewType m_a;
  xViewType m_x;
  yViewType m_y;
  ScalarType m_alpha;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedGemv(const ScalarType alpha, const AViewType &a, const xViewType &x, const yViewType &y)
      : m_a(a), m_x(x), m_y(y), m_alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto aa = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto xx = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto yy = Kokkos::subview(m_y, k, Kokkos::ALL());

    ScalarType beta(0);
    KokkosBlas::SerialGemv<ArgTrans, Algo::Gemv::Unblocked>::invoke(m_alpha, aa, xx, beta, yy);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGemv");
    const std::string name_value_type = TestUtils::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_x.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

/// \brief Implementation details of batched trmv analytical test
///        to confirm x:= alpha * op(A) * x
///        alpha = 1.5
///        4x4 matrix (upper)
///        U: [[1, -3, -2,  0],
///            [0,  3, -1, -2],
///            [0,  0,  9,  5],
///            [0,  0,  0, 27]]
///        x: [1, 2, 3, 4]
///        Ref:
///        (Upper/NoTrans/NoUnit) [-16.5, -7.5, 70.5, 162]
///        (Upper/NoTrans/Unit)   [-16.5, -13.5, 34.5, 6]
///        (Upper/Trans/NoUnit)   [1.5, 4.5, 34.5, 178.5]
///        (Upper/Trans/Unit)     [1.5, -1.5, -1.5, 22.5]
///
///        4x4 matrix (lower)
///        L: [[1,  0,  0,  0],
///            [3,  3,  0,  0],
///            [2,  1,  9,  0],
///            [0,  2, -5,  27]]
///        x: [1, 2, 3, 4]
///        Ref:
///        (Lower/NoTrans/NoUnit) [1.5,  13.5, 46.5, 145.5]
///        (Lower/NoTrans/Unit)   [1.5,   7.5, 10.5, -10.5]
///        (Lower/Trans/NoUnit)   [19.5, 25.5, 10.5, 162]
///        (Lower/Trans/Unit)     [19.5, 19.5, -25.5, 6]
///
/// \param Nb [in] Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_trmv_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo           = typename ParamTagType::uplo;
  using ArgTrans          = typename ParamTagType::trans;
  using ArgDiag           = typename ParamTagType::diag;

  const std::size_t BlkSize = 4;
  View3DType A("A", Nb, BlkSize, BlkSize);
  View2DType x("x", Nb, BlkSize), x_ref("x_ref", Nb, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, BlkSize, Nb * incx};
  StridedView2DType x_s("x_s", layout);

  // Only filling x2, A2 and deep_copy from its subview
  auto h_A     = Kokkos::create_mirror_view(A);
  auto h_x     = Kokkos::create_mirror_view(x);
  auto h_x_ref = Kokkos::create_mirror_view(x);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    if (std::same_as<ArgUplo, KokkosBatched::Uplo::Upper>) {
      h_A(ib, 0, 0) = 1;
      h_A(ib, 0, 1) = -3;
      h_A(ib, 0, 2) = -2;
      h_A(ib, 0, 3) = 0;
      h_A(ib, 1, 1) = 3;
      h_A(ib, 1, 2) = -1;
      h_A(ib, 1, 3) = -2;
      h_A(ib, 2, 2) = 9;
      h_A(ib, 2, 3) = 5;
      h_A(ib, 3, 3) = 27;
      if (std::same_as<ArgTrans, KokkosBatched::Trans::NoTranspose>) {
        if (std::same_as<ArgDiag, KokkosBatched::Diag::NonUnit>) {
          h_x_ref(ib, 0) = -16.5;
          h_x_ref(ib, 1) = -7.5;
          h_x_ref(ib, 2) = 70.5;
          h_x_ref(ib, 3) = 162;
        } else {
          h_x_ref(ib, 0) = -16.5;
          h_x_ref(ib, 1) = -13.5;
          h_x_ref(ib, 2) = 34.5;
          h_x_ref(ib, 3) = 6;
        }
      } else {
        if (std::same_as<ArgDiag, KokkosBatched::Diag::NonUnit>) {
          h_x_ref(ib, 0) = 1.5;
          h_x_ref(ib, 1) = 4.5;
          h_x_ref(ib, 2) = 34.5;
          h_x_ref(ib, 3) = 178.5;
        } else {
          h_x_ref(ib, 0) = 1.5;
          h_x_ref(ib, 1) = -1.5;
          h_x_ref(ib, 2) = -1.5;
          h_x_ref(ib, 3) = 22.5;
        }
      }
    } else {
      h_A(ib, 0, 0) = 1;
      h_A(ib, 1, 0) = 3;
      h_A(ib, 1, 1) = 3;
      h_A(ib, 2, 0) = 2;
      h_A(ib, 2, 1) = 1;
      h_A(ib, 2, 2) = 9;
      h_A(ib, 3, 0) = 0;
      h_A(ib, 3, 1) = 2;
      h_A(ib, 3, 2) = -5;
      h_A(ib, 3, 3) = 27;
      if (std::same_as<ArgTrans, KokkosBatched::Trans::NoTranspose>) {
        if (std::same_as<ArgDiag, KokkosBatched::Diag::NonUnit>) {
          h_x_ref(ib, 0) = 1.5;
          h_x_ref(ib, 1) = 13.5;
          h_x_ref(ib, 2) = 46.5;
          h_x_ref(ib, 3) = 145.5;
        } else {
          h_x_ref(ib, 0) = 1.5;
          h_x_ref(ib, 1) = 7.5;
          h_x_ref(ib, 2) = 10.5;
          h_x_ref(ib, 3) = -10.5;
        }
      } else {
        if (std::same_as<ArgDiag, KokkosBatched::Diag::NonUnit>) {
          h_x_ref(ib, 0) = 19.5;
          h_x_ref(ib, 1) = 25.5;
          h_x_ref(ib, 2) = 10.5;
          h_x_ref(ib, 3) = 162;
        } else {
          h_x_ref(ib, 0) = 19.5;
          h_x_ref(ib, 1) = 19.5;
          h_x_ref(ib, 2) = -25.5;
          h_x_ref(ib, 3) = 6;
        }
      }
    }

    for (std::size_t j = 0; j < BlkSize; j++) {
      h_x(ib, j) = static_cast<ScalarType>(j + 1);
    }
  }

  Kokkos::deep_copy(A, h_A);
  Kokkos::deep_copy(x, h_x);

  // Deep copy to strided views
  Kokkos::deep_copy(x_s, x);

  const ScalarType alpha = 1.5;

  auto info = Functor_BatchedTrmv<DeviceType, View3DType, View2DType, ScalarType, ParamTagType>(alpha, A, x).run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // With strided views
  info = Functor_BatchedTrmv<DeviceType, View3DType, StridedView2DType, ScalarType, ParamTagType>(alpha, A, x_s).run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  RealType eps = 1.0e1 * ats::epsilon();
  Kokkos::deep_copy(h_x, x);

  // Check if x:= alpha * op(A) * x
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }

  // Testing for strided views, reusing x
  Kokkos::deep_copy(x, x_s);
  Kokkos::deep_copy(h_x, x);

  // Check if x:= alpha * op(A) * x
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched trmv test
///
/// \param Nb [in] Batch size of matrices
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_trmv(const std::size_t Nb, const std::size_t BlkSize) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View2DType = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo    = typename ParamTagType::uplo;
  using ArgDiag    = typename ParamTagType::diag;

  View3DType A("A", Nb, BlkSize, BlkSize), A_ref("A_ref", Nb, BlkSize, BlkSize);
  View2DType x("x", Nb, BlkSize), x_ref("x_ref", Nb, BlkSize);

  // Create a random matrix A and x
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);

  // Upper or lower triangular part of A
  create_triangular_matrix<View3DType, View3DType, ArgUplo, ArgDiag>(A, A_ref);

  Kokkos::deep_copy(A, A_ref);

  const ScalarType alpha = 1.5;

  // Make reference with Beta = 0
  Functor_BatchedGemv<DeviceType, View3DType, View2DType, ScalarType, View2DType, ParamTagType>(alpha, A, x, x_ref)
      .run();

  auto info0 = Functor_BatchedTrmv<DeviceType, View3DType, View2DType, ScalarType, ParamTagType>(alpha, A, x).run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);

  RealType eps = 1.0e2 * ats::epsilon();

  auto h_x     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x);
  auto h_x_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x_ref);

  // Check if x:= alpha * op(A) * x
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

}  // namespace Trmv
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType>
int test_batched_trmv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Trmv::impl_test_batched_trmv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Trmv::impl_test_batched_trmv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Trmv::impl_test_batched_trmv<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Trmv::impl_test_batched_trmv<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Trmv::impl_test_batched_trmv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Trmv::impl_test_batched_trmv_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Trmv::impl_test_batched_trmv<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Trmv::impl_test_batched_trmv<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_serial_trmv_l_nt_u_float) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_nt_n_float) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_t_u_float) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_t_n_float) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Diag::NonUnit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_c_u_float) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_c_n_float) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_nt_u_float) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_nt_n_float) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_t_u_float) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_t_n_float) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Diag::NonUnit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_c_u_float) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_c_n_float) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_serial_trmv_l_nt_u_double) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_nt_n_double) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_t_u_double) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_t_n_double) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Diag::NonUnit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_c_u_double) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_c_n_double) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_nt_u_double) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_nt_n_double) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_t_u_double) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_t_n_double) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Diag::NonUnit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_c_u_double) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_c_n_double) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, double, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_serial_trmv_l_nt_u_fcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_nt_n_fcomplex) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_t_u_fcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_t_n_fcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Diag::NonUnit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_c_u_fcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_c_n_fcomplex) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_nt_u_fcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_nt_n_fcomplex) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_t_u_fcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_t_n_fcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Diag::NonUnit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_c_u_fcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_c_n_fcomplex) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_serial_trmv_l_nt_u_dcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_nt_n_dcomplex) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Lower, Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_t_u_dcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_t_n_dcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::Transpose, KokkosBatched::Diag::NonUnit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_c_u_dcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_l_c_n_dcomplex) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Lower, Trans::ConjTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_nt_u_dcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_nt_n_dcomplex) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Upper, Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_t_u_dcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_t_n_dcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::Transpose, KokkosBatched::Diag::NonUnit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_c_u_dcomplex) {
  using param_tag_type =
      ::Test::Trmv::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Diag::Unit, KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_serial_trmv_u_c_n_dcomplex) {
  using param_tag_type = ::Test::Trmv::ParamTag<Uplo::Upper, Trans::ConjTranspose, KokkosBatched::Diag::NonUnit,
                                                KokkosBatched::Mode::Serial>;
  test_batched_trmv<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
#endif
