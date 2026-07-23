// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Swap.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Swap {

template <typename DeviceType, typename XViewType, typename YViewType, typename ArgMode>
struct Functor_BatchedSwap {
  using execution_space = typename DeviceType::execution_space;
  using member_type     = typename Kokkos::TeamPolicy<execution_space>::member_type;
  XViewType m_x;
  YViewType m_y;

  Functor_BatchedSwap(const XViewType &x, const YViewType &y) : m_x(x), m_y(y) {}

  template <typename ViewType>
  KOKKOS_INLINE_FUNCTION auto slice(const ViewType &view, const int k) const {
    if constexpr (ViewType::rank() == 3) {
      return Kokkos::subview(view, k, Kokkos::ALL(), Kokkos::ALL());
    } else {
      return Kokkos::subview(view, k, Kokkos::ALL());
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(const member_type &member, int &info) const {
    const int k = member.league_rank();
    auto sub_x  = slice(m_x, k);
    auto sub_y  = slice(m_y, k);
    if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
      Kokkos::single(Kokkos::PerTeam(member), [&]() { info += KokkosBatched::SerialSwap::invoke(sub_x, sub_y); });
    } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Team>) {
      info += KokkosBatched::TeamSwap<member_type>::invoke(member, sub_x, sub_y);
    } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::TeamVector>) {
      info += KokkosBatched::TeamVectorSwap<member_type>::invoke(member, sub_x, sub_y);
    }
  }

  inline int run() {
    using value_type        = typename XViewType::non_const_value_type;
    std::string name_region = std::is_same_v<ArgMode, KokkosBatched::Mode::Serial> ? "KokkosBatched::Test::SerialSwap"
                              : std::is_same_v<ArgMode, KokkosBatched::Mode::Team>
                                  ? "KokkosBatched::Test::TeamSwap"
                                  : "KokkosBatched::Test::TeamVectorSwap";
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    const int league_size = m_x.extent_int(0);

    Kokkos::TeamPolicy<execution_space> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);

    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

/// \brief Implementation details of batched swap analytical test
/// x0 = [1, 3, 5], y0 = [2, 4, 6]
/// x1 = [[1, 3, 5], [7, 9, 11]], y1 = [[2, 4, 6], [8, 10, 12]]
///
/// \tparam DeviceType Kokkos device type
/// \tparam ScalarType Kokkos scalar type
/// \tparam LayoutType1 Kokkos layout type for the views
/// \tparam LayoutType2 Kokkos layout type for the strided views
/// \tparam ArgMode: one of Mode::Serial, Mode::Team, Mode::TeamVector
///
/// \param[in] Nb Batch size of vectors
template <typename DeviceType, typename ScalarType, typename LayoutType1, typename LayoutType2, typename ArgMode>
void impl_test_batched_swap_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using XView2DType       = Kokkos::View<ScalarType **, LayoutType1, DeviceType>;
  using YView2DType       = Kokkos::View<ScalarType **, LayoutType2, DeviceType>;
  using XView3DType       = Kokkos::View<ScalarType ***, LayoutType1, DeviceType>;
  using YView3DType       = Kokkos::View<ScalarType ***, LayoutType2, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using StridedView3DType = Kokkos::View<ScalarType ***, Kokkos::LayoutStride, DeviceType>;

  const std::size_t M = 2, N = 3;
  XView2DType x0("x0", Nb, N), x0_ref("x0_ref", Nb, N);
  YView2DType y0("y0", Nb, N), y0_ref("y0_ref", Nb, N);
  XView3DType x1("x1", Nb, M, N), x1_ref("x1_ref", Nb, M, N);
  YView3DType y1("y1", Nb, M, N), y1_ref("y1_ref", Nb, M, N);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout0{Nb, incx, N, Nb * incx};
  StridedView2DType x0_s("x0_s", layout0), y0_s("y0_s", layout0);

  Kokkos::LayoutStride layout1{Nb, incx, M, Nb * incx, N, Nb * incx * M};
  StridedView3DType x1_s("x1_s", layout1), y1_s("y1_s", layout1);

  auto h_x0 = Kokkos::create_mirror_view(x0);
  auto h_y0 = Kokkos::create_mirror_view(y0);
  auto h_x1 = Kokkos::create_mirror_view(x1);
  auto h_y1 = Kokkos::create_mirror_view(y1);

  constexpr bool is_complex = KokkosKernels::ArithTraits<ScalarType>::is_complex;

  for (std::size_t ib = 0; ib < Nb; ib++) {
    if constexpr (is_complex) {
      h_x0(ib, 0) = ScalarType(1, 7);
      h_x0(ib, 1) = ScalarType(3, 9);
      h_x0(ib, 2) = ScalarType(5, 11);

      h_y0(ib, 0) = ScalarType(2, 8);
      h_y0(ib, 1) = ScalarType(4, 10);
      h_y0(ib, 2) = ScalarType(6, 12);

      h_x1(ib, 0, 0) = ScalarType(1, 7);
      h_x1(ib, 0, 1) = ScalarType(3, 9);
      h_x1(ib, 0, 2) = ScalarType(5, 11);
      h_x1(ib, 1, 0) = ScalarType(13, 19);
      h_x1(ib, 1, 1) = ScalarType(15, 21);
      h_x1(ib, 1, 2) = ScalarType(17, 23);

      h_y1(ib, 0, 0) = ScalarType(2, 8);
      h_y1(ib, 0, 1) = ScalarType(4, 10);
      h_y1(ib, 0, 2) = ScalarType(6, 12);
      h_y1(ib, 1, 0) = ScalarType(14, 20);
      h_y1(ib, 1, 1) = ScalarType(16, 22);
      h_y1(ib, 1, 2) = ScalarType(18, 24);
    } else {
      h_x0(ib, 0) = ScalarType(1);
      h_x0(ib, 1) = ScalarType(3);
      h_x0(ib, 2) = ScalarType(5);

      h_y0(ib, 0) = ScalarType(2);
      h_y0(ib, 1) = ScalarType(4);
      h_y0(ib, 2) = ScalarType(6);

      h_x1(ib, 0, 0) = ScalarType(1);
      h_x1(ib, 0, 1) = ScalarType(3);
      h_x1(ib, 0, 2) = ScalarType(5);
      h_x1(ib, 1, 0) = ScalarType(13);
      h_x1(ib, 1, 1) = ScalarType(15);
      h_x1(ib, 1, 2) = ScalarType(17);

      h_y1(ib, 0, 0) = ScalarType(2);
      h_y1(ib, 0, 1) = ScalarType(4);
      h_y1(ib, 0, 2) = ScalarType(6);
      h_y1(ib, 1, 0) = ScalarType(14);
      h_y1(ib, 1, 1) = ScalarType(16);
      h_y1(ib, 1, 2) = ScalarType(18);
    }
  }

  Kokkos::deep_copy(x0, h_x0);
  Kokkos::deep_copy(y0, h_y0);
  Kokkos::deep_copy(x1, h_x1);
  Kokkos::deep_copy(y1, h_y1);

  // Deep copy to strided views
  Kokkos::deep_copy(x0_s, x0);
  Kokkos::deep_copy(y0_s, y0);
  Kokkos::deep_copy(x1_s, x1);
  Kokkos::deep_copy(y1_s, y1);

  // Reference results after swap
  Kokkos::deep_copy(x0_ref, y0);
  Kokkos::deep_copy(y0_ref, x0);
  Kokkos::deep_copy(x1_ref, y1);
  Kokkos::deep_copy(y1_ref, x1);

  auto info0 = Functor_BatchedSwap<DeviceType, XView2DType, YView2DType, ArgMode>(x0, y0).run();
  auto info1 = Functor_BatchedSwap<DeviceType, XView3DType, YView3DType, ArgMode>(x1, y1).run();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  // With strided views
  info0 = Functor_BatchedSwap<DeviceType, StridedView2DType, StridedView2DType, ArgMode>(x0_s, y0_s).run();
  info1 = Functor_BatchedSwap<DeviceType, StridedView3DType, StridedView3DType, ArgMode>(x1_s, y1_s).run();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  RealType eps = 1.0e1 * ats::epsilon();
  Kokkos::deep_copy(h_x0, x0);
  Kokkos::deep_copy(h_y0, y0);
  Kokkos::deep_copy(h_x1, x1);
  Kokkos::deep_copy(h_y1, y1);
  auto h_x0_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x0_ref);
  auto h_y0_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y0_ref);
  auto h_x1_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x1_ref);
  auto h_y1_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y1_ref);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < N; i++) {
      KK_EXPECT_NEAR(h_x0(ib, i), h_x0_ref(ib, i), eps);
      KK_EXPECT_NEAR(h_y0(ib, i), h_y0_ref(ib, i), eps);
    }
    for (std::size_t i = 0; i < M; i++) {
      for (std::size_t j = 0; j < N; j++) {
        KK_EXPECT_NEAR(h_x1(ib, i, j), h_x1_ref(ib, i, j), eps);
        KK_EXPECT_NEAR(h_y1(ib, i, j), h_y1_ref(ib, i, j), eps);
      }
    }
  }

  // Testing for strided views, reusing x and y
  Kokkos::deep_copy(x0, x0_s);
  Kokkos::deep_copy(y0, y0_s);
  Kokkos::deep_copy(x1, x1_s);
  Kokkos::deep_copy(y1, y1_s);
  Kokkos::deep_copy(h_x0, x0);
  Kokkos::deep_copy(h_y0, y0);
  Kokkos::deep_copy(h_x1, x1);
  Kokkos::deep_copy(h_y1, y1);
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < N; i++) {
      KK_EXPECT_NEAR(h_x0(ib, i), h_x0_ref(ib, i), eps);
      KK_EXPECT_NEAR(h_y0(ib, i), h_y0_ref(ib, i), eps);
    }
    for (std::size_t i = 0; i < M; i++) {
      for (std::size_t j = 0; j < N; j++) {
        KK_EXPECT_NEAR(h_x1(ib, i, j), h_x1_ref(ib, i, j), eps);
        KK_EXPECT_NEAR(h_y1(ib, i, j), h_y1_ref(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched swap test
///
/// \tparam DeviceType Kokkos device type
/// \tparam ScalarType Kokkos scalar type
/// \tparam LayoutType1 Kokkos layout type for the X views
/// \tparam LayoutType2 Kokkos layout type for the Y views
/// \tparam ArgMode: one of Mode::Serial, Mode::Team, Mode::TeamVector
///
/// \param[in] Nb Batch size of vectors
/// \param[in] N Length of the vector x
template <typename DeviceType, typename ScalarType, typename LayoutType1, typename LayoutType2, typename ArgMode>
void impl_test_batched_swap(const std::size_t Nb, const std::size_t M, const std::size_t N) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using XView2DType       = Kokkos::View<ScalarType **, LayoutType1, DeviceType>;
  using YView2DType       = Kokkos::View<ScalarType **, LayoutType2, DeviceType>;
  using XView3DType       = Kokkos::View<ScalarType ***, LayoutType1, DeviceType>;
  using YView3DType       = Kokkos::View<ScalarType ***, LayoutType2, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using StridedView3DType = Kokkos::View<ScalarType ***, Kokkos::LayoutStride, DeviceType>;

  XView2DType x0("x0", Nb, N), x0_ref("x0_ref", Nb, N);
  YView2DType y0("y0", Nb, N), y0_ref("y0_ref", Nb, N);
  XView3DType x1("x1", Nb, M, N), x1_ref("x1_ref", Nb, M, N);
  YView3DType y1("y1", Nb, M, N), y1_ref("y1_ref", Nb, M, N);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout0{Nb, incx, N, Nb * incx};
  StridedView2DType x0_s("x0_s", layout0), y0_s("y0_s", layout0);

  Kokkos::LayoutStride layout1{Nb, incx, M, Nb * incx, N, Nb * incx * M};
  StridedView3DType x1_s("x1_s", layout1), y1_s("y1_s", layout1);

  // Create random x
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(x0, rand_pool, randStart, randEnd);
  Kokkos::fill_random(y0, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x1, rand_pool, randStart, randEnd);
  Kokkos::fill_random(y1, rand_pool, randStart, randEnd);

  // Deep copy to strided views
  Kokkos::deep_copy(x0_s, x0);
  Kokkos::deep_copy(y0_s, y0);
  Kokkos::deep_copy(x1_s, x1);
  Kokkos::deep_copy(y1_s, y1);

  // Reference results after swap
  Kokkos::deep_copy(x0_ref, y0);
  Kokkos::deep_copy(y0_ref, x0);
  Kokkos::deep_copy(x1_ref, y1);
  Kokkos::deep_copy(y1_ref, x1);

  auto info0 = Functor_BatchedSwap<DeviceType, XView2DType, YView2DType, ArgMode>(x0, y0).run();
  auto info1 = Functor_BatchedSwap<DeviceType, XView3DType, YView3DType, ArgMode>(x1, y1).run();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  // With strided views
  info0 = Functor_BatchedSwap<DeviceType, StridedView2DType, StridedView2DType, ArgMode>(x0_s, y0_s).run();
  info1 = Functor_BatchedSwap<DeviceType, StridedView3DType, StridedView3DType, ArgMode>(x1_s, y1_s).run();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  RealType eps  = 1.0e1 * ats::epsilon();
  auto h_x0     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x0);
  auto h_y0     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y0);
  auto h_x1     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x1);
  auto h_y1     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y1);
  auto h_x0_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x0_ref);
  auto h_y0_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y0_ref);
  auto h_x1_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x1_ref);
  auto h_y1_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y1_ref);

  // Check if swap is correct
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < N; i++) {
      KK_EXPECT_NEAR(h_x0(ib, i), h_x0_ref(ib, i), eps);
      KK_EXPECT_NEAR(h_y0(ib, i), h_y0_ref(ib, i), eps);
    }
    for (std::size_t i = 0; i < M; i++) {
      for (std::size_t j = 0; j < N; j++) {
        KK_EXPECT_NEAR(h_x1(ib, i, j), h_x1_ref(ib, i, j), eps);
        KK_EXPECT_NEAR(h_y1(ib, i, j), h_y1_ref(ib, i, j), eps);
      }
    }
  }

  // Testing for strided views, reusing x and y
  Kokkos::deep_copy(x0, x0_s);
  Kokkos::deep_copy(y0, y0_s);
  Kokkos::deep_copy(x1, x1_s);
  Kokkos::deep_copy(y1, y1_s);
  Kokkos::deep_copy(h_x0, x0);
  Kokkos::deep_copy(h_y0, y0);
  Kokkos::deep_copy(h_x1, x1);
  Kokkos::deep_copy(h_y1, y1);
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < N; i++) {
      KK_EXPECT_NEAR(h_x0(ib, i), h_x0_ref(ib, i), eps);
      KK_EXPECT_NEAR(h_y0(ib, i), h_y0_ref(ib, i), eps);
    }
    for (std::size_t i = 0; i < M; i++) {
      for (std::size_t j = 0; j < N; j++) {
        KK_EXPECT_NEAR(h_x1(ib, i, j), h_x1_ref(ib, i, j), eps);
        KK_EXPECT_NEAR(h_y1(ib, i, j), h_y1_ref(ib, i, j), eps);
      }
    }
  }
}

}  // namespace Swap
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ArgMode>
int test_batched_swap() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Swap::impl_test_batched_swap_analytical<DeviceType, ScalarType, LayoutType, Kokkos::LayoutLeft, ArgMode>(1);
    Test::Swap::impl_test_batched_swap_analytical<DeviceType, ScalarType, LayoutType, Kokkos::LayoutLeft, ArgMode>(2);
    Test::Swap::impl_test_batched_swap_analytical<DeviceType, ScalarType, LayoutType, Kokkos::LayoutRight, ArgMode>(1);
    Test::Swap::impl_test_batched_swap_analytical<DeviceType, ScalarType, LayoutType, Kokkos::LayoutRight, ArgMode>(2);

    for (int ib = 0; ib < 5; ib++) {
      for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
          Test::Swap::impl_test_batched_swap<DeviceType, ScalarType, LayoutType, Kokkos::LayoutLeft, ArgMode>(ib, i, j);
          Test::Swap::impl_test_batched_swap<DeviceType, ScalarType, LayoutType, Kokkos::LayoutRight, ArgMode>(ib, i,
                                                                                                               j);
        }
      }
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Swap::impl_test_batched_swap_analytical<DeviceType, ScalarType, LayoutType, Kokkos::LayoutLeft, ArgMode>(1);
    Test::Swap::impl_test_batched_swap_analytical<DeviceType, ScalarType, LayoutType, Kokkos::LayoutLeft, ArgMode>(2);
    Test::Swap::impl_test_batched_swap_analytical<DeviceType, ScalarType, LayoutType, Kokkos::LayoutRight, ArgMode>(1);
    Test::Swap::impl_test_batched_swap_analytical<DeviceType, ScalarType, LayoutType, Kokkos::LayoutRight, ArgMode>(2);

    for (int ib = 0; ib < 5; ib++) {
      for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
          Test::Swap::impl_test_batched_swap<DeviceType, ScalarType, LayoutType, Kokkos::LayoutLeft, ArgMode>(ib, i, j);
          Test::Swap::impl_test_batched_swap<DeviceType, ScalarType, LayoutType, Kokkos::LayoutRight, ArgMode>(ib, i,
                                                                                                               j);
        }
      }
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_swap_float) {
  test_batched_swap<TestDevice, float, KokkosBatched::Mode::Serial>();
}
// Team
TEST_F(TestCategory, test_batched_team_swap_float) {
  test_batched_swap<TestDevice, float, KokkosBatched::Mode::Team>();
}
// TeamVector
TEST_F(TestCategory, test_batched_teamvector_swap_float) {
  test_batched_swap<TestDevice, float, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_swap_double) {
  test_batched_swap<TestDevice, double, KokkosBatched::Mode::Serial>();
}
// Team
TEST_F(TestCategory, test_batched_team_swap_double) {
  test_batched_swap<TestDevice, double, KokkosBatched::Mode::Team>();
}
// TeamVector
TEST_F(TestCategory, test_batched_teamvector_swap_double) {
  test_batched_swap<TestDevice, double, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_swap_fcomplex) {
  test_batched_swap<TestDevice, Kokkos::complex<float>, KokkosBatched::Mode::Serial>();
}
// Team
TEST_F(TestCategory, test_batched_team_swap_fcomplex) {
  test_batched_swap<TestDevice, Kokkos::complex<float>, KokkosBatched::Mode::Team>();
}
// TeamVector
TEST_F(TestCategory, test_batched_teamvector_swap_fcomplex) {
  test_batched_swap<TestDevice, Kokkos::complex<float>, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_swap_dcomplex) {
  test_batched_swap<TestDevice, Kokkos::complex<double>, KokkosBatched::Mode::Serial>();
}
// Team
TEST_F(TestCategory, test_batched_team_swap_dcomplex) {
  test_batched_swap<TestDevice, Kokkos::complex<double>, KokkosBatched::Mode::Team>();
}
// TeamVector
TEST_F(TestCategory, test_batched_teamvector_swap_dcomplex) {
  test_batched_swap<TestDevice, Kokkos::complex<double>, KokkosBatched::Mode::TeamVector>();
}
#endif
