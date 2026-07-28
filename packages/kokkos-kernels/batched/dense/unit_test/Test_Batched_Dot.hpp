// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Dot.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Dot {

template <typename DeviceType, typename ArgTrans, int Axis, typename XViewType, typename YViewType,
          typename NormViewType>
struct Functor_BatchedSerialDot {
  using execution_space = typename DeviceType::execution_space;
  XViewType m_x;
  YViewType m_y;
  NormViewType m_dot;

  Functor_BatchedSerialDot(const XViewType &x, const YViewType &y, const NormViewType &dot)
      : m_x(x), m_y(y), m_dot(dot) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    if constexpr (XViewType::rank() == 2) {
      auto sub_x   = Kokkos::subview(m_x, k, Kokkos::ALL());
      auto sub_y   = Kokkos::subview(m_y, k, Kokkos::ALL());
      auto sub_dot = Kokkos::subview(m_dot, k);

      info += KokkosBatched::SerialDot<ArgTrans, Axis>::invoke(sub_x, sub_y, sub_dot);
    } else {
      auto sub_x   = Kokkos::subview(m_x, k, Kokkos::ALL(), Kokkos::ALL());
      auto sub_y   = Kokkos::subview(m_y, k, Kokkos::ALL(), Kokkos::ALL());
      auto sub_dot = Kokkos::subview(m_dot, k, Kokkos::ALL());

      info += KokkosBatched::SerialDot<ArgTrans, Axis>::invoke(sub_x, sub_y, sub_dot);
    }
  }

  inline int run() {
    using value_type = typename XViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialDot");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_x.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename ArgTrans, int Axis, typename XViewType, typename YViewType,
          typename NormViewType, typename ArgMode>
struct Functor_BatchedTeamDot {
  using execution_space = typename DeviceType::execution_space;
  XViewType m_x;
  YViewType m_y;
  NormViewType m_dot;

  Functor_BatchedTeamDot(const XViewType &x, const YViewType &y, const NormViewType &dot)
      : m_x(x), m_y(y), m_dot(dot) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int k = member.league_rank();
    if constexpr (XViewType::rank() == 2) {
      auto sub_x   = Kokkos::subview(m_x, k, Kokkos::ALL());
      auto sub_y   = Kokkos::subview(m_y, k, Kokkos::ALL());
      auto sub_dot = Kokkos::subview(m_dot, k);

      if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Team>) {
        KokkosBatched::TeamDot<MemberType, ArgTrans, Axis>::invoke(member, sub_x, sub_y, sub_dot);
      } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::TeamVector>) {
        KokkosBatched::TeamVectorDot<MemberType, ArgTrans, Axis>::invoke(member, sub_x, sub_y, sub_dot);
      }
    } else {
      auto sub_x   = Kokkos::subview(m_x, k, Kokkos::ALL(), Kokkos::ALL());
      auto sub_y   = Kokkos::subview(m_y, k, Kokkos::ALL(), Kokkos::ALL());
      auto sub_dot = Kokkos::subview(m_dot, k, Kokkos::ALL());

      if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Team>) {
        KokkosBatched::TeamDot<MemberType, ArgTrans, Axis>::invoke(member, sub_x, sub_y, sub_dot);
      } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::TeamVector>) {
        KokkosBatched::TeamVectorDot<MemberType, ArgTrans, Axis>::invoke(member, sub_x, sub_y, sub_dot);
      }
    }
  }

  inline void run() {
    using value_type        = typename XViewType::non_const_value_type;
    std::string name_region = std::is_same_v<ArgMode, KokkosBatched::Mode::Team> ? "KokkosBatched::Test::TeamDot"
                                                                                 : "KokkosBatched::Test::TeamVectorDot";
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    const int league_size = m_x.extent_int(0);
    Kokkos::TeamPolicy<execution_space> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

/// \brief Implementation details of batched dot analytical test
///        to confirm the dot product is computed correctly
///        x0: [1, 2, 3]
///        y0: [4, 5, 6]
///        dot0: 32
///
///        x1: [[1, 2, 3],
///             [4, 5, 6]]
///        y1: [[4, 5, 6],
///             [7, 8, 9]]
///        dot1-ax0: [32, 50, 72]
///        dot1-ax1: [32, 77]
///
/// \param Nb [in] Batch size of vectors
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ArgTrans, typename ArgMode>
void impl_test_batched_dot_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View1DType        = Kokkos::View<ScalarType *, LayoutType, DeviceType>;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using StridedView3DType = Kokkos::View<ScalarType ***, Kokkos::LayoutStride, DeviceType>;

  const std::size_t m = 2, n = 3;
  View2DType x0("x0", Nb, n), y0("y0", Nb, n);
  View3DType x1("x1", Nb, m, n), y1("y1", Nb, m, n);
  View1DType dot0("dot0", Nb), ref_dot0("ref_dot0", Nb);
  View2DType dot1_ax0("dot1_ax0", Nb, n), dot1_ax1("dot1_ax1", Nb, m), ref_dot1_ax0("ref_dot1_ax0", Nb, n),
      ref_dot1_ax1("ref_dot1_ax1", Nb, m);

  // Testing incx/incy argument with strided views
  const std::size_t incx = 2, incy = 2;
  Kokkos::LayoutStride layout_x0{Nb, incx, n, Nb * incx}, layout_y0{Nb, incy, n, Nb * incy};
  Kokkos::LayoutStride layout_x1{Nb, incx, m, Nb * incx, n, Nb * incx * m},
      layout_y1{Nb, incy, m, Nb * incy, n, Nb * incy * m};
  StridedView2DType x0_s("x_s", layout_x0), y0_s("y_s", layout_y0);
  StridedView3DType x1_s("x_s", layout_x1), y1_s("y_s", layout_y1);
  View1DType dot0_s("dot_s", Nb);
  View2DType dot1_ax0_s("dot1_ax0_s", Nb, n), dot1_ax1_s("dot1_ax1_s", Nb, m);

  auto h_x0           = Kokkos::create_mirror_view(x0);
  auto h_y0           = Kokkos::create_mirror_view(y0);
  auto h_ref_dot0     = Kokkos::create_mirror_view(ref_dot0);
  auto h_x1           = Kokkos::create_mirror_view(x1);
  auto h_y1           = Kokkos::create_mirror_view(y1);
  auto h_ref_dot1_ax0 = Kokkos::create_mirror_view(ref_dot1_ax0);
  auto h_ref_dot1_ax1 = Kokkos::create_mirror_view(ref_dot1_ax1);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    h_x0(ib, 0)    = static_cast<ScalarType>(1);
    h_x0(ib, 1)    = static_cast<ScalarType>(2);
    h_x0(ib, 2)    = static_cast<ScalarType>(3);
    h_y0(ib, 0)    = static_cast<ScalarType>(4);
    h_y0(ib, 1)    = static_cast<ScalarType>(5);
    h_y0(ib, 2)    = static_cast<ScalarType>(6);
    h_ref_dot0(ib) = ScalarType(32);  // dot0 = 1*4 + 2*5 + 3*6 = 32

    h_x1(ib, 0, 0)        = static_cast<ScalarType>(1);
    h_x1(ib, 0, 1)        = static_cast<ScalarType>(2);
    h_x1(ib, 0, 2)        = static_cast<ScalarType>(3);
    h_x1(ib, 1, 0)        = static_cast<ScalarType>(4);
    h_x1(ib, 1, 1)        = static_cast<ScalarType>(5);
    h_x1(ib, 1, 2)        = static_cast<ScalarType>(6);
    h_y1(ib, 0, 0)        = static_cast<ScalarType>(4);
    h_y1(ib, 0, 1)        = static_cast<ScalarType>(5);
    h_y1(ib, 0, 2)        = static_cast<ScalarType>(6);
    h_y1(ib, 1, 0)        = static_cast<ScalarType>(7);
    h_y1(ib, 1, 1)        = static_cast<ScalarType>(8);
    h_y1(ib, 1, 2)        = static_cast<ScalarType>(9);
    h_ref_dot1_ax0(ib, 0) = ScalarType(32);   // dot1_ax0(0) = 1*4 + 4*7 = 32
    h_ref_dot1_ax0(ib, 1) = ScalarType(50);   // dot1_ax0(1) = 2*5 + 5*8 = 50
    h_ref_dot1_ax0(ib, 2) = ScalarType(72);   // dot1_ax0(2) = 3*6 + 6*9 = 72
    h_ref_dot1_ax1(ib, 0) = ScalarType(32);   // dot1_ax1(0) = 1*4 + 2*5 + 3*6 = 32
    h_ref_dot1_ax1(ib, 1) = ScalarType(122);  // dot1_ax1(1) = 4*7 + 5*8 + 6*9 = 122
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

  if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
    // 1D case
    auto info0 =
        Functor_BatchedSerialDot<DeviceType, ArgTrans, 0, View2DType, View2DType, View1DType>(x0, y0, dot0).run();
    EXPECT_EQ(info0, 0);

    // With strided views
    info0 = Functor_BatchedSerialDot<DeviceType, ArgTrans, 0, StridedView2DType, StridedView2DType, View1DType>(
                x0_s, y0_s, dot0_s)
                .run();
    EXPECT_EQ(info0, 0);

    // 2D case: axis 0
    auto info1_ax0 =
        Functor_BatchedSerialDot<DeviceType, ArgTrans, 0, View3DType, View3DType, View2DType>(x1, y1, dot1_ax0).run();
    EXPECT_EQ(info1_ax0, 0);

    // 2D case: axis 1
    auto info1_ax1 =
        Functor_BatchedSerialDot<DeviceType, ArgTrans, 1, View3DType, View3DType, View2DType>(x1, y1, dot1_ax1).run();
    EXPECT_EQ(info1_ax1, 0);

    // With strided views
    info1_ax0 = Functor_BatchedSerialDot<DeviceType, ArgTrans, 0, StridedView3DType, StridedView3DType, View2DType>(
                    x1_s, y1_s, dot1_ax0_s)
                    .run();
    EXPECT_EQ(info1_ax0, 0);

    info1_ax1 = Functor_BatchedSerialDot<DeviceType, ArgTrans, 1, StridedView3DType, StridedView3DType, View2DType>(
                    x1_s, y1_s, dot1_ax1_s)
                    .run();
    EXPECT_EQ(info1_ax1, 0);
  } else {
    // 1D case
    Functor_BatchedTeamDot<DeviceType, ArgTrans, 0, View2DType, View2DType, View1DType, ArgMode>(x0, y0, dot0).run();

    // With strided views
    Functor_BatchedTeamDot<DeviceType, ArgTrans, 0, StridedView2DType, StridedView2DType, View1DType, ArgMode>(
        x0_s, y0_s, dot0_s)
        .run();

    // 2D case: axis 0
    Functor_BatchedTeamDot<DeviceType, ArgTrans, 0, View3DType, View3DType, View2DType, ArgMode>(x1, y1, dot1_ax0)
        .run();

    // 2D case: axis 1
    Functor_BatchedTeamDot<DeviceType, ArgTrans, 1, View3DType, View3DType, View2DType, ArgMode>(x1, y1, dot1_ax1)
        .run();

    // With strided views
    Functor_BatchedTeamDot<DeviceType, ArgTrans, 0, StridedView3DType, StridedView3DType, View2DType, ArgMode>(
        x1_s, y1_s, dot1_ax0_s)
        .run();

    Functor_BatchedTeamDot<DeviceType, ArgTrans, 1, StridedView3DType, StridedView3DType, View2DType, ArgMode>(
        x1_s, y1_s, dot1_ax1_s)
        .run();
  }

  RealType eps      = 1.0e1 * ats::epsilon();
  auto h_dot0       = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dot0);
  auto h_dot0_s     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dot0_s);
  auto h_dot1_ax0   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dot1_ax0);
  auto h_dot1_ax1   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dot1_ax1);
  auto h_dot1_ax0_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dot1_ax0_s);
  auto h_dot1_ax1_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dot1_ax1_s);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    EXPECT_NEAR_KK(h_dot0(ib), h_ref_dot0(ib), eps);
    EXPECT_NEAR_KK(h_dot0_s(ib), h_ref_dot0(ib), eps);
    for (std::size_t j = 0; j < n; j++) {
      EXPECT_NEAR_KK(h_dot1_ax0(ib, j), h_ref_dot1_ax0(ib, j), eps);
      EXPECT_NEAR_KK(h_dot1_ax0_s(ib, j), h_ref_dot1_ax0(ib, j), eps);
    }
    for (std::size_t j = 0; j < m; j++) {
      EXPECT_NEAR_KK(h_dot1_ax1(ib, j), h_ref_dot1_ax1(ib, j), eps);
      EXPECT_NEAR_KK(h_dot1_ax1_s(ib, j), h_ref_dot1_ax1(ib, j), eps);
    }
  }
}

/// \brief Implementation details of batched dot test
///        Confirm dot = X^T * Y or X^H * Y
///
/// \param[in] Nb Batch size of vectors
/// \param[in] m  Number of rows
/// \param[in] n  Number of columns
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ArgTrans, typename ArgMode>
void impl_test_batched_dot(const std::size_t Nb, const std::size_t m, const std::size_t n) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View1DType = Kokkos::View<ScalarType *, LayoutType, DeviceType>;
  using View2DType = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  View2DType x0("x0", Nb, n), y0("y0", Nb, n);
  View3DType x1("x1", Nb, m, n), y1("y1", Nb, m, n);
  View1DType dot0("dot0", Nb), ref_dot0("ref_dot0", Nb);
  View2DType dot1_ax0("dot1_ax0", Nb, n), dot1_ax1("dot1_ax1", Nb, m), ref_dot1_ax0("ref_dot1_ax0", Nb, n),
      ref_dot1_ax1("ref_dot1_ax1", Nb, m);

  // Create random vectors x and y
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(x0, rand_pool, randStart, randEnd);
  Kokkos::fill_random(y0, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x1, rand_pool, randStart, randEnd);
  Kokkos::fill_random(y1, rand_pool, randStart, randEnd);

  using Op = std::conditional_t<std::is_same_v<ArgTrans, KokkosBatched::Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                KokkosBlas::Impl::OpID>;
  Op op;

  // Make a reference
  auto h_x0           = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x0);
  auto h_y0           = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y0);
  auto h_ref_dot0     = Kokkos::create_mirror_view(ref_dot0);
  auto h_x1           = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x1);
  auto h_y1           = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y1);
  auto h_ref_dot1_ax0 = Kokkos::create_mirror_view(ref_dot1_ax0);
  auto h_ref_dot1_ax1 = Kokkos::create_mirror_view(ref_dot1_ax1);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    ScalarType tmp_dot0 = 0;
    for (std::size_t j = 0; j < n; j++) {
      tmp_dot0 += op(h_x0(ib, j)) * h_y0(ib, j);
    }
    h_ref_dot0(ib) = tmp_dot0;
    for (std::size_t j = 0; j < n; j++) {
      ScalarType tmp_dot1 = 0;
      for (std::size_t i = 0; i < m; i++) {
        tmp_dot1 += op(h_x1(ib, i, j)) * h_y1(ib, i, j);
      }
      h_ref_dot1_ax0(ib, j) = tmp_dot1;
    }
    for (std::size_t j = 0; j < m; j++) {
      ScalarType tmp_dot1 = 0;
      for (std::size_t i = 0; i < n; i++) {
        tmp_dot1 += op(h_x1(ib, j, i)) * h_y1(ib, j, i);
      }
      h_ref_dot1_ax1(ib, j) = tmp_dot1;
    }
  }

  // Dot operation
  if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
    auto info0 =
        Functor_BatchedSerialDot<DeviceType, ArgTrans, 0, View2DType, View2DType, View1DType>(x0, y0, dot0).run();
    EXPECT_EQ(info0, 0);
    auto info1_ax0 =
        Functor_BatchedSerialDot<DeviceType, ArgTrans, 0, View3DType, View3DType, View2DType>(x1, y1, dot1_ax0).run();
    EXPECT_EQ(info1_ax0, 0);
    auto info1_ax1 =
        Functor_BatchedSerialDot<DeviceType, ArgTrans, 1, View3DType, View3DType, View2DType>(x1, y1, dot1_ax1).run();
    EXPECT_EQ(info1_ax1, 0);
  } else {
    Functor_BatchedTeamDot<DeviceType, ArgTrans, 0, View2DType, View2DType, View1DType, ArgMode>(x0, y0, dot0).run();
    Functor_BatchedTeamDot<DeviceType, ArgTrans, 0, View3DType, View3DType, View2DType, ArgMode>(x1, y1, dot1_ax0)
        .run();
    Functor_BatchedTeamDot<DeviceType, ArgTrans, 1, View3DType, View3DType, View2DType, ArgMode>(x1, y1, dot1_ax1)
        .run();
  }

  RealType eps    = 1.0e1 * ats::epsilon();
  auto h_dot0     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dot0);
  auto h_dot1_ax0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dot1_ax0);
  auto h_dot1_ax1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dot1_ax1);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    EXPECT_NEAR_KK(h_dot0(ib), h_ref_dot0(ib), eps);
    for (std::size_t j = 0; j < n; j++) {
      EXPECT_NEAR_KK(h_dot1_ax0(ib, j), h_ref_dot1_ax0(ib, j), eps);
    }
    for (std::size_t j = 0; j < m; j++) {
      EXPECT_NEAR_KK(h_dot1_ax1(ib, j), h_ref_dot1_ax1(ib, j), eps);
    }
  }
}

}  // namespace Dot
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ArgTrans, typename ArgMode>
int test_batched_dot() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Dot::impl_test_batched_dot_analytical<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(1);
    Test::Dot::impl_test_batched_dot_analytical<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(2);
    for (std::size_t i = 0; i < 3; i++) {
      for (std::size_t j = 0; j < 3; j++) {
        Test::Dot::impl_test_batched_dot<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(1, i, j);
        Test::Dot::impl_test_batched_dot<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(2, i, j);
      }
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Dot::impl_test_batched_dot_analytical<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(1);
    Test::Dot::impl_test_batched_dot_analytical<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(2);
    for (std::size_t i = 0; i < 3; i++) {
      for (std::size_t j = 0; j < 3; j++) {
        Test::Dot::impl_test_batched_dot<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(1, i, j);
        Test::Dot::impl_test_batched_dot<DeviceType, ScalarType, LayoutType, ArgTrans, ArgMode>(2, i, j);
      }
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_dot_t_float) {
  test_batched_dot<TestDevice, float, Trans::Transpose, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_dot_c_float) {
  test_batched_dot<TestDevice, float, Trans::ConjTranspose, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_dot_t_float) {
  test_batched_dot<TestDevice, float, Trans::Transpose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_dot_c_float) {
  test_batched_dot<TestDevice, float, Trans::ConjTranspose, KokkosBatched::Mode::Team>();
}

// TeamVector
TEST_F(TestCategory, test_batched_team_vector_dot_t_float) {
  test_batched_dot<TestDevice, float, Trans::Transpose, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_team_vector_dot_c_float) {
  test_batched_dot<TestDevice, float, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_dot_t_double) {
  test_batched_dot<TestDevice, double, Trans::Transpose, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_dot_c_double) {
  test_batched_dot<TestDevice, double, Trans::ConjTranspose, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_dot_t_double) {
  test_batched_dot<TestDevice, double, Trans::Transpose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_dot_c_double) {
  test_batched_dot<TestDevice, double, Trans::ConjTranspose, KokkosBatched::Mode::Team>();
}

// TeamVector
TEST_F(TestCategory, test_batched_team_vector_dot_t_double) {
  test_batched_dot<TestDevice, double, Trans::Transpose, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_team_vector_dot_c_double) {
  test_batched_dot<TestDevice, double, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_dot_t_fcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<float>, Trans::Transpose, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_dot_c_fcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<float>, Trans::ConjTranspose, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_dot_t_fcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<float>, Trans::Transpose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_dot_c_fcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<float>, Trans::ConjTranspose, KokkosBatched::Mode::Team>();
}

// TeamVector
TEST_F(TestCategory, test_batched_team_vector_dot_t_fcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<float>, Trans::Transpose, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_team_vector_dot_c_fcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<float>, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_dot_t_dcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<double>, Trans::Transpose, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_dot_c_dcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<double>, Trans::ConjTranspose, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_dot_t_dcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<double>, Trans::Transpose, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_dot_c_dcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<double>, Trans::ConjTranspose, KokkosBatched::Mode::Team>();
}

// TeamVector
TEST_F(TestCategory, test_batched_team_vector_dot_t_dcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<double>, Trans::Transpose, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_team_vector_dot_c_dcomplex) {
  test_batched_dot<TestDevice, Kokkos::complex<double>, Trans::ConjTranspose, KokkosBatched::Mode::TeamVector>();
}
#endif
