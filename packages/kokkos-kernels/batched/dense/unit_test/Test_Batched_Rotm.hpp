// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Rotm.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Rotm {

template <typename DeviceType, typename XViewType, typename YViewType, typename ParamViewType>
struct Functor_BatchedSerialRotm {
  using execution_space = typename DeviceType::execution_space;
  XViewType m_x;
  YViewType m_y;
  ParamViewType m_param;

  Functor_BatchedSerialRotm(const XViewType &x, const YViewType &y, const ParamViewType &param)
      : m_x(x), m_y(y), m_param(param) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    auto sub_x     = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto sub_y     = Kokkos::subview(m_y, k, Kokkos::ALL());
    auto sub_param = Kokkos::subview(m_param, k, Kokkos::ALL());

    info += KokkosBatched::SerialRotm::invoke(sub_x, sub_y, sub_param);
  }

  inline int run() {
    using value_type = typename XViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialRotm");
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

template <typename DeviceType, typename XViewType, typename YViewType, typename ParamViewType, typename ArgMode>
struct Functor_BatchedTeamRotm {
  using execution_space = typename DeviceType::execution_space;
  XViewType m_x;
  YViewType m_y;
  ParamViewType m_param;

  Functor_BatchedTeamRotm(const XViewType &x, const YViewType &y, const ParamViewType &param)
      : m_x(x), m_y(y), m_param(param) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int k    = member.league_rank();
    auto sub_x     = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto sub_y     = Kokkos::subview(m_y, k, Kokkos::ALL());
    auto sub_param = Kokkos::subview(m_param, k, Kokkos::ALL());
    if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Team>) {
      KokkosBatched::TeamRotm<MemberType>::invoke(member, sub_x, sub_y, sub_param);
    } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::TeamVector>) {
      KokkosBatched::TeamVectorRotm<MemberType>::invoke(member, sub_x, sub_y, sub_param);
    }
  }

  inline void run() {
    using value_type                  = typename XViewType::non_const_value_type;
    std::string name_region           = std::is_same_v<ArgMode, KokkosBatched::Mode::Team>
                                            ? "KokkosBatched::Test::TeamRotm"
                                            : "KokkosBatched::Test::TeamVectorRotm";
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    const int league_size = m_x.extent_int(0);
    Kokkos::TeamPolicy<execution_space> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

/// \brief Implementation details of batched rotm analytical test
/// Applies the modified Givens rotation to vectors x and y:
///  x(i) := h11*x(i) + h12*y(i)
///  y(i) := h21*x(i) + h22*y(i)
///
/// The matrix H is given by
/// flag == -1.0  flag ==  0.0  flag ==  1.0  flag == -2.0
/// [[h11, h12],  [[1, h12],    [[h11, 1],    [[1, 0],
///  [h21, h22]]   [h21, 1]]     [[-1, h12]]   [[0, 1]]
///
///  x: [1, 2, 3]
///  y: [4, 5, 6]
///  param: [h11, h21, h12, h22] = [2, 0.5, -1, 3]
///
/// flag == -1.0        flag ==  0.0   flag ==  1.0   flag == -2.0
/// x: [-2, -1, 0]      [-3, -3, -3]   [6, 9, 12]     [1, 2, 3]
/// y: [12.5, 16, 19.5] [4.5, 6, 7.5]  [11, 13, 15]   [4, 5, 6]
///
/// \tparam DeviceType Kokkos device type
/// \tparam ScalarType Kokkos scalar type
/// \tparam LayoutType Kokkos layout type for the views
/// \tparam Flag: one of -1, 0, 1, or -2, specifying the values of h11, h21, h12, and h22 as described above
/// \tparam ArgMode: one of Mode::Serial, Mode::Team, Mode::TeamVector
///
/// \param[in] Nb Batch size of vectors
template <typename DeviceType, typename ScalarType, typename LayoutType, int Flag, typename ArgMode>
void impl_test_batched_rotm_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;

  const std::size_t N = 3;
  View2DType x("x", Nb, N), y("y", Nb, N);
  View2DType param("param", Nb, 5);
  View2DType x_ref("x_ref", Nb, N), y_ref("y_ref", Nb, N);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, N, Nb * incx};
  StridedView2DType x_s("x_s", layout), y_s("y_s", layout);

  auto h_x     = Kokkos::create_mirror_view(x);
  auto h_y     = Kokkos::create_mirror_view(y);
  auto h_x_ref = Kokkos::create_mirror_view(x_ref);
  auto h_y_ref = Kokkos::create_mirror_view(y_ref);
  auto h_param = Kokkos::create_mirror_view(param);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    h_param(ib, 0) = static_cast<ScalarType>(Flag);
    h_param(ib, 1) = static_cast<ScalarType>(2);
    h_param(ib, 2) = static_cast<ScalarType>(0.5);
    h_param(ib, 3) = static_cast<ScalarType>(-1);
    h_param(ib, 4) = static_cast<ScalarType>(3);

    h_x(ib, 0) = ScalarType(1);
    h_x(ib, 1) = ScalarType(2);
    h_x(ib, 2) = ScalarType(3);

    h_y(ib, 0) = ScalarType(4);
    h_y(ib, 1) = ScalarType(5);
    h_y(ib, 2) = ScalarType(6);
    if constexpr (Flag == -1) {
      h_x_ref(ib, 0) = ScalarType(-2);
      h_x_ref(ib, 1) = ScalarType(-1);
      h_x_ref(ib, 2) = ScalarType(0);

      h_y_ref(ib, 0) = ScalarType(12.5);
      h_y_ref(ib, 1) = ScalarType(16);
      h_y_ref(ib, 2) = ScalarType(19.5);
    } else if constexpr (Flag == 0) {
      h_x_ref(ib, 0) = ScalarType(-3);
      h_x_ref(ib, 1) = ScalarType(-3);
      h_x_ref(ib, 2) = ScalarType(-3);

      h_y_ref(ib, 0) = ScalarType(4.5);
      h_y_ref(ib, 1) = ScalarType(6);
      h_y_ref(ib, 2) = ScalarType(7.5);
    } else if constexpr (Flag == 1) {
      h_x_ref(ib, 0) = ScalarType(6);
      h_x_ref(ib, 1) = ScalarType(9);
      h_x_ref(ib, 2) = ScalarType(12);

      h_y_ref(ib, 0) = ScalarType(11);
      h_y_ref(ib, 1) = ScalarType(13);
      h_y_ref(ib, 2) = ScalarType(15);
    } else if constexpr (Flag == -2) {
      h_x_ref(ib, 0) = ScalarType(1);
      h_x_ref(ib, 1) = ScalarType(2);
      h_x_ref(ib, 2) = ScalarType(3);

      h_y_ref(ib, 0) = ScalarType(4);
      h_y_ref(ib, 1) = ScalarType(5);
      h_y_ref(ib, 2) = ScalarType(6);
    }
  }

  Kokkos::deep_copy(x, h_x);
  Kokkos::deep_copy(y, h_y);
  Kokkos::deep_copy(param, h_param);

  // Deep copy to strided views
  Kokkos::deep_copy(x_s, x);
  Kokkos::deep_copy(y_s, y);

  if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
    auto info = Functor_BatchedSerialRotm<DeviceType, View2DType, View2DType, View2DType>(x, y, param).run();
    EXPECT_EQ(info, 0);
  } else {
    Functor_BatchedTeamRotm<DeviceType, View2DType, View2DType, View2DType, ArgMode>(x, y, param).run();
  }

  // With strided views
  if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
    auto info =
        Functor_BatchedSerialRotm<DeviceType, StridedView2DType, StridedView2DType, View2DType>(x_s, y_s, param).run();
    EXPECT_EQ(info, 0);
  } else {
    Functor_BatchedTeamRotm<DeviceType, StridedView2DType, StridedView2DType, View2DType, ArgMode>(x_s, y_s, param)
        .run();
  }

  RealType eps = 1.0e1 * ats::epsilon();
  Kokkos::deep_copy(h_x, x);
  Kokkos::deep_copy(h_y, y);

  // Check if x(i) := h11*x(i) + h12*y(i)
  //          y(i) := h21*x(i) + h22*y(i)
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < N; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
      EXPECT_NEAR_KK(h_y(ib, i), h_y_ref(ib, i), eps);
    }
  }

  // Testing for strided views x_s and y_s, reusing x and y
  Kokkos::deep_copy(x, x_s);
  Kokkos::deep_copy(y, y_s);
  Kokkos::deep_copy(h_x, x);
  Kokkos::deep_copy(h_y, y);
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < N; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
      EXPECT_NEAR_KK(h_y(ib, i), h_y_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched rotm test
///
/// \tparam DeviceType Kokkos device type
/// \tparam ScalarType Kokkos scalar type
/// \tparam LayoutType Kokkos layout type for the views
/// \tparam Flag: one of -1, 0, 1, or -2, specifying the values of h11, h21, h12, and h22 as described above
/// \tparam ArgMode: one of Mode::Serial, Mode::Team, Mode::TeamVector
///
/// \param[in] Nb Batch size of vectors
/// \param[in] N Length of vectors x and y
template <typename DeviceType, typename ScalarType, typename LayoutType, int Flag, typename ArgMode>
void impl_test_batched_rotm(const std::size_t Nb, const std::size_t N) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;

  View2DType x("x", Nb, N), y("y", Nb, N), param("param", Nb, 5);
  View2DType x_ref("x_ref", Nb, N), y_ref("y_ref", Nb, N);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, N, Nb * incx};
  StridedView2DType x_s("x_s", layout), y_s("y_s", layout);

  // Create random x and y
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);
  Kokkos::fill_random(y, rand_pool, randStart, randEnd);
  Kokkos::fill_random(param, rand_pool, randStart, randEnd);

  // Set the flag value in param
  auto sub_param = Kokkos::subview(param, Kokkos::ALL(), 0);
  Kokkos::deep_copy(sub_param, static_cast<ScalarType>(Flag));

  // Save copies for reference
  Kokkos::deep_copy(x_ref, x);
  Kokkos::deep_copy(y_ref, y);

  // Deep copy to strided views
  Kokkos::deep_copy(x_s, x);
  Kokkos::deep_copy(y_s, y);

  // Run rotm on (x, y)
  if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
    auto info = Functor_BatchedSerialRotm<DeviceType, View2DType, View2DType, View2DType>(x, y, param).run();
    EXPECT_EQ(info, 0);
  } else {
    Functor_BatchedTeamRotm<DeviceType, View2DType, View2DType, View2DType, ArgMode>(x, y, param).run();
  }

  // With strided views
  if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
    auto info =
        Functor_BatchedSerialRotm<DeviceType, StridedView2DType, StridedView2DType, View2DType>(x_s, y_s, param).run();
    EXPECT_EQ(info, 0);
  } else {
    Functor_BatchedTeamRotm<DeviceType, StridedView2DType, StridedView2DType, View2DType, ArgMode>(x_s, y_s, param)
        .run();
  }

  // Make a reference at host
  auto h_x_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x_ref);
  auto h_y_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y_ref);
  auto h_param = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, param);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < N; i++) {
      // No modification is needed for flag == -2, as the reference is the same as the input
      if constexpr (Flag == -1) {
        auto temp      = h_param(ib, 1) * h_x_ref(ib, i) + h_param(ib, 3) * h_y_ref(ib, i);
        h_y_ref(ib, i) = h_param(ib, 2) * h_x_ref(ib, i) + h_param(ib, 4) * h_y_ref(ib, i);
        h_x_ref(ib, i) = temp;
      } else if constexpr (Flag == 0) {
        auto temp      = h_x_ref(ib, i) + h_param(ib, 3) * h_y_ref(ib, i);
        h_y_ref(ib, i) = h_param(ib, 2) * h_x_ref(ib, i) + h_y_ref(ib, i);
        h_x_ref(ib, i) = temp;
      } else if constexpr (Flag == 1) {
        auto temp      = h_param(ib, 1) * h_x_ref(ib, i) + h_y_ref(ib, i);
        h_y_ref(ib, i) = -1.0 * h_x_ref(ib, i) + h_param(ib, 4) * h_y_ref(ib, i);
        h_x_ref(ib, i) = temp;
      }
    }
  }

  RealType eps = 1.0e1 * ats::epsilon();
  auto h_x     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x);
  auto h_y     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, y);

  // Check if x(i) := h11*x(i) + h12*y(i)
  //          y(i) := h21*x(i) + h22*y(i)
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < N; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
      EXPECT_NEAR_KK(h_y(ib, i), h_y_ref(ib, i), eps);
    }
  }

  // Testing for strided views x_s and y_s, reusing x and y
  Kokkos::deep_copy(x, x_s);
  Kokkos::deep_copy(y, y_s);
  Kokkos::deep_copy(h_x, x);
  Kokkos::deep_copy(h_y, y);
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < N; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
      EXPECT_NEAR_KK(h_y(ib, i), h_y_ref(ib, i), eps);
    }
  }
}

}  // namespace Rotm
}  // namespace Test

template <typename DeviceType, typename ScalarType, int Flag, typename ArgMode>
int test_batched_rotm() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Rotm::impl_test_batched_rotm_analytical<DeviceType, ScalarType, LayoutType, Flag, ArgMode>(1);
    Test::Rotm::impl_test_batched_rotm_analytical<DeviceType, ScalarType, LayoutType, Flag, ArgMode>(2);
    for (int i = 0; i < 10; i++) {
      Test::Rotm::impl_test_batched_rotm<DeviceType, ScalarType, LayoutType, Flag, ArgMode>(1, i);
      Test::Rotm::impl_test_batched_rotm<DeviceType, ScalarType, LayoutType, Flag, ArgMode>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Rotm::impl_test_batched_rotm_analytical<DeviceType, ScalarType, LayoutType, Flag, ArgMode>(1);
    Test::Rotm::impl_test_batched_rotm_analytical<DeviceType, ScalarType, LayoutType, Flag, ArgMode>(2);
    for (int i = 0; i < 10; i++) {
      Test::Rotm::impl_test_batched_rotm<DeviceType, ScalarType, LayoutType, Flag, ArgMode>(1, i);
      Test::Rotm::impl_test_batched_rotm<DeviceType, ScalarType, LayoutType, Flag, ArgMode>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_rotm_float) {
  test_batched_rotm<TestDevice, float, -1, KokkosBatched::Mode::Serial>();
  test_batched_rotm<TestDevice, float, 0, KokkosBatched::Mode::Serial>();
  test_batched_rotm<TestDevice, float, 1, KokkosBatched::Mode::Serial>();
  test_batched_rotm<TestDevice, float, -2, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_rotm_float) {
  test_batched_rotm<TestDevice, float, -1, KokkosBatched::Mode::Team>();
  test_batched_rotm<TestDevice, float, 0, KokkosBatched::Mode::Team>();
  test_batched_rotm<TestDevice, float, 1, KokkosBatched::Mode::Team>();
  test_batched_rotm<TestDevice, float, -2, KokkosBatched::Mode::Team>();
}

// TeamVector
TEST_F(TestCategory, test_batched_teamvector_rotm_float) {
  test_batched_rotm<TestDevice, float, -1, KokkosBatched::Mode::TeamVector>();
  test_batched_rotm<TestDevice, float, 0, KokkosBatched::Mode::TeamVector>();
  test_batched_rotm<TestDevice, float, 1, KokkosBatched::Mode::TeamVector>();
  test_batched_rotm<TestDevice, float, -2, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_rotm_double) {
  test_batched_rotm<TestDevice, double, -1, KokkosBatched::Mode::Serial>();
  test_batched_rotm<TestDevice, double, 0, KokkosBatched::Mode::Serial>();
  test_batched_rotm<TestDevice, double, 1, KokkosBatched::Mode::Serial>();
  test_batched_rotm<TestDevice, double, -2, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_rotm_double) {
  test_batched_rotm<TestDevice, double, -1, KokkosBatched::Mode::Team>();
  test_batched_rotm<TestDevice, double, 0, KokkosBatched::Mode::Team>();
  test_batched_rotm<TestDevice, double, 1, KokkosBatched::Mode::Team>();
  test_batched_rotm<TestDevice, double, -2, KokkosBatched::Mode::Team>();
}

// TeamVector
TEST_F(TestCategory, test_batched_teamvector_rotm_double) {
  test_batched_rotm<TestDevice, double, -1, KokkosBatched::Mode::TeamVector>();
  test_batched_rotm<TestDevice, double, 0, KokkosBatched::Mode::TeamVector>();
  test_batched_rotm<TestDevice, double, 1, KokkosBatched::Mode::TeamVector>();
  test_batched_rotm<TestDevice, double, -2, KokkosBatched::Mode::TeamVector>();
}
#endif
