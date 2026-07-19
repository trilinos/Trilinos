// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Nrm.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Nrm {

template <typename DeviceType, typename XViewType, typename NormViewType, typename ArgNorm, typename ArgMode>
struct Functor_BatchedNrm {
  using execution_space = typename DeviceType::execution_space;
  using member_type     = typename Kokkos::TeamPolicy<execution_space>::member_type;
  XViewType m_x;
  NormViewType m_norm;

  Functor_BatchedNrm(const XViewType &x, const NormViewType &norm) : m_x(x), m_norm(norm) {}

  KOKKOS_INLINE_FUNCTION void operator()(const member_type &member, int &info) const {
    const int k   = member.league_rank();
    auto sub_x    = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto sub_norm = Kokkos::subview(m_norm, k);
    if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Serial>) {
      info += KokkosBatched::SerialNrm<ArgNorm>::invoke(sub_x, sub_norm);
    } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::Team>) {
      info += KokkosBatched::TeamNrm<member_type, ArgNorm>::invoke(member, sub_x, sub_norm);
    } else if constexpr (std::is_same_v<ArgMode, KokkosBatched::Mode::TeamVector>) {
      info += KokkosBatched::TeamVectorNrm<member_type, ArgNorm>::invoke(member, sub_x, sub_norm);
    }
  }

  inline int run() {
    using value_type        = typename XViewType::non_const_value_type;
    std::string name_region = std::is_same_v<ArgMode, KokkosBatched::Mode::Serial> ? "KokkosBatched::Test::SerialNrm"
                              : std::is_same_v<ArgMode, KokkosBatched::Mode::Team>
                                  ? "KokkosBatched::Test::TeamNrm"
                                  : "KokkosBatched::Test::TeamVectorNrm";
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

/// \brief Implementation details of batched nrm analytical test
/// x = [1, 3, 5]
/// z = [1 + 2j, 3 + 4j, 5 + 6j]
/// L1 norm: x: 9, z: 21
/// L2 norm: x: sqrt(1^2 + 3^2 + 5^2) = sqrt(35), z: sqrt((1^2 + 2^2) + (3^2 + 4^2) + (5^2 + 6^2)) = sqrt(91)
/// Linf norm: x: 5, z: sqrt(5^2 + 6^2) = sqrt(61)
///
/// \tparam DeviceType Kokkos device type
/// \tparam ScalarType Kokkos scalar type
/// \tparam LayoutType Kokkos layout type for the views
/// \tparam ArgNorm: one of L1, L2, Linf, ScaledL2
/// \tparam ArgMode: one of Mode::Serial, Mode::Team, Mode::TeamVector
///
/// \param[in] Nb Batch size of vectors
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ArgNorm, typename ArgMode>
void impl_test_batched_nrm_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using NormViewType      = Kokkos::View<RealType *, LayoutType, DeviceType>;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;

  const std::size_t N = 3;
  View2DType x("x", Nb, N);
  NormViewType norm("norm", Nb), norm_s("norm_s", Nb), norm_ref("norm_ref", Nb);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, N, Nb * incx};
  StridedView2DType x_s("x_s", layout);

  auto h_x        = Kokkos::create_mirror_view(x);
  auto h_norm_ref = Kokkos::create_mirror_view(norm_ref);

  constexpr bool is_complex = KokkosKernels::ArithTraits<ScalarType>::is_complex;

  for (std::size_t ib = 0; ib < Nb; ib++) {
    if constexpr (is_complex) {
      h_x(ib, 0) = ScalarType(1, 2);
      h_x(ib, 1) = ScalarType(3, 4);
      h_x(ib, 2) = ScalarType(5, 6);

      // L1 norm: 21, L2 norm: sqrt(91), Linf norm: sqrt(61)
      h_norm_ref(ib) = std::is_same_v<ArgNorm, Norm::L1>         ? RealType(21)
                       : std::is_same_v<ArgNorm, Norm::L2>       ? Kokkos::sqrt(RealType(91))
                       : std::is_same_v<ArgNorm, Norm::LInf>     ? Kokkos::sqrt(RealType(61))
                       : std::is_same_v<ArgNorm, Norm::ScaledL2> ? Kokkos::sqrt(RealType(91))
                                                                 : RealType(-1);
    } else {
      h_x(ib, 0) = ScalarType(1);
      h_x(ib, 1) = ScalarType(3);
      h_x(ib, 2) = ScalarType(5);

      // L1 norm: 9, L2 norm: sqrt(35), Linf norm: 5
      h_norm_ref(ib) = std::is_same_v<ArgNorm, Norm::L1>         ? RealType(9)
                       : std::is_same_v<ArgNorm, Norm::L2>       ? ats::sqrt(RealType(35))
                       : std::is_same_v<ArgNorm, Norm::LInf>     ? RealType(5)
                       : std::is_same_v<ArgNorm, Norm::ScaledL2> ? ats::sqrt(RealType(35))
                                                                 : RealType(-1);
    }
  }

  Kokkos::deep_copy(x, h_x);

  // Deep copy to strided views
  Kokkos::deep_copy(x_s, x);

  auto info = Functor_BatchedNrm<DeviceType, View2DType, NormViewType, ArgNorm, ArgMode>(x, norm).run();
  EXPECT_EQ(info, 0);

  // With strided views
  info = Functor_BatchedNrm<DeviceType, StridedView2DType, NormViewType, ArgNorm, ArgMode>(x_s, norm_s).run();
  EXPECT_EQ(info, 0);

  RealType eps  = 1.0e1 * ats::epsilon();
  auto h_norm   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, norm);
  auto h_norm_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, norm_s);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    KK_EXPECT_NEAR(h_norm(ib), h_norm_ref(ib), eps);
    KK_EXPECT_NEAR(h_norm_s(ib), h_norm_ref(ib), eps);
  }
}

/// \brief Implementation details of batched nrm overflow test
/// 3.4e38 and 1e308 are the largest representable finite values for fp32 and fp64 respectively.
/// For fp32, where x^2 is expected to overflow
/// x = [1e20, 1e20, 1e20]
/// z = [1e20 + 1e20j, 1e20 + 1e20j, 1e20 + 1e20j]
/// L2 norm: x: sqrt(3) * 1e20, z: sqrt(6) * 1e20
///
/// For fp64, where x^2 is expected to overflow
/// x = [1e200, 1e200, 1e200]
/// z = [1e200 + 1e200j, 1e200 + 1e200j, 1e200 + 1e200j]
/// L2 norm: x: sqrt(3) * 1e200, z: sqrt(6) * 1e200
///
/// \tparam DeviceType Kokkos device type
/// \tparam ScalarType Kokkos scalar type
/// \tparam LayoutType Kokkos layout type for the views
/// \tparam ArgMode: one of Mode::Serial, Mode::Team, Mode::TeamVector
///
/// \param[in] Nb Batch size of vectors
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ArgMode>
void impl_test_batched_nrm_overflow(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using NormViewType      = Kokkos::View<RealType *, LayoutType, DeviceType>;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;

  const std::size_t N = 3;
  View2DType x("x", Nb, N);
  NormViewType norm("norm", Nb), norm_s("norm_s", Nb), norm_ref("norm_ref", Nb);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, N, Nb * incx};
  StridedView2DType x_s("x_s", layout);

  auto h_x        = Kokkos::create_mirror_view(x);
  auto h_norm_ref = Kokkos::create_mirror_view(norm_ref);

  constexpr bool is_complex = KokkosKernels::ArithTraits<ScalarType>::is_complex;
  const RealType large_val  = std::is_same_v<RealType, float> ? RealType(1e20) : RealType(1e200);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    if constexpr (is_complex) {
      h_x(ib, 0) = ScalarType(large_val, large_val);
      h_x(ib, 1) = ScalarType(large_val, large_val);
      h_x(ib, 2) = ScalarType(large_val, large_val);

      // L2 norm: sqrt(6) * large_val
      h_norm_ref(ib) = Kokkos::sqrt(RealType(6)) * large_val;
    } else {
      h_x(ib, 0) = ScalarType(large_val);
      h_x(ib, 1) = ScalarType(large_val);
      h_x(ib, 2) = ScalarType(large_val);

      // L2 norm: sqrt(3) * large_val
      h_norm_ref(ib) = Kokkos::sqrt(RealType(3)) * large_val;
    }
  }

  Kokkos::deep_copy(x, h_x);

  // Deep copy to strided views
  Kokkos::deep_copy(x_s, x);

  using ArgNorm = KokkosBatched::Norm::ScaledL2;
  auto info     = Functor_BatchedNrm<DeviceType, View2DType, NormViewType, ArgNorm, ArgMode>(x, norm).run();
  EXPECT_EQ(info, 0);

  // With strided views
  info = Functor_BatchedNrm<DeviceType, StridedView2DType, NormViewType, ArgNorm, ArgMode>(x_s, norm_s).run();
  EXPECT_EQ(info, 0);

  RealType eps  = 1.0e1 * ats::epsilon();
  auto h_norm   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, norm);
  auto h_norm_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, norm_s);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    KK_EXPECT_NEAR(h_norm(ib), h_norm_ref(ib), eps);
    KK_EXPECT_NEAR(h_norm_s(ib), h_norm_ref(ib), eps);
  }
}

/// \brief Implementation details of batched nrm test
///
/// \tparam DeviceType Kokkos device type
/// \tparam ScalarType Kokkos scalar type
/// \tparam LayoutType Kokkos layout type for the views
/// \tparam ArgNorm: one of L1, L2, Linf, ScaledL2
/// \tparam ArgMode: one of Mode::Serial, Mode::Team, Mode::TeamVector
///
/// \param[in] Nb Batch size of vectors
/// \param[in] N Length of the vector x
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ArgNorm, typename ArgMode>
void impl_test_batched_nrm(const std::size_t Nb, const std::size_t N) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using NormViewType      = Kokkos::View<RealType *, LayoutType, DeviceType>;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;

  View2DType x("x", Nb, N);
  NormViewType norm("norm", Nb), norm_s("norm_s", Nb), norm_ref("norm_ref", Nb);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, N, Nb * incx};
  StridedView2DType x_s("x_s", layout);

  // Create random x
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);

  // Deep copy to strided views
  Kokkos::deep_copy(x_s, x);

  auto info = Functor_BatchedNrm<DeviceType, View2DType, NormViewType, ArgNorm, ArgMode>(x, norm).run();
  EXPECT_EQ(info, 0);

  // With strided views
  info = Functor_BatchedNrm<DeviceType, StridedView2DType, NormViewType, ArgNorm, ArgMode>(x_s, norm_s).run();
  EXPECT_EQ(info, 0);

  auto asum = [](const ScalarType &val) {
    if constexpr (ats::is_complex) {
      return ats::abs(val.real()) + ats::abs(val.imag());
    } else {
      return ats::abs(val);
    }
  };

  // Make a reference at host
  auto h_x        = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, x);
  auto h_norm_ref = Kokkos::create_mirror_view(norm_ref);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    RealType norm_tmp = 0;
    if constexpr (std::is_same_v<ArgNorm, Norm::L1>) {
      for (std::size_t i = 0; i < N; i++) {
        norm_tmp += asum(h_x(ib, i));
      }
    } else if constexpr (std::is_same_v<ArgNorm, Norm::L2> || std::is_same_v<ArgNorm, Norm::ScaledL2>) {
      for (std::size_t i = 0; i < N; i++) {
        const auto abs_val = ats::abs(h_x(ib, i));
        norm_tmp += abs_val * abs_val;
      }
      norm_tmp = Kokkos::sqrt(norm_tmp);
    } else if constexpr (std::is_same_v<ArgNorm, Norm::LInf>) {
      for (std::size_t i = 0; i < N; i++) {
        const auto abs_val = ats::abs(h_x(ib, i));
        if (abs_val > norm_tmp) norm_tmp = abs_val;
      }
    }
    h_norm_ref(ib) = norm_tmp;
  }

  RealType eps  = 1.0e1 * ats::epsilon();
  auto h_norm   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, norm);
  auto h_norm_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, norm_s);

  // Check if norm is correct
  for (std::size_t ib = 0; ib < Nb; ib++) {
    KK_EXPECT_NEAR(h_norm(ib), h_norm_ref(ib), eps);
    KK_EXPECT_NEAR(h_norm_s(ib), h_norm_ref(ib), eps);
  }
}

}  // namespace Nrm
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ArgNorm, typename ArgMode>
int test_batched_nrm() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Nrm::impl_test_batched_nrm_analytical<DeviceType, ScalarType, LayoutType, ArgNorm, ArgMode>(1);
    Test::Nrm::impl_test_batched_nrm_analytical<DeviceType, ScalarType, LayoutType, ArgNorm, ArgMode>(2);
    for (int i = 0; i < 5; i++) {
      Test::Nrm::impl_test_batched_nrm<DeviceType, ScalarType, LayoutType, ArgNorm, ArgMode>(1, i);
      Test::Nrm::impl_test_batched_nrm<DeviceType, ScalarType, LayoutType, ArgNorm, ArgMode>(2, i);
    }

    if constexpr (std::is_same_v<ArgNorm, KokkosBatched::Norm::ScaledL2>) {
      // Additional test for large values to check for overflow handling in ScaledL2 norm
      Test::Nrm::impl_test_batched_nrm_overflow<DeviceType, ScalarType, LayoutType, ArgMode>(3);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Nrm::impl_test_batched_nrm_analytical<DeviceType, ScalarType, LayoutType, ArgNorm, ArgMode>(1);
    Test::Nrm::impl_test_batched_nrm_analytical<DeviceType, ScalarType, LayoutType, ArgNorm, ArgMode>(2);
    for (int i = 0; i < 5; i++) {
      Test::Nrm::impl_test_batched_nrm<DeviceType, ScalarType, LayoutType, ArgNorm, ArgMode>(1, i);
      Test::Nrm::impl_test_batched_nrm<DeviceType, ScalarType, LayoutType, ArgNorm, ArgMode>(2, i);
    }

    if constexpr (std::is_same_v<ArgNorm, KokkosBatched::Norm::ScaledL2>) {
      // Additional test for large values to check for overflow handling in ScaledL2 norm
      Test::Nrm::impl_test_batched_nrm_overflow<DeviceType, ScalarType, LayoutType, ArgMode>(3);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_nrm_l1_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::L1, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_nrm_l2_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::L2, KokkosBatched::Mode::Serial>();
}

TEST_F(TestCategory, test_batched_serial_nrm_linf_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::LInf, KokkosBatched::Mode::Serial>();
}

TEST_F(TestCategory, test_batched_serial_nrm_scaled_l2_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::ScaledL2, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_nrm_l1_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::L1, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_l2_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::L2, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_linf_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::LInf, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_scaled_l2_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::ScaledL2, KokkosBatched::Mode::Team>();
}

// TeamVector
TEST_F(TestCategory, test_batched_teamvector_nrm_l1_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::L1, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_l2_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::L2, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_linf_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::LInf, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_scaled_l2_float) {
  test_batched_nrm<TestDevice, float, KokkosBatched::Norm::ScaledL2, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_nrm_l1_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::L1, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_nrm_l2_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::L2, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_nrm_linf_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::LInf, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_nrm_scaled_l2_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::ScaledL2, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_nrm_l1_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::L1, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_l2_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::L2, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_linf_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::LInf, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_scaled_l2_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::ScaledL2, KokkosBatched::Mode::Team>();
}

// TeamVector
TEST_F(TestCategory, test_batched_teamvector_nrm_l1_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::L1, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_l2_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::L2, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_linf_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::LInf, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_scaled_l2_double) {
  test_batched_nrm<TestDevice, double, KokkosBatched::Norm::ScaledL2, KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
// Serial
TEST_F(TestCategory, test_batched_serial_nrm_l1_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::L1, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_nrm_l2_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::L2, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_nrm_linf_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::LInf, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_nrm_scaled_l2_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::ScaledL2, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_nrm_l1_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::L1, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_l2_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::L2, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_linf_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::LInf, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_scaled_l2_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::ScaledL2, KokkosBatched::Mode::Team>();
}

// TeamVector
TEST_F(TestCategory, test_batched_teamvector_nrm_l1_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::L1, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_l2_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::L2, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_linf_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::LInf, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_scaled_l2_fcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<float>, KokkosBatched::Norm::ScaledL2,
                   KokkosBatched::Mode::TeamVector>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
// Serial
TEST_F(TestCategory, test_batched_serial_nrm_l1_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::L1, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_nrm_l2_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::L2, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_nrm_linf_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::LInf, KokkosBatched::Mode::Serial>();
}
TEST_F(TestCategory, test_batched_serial_nrm_scaled_l2_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::ScaledL2, KokkosBatched::Mode::Serial>();
}

// Team
TEST_F(TestCategory, test_batched_team_nrm_l1_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::L1, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_l2_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::L2, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_linf_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::LInf, KokkosBatched::Mode::Team>();
}
TEST_F(TestCategory, test_batched_team_nrm_scaled_l2_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::ScaledL2, KokkosBatched::Mode::Team>();
}

// TeamVector
TEST_F(TestCategory, test_batched_teamvector_nrm_l1_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::L1, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_l2_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::L2, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_linf_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::LInf, KokkosBatched::Mode::TeamVector>();
}
TEST_F(TestCategory, test_batched_teamvector_nrm_scaled_l2_dcomplex) {
  test_batched_nrm<TestDevice, Kokkos::complex<double>, KokkosBatched::Norm::ScaledL2,
                   KokkosBatched::Mode::TeamVector>();
}
#endif
