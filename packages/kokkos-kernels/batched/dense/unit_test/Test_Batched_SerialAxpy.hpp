// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Kim Liegeois (knliege@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Axpy.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace Axpy {

template <typename DeviceType, typename ViewType, typename alphaViewType, bool KeepDim>
struct Functor_TestBatchedSerialAxpy {
  using execution_space = typename DeviceType::execution_space;
  const alphaViewType m_alpha;
  const ViewType m_X;
  const ViewType m_Y;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedSerialAxpy(const alphaViewType &alpha, const ViewType &X, const ViewType &Y)
      : m_alpha(alpha), m_X(X), m_Y(Y) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    if constexpr (KeepDim) {
      auto sub_alpha = Kokkos::subview(m_alpha, Kokkos::make_pair(k, k + 1));
      auto sub_x     = Kokkos::subview(m_X, Kokkos::make_pair(k, k + 1), Kokkos::ALL);
      auto sub_y     = Kokkos::subview(m_Y, Kokkos::make_pair(k, k + 1), Kokkos::ALL);
      KokkosBatched::SerialAxpy::invoke(sub_alpha, sub_x, sub_y);
    } else {
      auto alpha = m_alpha(k);
      auto sub_x = Kokkos::subview(m_X, k, Kokkos::ALL);
      auto sub_y = Kokkos::subview(m_Y, k, Kokkos::ALL);
      KokkosBatched::SerialAxpy::invoke(alpha, sub_x, sub_y);
    }
  }

  inline void run() {
    using value_type = typename ViewType::value_type;
    std::string name_region("KokkosBatched::Test::SerialAxpy");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_X.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

template <typename DeviceType, typename ViewType, typename alphaViewType>
void impl_test_batched_axpy(const int N, const int BlkSize) {
  using value_type             = typename ViewType::value_type;
  using const_value_type       = typename ViewType::const_value_type;
  using alpha_const_value_type = typename alphaViewType::const_value_type;
  using ats                    = KokkosKernels::ArithTraits<value_type>;

  ViewType X0("x0", N, BlkSize), X1("x1", N, BlkSize), Y0("y0", N, BlkSize), Y1("y1", N, BlkSize), Y2("y2", N, BlkSize);

  alphaViewType alpha("alpha", N);

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  Kokkos::fill_random(X0, random, const_value_type(1.0));
  Kokkos::fill_random(Y0, random, const_value_type(1.0));
  Kokkos::fill_random(alpha, random, alpha_const_value_type(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(X1, X0);
  Kokkos::deep_copy(Y1, Y0);
  Kokkos::deep_copy(Y2, Y0);

  /// test body
  auto alpha_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, alpha);
  auto X0_host    = Kokkos::create_mirror_view(X0);
  auto Y0_host    = Kokkos::create_mirror_view(Y0);

  Kokkos::deep_copy(X0_host, X0);
  Kokkos::deep_copy(Y0_host, Y0);

  for (int l = 0; l < N; ++l)
    for (int i = 0; i < BlkSize; ++i) Y0_host(l, i) += alpha_host(l) * X0_host(l, i);

  Functor_TestBatchedSerialAxpy<DeviceType, ViewType, alphaViewType, true>(alpha, X1, Y1).run();
  Functor_TestBatchedSerialAxpy<DeviceType, ViewType, alphaViewType, false>(alpha, X1, Y2).run();

  Kokkos::fence();

  /// for comparison send it to host
  auto Y1_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Y1);
  auto Y2_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Y2);

  /// check c0 = c1 ; this eps is about 10^-14
  using mag_type = typename ats::mag_type;
  mag_type sum(1), diff1(0), diff2(0);
  const mag_type eps = 1.0e3 * ats::epsilon();

  for (int l = 0; l < N; ++l)
    for (int i = 0; i < BlkSize; ++i) {
      sum += ats::abs(Y0_host(l, i));
      diff1 += ats::abs(Y0_host(l, i) - Y1_host(l, i));
      diff2 += ats::abs(Y0_host(l, i) - Y2_host(l, i));
    }
  EXPECT_NEAR_KK(diff1 / sum, 0, eps);
  EXPECT_NEAR_KK(diff2 / sum, 0, eps);
}
}  // namespace Axpy
}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType>
int test_batched_axpy() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using ViewType      = Kokkos::View<ValueType **, Kokkos::LayoutLeft, DeviceType>;
    using alphaViewType = Kokkos::View<ScalarType *, Kokkos::LayoutLeft, DeviceType>;

    for (int i = 3; i < 10; ++i) {
      Test::Axpy::impl_test_batched_axpy<DeviceType, ViewType, alphaViewType>(1024, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using ViewType      = Kokkos::View<ValueType **, Kokkos::LayoutRight, DeviceType>;
    using alphaViewType = Kokkos::View<ScalarType *, Kokkos::LayoutRight, DeviceType>;

    for (int i = 3; i < 10; ++i) {
      Test::Axpy::impl_test_batched_axpy<DeviceType, ViewType, alphaViewType>(1024, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_axpy_nt_float_float) { test_batched_axpy<TestDevice, float, float>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_axpy_nt_double_double) { test_batched_axpy<TestDevice, double, double>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_axpy_nt_fcomplex_fcomplex) {
  test_batched_axpy<TestDevice, Kokkos::complex<float>, Kokkos::complex<float>>();
}

TEST_F(TestCategory, batched_scalar_serial_axpy_nt_fcomplex_float) {
  test_batched_axpy<TestDevice, Kokkos::complex<float>, float>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_axpy_nt_dcomplex_dcomplex) {
  test_batched_axpy<TestDevice, Kokkos::complex<double>, Kokkos::complex<double>>();
}

TEST_F(TestCategory, batched_scalar_serial_axpy_nt_dcomplex_double) {
  test_batched_axpy<TestDevice, Kokkos::complex<double>, double>();
}
#endif
