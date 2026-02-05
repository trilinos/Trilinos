// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Lacgv.hpp>

namespace Test {
namespace Lacgv {

template <typename DeviceType, typename XViewType>
struct Functor_BatchedSerialLacgv {
  using execution_space = typename DeviceType::execution_space;
  XViewType m_x;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialLacgv(const XViewType &x) : m_x(x) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    auto sub_x = Kokkos::subview(m_x, k, Kokkos::ALL());
    info += KokkosBatched::SerialLacgv::invoke(sub_x);
  }

  inline int run() {
    using value_type = typename XViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialLacgv");
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

/// \brief Implementation details of batched lacgv analytical test
///
/// \param Nb [in] Batch size of vectors
///        1D complex vector
///        x: [1 + 1j, -3+2j, -2-2j,  0+1j]
///        conj(x): [1 - 1j, -3-2j, -2+2j,  0-1j]
template <typename DeviceType, typename ScalarType, typename LayoutType>
void impl_test_batched_lacgv_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;

  const std::size_t BlkSize = 4;
  View2DType x("x", Nb, BlkSize), x_ref("x_ref", Nb, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, BlkSize, Nb * incx};
  StridedView2DType x_s("x_s", layout);

  // Initialize a vector x
  auto h_x     = Kokkos::create_mirror_view(x);
  auto h_x_ref = Kokkos::create_mirror_view(x_ref);

  constexpr bool is_complex = KokkosKernels::ArithTraits<ScalarType>::is_complex;

  for (std::size_t ib = 0; ib < Nb; ib++) {
    if constexpr (is_complex) {
      h_x(ib, 0) = ScalarType(1.0, 1.0);
      h_x(ib, 1) = ScalarType(-3.0, 2.0);
      h_x(ib, 2) = ScalarType(-2.0, -2.0);
      h_x(ib, 3) = ScalarType(0.0, 1.0);

      h_x_ref(ib, 0) = ScalarType(1.0, -1.0);
      h_x_ref(ib, 1) = ScalarType(-3.0, -2.0);
      h_x_ref(ib, 2) = ScalarType(-2.0, 2.0);
      h_x_ref(ib, 3) = ScalarType(0.0, -1.0);
    } else {
      h_x(ib, 0) = 1.0;
      h_x(ib, 1) = -3.0;
      h_x(ib, 2) = -2.0;
      h_x(ib, 3) = 0.0;

      h_x_ref(ib, 0) = 1.0;
      h_x_ref(ib, 1) = -3.0;
      h_x_ref(ib, 2) = -2.0;
      h_x_ref(ib, 3) = 0.0;
    }
  }

  Kokkos::deep_copy(x, h_x);
  Kokkos::deep_copy(x_s, x);

  auto info0 = Functor_BatchedSerialLacgv<DeviceType, View2DType>(x).run();

  // With strided Views
  auto info1 = Functor_BatchedSerialLacgv<DeviceType, StridedView2DType>(x_s).run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  RealType eps = 1.0e1 * ats::epsilon();
  Kokkos::deep_copy(h_x, x);
  // Check if x is conjugated
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }

  // Reuse x to compare x_s and x_ref
  Kokkos::deep_copy(x, x_s);
  Kokkos::deep_copy(h_x, x);
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched lacgv test
///
/// \param Nb [in] Batch size of vectors
/// \param BlkSize [in] Length of vector X
template <typename DeviceType, typename ScalarType, typename LayoutType>
void impl_test_batched_lacgv(const std::size_t Nb, const std::size_t BlkSize) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;

  View2DType x("x", Nb, BlkSize), x_ref("x_ref", Nb, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, BlkSize, Nb * incx};
  StridedView2DType x_s("x_s", layout);

  // Create a random vector x
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(x_ref, x);

  // Deep copy to strided views
  Kokkos::deep_copy(x_s, x);

  auto info0 = Functor_BatchedSerialLacgv<DeviceType, View2DType>(x).run();

  // With strided Views
  auto info1 = Functor_BatchedSerialLacgv<DeviceType, StridedView2DType>(x_s).run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  // Make a reference at host
  auto h_x_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x_ref);
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      h_x_ref(ib, i) = ats::conj(h_x_ref(ib, i));
    }
  }

  RealType eps = 1.0e1 * ats::epsilon();
  auto h_x     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
  // Check if x is conjugated
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }

  // Reuse x to compare x_s and x_ref
  Kokkos::deep_copy(x, x_s);
  Kokkos::deep_copy(h_x, x);
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

}  // namespace Lacgv
}  // namespace Test

template <typename DeviceType, typename ScalarType>
int test_batched_lacgv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Lacgv::impl_test_batched_lacgv_analytical<DeviceType, ScalarType, LayoutType>(1);
    Test::Lacgv::impl_test_batched_lacgv_analytical<DeviceType, ScalarType, LayoutType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Lacgv::impl_test_batched_lacgv<DeviceType, ScalarType, LayoutType>(1, i);
      Test::Lacgv::impl_test_batched_lacgv<DeviceType, ScalarType, LayoutType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Lacgv::impl_test_batched_lacgv_analytical<DeviceType, ScalarType, LayoutType>(1);
    Test::Lacgv::impl_test_batched_lacgv_analytical<DeviceType, ScalarType, LayoutType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Lacgv::impl_test_batched_lacgv<DeviceType, ScalarType, LayoutType>(1, i);
      Test::Lacgv::impl_test_batched_lacgv<DeviceType, ScalarType, LayoutType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_lacgv_float) { test_batched_lacgv<TestDevice, float>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_lacgv_double) { test_batched_lacgv<TestDevice, double>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_lacgv_fcomplex) { test_batched_lacgv<TestDevice, Kokkos::complex<float>>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_lacgv_dcomplex) { test_batched_lacgv<TestDevice, Kokkos::complex<double>>(); }
#endif
