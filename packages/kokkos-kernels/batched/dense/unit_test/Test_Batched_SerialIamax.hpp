// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Iamax.hpp"

using namespace KokkosBatched;

namespace Test {
namespace Iamax {

template <typename DeviceType, typename XViewType, typename RViewType>
struct Functor_BatchedSerialIamax {
  using execution_space = typename DeviceType::execution_space;
  XViewType m_x;
  RViewType m_r;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialIamax(const XViewType &x, const RViewType &r) : m_x(x), m_r(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto sub_x = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto iamax = KokkosBatched::SerialIamax::invoke(sub_x);
    m_r(k)     = static_cast<int>(iamax);
  }

  inline void run() {
    using value_type = typename XViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialIamax");
    std::string name_value_type = Test::value_type_name<value_type>();
    std::string name            = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_x.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

/// \brief Implementation details of batched iamax analytical test
///        A0: [1, 2, 0] -> 1
///        A1: [-5, 4, 3] -> 0
///        A2: [0, 0, 0] -> 0
///        A3: [0, -1, -1] -> 1
///
/// \param N [in] Batch size of A
template <typename DeviceType, typename ScalarType, typename LayoutType>
void impl_test_batched_iamax_analytical(const std::size_t N) {
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using MaxView1DType     = Kokkos::View<int *, LayoutType, DeviceType>;

  View2DType A0("A0", N, 3), A1("A1", N, 3), A2("A2", N, 3), A3("A3", N, 3);
  MaxView1DType iamax0("iamax0", N), iamax_ref0("iamax_ref0", N), iamax1("iamax1", N), iamax_ref1("iamax_ref1", N),
      iamax2("iamax2", N), iamax_ref2("iamax_ref2", N), iamax3("iamax3", N), iamax_ref3("iamax_ref3", N);

  // Testing incx argument with strided views
  constexpr std::size_t incx = 2;
  Kokkos::LayoutStride layout{N, incx, 3, N * incx};
  StridedView2DType A0_s("A0_s", layout), A1_s("A1_s", layout), A2_s("A2_s", layout), A3_s("A3_s", layout);
  MaxView1DType iamax_s0("iamax_s0", N), iamax_s1("iamax_s1", N), iamax_s2("iamax_s2", N), iamax_s3("iamax_s3", N);

  // Initialize A0, A1, A2, A3
  auto h_A0 = Kokkos::create_mirror_view(A0);
  auto h_A1 = Kokkos::create_mirror_view(A1);
  auto h_A2 = Kokkos::create_mirror_view(A2);
  auto h_A3 = Kokkos::create_mirror_view(A3);

  auto h_iamax_ref0 = Kokkos::create_mirror_view(iamax_ref0);
  auto h_iamax_ref1 = Kokkos::create_mirror_view(iamax_ref1);
  auto h_iamax_ref2 = Kokkos::create_mirror_view(iamax_ref2);
  auto h_iamax_ref3 = Kokkos::create_mirror_view(iamax_ref3);
  for (std::size_t k = 0; k < N; k++) {
    h_A0(k, 0) = 1;
    h_A0(k, 1) = 2;
    h_A0(k, 2) = 0;

    h_A1(k, 0) = -5;
    h_A1(k, 1) = 4;
    h_A1(k, 2) = 3;

    h_A2(k, 0) = 0;
    h_A2(k, 1) = 0;
    h_A2(k, 2) = 0;

    h_A3(k, 0) = 0;
    h_A3(k, 1) = -1;
    h_A3(k, 2) = -1;

    h_iamax_ref0(k) = 1;
    h_iamax_ref1(k) = 0;
    h_iamax_ref2(k) = 0;
    h_iamax_ref3(k) = 1;
  }
  Kokkos::deep_copy(A0, h_A0);
  Kokkos::deep_copy(A1, h_A1);
  Kokkos::deep_copy(A2, h_A2);
  Kokkos::deep_copy(A3, h_A3);

  // Strided view can be copied only on the same device
  Kokkos::deep_copy(A0_s, A0);
  Kokkos::deep_copy(A1_s, A1);
  Kokkos::deep_copy(A2_s, A2);
  Kokkos::deep_copy(A3_s, A3);

  Functor_BatchedSerialIamax<DeviceType, View2DType, MaxView1DType>(A0, iamax0).run();
  Functor_BatchedSerialIamax<DeviceType, View2DType, MaxView1DType>(A1, iamax1).run();
  Functor_BatchedSerialIamax<DeviceType, View2DType, MaxView1DType>(A2, iamax2).run();
  Functor_BatchedSerialIamax<DeviceType, View2DType, MaxView1DType>(A3, iamax3).run();

  // For strided views
  Functor_BatchedSerialIamax<DeviceType, StridedView2DType, MaxView1DType>(A0_s, iamax_s0).run();
  Functor_BatchedSerialIamax<DeviceType, StridedView2DType, MaxView1DType>(A1_s, iamax_s1).run();
  Functor_BatchedSerialIamax<DeviceType, StridedView2DType, MaxView1DType>(A2_s, iamax_s2).run();
  Functor_BatchedSerialIamax<DeviceType, StridedView2DType, MaxView1DType>(A3_s, iamax_s3).run();

  Kokkos::fence();

  // Copy to host for comparison
  auto h_iamax0   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iamax0);
  auto h_iamax1   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iamax1);
  auto h_iamax2   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iamax2);
  auto h_iamax3   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iamax3);
  auto h_iamax_s0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iamax_s0);
  auto h_iamax_s1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iamax_s1);
  auto h_iamax_s2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iamax_s2);
  auto h_iamax_s3 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iamax_s3);

  // Check if max index is correct
  for (std::size_t k = 0; k < N; k++) {
    EXPECT_EQ(h_iamax0(k), h_iamax_ref0(k));
    EXPECT_EQ(h_iamax1(k), h_iamax_ref1(k));
    EXPECT_EQ(h_iamax2(k), h_iamax_ref2(k));
    EXPECT_EQ(h_iamax3(k), h_iamax_ref3(k));
    EXPECT_EQ(h_iamax_s0(k), h_iamax_ref0(k));
    EXPECT_EQ(h_iamax_s1(k), h_iamax_ref1(k));
    EXPECT_EQ(h_iamax_s2(k), h_iamax_ref2(k));
    EXPECT_EQ(h_iamax_s3(k), h_iamax_ref3(k));
  }
}

/// \brief Implementation details of batched pbtrs test
///
/// \param N [in] Batch size of RHS (banded matrix can also be batched matrix)
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType>
void impl_test_batched_iamax(const std::size_t N, const std::size_t BlkSize) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using MaxView1DType     = Kokkos::View<int *, LayoutType, DeviceType>;

  View2DType A("A", N, BlkSize);
  MaxView1DType iamax("iamax", N), iamax_ref("iamax_ref", N);

  // Testing incx argument with strided views
  constexpr std::size_t incx = 2;
  Kokkos::LayoutStride layout{N, incx, BlkSize, N * incx};
  StridedView2DType A_s("A_s", layout);
  MaxView1DType iamax_s("iamax_s", N);

  // Initialize A with random values
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);

  // Strided view can be copied only on the same device
  Kokkos::deep_copy(A_s, A);

  Functor_BatchedSerialIamax<DeviceType, View2DType, MaxView1DType>(A, iamax).run();

  // For strided views
  Functor_BatchedSerialIamax<DeviceType, StridedView2DType, MaxView1DType>(A_s, iamax_s).run();

  Kokkos::fence();

  // Reference
  auto h_iamax_ref = Kokkos::create_mirror_view(iamax_ref);
  if (BlkSize == 0) {
    // As well as blas, we store 0 (0 in Fortran) for empty matrix
    for (std::size_t k = 0; k < N; k++) {
      h_iamax_ref(k) = 0;
    }
  } else {
    auto h_A = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
    for (std::size_t k = 0; k < N; k++) {
      RealType amax = Kokkos::abs(h_A(k, 0));
      int iamax_tmp = 0;
      for (std::size_t i = 1; i < BlkSize; i++) {
        const RealType abs_A_i = Kokkos::abs(h_A(k, i));
        if (abs_A_i > amax) {
          amax      = abs_A_i;
          iamax_tmp = static_cast<int>(i);
        }
      }
      h_iamax_ref(k) = iamax_tmp;
    }
  }

  // Copy to host for comparison
  auto h_iamax   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iamax);
  auto h_iamax_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iamax_s);

  // Check if max index is correct
  for (std::size_t k = 0; k < N; k++) {
    EXPECT_EQ(h_iamax(k), h_iamax_ref(k));
    EXPECT_EQ(h_iamax_s(k), h_iamax_ref(k));
  }
}

}  // namespace Iamax
}  // namespace Test

template <typename DeviceType, typename ScalarType>
int test_batched_iamax() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Iamax::impl_test_batched_iamax_analytical<DeviceType, ScalarType, LayoutType>(1);
    Test::Iamax::impl_test_batched_iamax_analytical<DeviceType, ScalarType, LayoutType>(2);
    for (std::size_t i = 0; i < 10; i++) {
      Test::Iamax::impl_test_batched_iamax<DeviceType, ScalarType, LayoutType>(1, i);
      Test::Iamax::impl_test_batched_iamax<DeviceType, ScalarType, LayoutType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Iamax::impl_test_batched_iamax_analytical<DeviceType, ScalarType, LayoutType>(1);
    Test::Iamax::impl_test_batched_iamax_analytical<DeviceType, ScalarType, LayoutType>(2);
    for (std::size_t i = 0; i < 10; i++) {
      Test::Iamax::impl_test_batched_iamax<DeviceType, ScalarType, LayoutType>(1, i);
      Test::Iamax::impl_test_batched_iamax<DeviceType, ScalarType, LayoutType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_iamax_float) { test_batched_iamax<TestDevice, float>(); }
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_iamax_double) { test_batched_iamax<TestDevice, double>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_iamax_fcomplex) { test_batched_iamax<TestDevice, Kokkos::complex<float>>(); }
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_iamax_dcomplex) { test_batched_iamax<TestDevice, Kokkos::complex<double>>(); }
#endif
