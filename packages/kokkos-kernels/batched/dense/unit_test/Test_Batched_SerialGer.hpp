// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Ger.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Ger {

template <typename T>
struct ParamTag {
  using trans = T;
};

template <typename DeviceType, typename XViewType, typename YViewType, typename AViewType, typename ScalarType,
          typename ParamTagType>
struct Functor_BatchedSerialGer {
  using execution_space = typename DeviceType::execution_space;
  XViewType m_x;
  YViewType m_y;
  AViewType m_A;
  ScalarType m_alpha;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGer(const ScalarType alpha, const XViewType &x, const YViewType &y, const AViewType &A)
      : m_x(x), m_y(y), m_A(A), m_alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    auto sub_x = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto sub_y = Kokkos::subview(m_y, k, Kokkos::ALL());
    auto sub_A = Kokkos::subview(m_A, k, Kokkos::ALL(), Kokkos::ALL());

    info += KokkosBatched::SerialGer<typename ParamTagType::trans>::invoke(m_alpha, sub_x, sub_y, sub_A);
  }

  inline int run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGer");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_A.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

/// \brief Implementation details of batched ger analytical test
///        3x4 matrix
///        A: [[1, -3, -2,  0],
///            [-1, 1, -3, -2],
///            [2, -1,  1, -3]]
///        x0: [1, 2, 3]
///        y:  [0, 1, 2, 3]
///        Ref: [[ 1.,  -1.5,  1.,   4.5],
///              [-1.,   4.,   3.,   7. ],
///              [ 2.,   3.5, 10.,  10.5],]
///
///        4x4 matrix
///        A: [[1, -3, -2,  0],
///            [-1, 1, -3, -2],
///            [2, -1,  1, -3],
///            [0,  2, -1,  1]]
///        x1: [1, 2, 3, 4]
///        y:  [0, 1, 2, 3]
///        Ref: [[ 1.,  -1.5,  1.,   4.5],
///              [-1.,   4.,   3.,   7. ],
///              [ 2.,   3.5, 10.,  10.5],
///              [ 0.,   8.,  11.,  19. ]]
///
///        5x4 matrix
///        A: [[1, -3, -2,  0],
///            [-1, 1, -3, -2],
///            [2, -1,  1, -3],
///            [0,  2, -1,  1],
///            [0,  0,  2, -1]]
///        x1: [1, 2, 3, 4, 5]
///        y:  [0, 1, 2, 3]
///        Ref: [[ 1.,  -1.5,  1.,   4.5],
///              [-1.,   4.,   3.,   7. ],
///              [ 2.,   3.5, 10.,  10.5],
///              [ 0.,   8.,  11.,  19. ],
///              [ 0.,   7.5, 17.,  21.5]]
///
/// \param Nb [in] Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_ger_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  const std::size_t BlkSize = 4, BlkSize_s = 3, BlkSize_l = 5;
  View3DType A0("A0", Nb, BlkSize_s, BlkSize), A0_s("A0_s", Nb, BlkSize_s, BlkSize),
      A0_ref("A0_ref", Nb, BlkSize_s, BlkSize);
  View3DType A1_s("A1_s", Nb, BlkSize, BlkSize), A1("A1", Nb, BlkSize, BlkSize), A1_ref("A1_ref", Nb, BlkSize, BlkSize);
  View3DType A2("A2", Nb, BlkSize_l, BlkSize), A2_s("A2_s", Nb, BlkSize_l, BlkSize),
      A2_ref("A2_ref", Nb, BlkSize_l, BlkSize);

  View2DType x0("x0", Nb, BlkSize_s), x1("x1", Nb, BlkSize), x2("x2", Nb, BlkSize_l), y("y", Nb, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout0{Nb, incx, BlkSize_s, Nb * incx}, layout1{Nb, incx, BlkSize, Nb * incx},
      layout2{Nb, incx, BlkSize_l, Nb * incx};
  StridedView2DType x0_s("x0_s", layout0), x1_s("x1_s", layout1), x2_s("x2_s", layout2), y_s("y_s", layout1);

  // Only filling x2, A2 and deep_copy from its subview
  auto h_A2     = Kokkos::create_mirror_view(A2);
  auto h_A2_ref = Kokkos::create_mirror_view(A2_ref);
  auto h_x2     = Kokkos::create_mirror_view(x2);
  auto h_y      = Kokkos::create_mirror_view(y);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    h_A2(ib, 0, 0) = 1;
    h_A2(ib, 0, 1) = -3;
    h_A2(ib, 0, 2) = -2;
    h_A2(ib, 0, 3) = 0;
    h_A2(ib, 1, 0) = -1;
    h_A2(ib, 1, 1) = 1;
    h_A2(ib, 1, 2) = -3;
    h_A2(ib, 1, 3) = -2;
    h_A2(ib, 2, 0) = 2;
    h_A2(ib, 2, 1) = -1;
    h_A2(ib, 2, 2) = 1;
    h_A2(ib, 2, 3) = -3;
    h_A2(ib, 3, 0) = 0;
    h_A2(ib, 3, 1) = 2;
    h_A2(ib, 3, 2) = -1;
    h_A2(ib, 3, 3) = 1;
    h_A2(ib, 4, 2) = 2;
    h_A2(ib, 4, 3) = -1;

    h_A2_ref(ib, 0, 0) = 1;
    h_A2_ref(ib, 0, 1) = -1.5;
    h_A2_ref(ib, 0, 2) = 1;
    h_A2_ref(ib, 0, 3) = 4.5;
    h_A2_ref(ib, 1, 0) = -1;
    h_A2_ref(ib, 1, 1) = 4;
    h_A2_ref(ib, 1, 2) = 3;
    h_A2_ref(ib, 1, 3) = 7;
    h_A2_ref(ib, 2, 0) = 2;
    h_A2_ref(ib, 2, 1) = 3.5;
    h_A2_ref(ib, 2, 2) = 10;
    h_A2_ref(ib, 2, 3) = 10.5;
    h_A2_ref(ib, 3, 0) = 0;
    h_A2_ref(ib, 3, 1) = 8;
    h_A2_ref(ib, 3, 2) = 11;
    h_A2_ref(ib, 3, 3) = 19;
    h_A2_ref(ib, 4, 1) = 7.5;
    h_A2_ref(ib, 4, 2) = 17;
    h_A2_ref(ib, 4, 3) = 21.5;

    for (std::size_t i = 0; i < BlkSize_l; i++) {
      h_x2(ib, i) = i + 1;
    }

    for (std::size_t j = 0; j < BlkSize; j++) {
      h_y(ib, j) = j;
    }
  }

  Kokkos::deep_copy(A2, h_A2);
  Kokkos::deep_copy(x2, h_x2);
  Kokkos::deep_copy(y, h_y);

  auto A2_m3       = Kokkos::subview(A2, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize_s), Kokkos::ALL);
  auto A2_m4       = Kokkos::subview(A2, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize), Kokkos::ALL);
  auto h_A2_ref_m3 = Kokkos::subview(h_A2_ref, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize_s), Kokkos::ALL);
  auto h_A2_ref_m4 = Kokkos::subview(h_A2_ref, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize), Kokkos::ALL);

  Kokkos::deep_copy(A0, A2_m3);  // Extract 3x4 matrix
  Kokkos::deep_copy(A1, A2_m4);  // Extract 4x4 matrix

  auto h_A0_ref = Kokkos::create_mirror_view(A0_ref);
  auto h_A1_ref = Kokkos::create_mirror_view(A1_ref);
  Kokkos::deep_copy(h_A0_ref, h_A2_ref_m3);  // Extract 3x4 matrix
  Kokkos::deep_copy(h_A1_ref, h_A2_ref_m4);  // Extract 4x4 matrix

  auto x2_m3 = Kokkos::subview(x2, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize_s));
  auto x2_m4 = Kokkos::subview(x2, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize));
  Kokkos::deep_copy(x0, x2_m3);
  Kokkos::deep_copy(x1, x2_m4);

  // Deep copy to strided views
  Kokkos::deep_copy(A0_s, A0);
  Kokkos::deep_copy(A1_s, A1);
  Kokkos::deep_copy(A2_s, A2);
  Kokkos::deep_copy(x0_s, x0);
  Kokkos::deep_copy(x1_s, x1);
  Kokkos::deep_copy(x2_s, x2);
  Kokkos::deep_copy(y_s, y);

  const ScalarType alpha = 1.5;

  auto info0 = Functor_BatchedSerialGer<DeviceType, View2DType, View2DType, View3DType, ScalarType, ParamTagType>(
                   alpha, x0, y, A0)
                   .run();
  auto info1 = Functor_BatchedSerialGer<DeviceType, View2DType, View2DType, View3DType, ScalarType, ParamTagType>(
                   alpha, x1, y, A1)
                   .run();
  auto info2 = Functor_BatchedSerialGer<DeviceType, View2DType, View2DType, View3DType, ScalarType, ParamTagType>(
                   alpha, x2, y, A2)
                   .run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);
  EXPECT_EQ(info2, 0);

  // For strided views
  info0 =
      Functor_BatchedSerialGer<DeviceType, StridedView2DType, StridedView2DType, View3DType, ScalarType, ParamTagType>(
          alpha, x0_s, y_s, A0_s)
          .run();
  info1 =
      Functor_BatchedSerialGer<DeviceType, StridedView2DType, StridedView2DType, View3DType, ScalarType, ParamTagType>(
          alpha, x1_s, y_s, A1_s)
          .run();
  info2 =
      Functor_BatchedSerialGer<DeviceType, StridedView2DType, StridedView2DType, View3DType, ScalarType, ParamTagType>(
          alpha, x2_s, y_s, A2_s)
          .run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);
  EXPECT_EQ(info2, 0);

  RealType eps = 1.0e1 * ats::epsilon();
  Kokkos::deep_copy(h_A2, A2);
  auto h_A0   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A0);
  auto h_A1   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A1);
  auto h_A0_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A0_s);
  auto h_A1_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A1_s);
  auto h_A2_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A2_s);

  // Check if A:= alpha * x * y**T + A
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize_s; i++) {
      for (std::size_t j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_A0(ib, i, j), h_A0_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_A0_s(ib, i, j), h_A0_ref(ib, i, j), eps);
      }
    }
    for (std::size_t i = 0; i < BlkSize; i++) {
      for (std::size_t j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_A1(ib, i, j), h_A1_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_A1_s(ib, i, j), h_A1_ref(ib, i, j), eps);
      }
    }
    for (std::size_t i = 0; i < BlkSize_l; i++) {
      for (std::size_t j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_A2(ib, i, j), h_A2_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_A2_s(ib, i, j), h_A2_ref(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched ger test
///
/// \param N [in] Batch size of matrices
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_ger(const std::size_t N, const std::size_t BlkSize) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  View3DType A("A", N, BlkSize, BlkSize), A0("A0", N, BlkSize, BlkSize), A_s("A_s", N, BlkSize, BlkSize),
      A0_s("A0_s", N, BlkSize, BlkSize), A_ref("A_ref", N, BlkSize, BlkSize), A0_ref("A0_ref", N, BlkSize, BlkSize);

  View2DType x("x", N, BlkSize), y("y", N, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{N, incx, BlkSize, N * incx};
  StridedView2DType x_s("x_s", layout), y_s("y_s", layout);

  // Create a random matrix A and make it Positive Definite Symmetric
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  // Initialize A, x, y with random matrix
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);
  Kokkos::fill_random(y, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(A_ref, A);

  Kokkos::deep_copy(A_s, A);
  Kokkos::deep_copy(x_s, x);
  Kokkos::deep_copy(y_s, y);

  // When A0 is zero
  const ScalarType alpha = 1.5;
  auto info0 = Functor_BatchedSerialGer<DeviceType, View2DType, View2DType, View3DType, ScalarType, ParamTagType>(
                   alpha, x, y, A0)
                   .run();

  // When A is a random matrix
  auto info1 =
      Functor_BatchedSerialGer<DeviceType, View2DType, View2DType, View3DType, ScalarType, ParamTagType>(alpha, x, y, A)
          .run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  // With strided Views
  info0 =
      Functor_BatchedSerialGer<DeviceType, StridedView2DType, StridedView2DType, View3DType, ScalarType, ParamTagType>(
          alpha, x_s, y_s, A0_s)
          .run();

  // When A is a random matrix
  info1 =
      Functor_BatchedSerialGer<DeviceType, StridedView2DType, StridedView2DType, View3DType, ScalarType, ParamTagType>(
          alpha, x_s, y_s, A_s)
          .run();

  // Make a reference at host
  auto h_x      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
  auto h_y      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  auto h_A_ref  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_ref);
  auto h_A0_ref = Kokkos::create_mirror_view(Kokkos::HostSpace(), A0_ref);

  const bool is_conj = std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>;
  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t j = 0; j < BlkSize; j++) {
      if (h_y(ib, j) != 0) {
        auto temp = is_conj ? alpha * KokkosKernels::ArithTraits<ScalarType>::conj(h_y(ib, j)) : alpha * h_y(ib, j);
        for (std::size_t i = 0; i < BlkSize; i++) {
          h_A_ref(ib, i, j)  = h_A_ref(ib, i, j) + h_x(ib, i) * temp;
          h_A0_ref(ib, i, j) = h_x(ib, i) * temp;
        }
      }
    }
  }

  RealType eps = 1.0e1 * ats::epsilon();

  auto h_A    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  auto h_A0   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A0);
  auto h_A_s  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_s);
  auto h_A0_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A0_s);

  // Check if A:= alpha * x * y**T + A or A:= alpha * x * y**H + A
  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      for (std::size_t j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_A(ib, i, j), h_A_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_A0(ib, i, j), h_A0_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_A_s(ib, i, j), h_A_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_A0_s(ib, i, j), h_A0_ref(ib, i, j), eps);
      }
    }
  }
}

}  // namespace Ger
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType>
int test_batched_ger() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Ger::impl_test_batched_ger_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Ger::impl_test_batched_ger_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Ger::impl_test_batched_ger<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Ger::impl_test_batched_ger<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Ger::impl_test_batched_ger_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Ger::impl_test_batched_ger_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Ger::impl_test_batched_ger<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Ger::impl_test_batched_ger<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_ger_t_float) {
  using param_tag_type = ::Test::Ger::ParamTag<Trans::Transpose>;
  test_batched_ger<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_ger_c_float) {
  using param_tag_type = ::Test::Ger::ParamTag<Trans::ConjTranspose>;
  test_batched_ger<TestDevice, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_ger_t_double) {
  using param_tag_type = ::Test::Ger::ParamTag<Trans::Transpose>;
  test_batched_ger<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_ger_c_double) {
  using param_tag_type = ::Test::Ger::ParamTag<Trans::ConjTranspose>;
  test_batched_ger<TestDevice, double, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_ger_t_fcomplex) {
  using param_tag_type = ::Test::Ger::ParamTag<Trans::Transpose>;
  test_batched_ger<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_ger_c_fcomplex) {
  using param_tag_type = ::Test::Ger::ParamTag<Trans::ConjTranspose>;
  test_batched_ger<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_ger_t_dcomplex) {
  using param_tag_type = ::Test::Ger::ParamTag<Trans::Transpose>;
  test_batched_ger<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_ger_c_dcomplex) {
  using param_tag_type = ::Test::Ger::ParamTag<Trans::ConjTranspose>;
  test_batched_ger<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
#endif
