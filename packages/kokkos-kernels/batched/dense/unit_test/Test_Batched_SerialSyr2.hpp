// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Syr2.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Syr2 {

template <typename U, typename T>
struct ParamTag {
  using uplo  = U;
  using trans = T;
};

template <typename DeviceType, typename XViewType, typename YViewType, typename AViewType, typename ScalarType,
          typename ParamTagType>
struct Functor_BatchedSerialSyr2 {
  using execution_space = typename DeviceType::execution_space;
  XViewType m_x;
  YViewType m_y;
  AViewType m_A;
  ScalarType m_alpha;

  Functor_BatchedSerialSyr2(const ScalarType alpha, const XViewType &x, const YViewType &y, const AViewType &A)
      : m_x(x), m_y(y), m_A(A), m_alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    auto sub_x = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto sub_y = Kokkos::subview(m_y, k, Kokkos::ALL());
    auto sub_A = Kokkos::subview(m_A, k, Kokkos::ALL(), Kokkos::ALL());

    info += KokkosBatched::SerialSyr2<typename ParamTagType::uplo, typename ParamTagType::trans>::invoke(m_alpha, sub_x,
                                                                                                         sub_y, sub_A);
  }

  inline int run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialSyr2");
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

/// \brief Implementation details of batched syr2 analytical test
///        to confirm A:= alpha*x*y**T + alpha*y*x**T + A is computed correctly
///        alpha = 1.5
///        4x4 matrix (upper)
///        U: [[1, -3, -2,  0],
///            [0,  1, -3, -2],
///            [0,  0,  1, -3],
///            [0,  0,  0,  1]]
///        x: [1, 2, 3, 4]
///        y: [4, 3, 2, 1]
///        Ref: [[13.,  13.5,  19.,  25.5],
///              [ 0.,  19.,  16.5,  19., ],
///              [ 0.,   0.,  19.,  13.5, ],
///              [ 0.,   0.,   0.,  13.,  ]]
///
///        4x4 matrix (lower)
///        L: [[1,  0,  0,  0],
///            [-1, 1,  0,  0],
///            [2, -1,  1,  0],
///            [0,  2, -1,  1]]
///        x: [1, 2, 3, 4]
///        y: [4, 3, 2, 1]
///        Ref: [[13.,   0.,   0.,   0., ],
///              [15.5, 19.,   0.,   0., ],
///              [23.,  18.5, 19.,   0., ],
///              [25.5, 23.,  15.5, 13., ]]
///
/// \param Nb [in] Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_syr2_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo           = typename ParamTagType::uplo;

  const std::size_t BlkSize = 4;
  View3DType A("A", Nb, BlkSize, BlkSize), A_s("A_s", Nb, BlkSize, BlkSize), A_ref("A_ref", Nb, BlkSize, BlkSize);
  View2DType x("x", Nb, BlkSize), y("y", Nb, BlkSize);

  const std::size_t incx = 2, incy = 2;
  // Testing incx/incy argument with strided views
  Kokkos::LayoutStride layout_x{Nb, incx, BlkSize, Nb * incx};
  Kokkos::LayoutStride layout_y{Nb, incy, BlkSize, Nb * incy};
  StridedView2DType x_s("x_s", layout_x);
  StridedView2DType y_s("y_s", layout_y);

  auto h_A     = Kokkos::create_mirror_view(A);
  auto h_A_ref = Kokkos::create_mirror_view(A_ref);
  auto h_x     = Kokkos::create_mirror_view(x);
  auto h_y     = Kokkos::create_mirror_view(y);

  for (std::size_t ib = 0; ib < Nb; ib++) {
    h_A(ib, 0, 0) = 1;
    h_A(ib, 0, 1) = -3;
    h_A(ib, 0, 2) = -2;
    h_A(ib, 0, 3) = 0;
    h_A(ib, 1, 0) = -1;
    h_A(ib, 1, 1) = 1;
    h_A(ib, 1, 2) = -3;
    h_A(ib, 1, 3) = -2;
    h_A(ib, 2, 0) = 2;
    h_A(ib, 2, 1) = -1;
    h_A(ib, 2, 2) = 1;
    h_A(ib, 2, 3) = -3;
    h_A(ib, 3, 0) = 0;
    h_A(ib, 3, 1) = 2;
    h_A(ib, 3, 2) = -1;
    h_A(ib, 3, 3) = 1;

    if (std::is_same_v<ArgUplo, KokkosBatched::Uplo::Upper>) {
      h_A_ref(ib, 0, 0) = 13;
      h_A_ref(ib, 0, 1) = 13.5;
      h_A_ref(ib, 0, 2) = 19;
      h_A_ref(ib, 0, 3) = 25.5;
      h_A_ref(ib, 1, 0) = 0;
      h_A_ref(ib, 1, 1) = 19;
      h_A_ref(ib, 1, 2) = 16.5;
      h_A_ref(ib, 1, 3) = 19;
      h_A_ref(ib, 2, 0) = 0;
      h_A_ref(ib, 2, 1) = 0;
      h_A_ref(ib, 2, 2) = 19;
      h_A_ref(ib, 2, 3) = 13.5;
      h_A_ref(ib, 3, 0) = 0;
      h_A_ref(ib, 3, 1) = 0;
      h_A_ref(ib, 3, 2) = 0;
      h_A_ref(ib, 3, 3) = 13;
    } else {
      h_A_ref(ib, 0, 0) = 13;
      h_A_ref(ib, 0, 1) = 0;
      h_A_ref(ib, 0, 2) = 0;
      h_A_ref(ib, 0, 3) = 0;
      h_A_ref(ib, 1, 0) = 15.5;
      h_A_ref(ib, 1, 1) = 19;
      h_A_ref(ib, 1, 2) = 0;
      h_A_ref(ib, 1, 3) = 0;
      h_A_ref(ib, 2, 0) = 23;
      h_A_ref(ib, 2, 1) = 18.5;
      h_A_ref(ib, 2, 2) = 19;
      h_A_ref(ib, 2, 3) = 0;
      h_A_ref(ib, 3, 0) = 25.5;
      h_A_ref(ib, 3, 1) = 23;
      h_A_ref(ib, 3, 2) = 15.5;
      h_A_ref(ib, 3, 3) = 13;
    }

    for (std::size_t j = 0; j < BlkSize; j++) {
      h_x(ib, j) = static_cast<ScalarType>(j + 1);
      h_y(ib, j) = static_cast<ScalarType>(BlkSize - j);
    }
  }

  Kokkos::deep_copy(A, h_A);
  Kokkos::deep_copy(x, h_x);
  Kokkos::deep_copy(y, h_y);

  // Upper or lower diagonal part of A into A_s
  create_triangular_matrix<View3DType, View3DType, ArgUplo, KokkosBatched::Diag::NonUnit>(A, A_s);

  Kokkos::deep_copy(A, A_s);

  // Deep copy to strided views
  Kokkos::deep_copy(x_s, x);
  Kokkos::deep_copy(y_s, y);

  const ScalarType alpha = 1.5;

  auto info = Functor_BatchedSerialSyr2<DeviceType, View2DType, View2DType, View3DType, ScalarType, ParamTagType>(
                  alpha, x, y, A)
                  .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // With strided views
  info =
      Functor_BatchedSerialSyr2<DeviceType, StridedView2DType, StridedView2DType, View3DType, ScalarType, ParamTagType>(
          alpha, x_s, y_s, A_s)
          .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  RealType eps = 1.0e1 * ats::epsilon();
  Kokkos::deep_copy(h_A, A);
  auto h_A_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_s);

  // Check if A:= alpha * x * y**T + alpha * y * x**T + A
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      for (std::size_t j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_A(ib, i, j), h_A_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_A_s(ib, i, j), h_A_ref(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched syr2 test
///
/// \param Nb [in] Batch size of matrices
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_syr2(const std::size_t Nb, const std::size_t BlkSize) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo           = typename ParamTagType::uplo;

  View3DType A("A", Nb, BlkSize, BlkSize), A0("A0", Nb, BlkSize, BlkSize), A_s("A_s", Nb, BlkSize, BlkSize),
      A0_s("A0_s", Nb, BlkSize, BlkSize), A_ref("A_ref", Nb, BlkSize, BlkSize), A0_ref("A0_ref", Nb, BlkSize, BlkSize);
  View2DType x("x", Nb, BlkSize), y("y", Nb, BlkSize);

  const std::size_t incx = 2, incy = 2;
  // Testing incx/incy argument with strided views
  Kokkos::LayoutStride layout_x{Nb, incx, BlkSize, Nb * incx};
  Kokkos::LayoutStride layout_y{Nb, incy, BlkSize, Nb * incy};
  StridedView2DType x_s("x_s", layout_x);
  StridedView2DType y_s("y_s", layout_y);

  // Create a random matrix A, x, and y
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);
  Kokkos::fill_random(y, rand_pool, randStart, randEnd);

  // Upper or lower triangular part of A
  create_triangular_matrix<View3DType, View3DType, ArgUplo, KokkosBatched::Diag::NonUnit>(A, A_ref);

  Kokkos::deep_copy(A, A_ref);

  // Deep copy to strided views
  Kokkos::deep_copy(A_s, A);
  Kokkos::deep_copy(x_s, x);
  Kokkos::deep_copy(y_s, y);

  // When A0 is zero
  const ScalarType alpha = 1.5;
  auto info0 = Functor_BatchedSerialSyr2<DeviceType, View2DType, View2DType, View3DType, ScalarType, ParamTagType>(
                   alpha, x, y, A0)
                   .run();

  // When A is a random matrix
  auto info1 = Functor_BatchedSerialSyr2<DeviceType, View2DType, View2DType, View3DType, ScalarType, ParamTagType>(
                   alpha, x, y, A)
                   .run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  // With strided Views
  info0 =
      Functor_BatchedSerialSyr2<DeviceType, StridedView2DType, StridedView2DType, View3DType, ScalarType, ParamTagType>(
          alpha, x_s, y_s, A0_s)
          .run();

  // When A is a random matrix
  info1 =
      Functor_BatchedSerialSyr2<DeviceType, StridedView2DType, StridedView2DType, View3DType, ScalarType, ParamTagType>(
          alpha, x_s, y_s, A_s)
          .run();

  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  // Make a reference at host
  auto h_x      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
  auto h_y      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  auto h_A_ref  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_ref);
  auto h_A0_ref = Kokkos::create_mirror_view(Kokkos::HostSpace(), A0_ref);

  // Note: ConjTranspose corresponds to {c,z}her2 for Hermitian matrix
  const bool is_conj = std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>;
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t j = 0; j < BlkSize; j++) {
      if (h_x(ib, j) != ScalarType(0) || h_y(ib, j) != ScalarType(0)) {
        auto temp1 = is_conj ? alpha * KokkosKernels::ArithTraits<ScalarType>::conj(h_y(ib, j)) : alpha * h_y(ib, j);
        auto temp2 = is_conj ? KokkosKernels::ArithTraits<ScalarType>::conj(alpha * h_x(ib, j)) : alpha * h_x(ib, j);

        if (std::is_same_v<ArgUplo, KokkosBatched::Uplo::Upper>) {
          for (std::size_t i = 0; i < j + 1; i++) {
            h_A_ref(ib, i, j)  = h_A_ref(ib, i, j) + h_x(ib, i) * temp1 + h_y(ib, i) * temp2;
            h_A0_ref(ib, i, j) = h_x(ib, i) * temp1 + h_y(ib, i) * temp2;
          }
        } else {
          for (std::size_t i = j; i < BlkSize; i++) {
            h_A_ref(ib, i, j)  = h_A_ref(ib, i, j) + h_x(ib, i) * temp1 + h_y(ib, i) * temp2;
            h_A0_ref(ib, i, j) = h_x(ib, i) * temp1 + h_y(ib, i) * temp2;
          }
        }
        h_A_ref(ib, j, j) =
            is_conj ? KokkosKernels::ArithTraits<ScalarType>::real(h_A_ref(ib, j, j)) : h_A_ref(ib, j, j);
        h_A0_ref(ib, j, j) =
            is_conj ? KokkosKernels::ArithTraits<ScalarType>::real(h_A0_ref(ib, j, j)) : h_A0_ref(ib, j, j);
      }
    }
  }

  RealType eps = 1.0e1 * ats::epsilon();

  auto h_A    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  auto h_A0   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A0);
  auto h_A_s  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_s);
  auto h_A0_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A0_s);

  // Check if A:= alpha * x * y**T + alpha * y * x**T + A or A:= alpha * x * y**H + conjg(alpha) * y * x**H + A
  for (std::size_t ib = 0; ib < Nb; ib++) {
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

}  // namespace Syr2
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType>
int test_batched_syr2() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Syr2::impl_test_batched_syr2_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Syr2::impl_test_batched_syr2_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Syr2::impl_test_batched_syr2<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Syr2::impl_test_batched_syr2<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Syr2::impl_test_batched_syr2_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Syr2::impl_test_batched_syr2_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Syr2::impl_test_batched_syr2<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Syr2::impl_test_batched_syr2<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_syr2_l_t_float) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Lower, Trans::Transpose>;
  test_batched_syr2<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_l_c_float) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Lower, Trans::ConjTranspose>;
  test_batched_syr2<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_u_t_float) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Upper, Trans::Transpose>;
  test_batched_syr2<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_u_c_float) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Upper, Trans::ConjTranspose>;
  test_batched_syr2<TestDevice, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_syr2_l_t_double) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Lower, Trans::Transpose>;
  test_batched_syr2<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_l_c_double) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Lower, Trans::ConjTranspose>;
  test_batched_syr2<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_u_t_double) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Upper, Trans::Transpose>;
  test_batched_syr2<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_u_c_double) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Upper, Trans::ConjTranspose>;
  test_batched_syr2<TestDevice, double, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_syr2_l_t_fcomplex) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Lower, Trans::Transpose>;
  test_batched_syr2<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_l_c_fcomplex) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Lower, Trans::ConjTranspose>;
  test_batched_syr2<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_u_t_fcomplex) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Upper, Trans::Transpose>;
  test_batched_syr2<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_u_c_fcomplex) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Upper, Trans::ConjTranspose>;
  test_batched_syr2<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_syr2_l_t_dcomplex) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Lower, Trans::Transpose>;
  test_batched_syr2<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_l_c_dcomplex) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Lower, Trans::ConjTranspose>;
  test_batched_syr2<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_u_t_dcomplex) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Upper, Trans::Transpose>;
  test_batched_syr2<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr2_u_c_dcomplex) {
  using param_tag_type = ::Test::Syr2::ParamTag<Uplo::Upper, Trans::ConjTranspose>;
  test_batched_syr2<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
#endif
