// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Syr.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Syr {

template <typename U, typename T>
struct ParamTag {
  using uplo  = U;
  using trans = T;
};

template <typename DeviceType, typename XViewType, typename AViewType, typename ScalarType, typename ParamTagType>
struct Functor_BatchedSerialSyr {
  using execution_space = typename DeviceType::execution_space;
  XViewType m_x;
  AViewType m_A;
  ScalarType m_alpha;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialSyr(const ScalarType alpha, const XViewType &x, const AViewType &A)
      : m_x(x), m_A(A), m_alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    auto sub_x = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto sub_A = Kokkos::subview(m_A, k, Kokkos::ALL(), Kokkos::ALL());

    info += KokkosBatched::SerialSyr<typename ParamTagType::uplo, typename ParamTagType::trans>::invoke(m_alpha, sub_x,
                                                                                                        sub_A);
  }

  inline int run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialSyr");
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

/// \brief Implementation details of batched syr analytical test
///        to confirm A:= x*x**T + A is computed correctly
/// \param Nb [in] Batch size
///        alpha = 1.5
///        4x4 matrix (upper)
///        U: [[1, -3, -2,  0],
///            [0,  1, -3, -2],
///            [0,  0,  1, -3],
///            [0,  0,  0,  1]]
///        x: [1, 2, 3, 4]
///        Ref: [[ 2.5,  0.,   2.5,  6., ],
///              [ 0.,   7.,   6.,  10., ],
///              [ 0.,   0.,  14.5, 15., ],
///              [ 0.,   0.,   0.,  25., ]]
///
///        4x4 matrix (lower)
///        L: [[1,  0,  0,  0],
///            [-1, 1,  0,  0],
///            [2, -1,  1,  0],
///            [0,  2, -1,  1]]
///        x: [1, 2, 3, 4]
///        Ref: [[ 2.5,  0.,   0.,   0., ],
///              [ 2.,   7.,   0.,   0., ],
///              [ 6.5,  8.,  14.5,  0., ],
///              [ 6.,  14.,  17.,  25., ]]
///
/// \param Nb [in] Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_syr_analytical(const std::size_t Nb) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo           = typename ParamTagType::uplo;

  const std::size_t BlkSize = 4;
  View3DType A("A", Nb, BlkSize, BlkSize), A_s("A_s", Nb, BlkSize, BlkSize), A_ref("A_ref", Nb, BlkSize, BlkSize);
  View2DType x("x", Nb, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, BlkSize, Nb * incx};
  StridedView2DType x_s("x_s", layout);

  // Only filling x2, A2 and deep_copy from its subview
  auto h_A     = Kokkos::create_mirror_view(A);
  auto h_A_ref = Kokkos::create_mirror_view(A_ref);
  auto h_x     = Kokkos::create_mirror_view(x);

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
      h_A_ref(ib, 0, 0) = 2.5;
      h_A_ref(ib, 0, 1) = 0;
      h_A_ref(ib, 0, 2) = 2.5;
      h_A_ref(ib, 0, 3) = 6;
      h_A_ref(ib, 1, 0) = 0;
      h_A_ref(ib, 1, 1) = 7;
      h_A_ref(ib, 1, 2) = 6;
      h_A_ref(ib, 1, 3) = 10;
      h_A_ref(ib, 2, 0) = 0;
      h_A_ref(ib, 2, 1) = 0;
      h_A_ref(ib, 2, 2) = 14.5;
      h_A_ref(ib, 2, 3) = 15.;
      h_A_ref(ib, 3, 0) = 0;
      h_A_ref(ib, 3, 1) = 0;
      h_A_ref(ib, 3, 2) = 0;
      h_A_ref(ib, 3, 3) = 25.;
    } else {
      h_A_ref(ib, 0, 0) = 2.5;
      h_A_ref(ib, 0, 1) = 0;
      h_A_ref(ib, 0, 2) = 0;
      h_A_ref(ib, 0, 3) = 0;
      h_A_ref(ib, 1, 0) = 2;
      h_A_ref(ib, 1, 1) = 7;
      h_A_ref(ib, 1, 2) = 0;
      h_A_ref(ib, 1, 3) = 0;
      h_A_ref(ib, 2, 0) = 6.5;
      h_A_ref(ib, 2, 1) = 8;
      h_A_ref(ib, 2, 2) = 14.5;
      h_A_ref(ib, 2, 3) = 0;
      h_A_ref(ib, 3, 0) = 6;
      h_A_ref(ib, 3, 1) = 14;
      h_A_ref(ib, 3, 2) = 17;
      h_A_ref(ib, 3, 3) = 25;
    }

    for (std::size_t j = 0; j < BlkSize; j++) {
      h_x(ib, j) = static_cast<ScalarType>(j + 1);
    }
  }

  Kokkos::deep_copy(A, h_A);
  Kokkos::deep_copy(x, h_x);

  // Upper or lower diagnoal part of A into A_s
  create_triangular_matrix<View3DType, View3DType, ArgUplo, KokkosBatched::Diag::NonUnit>(A, A_s);

  Kokkos::deep_copy(A, A_s);

  // Deep copy to strided views
  Kokkos::deep_copy(x_s, x);

  const ScalarType alpha = 1.5;

  auto info = Functor_BatchedSerialSyr<DeviceType, View2DType, View3DType, ScalarType, ParamTagType>(alpha, x, A).run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // With strided views
  info = Functor_BatchedSerialSyr<DeviceType, StridedView2DType, View3DType, ScalarType, ParamTagType>(alpha, x_s, A_s)
             .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  RealType eps = 1.0e1 * ats::epsilon();
  Kokkos::deep_copy(h_A, A);
  auto h_A_s = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_s);

  // Check if A:= alpha * x * y**T + A
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      for (std::size_t j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_A(ib, i, j), h_A_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_A_s(ib, i, j), h_A_ref(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched syr test
///
/// \param Nb [in] Batch size of matrices
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType>
void impl_test_batched_syr(const std::size_t Nb, const std::size_t BlkSize) {
  using ats               = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType          = typename ats::mag_type;
  using View2DType        = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ScalarType **, Kokkos::LayoutStride, DeviceType>;
  using View3DType        = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using ArgUplo           = typename ParamTagType::uplo;

  View3DType A("A", Nb, BlkSize, BlkSize), A0("A0", Nb, BlkSize, BlkSize), A_s("A_s", Nb, BlkSize, BlkSize),
      A0_s("A0_s", Nb, BlkSize, BlkSize), A_ref("A_ref", Nb, BlkSize, BlkSize), A0_ref("A0_ref", Nb, BlkSize, BlkSize);
  View2DType x("x", Nb, BlkSize);

  const std::size_t incx = 2;
  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{Nb, incx, BlkSize, Nb * incx};
  StridedView2DType x_s("x_s", layout);

  // Create a random matrix A and x
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);

  // Upper or lower triangular part of A
  create_triangular_matrix<View3DType, View3DType, ArgUplo, KokkosBatched::Diag::NonUnit>(A, A_ref);

  Kokkos::deep_copy(A, A_ref);

  // Deep copy to strided views
  Kokkos::deep_copy(A_s, A);
  Kokkos::deep_copy(x_s, x);

  // When A0 is zero
  const ScalarType alpha = 1.5;
  auto info0 =
      Functor_BatchedSerialSyr<DeviceType, View2DType, View3DType, ScalarType, ParamTagType>(alpha, x, A0).run();

  // When A is a random matrix
  auto info1 =
      Functor_BatchedSerialSyr<DeviceType, View2DType, View3DType, ScalarType, ParamTagType>(alpha, x, A).run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  // With strided Views
  info0 =
      Functor_BatchedSerialSyr<DeviceType, StridedView2DType, View3DType, ScalarType, ParamTagType>(alpha, x_s, A0_s)
          .run();

  // When A is a random matrix
  info1 = Functor_BatchedSerialSyr<DeviceType, StridedView2DType, View3DType, ScalarType, ParamTagType>(alpha, x_s, A_s)
              .run();

  // Make a reference at host
  auto h_x      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
  auto h_A_ref  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_ref);
  auto h_A0_ref = Kokkos::create_mirror_view(Kokkos::HostSpace(), A0_ref);

  // Note: ConjTranspose corresponds to {c,z}her for Hermitian matrix
  const bool is_conj = std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>;
  for (std::size_t ib = 0; ib < Nb; ib++) {
    for (std::size_t j = 0; j < BlkSize; j++) {
      if (h_x(ib, j) != 0) {
        auto temp = is_conj ? alpha * KokkosKernels::ArithTraits<ScalarType>::conj(h_x(ib, j)) : alpha * h_x(ib, j);

        if (std::is_same_v<ArgUplo, KokkosBatched::Uplo::Upper>) {
          for (std::size_t i = 0; i < j + 1; i++) {
            h_A_ref(ib, i, j)  = h_A_ref(ib, i, j) + h_x(ib, i) * temp;
            h_A0_ref(ib, i, j) = h_x(ib, i) * temp;
          }
        } else {
          for (std::size_t i = j; i < BlkSize; i++) {
            h_A_ref(ib, i, j)  = h_A_ref(ib, i, j) + h_x(ib, i) * temp;
            h_A0_ref(ib, i, j) = h_x(ib, i) * temp;
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

  // Check if A:= alpha * x * y**T + A or A:= alpha * x * y**H + A
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

}  // namespace Syr
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType>
int test_batched_syr() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Syr::impl_test_batched_syr_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Syr::impl_test_batched_syr_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Syr::impl_test_batched_syr<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Syr::impl_test_batched_syr<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Syr::impl_test_batched_syr_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(1);
    Test::Syr::impl_test_batched_syr_analytical<DeviceType, ScalarType, LayoutType, ParamTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Syr::impl_test_batched_syr<DeviceType, ScalarType, LayoutType, ParamTagType>(1, i);
      Test::Syr::impl_test_batched_syr<DeviceType, ScalarType, LayoutType, ParamTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_syr_l_t_float) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Lower, Trans::Transpose>;
  test_batched_syr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_l_c_float) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Lower, Trans::ConjTranspose>;
  test_batched_syr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_u_t_float) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Upper, Trans::Transpose>;
  test_batched_syr<TestDevice, float, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_u_c_float) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Upper, Trans::ConjTranspose>;
  test_batched_syr<TestDevice, float, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_syr_l_t_double) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Lower, Trans::Transpose>;
  test_batched_syr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_l_c_double) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Lower, Trans::ConjTranspose>;
  test_batched_syr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_u_t_double) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Upper, Trans::Transpose>;
  test_batched_syr<TestDevice, double, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_u_c_double) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Upper, Trans::ConjTranspose>;
  test_batched_syr<TestDevice, double, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_syr_l_t_fcomplex) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Lower, Trans::Transpose>;
  test_batched_syr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_l_c_fcomplex) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Lower, Trans::ConjTranspose>;
  test_batched_syr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_u_t_fcomplex) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Upper, Trans::Transpose>;
  test_batched_syr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_u_c_fcomplex) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Upper, Trans::ConjTranspose>;
  test_batched_syr<TestDevice, Kokkos::complex<float>, param_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_syr_l_t_dcomplex) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Lower, Trans::Transpose>;
  test_batched_syr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_l_c_dcomplex) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Lower, Trans::ConjTranspose>;
  test_batched_syr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_u_t_dcomplex) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Upper, Trans::Transpose>;
  test_batched_syr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
TEST_F(TestCategory, test_batched_syr_u_c_dcomplex) {
  using param_tag_type = ::Test::Syr::ParamTag<Uplo::Upper, Trans::ConjTranspose>;
  test_batched_syr<TestDevice, Kokkos::complex<double>, param_tag_type>();
}
#endif
