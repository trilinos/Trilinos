// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "Kokkos_DynRankView.hpp"

#include "KokkosBatched_Util.hpp"
#include "KokkosBlas2_gemv.hpp"
#include "KokkosBatched_Trsv_Decl.hpp"
#include "Test_Batched_DenseUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace Trsv {

template <typename U, typename T, typename D>
struct ParamTag {
  using uplo  = U;
  using trans = T;
  using diag  = D;
};

template <typename DeviceType, typename AViewType, typename BViewType, typename ScalarType, typename ParamTagType,
          typename AlgoTagType>
struct Functor_BatchedSerialTrsv {
  using execution_space = typename DeviceType::execution_space;
  AViewType m_a;
  BViewType m_b;

  ScalarType m_alpha;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialTrsv(const ScalarType alpha, const AViewType &a, const BViewType &b)
      : m_a(a), m_b(b), m_alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k, int &info) const {
    auto aa = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(m_b, k, Kokkos::ALL());

    info += SerialTrsv<typename ParamTagType::uplo, typename ParamTagType::trans, typename ParamTagType::diag,
                       AlgoTagType>::invoke(m_alpha, aa, bb);
  }

  inline int run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialTrsv");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, m_b.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename ScalarType, typename AViewType, typename xViewType, typename yViewType,
          typename ParamTagType>
struct Functor_BatchedSerialGemv {
  using execution_space = typename DeviceType::execution_space;
  AViewType m_a;
  xViewType m_x;
  yViewType m_y;
  ScalarType m_alpha, m_beta;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGemv(const ScalarType alpha, const AViewType &a, const xViewType &x, const ScalarType beta,
                            const yViewType &y)
      : m_a(a), m_x(x), m_y(y), m_alpha(alpha), m_beta(beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k) const {
    auto aa = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto xx = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto yy = Kokkos::subview(m_y, k, Kokkos::ALL());

    KokkosBlas::SerialGemv<typename ParamTagType::trans, Algo::Gemv::Unblocked>::invoke(m_alpha, aa, xx, m_beta, yy);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGemv");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, m_x.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename ScalarType, typename ValueType, typename LayoutType, typename ParamTagType,
          typename AlgoTagType>
void impl_test_batched_trsv_blocking(const int N, const int BlkSize) {
  using ats        = KokkosKernels::ArithTraits<ValueType>;
  using View2DType = Kokkos::View<ValueType **, LayoutType, DeviceType>;
  using View3DType = Kokkos::View<ValueType ***, LayoutType, DeviceType>;

  /// randomized input testing views
  ScalarType alpha(1.5);

  View3DType a0("a0", N, BlkSize, BlkSize), a1("a1", N, BlkSize, BlkSize);
  View2DType b0("b0", N, BlkSize), b1("b1", N, BlkSize);

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  Kokkos::fill_random(a0, random, ValueType(1.0));
  Kokkos::fill_random(b0, random, ValueType(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(b0, 1.0);

  Kokkos::deep_copy(a1, a0);
  Kokkos::deep_copy(b1, b0);

  auto info0 =
      Functor_BatchedSerialTrsv<DeviceType, View3DType, View2DType, ScalarType, ParamTagType, Algo::Trsv::Blocked>(
          alpha, a0, b0)
          .run();
  auto info1 = Functor_BatchedSerialTrsv<DeviceType, View3DType, View2DType, ScalarType, ParamTagType, AlgoTagType>(
                   alpha, a1, b1)
                   .run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);

  /// for comparison send it to host
  auto b0_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b0);
  auto b1_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b1);

  /// this eps is about 10^-14
  using mag_type = typename ats::mag_type;
  mag_type sum(1), diff(0);
  const mag_type eps = 1.0e3 * ats::epsilon();

  /// check b0 = b1 ;
  for (int k = 0; k < N; ++k)
    for (int i = 0; i < BlkSize; ++i) {
      sum += ats::abs(b0_host(k, i));
      diff += ats::abs(b0_host(k, i) - b1_host(k, i));
    }
  EXPECT_NEAR(diff / sum, 0.0, eps);
}

/// \brief Implementation details of batched trsv analytical test
/// Confirm A * x = b, where
///        A: [[1, 1, 1],
///            [1, 2, 2],
///            [1, 2, 3]]
///        b: [1, 1, 1]
///
/// Upper and Non-transpose
/// x0 + x1 +  x2 = 1
///    +2x1 + 2x2 = 1
///           3x2 = 1
/// x = [1/2, 1/6, 1/3]
///
/// Upper and Transpose
/// x0            = 1
/// x0 +2x1       = 1
/// x0 +2x1 + 3x2 = 1
/// x = [1, 0, 0]
///
/// Lower, Non-transpose
/// x0            = 1
/// x0 +2x1       = 1
/// x0 +2x1 + 3x2 = 1
/// x = [1, 0, 0]
///
/// Lower, Transpose
/// x0 + x1 +  x2 = 1
///    +2x1 + 2x2 = 1
///           3x2 = 1
/// x = [1/2, 1/6, 1/3]
/// \param N [in] Batch size of matrices and RHS
template <typename DeviceType, typename ScalarType, typename ValueType, typename LayoutType, typename ParamTagType,
          typename AlgoTagType>
void impl_test_batched_trsv_analytical(const std::size_t N) {
  using ats      = typename KokkosKernels::ArithTraits<ValueType>;
  using RealType = typename ats::mag_type;

  using View2DType        = Kokkos::View<ValueType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ValueType **, Kokkos::LayoutStride, DeviceType>;
  using View3DType        = Kokkos::View<ValueType ***, LayoutType, DeviceType>;

  constexpr std::size_t BlkSize = 3, incx = 2;
  View3DType A("A", N, BlkSize, BlkSize);
  View2DType x("x", N, BlkSize), x_ref("x_ref", N, BlkSize);  // Solutions

  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{N, incx, BlkSize, N * incx};
  StridedView2DType x_s("x_s", layout);  // Solutions

  auto h_A     = Kokkos::create_mirror_view(A);
  auto h_x_ref = Kokkos::create_mirror_view(x_ref);
  for (std::size_t ib = 0; ib < N; ib++) {
    h_A(ib, 0, 0) = 1.0;
    h_A(ib, 0, 1) = 1.0;
    h_A(ib, 0, 2) = 1.0;
    h_A(ib, 1, 0) = 1.0;
    h_A(ib, 1, 1) = 2.0;
    h_A(ib, 1, 2) = 2.0;
    h_A(ib, 2, 0) = 1.0;
    h_A(ib, 2, 1) = 2.0;
    h_A(ib, 2, 2) = 3.0;

    if (std::is_same_v<typename ParamTagType::uplo, KokkosBatched::Uplo::Upper>) {
      if (std::is_same_v<typename ParamTagType::trans, Trans::NoTranspose>) {
        if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
          h_x_ref(ib, 0) = 1.0 / 2.0;
          h_x_ref(ib, 1) = 1.0 / 6.0;
          h_x_ref(ib, 2) = 1.0 / 3.0;
        } else {
          h_x_ref(ib, 0) = 1.0;
          h_x_ref(ib, 1) = -1.0;
          h_x_ref(ib, 2) = 1.0;
        }
      } else {
        // Diag::NonUnit does not matter
        h_x_ref(ib, 0) = 1.0;
        h_x_ref(ib, 1) = 0.0;
        h_x_ref(ib, 2) = 0.0;
      }
    } else {
      if (std::is_same_v<typename ParamTagType::trans, Trans::NoTranspose>) {
        // Diag::NonUnit does not matter
        h_x_ref(ib, 0) = 1.0;
        h_x_ref(ib, 1) = 0.0;
        h_x_ref(ib, 2) = 0.0;
      } else {
        if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
          h_x_ref(ib, 0) = 1.0 / 2.0;
          h_x_ref(ib, 1) = 1.0 / 6.0;
          h_x_ref(ib, 2) = 1.0 / 3.0;
        } else {
          h_x_ref(ib, 0) = 1.0;
          h_x_ref(ib, 1) = -1.0;
          h_x_ref(ib, 2) = 1.0;
        }
      }
    }
  }

  Kokkos::deep_copy(A, h_A);

  // Set RHS as [1.0, 1.0, 1.0]
  Kokkos::deep_copy(x, 1.0);
  Kokkos::deep_copy(x_s, x);

  // trsv to solve U * x = b or L * x = b
  auto info =
      Functor_BatchedSerialTrsv<DeviceType, View3DType, View2DType, ScalarType, ParamTagType, Algo::Trsv::Unblocked>(
          1.0, A, x)
          .run();
  auto info_s = Functor_BatchedSerialTrsv<DeviceType, View3DType, StridedView2DType, ScalarType, ParamTagType,
                                          Algo::Trsv::Unblocked>(1.0, A, x_s)
                    .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);
  EXPECT_EQ(info_s, 0);

  // Check x = x_ref
  RealType eps = 1.0e1 * ats::epsilon();
  auto h_x     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }

  // Testing for strided views, reusing x
  Kokkos::deep_copy(x, x_s);
  Kokkos::deep_copy(h_x, x);
  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched trsv test
/// \param N [in] Batch size of matrices and RHS
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename ValueType, typename LayoutType, typename ParamTagType,
          typename AlgoTagType>
void impl_test_batched_trsv(const std::size_t N, const std::size_t BlkSize) {
  using ats      = typename KokkosKernels::ArithTraits<ValueType>;
  using RealType = typename ats::mag_type;

  using View2DType        = Kokkos::View<ValueType **, LayoutType, DeviceType>;
  using StridedView2DType = Kokkos::View<ValueType **, Kokkos::LayoutStride, DeviceType>;
  using View3DType        = Kokkos::View<ValueType ***, LayoutType, DeviceType>;
  using DynViewType       = Kokkos::DynRankView<ValueType, LayoutType, DeviceType>;

  constexpr std::size_t incx = 2;
  View3DType A("A", N, BlkSize, BlkSize), AT("AT", N, BlkSize, BlkSize);
  View2DType x("x", N, BlkSize), y("y", N, BlkSize), x_ref("x_ref", N, BlkSize);  // Solutions

  // Testing incx argument with strided views
  Kokkos::LayoutStride layout{N, incx, BlkSize, N * incx};
  StridedView2DType x_s("x_s", layout), y_s("y_s", layout);  // Solutions

  // Testing DynViewType
  DynViewType A_d("A_d", N, BlkSize, BlkSize);
  DynViewType x_d("x_d", N, BlkSize), y_d("y_d", N, BlkSize);

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);
  Kokkos::deep_copy(x_ref, x);  // Keep reference solution
  Kokkos::deep_copy(x_s, x);
  Kokkos::deep_copy(x_d, x);
  Kokkos::deep_copy(A_d, A);

  // Create triangluar matrix
  create_triangular_matrix<View3DType, View3DType, typename ParamTagType::uplo, typename ParamTagType::diag>(A, AT);

  // trsv to solve U * x = b or L * x = b
  auto info =
      Functor_BatchedSerialTrsv<DeviceType, View3DType, View2DType, ScalarType, ParamTagType, Algo::Trsv::Unblocked>(
          1.0, A, x)
          .run();

  auto info_s = Functor_BatchedSerialTrsv<DeviceType, View3DType, StridedView2DType, ScalarType, ParamTagType,
                                          Algo::Trsv::Unblocked>(1.0, A, x_s)
                    .run();

  auto info_d =
      Functor_BatchedSerialTrsv<DeviceType, DynViewType, DynViewType, ScalarType, ParamTagType, Algo::Trsv::Unblocked>(
          1.0, A_d, x_d)
          .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);
  EXPECT_EQ(info_s, 0);
  EXPECT_EQ(info_d, 0);

  // Compute A * x by gemv
  // Gemv to compute A*x, this should be identical to b
  Functor_BatchedSerialGemv<DeviceType, ScalarType, View3DType, View2DType, View2DType, ParamTagType>(1.0, AT, x, 0.0,
                                                                                                      y)
      .run();

  // Gemv to compute A*x, this should be identical to b
  Functor_BatchedSerialGemv<DeviceType, ScalarType, View3DType, StridedView2DType, StridedView2DType, ParamTagType>(
      1.0, AT, x_s, 0.0, y_s)
      .run();

  // Gemv to compute A*x, this should be identical to b
  Functor_BatchedSerialGemv<DeviceType, ScalarType, View3DType, DynViewType, DynViewType, ParamTagType>(1.0, AT, x_d,
                                                                                                        0.0, y_d)
      .run();

  // Check A*x = x_ref
  RealType eps = 1.0e3 * ats::epsilon();
  auto h_y     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  auto h_x_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x_ref);
  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_y(ib, i), h_x_ref(ib, i), eps);
    }
  }

  // Testing for strided views, reusing y
  Kokkos::deep_copy(y, y_s);
  Kokkos::deep_copy(h_y, y);
  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_y(ib, i), h_x_ref(ib, i), eps);
    }
  }

  // Testing for dynamic views, reusing y
  auto h_y_d = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y_d);
  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_y_d(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

}  // namespace Trsv
}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_trsv() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Trsv::impl_test_batched_trsv_analytical<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(1);
    Test::Trsv::impl_test_batched_trsv_analytical<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(2);

    // FIXME: ConjTranspose with blocking is not implemented yet
    if constexpr (!std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>) {
      Test::Trsv::impl_test_batched_trsv_blocking<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(0, 10);
    }

    for (int i = 0; i < 10; ++i) {
      // FIXME: ConjTranspose with blocking is not implemented yet
      if constexpr (!std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>) {
        Test::Trsv::impl_test_batched_trsv_blocking<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                    AlgoTagType>(1, i);
      }
      Test::Trsv::impl_test_batched_trsv<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType, AlgoTagType>(1,
                                                                                                                   i);
      Test::Trsv::impl_test_batched_trsv<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType, AlgoTagType>(2,
                                                                                                                   i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Trsv::impl_test_batched_trsv_analytical<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(1);
    Test::Trsv::impl_test_batched_trsv_analytical<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(2);

    // FIXME: ConjTranspose with blocking is not implemented yet
    if constexpr (!std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>) {
      Test::Trsv::impl_test_batched_trsv_blocking<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(0, 10);
    }

    for (int i = 0; i < 10; ++i) {
      // FIXME: ConjTranspose with blocking is not implemented yet
      if constexpr (!std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>) {
        Test::Trsv::impl_test_batched_trsv_blocking<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                    AlgoTagType>(1, i);
      }
      Test::Trsv::impl_test_batched_trsv<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType, AlgoTagType>(1,
                                                                                                                   i);
      Test::Trsv::impl_test_batched_trsv<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType, AlgoTagType>(2,
                                                                                                                   i);
    }
  }
#endif

  return 0;
}
