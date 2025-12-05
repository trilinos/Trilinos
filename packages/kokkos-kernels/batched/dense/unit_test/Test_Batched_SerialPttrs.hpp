// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Pttrf.hpp"
#include "KokkosBatched_Pttrs.hpp"
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Pttrs {

template <typename U>
struct ParamTag {
  using uplo = U;
};

template <typename DeviceType, typename DViewType, typename EViewType, typename AlgoTagType>
struct Functor_BatchedSerialPttrf {
  using execution_space = typename DeviceType::execution_space;
  DViewType m_d;
  EViewType m_e;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialPttrf(const DViewType &d, const EViewType &e) : m_d(d), m_e(e) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto sub_d = Kokkos::subview(m_d, k, Kokkos::ALL());
    auto sub_e = Kokkos::subview(m_e, k, Kokkos::ALL());

    KokkosBatched::SerialPttrf<AlgoTagType>::invoke(sub_d, sub_e);
  }

  inline void run() {
    using value_type = typename EViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPttrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_d.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename DViewType, typename EViewType, typename BViewType, typename ParamTagType,
          typename AlgoTagType>
struct Functor_BatchedSerialPttrs {
  using execution_space = typename DeviceType::execution_space;
  DViewType m_d;
  EViewType m_e;
  BViewType m_b;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialPttrs(const DViewType &d, const EViewType &e, const BViewType &b) : m_d(d), m_e(e), m_b(b) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k, int &info) const {
    auto sub_d = Kokkos::subview(m_d, k, Kokkos::ALL());
    auto sub_e = Kokkos::subview(m_e, k, Kokkos::ALL());
    auto sub_b = Kokkos::subview(m_b, k, Kokkos::ALL());

    info += KokkosBatched::SerialPttrs<typename ParamTagType::uplo, AlgoTagType>::invoke(sub_d, sub_e, sub_b);
  }

  inline int run() {
    using value_type = typename BViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPttrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, m_d.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
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
  void operator()(const int k) const {
    auto sub_a = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_x = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto sub_y = Kokkos::subview(m_y, k, Kokkos::ALL());

    KokkosBlas::SerialGemv<Trans::NoTranspose, Algo::Gemv::Unblocked>::invoke(m_alpha, sub_a, sub_x, m_beta, sub_y);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPttrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename ScalarType, typename AViewType, typename BViewType, typename CViewType,
          typename ArgTransB>
struct Functor_BatchedSerialGemm {
  using execution_space = typename DeviceType::execution_space;
  AViewType m_a;
  BViewType m_b;
  CViewType m_c;
  ScalarType m_alpha, m_beta;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGemm(const ScalarType alpha, const AViewType &a, const BViewType &b, const ScalarType beta,
                            const CViewType &c)
      : m_a(a), m_b(b), m_c(c), m_alpha(alpha), m_beta(beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto sub_a = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_b = Kokkos::subview(m_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_c = Kokkos::subview(m_c, k, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialGemm<Trans::NoTranspose, ArgTransB, Algo::Gemm::Unblocked>::invoke(m_alpha, sub_a, sub_b,
                                                                                            m_beta, sub_c);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPttrf");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

/// \brief Implementation details of batched pttrs test
/// Confirm A * x = b, where
///   A: [[4, 1, 0],
///       [1, 4, 1],
///       [0, 1, 4]]
///   b: [1, 1, 1]
///   x: [3/14, 1/7, 3/14]
///
/// This corresponds to the following system of equations:
///        4 x0 +   x1        = 1
///          x0 + 4 x1 +   x2 = 1
///                 x1 + 4 x2 = 1
///
/// We confirm this with the factorized tridiagonal matrix L
/// For packed storage,
///   D: [4, 4-1*(1/4), 4-1*(1/(4-1*(1/4)))]
///   L: [1/4, 1/(4-1*(1/4))]
///
/// \param N [in] Batch size of matrix A and RHS
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_pttrs_analytical(const int N) {
  using ats            = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;

  constexpr int BlkSize = 3;
  RealView2DType d(Kokkos::view_alloc("d", Kokkos::WithoutInitializing), N,
                   BlkSize);  // Diagonal components
  View2DType e(Kokkos::view_alloc("e", Kokkos::WithoutInitializing), N,
               BlkSize - 1);                                  // Upper and lower diagonal components (identical)
  View2DType x("x", N, BlkSize), x_ref("x_ref", N, BlkSize);  // Solutions

  auto h_d     = Kokkos::create_mirror_view(d);
  auto h_e     = Kokkos::create_mirror_view(e);
  auto h_x_ref = Kokkos::create_mirror_view(x_ref);

  RealType d0   = 4.0;
  ScalarType e0 = ScalarType(1.0);

  for (int ib = 0; ib < N; ib++) {
    h_d(ib, 0) = d0;

    h_e(ib, 0) = e0 / h_d(ib, 0);
    h_d(ib, 1) = d0 - ats::real(h_e(ib, 0)) * ats::real(e0);

    h_e(ib, 1) = e0 / h_d(ib, 1);
    h_d(ib, 2) = d0 - ats::real(h_e(ib, 1)) * ats::real(e0);

    h_x_ref(ib, 0) = ScalarType(3.0 / 14.0);
    h_x_ref(ib, 1) = ScalarType(1.0 / 7.0);
    h_x_ref(ib, 2) = ScalarType(3.0 / 14.0);
  }

  Kokkos::fence();

  Kokkos::deep_copy(x, ScalarType(1.0));
  Kokkos::deep_copy(d, h_d);
  Kokkos::deep_copy(e, h_e);

  // pttrs (Note, d and e must be factorized by pttrf)
  auto info =
      Functor_BatchedSerialPttrs<DeviceType, RealView2DType, View2DType, View2DType, ParamTagType, AlgoTagType>(d, e, x)
          .run();

  Kokkos::fence();

  EXPECT_EQ(info, 0);

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  auto h_x = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);

  // Check x = [3/14, 1/7, 3/14]
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched pttrs test
///
/// \param N [in] Batch size of matrix A and RHS
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_pttrs(const int N, const int BlkSize) {
  using ats            = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType     = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  View3DType A("A", N, BlkSize, BlkSize), EL("EL", N, BlkSize, BlkSize), EU("EU", N, BlkSize, BlkSize),
      D("D", N, BlkSize, BlkSize), I("I", N, BlkSize, BlkSize);
  RealView2DType d(Kokkos::view_alloc("d", Kokkos::WithoutInitializing), N, BlkSize),  // Diagonal components
      ones(Kokkos::view_alloc("ones", Kokkos::WithoutInitializing), N, BlkSize);
  View2DType e_upper("e_upper", N, BlkSize - 1),
      e_lower("e_lower", N, BlkSize - 1);                                         // upper and lower diagonal components
  View2DType x("x", N, BlkSize), b_ref("x_ref", N, BlkSize), b("b", N, BlkSize);  // Solutions

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  RealType realRandStart, realRandEnd;
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, realRandStart, realRandEnd);
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);

  // Add BlkSize to ensure positive definiteness
  Kokkos::fill_random(d, rand_pool, realRandStart + BlkSize, realRandEnd + BlkSize);
  Kokkos::fill_random(e_upper, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);

  auto h_e_upper = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), e_upper);
  auto h_e_lower = Kokkos::create_mirror_view(e_lower);

  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize - 1; i++) {
      // Fill the lower diagonal with conjugate of the upper diagonal
      h_e_lower(ib, i) = KokkosKernels::ArithTraits<ScalarType>::conj(h_e_upper(ib, i));
    }
  }

  Kokkos::deep_copy(e_lower, h_e_lower);
  Kokkos::deep_copy(b_ref, x);
  Kokkos::deep_copy(ones, RealType(1.0));

  // Reconstruct Tridiagonal matrix A
  // A = D + EL + EU
  create_diagonal_matrix(e_lower, EL, -1);
  create_diagonal_matrix(e_upper, EU, 1);
  create_diagonal_matrix(d, D);
  create_diagonal_matrix(ones, I);

  // Matrix matrix addition by Gemm
  // D + EU by D * I + EU (result stored in EU)
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, D, I,
                                                                                                            1.0, EU)
      .run();

  // Copy EL to A
  Kokkos::deep_copy(A, EL);

  // EU + EL by EU * I + A (result stored in A)
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, EU, I,
                                                                                                            1.0, A)
      .run();

  Kokkos::fence();

  int info = 0;
  if (std::is_same_v<typename ParamTagType::uplo, KokkosBatched::Uplo::Upper>) {
    // Factorize matrix A -> U**H * D * U
    // d and e are updated by pttrf
    Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e_upper).run();

    // pttrs (Note, d and e must be factorized by pttrf)
    info = Functor_BatchedSerialPttrs<DeviceType, RealView2DType, View2DType, View2DType, ParamTagType, AlgoTagType>(
               d, e_upper, x)
               .run();
  } else {
    // Factorize matrix A -> L * D * L**H
    // d and e are updated by pttrf
    Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e_lower).run();

    // pttrs (Note, d and e must be factorized by pttrf)
    info = Functor_BatchedSerialPttrs<DeviceType, RealView2DType, View2DType, View2DType, ParamTagType, AlgoTagType>(
               d, e_lower, x)
               .run();
  }

  Kokkos::fence();

  EXPECT_EQ(info, 0);

  // Gemv to compute b = A * x, this should be identical to b_ref
  Functor_BatchedSerialGemv<DeviceType, ScalarType, View3DType, View2DType, View2DType>(1.0, A, x, 0.0, b).run();

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  auto h_b     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b);
  auto h_b_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b_ref);

  // Check A * x = b_ref
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_b(ib, i), h_b_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched pttrs test for early return
///        BlkSize must be 0 or 1
///
/// \param N [in] Batch size of matrix A and RHS
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_pttrs_quick_return(const int N, const int BlkSize) {
  using ats            = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;

  if (BlkSize > 1) return;

  const int BlkSize_minus_1 = BlkSize > 0 ? BlkSize - 1 : 0;

  RealView2DType d(Kokkos::view_alloc("d", Kokkos::WithoutInitializing), N,
                   BlkSize);  // Diagonal components
  View2DType e(Kokkos::view_alloc("e", Kokkos::WithoutInitializing), N,
               BlkSize_minus_1);  // lower diagonal components
  View2DType x("x", N, BlkSize);  // Solutions

  const RealType reference_value = 4.0;

  Kokkos::deep_copy(d, reference_value);
  Kokkos::deep_copy(e, ScalarType(1.0));
  Kokkos::deep_copy(x, ScalarType(1.0));

  // pttrs (Note, d and e must be factorized by pttrf)
  // For 1x1 case, pttrf does not change d and e
  auto info =
      Functor_BatchedSerialPttrs<DeviceType, RealView2DType, View2DType, View2DType, ParamTagType, AlgoTagType>(d, e, x)
          .run();

  Kokkos::fence();

  EXPECT_EQ(info, 0);

  RealType eps = ats::epsilon();

  auto h_x = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);

  // Check x = x_ref
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), ScalarType(1.0 / reference_value), eps);
    }
  }
}

}  // namespace Pttrs
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_pttrs() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Pttrs::impl_test_batched_pttrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Pttrs::impl_test_batched_pttrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 2; i++) {
      Test::Pttrs::impl_test_batched_pttrs_quick_return<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
          1, i);
      Test::Pttrs::impl_test_batched_pttrs_quick_return<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
          2, i);
    }
    for (int i = 2; i < 10; i++) {
      Test::Pttrs::impl_test_batched_pttrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, i);
      Test::Pttrs::impl_test_batched_pttrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Pttrs::impl_test_batched_pttrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Pttrs::impl_test_batched_pttrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 2; i++) {
      Test::Pttrs::impl_test_batched_pttrs_quick_return<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
          1, i);
      Test::Pttrs::impl_test_batched_pttrs_quick_return<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
          2, i);
    }
    for (int i = 2; i < 10; i++) {
      Test::Pttrs::impl_test_batched_pttrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, i);
      Test::Pttrs::impl_test_batched_pttrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_pttrs_l_float) {
  using algo_tag_type  = typename Algo::Pttrs::Unblocked;
  using param_tag_type = ::Test::Pttrs::ParamTag<Uplo::Lower>;

  test_batched_pttrs<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pttrs_u_float) {
  using algo_tag_type  = typename Algo::Pttrs::Unblocked;
  using param_tag_type = ::Test::Pttrs::ParamTag<Uplo::Upper>;
  test_batched_pttrs<TestDevice, float, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_pttrs_l_double) {
  using algo_tag_type  = typename Algo::Pttrs::Unblocked;
  using param_tag_type = ::Test::Pttrs::ParamTag<Uplo::Lower>;

  test_batched_pttrs<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pttrs_u_double) {
  using algo_tag_type  = typename Algo::Pttrs::Unblocked;
  using param_tag_type = ::Test::Pttrs::ParamTag<Uplo::Upper>;
  test_batched_pttrs<TestDevice, double, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_pttrs_l_fcomplex) {
  using algo_tag_type  = typename Algo::Pttrs::Unblocked;
  using param_tag_type = ::Test::Pttrs::ParamTag<Uplo::Lower>;

  test_batched_pttrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pttrs_u_fcomplex) {
  using algo_tag_type  = typename Algo::Pttrs::Unblocked;
  using param_tag_type = ::Test::Pttrs::ParamTag<Uplo::Upper>;
  test_batched_pttrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_pttrs_l_dcomplex) {
  using algo_tag_type  = typename Algo::Pttrs::Unblocked;
  using param_tag_type = ::Test::Pttrs::ParamTag<Uplo::Lower>;

  test_batched_pttrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pttrs_u_dcomplex) {
  using algo_tag_type  = typename Algo::Pttrs::Unblocked;
  using param_tag_type = ::Test::Pttrs::ParamTag<Uplo::Upper>;
  test_batched_pttrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
#endif
