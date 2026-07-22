// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Pttrf.hpp"
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Pttrf {

template <typename DeviceType, typename DViewType, typename EViewType, typename AlgoTagType>
struct Functor_BatchedSerialPttrf {
  using execution_space = typename DeviceType::execution_space;
  DViewType m_d;
  EViewType m_e;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialPttrf(const DViewType &d, const EViewType &e) : m_d(d), m_e(e) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    auto sub_d = Kokkos::subview(m_d, k, Kokkos::ALL());
    auto sub_e = Kokkos::subview(m_e, k, Kokkos::ALL());

    info += KokkosBatched::SerialPttrf<AlgoTagType>::invoke(sub_d, sub_e);
  }

  inline int run() {
    using value_type = typename EViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPttrf");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_d.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
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

/// \brief Implementation details of batched pttrf test for random matrix
///
/// \param N [in] Batch size of matrix A
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename AlgoTagType>
void impl_test_batched_pttrf(const int N, const int BlkSize) {
  using ats            = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType     = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  View3DType A("A", N, BlkSize, BlkSize), A_reconst("A_reconst", N, BlkSize, BlkSize);
  View3DType EL("EL", N, BlkSize, BlkSize), EU("EU", N, BlkSize, BlkSize), D("D", N, BlkSize, BlkSize),
      LD("LD", N, BlkSize, BlkSize), L("L", N, BlkSize, BlkSize), I("I", N, BlkSize, BlkSize);
  RealView2DType d("d", N, BlkSize),  // Diagonal components
      ones(Kokkos::view_alloc("ones", Kokkos::WithoutInitializing), N, BlkSize);
  View2DType e_upper("e_upper", N, BlkSize - 1), e_lower("e_lower", N,
                                                         BlkSize - 1);  // upper and lower diagonal components

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  RealType realRandStart, realRandEnd;
  ScalarType randStart, randEnd;

  KokkosKernels::Impl::getRandomBounds(1.0, realRandStart, realRandEnd);
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);

  // Add BlkSize to ensure positive definiteness
  Kokkos::fill_random(d, rand_pool, realRandStart + BlkSize, realRandEnd + BlkSize);
  Kokkos::fill_random(e_upper, rand_pool, randStart, randEnd);

  auto h_e_upper = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), e_upper);
  auto h_e_lower = Kokkos::create_mirror_view(e_lower);

  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize - 1; i++) {
      // Fill the lower diagonal with conjugate of the upper diagonal
      h_e_lower(ib, i) = KokkosKernels::ArithTraits<ScalarType>::conj(h_e_upper(ib, i));
    }
  }

  Kokkos::deep_copy(e_lower, h_e_lower);
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

  // Factorize matrix A -> L * D * L**H
  // d and e are updated by pttrf
  auto info = Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e_lower).run();

  Kokkos::fence();

  EXPECT_EQ(info, 0);

  // Reconstruct L and D from factorized matrix
  create_diagonal_matrix(e_lower, EL, -1);
  create_diagonal_matrix(d, D);

  // Copy I to L
  Kokkos::deep_copy(L, I);

  // EL + I by EL * I + L (result stored in L)
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, EL, I,
                                                                                                            1.0, L)
      .run();

  // Reconstruct A by L*D*L**H
  // Gemm to compute L*D -> LD
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose>(1.0, L, D,
                                                                                                            0.0, LD)
      .run();

  // Gemm to compute (L*D)*L**H -> A_reconst
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::ConjTranspose>(
      1.0, LD, L, 0.0, A_reconst)
      .run();
  Kokkos::fence();

  RealType eps = 1.0e3 * ats::epsilon();

  // Check A = L*D*L**H
  auto h_A         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  auto h_A_reconst = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_reconst);
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      for (int j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_A_reconst(ib, i, j), h_A(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched pttrf test for early return
///        BlkSize must be 0 or 1
///
/// \param N [in] Batch size of matrix A
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename AlgoTagType>
void impl_test_batched_pttrf_quick_return(const int N, const int BlkSize) {
  using ats            = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;

  if (BlkSize > 1) return;

  const int BlkSize_minus_1 = BlkSize > 0 ? BlkSize - 1 : 0;

  RealView2DType d("d", N, BlkSize), d2("d2", N, BlkSize);  // Diagonal components
  View2DType e("e", N, BlkSize_minus_1);                    // lower diagonal components

  const RealType reference_value = 4.0;

  Kokkos::deep_copy(d, reference_value);
  Kokkos::deep_copy(d2, -reference_value);
  Kokkos::deep_copy(e, ScalarType(1.0));

  // Factorize matrix A -> L * D * L**H
  // d and e are updated by pttrf
  // Early return if BlkSize is 0 or 1
  auto info = Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e).run();

  // For negative values, info should be 1 for BlkSize = 1
  auto info2 = Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d2, e).run();

  Kokkos::fence();

  int expected_info2 = BlkSize == 0 ? 0 : N;
  EXPECT_EQ(info, 0);
  EXPECT_EQ(info2, expected_info2);

  RealType eps = 1.0e3 * ats::epsilon();

  // Check if d is unchanged
  auto h_d  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d);
  auto h_d2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d2);
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_d(ib, i), reference_value, eps);
      EXPECT_NEAR_KK(h_d2(ib, i), -reference_value, eps);
    }
  }
}

/// \brief Implementation details of batched pttrf test
/// Confirm that A = L*D*L**T (or L*D*L**H) where
/// for conventional storage,
/// A: [[4, 1, 0],
///     [1, 4, 1],
///     [0, 1, 4]]
/// D: [[4, 0, 0],
///     [0, 4-1*(1/4), 0],
///     [0, 0, 4-1*(1/(4-1*(1/4)))]]
/// L: [[0, 0, 0],
///     [1/4, 0, 0],
///     [0, 1/(4-1*(1/4)), 0]]
///
/// for packed storage for tridiagonal matrix,
/// AD: [4, 4, 4], AE: [1, 1]
///
/// D: [4, 4-1*(1/4), 4-1*(1/(4-1*(1/4)))]
/// L: [1/4, 1/(4-1*(1/4))]
///
/// \param N [in] Batch size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename AlgoTagType>
void impl_test_batched_pttrf_analytical(const int N) {
  using ats            = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View2DType     = Kokkos::View<ScalarType **, LayoutType, DeviceType>;

  const int BlkSize = 3;

  // Diagonal components
  RealView2DType d(Kokkos::view_alloc("d", Kokkos::WithoutInitializing), N, BlkSize), d_ref("d_ref", N, BlkSize);
  // Upper and lower diagonal components (identical)
  View2DType e(Kokkos::view_alloc("e", Kokkos::WithoutInitializing), N, BlkSize - 1), e_ref("e_ref", N, BlkSize - 1);

  Kokkos::deep_copy(d, RealType(4.0));
  Kokkos::deep_copy(e, ScalarType(1.0));

  // Make a reference using the analytical formula
  auto h_d     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d);
  auto h_e     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), e);
  auto h_d_ref = Kokkos::create_mirror_view(d_ref);
  auto h_e_ref = Kokkos::create_mirror_view(e_ref);

  for (int ib = 0; ib < N; ib++) {
    h_d_ref(ib, 0) = h_d(ib, 0);

    h_e_ref(ib, 0) = h_e(ib, 0) / h_d(ib, 0);
    h_d_ref(ib, 1) = h_d(ib, 1) - ats::real(h_e_ref(ib, 0)) * ats::real(h_e(ib, 0));

    h_e_ref(ib, 1) = h_e(ib, 1) / h_d_ref(ib, 1);
    h_d_ref(ib, 2) = h_d(ib, 2) - ats::real(h_e_ref(ib, 1)) * ats::real(h_e(ib, 1));
  }

  // Factorize matrix A -> L * D * L**T
  // d and e are updated by pttrf
  auto info = Functor_BatchedSerialPttrf<DeviceType, RealView2DType, View2DType, AlgoTagType>(d, e).run();
  Kokkos::fence();

  EXPECT_EQ(info, 0);

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  // Check d and e are equal to the analytical values
  Kokkos::deep_copy(h_d, d);
  Kokkos::deep_copy(h_e, e);
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_d(ib, i), h_d_ref(ib, i), eps);
    }
    for (int i = 0; i < BlkSize - 1; i++) {
      EXPECT_NEAR_KK(h_e(ib, i), h_e_ref(ib, i), eps);
    }
  }
}

}  // namespace Pttrf
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename AlgoTagType>
int test_batched_pttrf() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    for (int i = 0; i < 2; i++) {
      Test::Pttrf::impl_test_batched_pttrf_quick_return<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Pttrf::impl_test_batched_pttrf_quick_return<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
    Test::Pttrf::impl_test_batched_pttrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(1);
    Test::Pttrf::impl_test_batched_pttrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(2);
    for (int i = 2; i < 10; i++) {
      Test::Pttrf::impl_test_batched_pttrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Pttrf::impl_test_batched_pttrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    for (int i = 0; i < 2; i++) {
      Test::Pttrf::impl_test_batched_pttrf_quick_return<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Pttrf::impl_test_batched_pttrf_quick_return<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
    Test::Pttrf::impl_test_batched_pttrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(1);
    Test::Pttrf::impl_test_batched_pttrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(2);
    for (int i = 2; i < 10; i++) {
      Test::Pttrf::impl_test_batched_pttrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Pttrf::impl_test_batched_pttrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_pttrf_float) {
  using algo_tag_type = typename Algo::Pttrf::Unblocked;

  test_batched_pttrf<TestDevice, float, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_pttrf_double) {
  using algo_tag_type = typename Algo::Pttrf::Unblocked;

  test_batched_pttrf<TestDevice, double, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_pttrf_fcomplex) {
  using algo_tag_type = typename Algo::Pttrf::Unblocked;

  test_batched_pttrf<TestDevice, Kokkos::complex<float>, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_pttrf_dcomplex) {
  using algo_tag_type = typename Algo::Pttrf::Unblocked;

  test_batched_pttrf<TestDevice, Kokkos::complex<double>, algo_tag_type>();
}
#endif
