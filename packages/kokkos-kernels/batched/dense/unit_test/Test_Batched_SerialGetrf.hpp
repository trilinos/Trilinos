// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Getrf.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Getrf {

template <typename DeviceType, typename AViewType, typename PivViewType, typename AlgoTagType>
struct Functor_BatchedSerialGetrf {
  using execution_space = typename DeviceType::execution_space;
  AViewType m_a;
  PivViewType m_ipiv;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGetrf(const AViewType &a, const PivViewType &ipiv) : m_a(a), m_ipiv(ipiv) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    auto sub_a    = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_ipiv = Kokkos::subview(m_ipiv, k, Kokkos::ALL());

    info += KokkosBatched::SerialGetrf<AlgoTagType>::invoke(sub_a, sub_ipiv);
  }

  inline int run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGetrf");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_a.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
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

    KokkosBatched::SerialGemm<Trans::NoTranspose, Trans::NoTranspose, Algo::Gemm::Unblocked>::invoke(
        m_alpha, sub_a, sub_b, m_beta, sub_c);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGetrf");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

/// \brief Implementation details of batched getrf test
///        LU factorization with partial pivoting
///        4x4 matrix
///        A = [[1. 0. 0. 0.]
///             [0. 1. 0. 0.]
///             [0. 0. 1. 0.]
///             [0. 0. 0. 1.]]
///        LU = [[1. 0. 0. 0.]
///              [0. 1. 0. 0.]
///              [0. 0. 1. 0.]
///              [0. 0. 0. 1.]]
///        piv = [0 1 2 3]
///
///        3x4 matrix
///        A1 = [[1. 0. 0. 0.]
///              [0. 1. 0. 0.]
///              [0. 0. 1. 0.]]
///        LU1 = [[1. 0. 0. 0.]
///               [0. 1. 0. 0.]
///               [0. 0. 1. 0.]]
///        piv1 = [0 1 2]
///
///        4x3 matrix
///        A2 = [[1. 0. 0.]
///              [0. 1. 0.]
///              [0. 0. 1.]
///              [0. 0. 0.]]
///        LU2 = [[1. 0. 0.]
///               [0. 1. 0.]
///               [0. 0. 1.]
///               [0. 0. 0.]]
///        piv2 = [0 1 2]
///        3x3 more general matrix
///        which satisfies PA = LU
///        P = [[0 0 1]
///             [1 0 0]
///             [0 1 0]]
///        A = [[1  2  3]
///             [2 -4  6]
///             [3 -9 -3]]
///        L = [[1   0   0]
///             [1/3 1   0]
///             [2/3 2/5 1]]
///        U = [[-3  -9   -3]
///             [ 0   5    4]
///             [ 0   0 32/5]]
///        Note P is obtained by piv = [2 2 2]
///        We compare the non-diagnoal elements of L only, which is
///        NL = [[0   0   0]
///              [1/3 0   0]
///              [2/3 2/5 0]]
/// \param Nb [in] Batch size of matrices
template <typename DeviceType, typename ScalarType, typename LayoutType, typename AlgoTagType>
void impl_test_batched_getrf_analytical(const int Nb) {
  using ats           = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType      = typename ats::mag_type;
  using View3DType    = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using PivView2DType = Kokkos::View<int **, LayoutType, DeviceType>;

  constexpr int M = 4, N = 3;
  View3DType A0("A0", Nb, M, M), LU0("LU0", Nb, M, M);
  PivView2DType ipiv0("ipiv0", Nb, M), ipiv0_ref("ipiv0_ref", Nb, M);

  // Non-square matrix
  View3DType A1("A1", Nb, N, M), LU1("LU1", Nb, N, M);
  PivView2DType ipiv1("ipiv1", Nb, N), ipiv1_ref("ipiv1_ref", Nb, N);

  View3DType A2("A2", Nb, M, N), LU2("LU2", Nb, M, N);
  PivView2DType ipiv2("ipiv2", Nb, N), ipiv2_ref("ipiv1_ref", Nb, N);

  // Complicated matrix
  View3DType A3("A3", Nb, N, N), LU3("LU3", Nb, N, N), L3("L3", Nb, N, N), U3("U3", Nb, N, N),
      L3_ref("L3_ref", Nb, N, N), U3_ref("U3_ref", Nb, N, N);
  PivView2DType ipiv3("ipiv3", Nb, N), ipiv3_ref("ipiv3_ref", Nb, N);

  auto h_A0        = Kokkos::create_mirror_view(A0);
  auto h_A1        = Kokkos::create_mirror_view(A1);
  auto h_A2        = Kokkos::create_mirror_view(A2);
  auto h_A3        = Kokkos::create_mirror_view(A3);
  auto h_L3_ref    = Kokkos::create_mirror_view(L3_ref);
  auto h_U3_ref    = Kokkos::create_mirror_view(U3_ref);
  auto h_ipiv0_ref = Kokkos::create_mirror_view(ipiv0_ref);
  auto h_ipiv1_ref = Kokkos::create_mirror_view(ipiv1_ref);
  auto h_ipiv2_ref = Kokkos::create_mirror_view(ipiv2_ref);
  auto h_ipiv3_ref = Kokkos::create_mirror_view(ipiv3_ref);
  for (int ib = 0; ib < Nb; ib++) {
    for (int i = 0; i < M; i++) {
      h_ipiv0_ref(ib, i) = i;
      for (int j = 0; j < M; j++) {
        h_A0(ib, i, j) = i == j ? 1.0 : 0.0;
      }
    }

    for (int i = 0; i < N; i++) {
      h_ipiv1_ref(ib, i) = i;
      h_ipiv2_ref(ib, i) = i;
      for (int j = 0; j < M; j++) {
        h_A1(ib, i, j) = i == j ? 1.0 : 0.0;
        h_A2(ib, j, i) = i == j ? 1.0 : 0.0;
      }
    }

    h_A3(ib, 0, 0) = 1.0;
    h_A3(ib, 0, 1) = 2.0;
    h_A3(ib, 0, 2) = 3.0;
    h_A3(ib, 1, 0) = 2.0;
    h_A3(ib, 1, 1) = -4.0;
    h_A3(ib, 1, 2) = 6.0;
    h_A3(ib, 2, 0) = 3.0;
    h_A3(ib, 2, 1) = -9.0;
    h_A3(ib, 2, 2) = -3.0;

    h_L3_ref(ib, 0, 0) = 0.0;
    h_L3_ref(ib, 0, 1) = 0.0;
    h_L3_ref(ib, 0, 2) = 0.0;
    h_L3_ref(ib, 1, 0) = 1.0 / 3.0;
    h_L3_ref(ib, 1, 1) = 0.0;
    h_L3_ref(ib, 1, 2) = 0.0;
    h_L3_ref(ib, 2, 0) = 2.0 / 3.0;
    h_L3_ref(ib, 2, 1) = 2.0 / 5.0;
    h_L3_ref(ib, 2, 2) = 0.0;

    h_U3_ref(ib, 0, 0) = 3.0;
    h_U3_ref(ib, 0, 1) = -9.0;
    h_U3_ref(ib, 0, 2) = -3.0;
    h_U3_ref(ib, 1, 0) = 0.0;
    h_U3_ref(ib, 1, 1) = 5.0;
    h_U3_ref(ib, 1, 2) = 4.0;
    h_U3_ref(ib, 2, 0) = 0.0;
    h_U3_ref(ib, 2, 1) = 0.0;
    h_U3_ref(ib, 2, 2) = 32.0 / 5.0;

    h_ipiv3_ref(ib, 0) = 2;
    h_ipiv3_ref(ib, 1) = 2;
    h_ipiv3_ref(ib, 2) = 2;
  }

  Kokkos::deep_copy(A0, h_A0);
  Kokkos::deep_copy(A1, h_A1);
  Kokkos::deep_copy(A2, h_A2);
  Kokkos::deep_copy(A3, h_A3);
  Kokkos::deep_copy(LU0, A0);
  Kokkos::deep_copy(LU1, A1);
  Kokkos::deep_copy(LU2, A2);
  Kokkos::deep_copy(LU3, A3);

  // getrf to factorize matrix A = P * L * U
  auto info0 = Functor_BatchedSerialGetrf<DeviceType, View3DType, PivView2DType, AlgoTagType>(LU0, ipiv0).run();
  auto info1 = Functor_BatchedSerialGetrf<DeviceType, View3DType, PivView2DType, AlgoTagType>(LU1, ipiv1).run();
  auto info2 = Functor_BatchedSerialGetrf<DeviceType, View3DType, PivView2DType, AlgoTagType>(LU2, ipiv2).run();
  auto info3 = Functor_BatchedSerialGetrf<DeviceType, View3DType, PivView2DType, AlgoTagType>(LU3, ipiv3).run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);
  EXPECT_EQ(info2, 0);
  EXPECT_EQ(info3, 0);

  auto h_ipiv0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv0);
  auto h_ipiv1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv1);
  auto h_ipiv2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv2);
  auto h_ipiv3 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv3);

  for (int ib = 0; ib < Nb; ib++) {
    // Check if piv0 = [0 1 2 3]
    for (int i = 0; i < M; i++) {
      EXPECT_EQ(h_ipiv0(ib, i), h_ipiv0_ref(ib, i));
    }
    // Check if piv1 = [0 1 2] and piv2 = [0 1 2]
    for (int i = 0; i < N; i++) {
      EXPECT_EQ(h_ipiv1(ib, i), h_ipiv1_ref(ib, i));
      EXPECT_EQ(h_ipiv2(ib, i), h_ipiv2_ref(ib, i));
    }
    // Check if piv3 = [2 2 2]
    for (int i = 0; i < N; i++) {
      EXPECT_EQ(h_ipiv3(ib, i), h_ipiv3_ref(ib, i));
    }
  }

  // Reconstruct L and U from Factorized matrix A
  // Copy non-diagonal lower triangular components to NL
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Lower, KokkosBatched::Diag::NonUnit>(LU3, L3,
                                                                                                             -1);

  // Copy upper triangular components to U
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Upper, KokkosBatched::Diag::NonUnit>(LU3, U3);

  RealType eps = 1.0e1 * ats::epsilon();
  auto h_LU0   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LU0);
  auto h_LU1   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LU1);
  auto h_LU2   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LU2);

  // Check if LU = A (permuted)
  for (int ib = 0; ib < Nb; ib++) {
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < M; j++) {
        EXPECT_NEAR_KK(h_LU0(ib, i, j), h_A0(ib, i, j), eps);
      }
    }
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
        EXPECT_NEAR_KK(h_LU1(ib, i, j), h_A1(ib, i, j), eps);
        EXPECT_NEAR_KK(h_LU2(ib, j, i), h_A2(ib, j, i), eps);
      }
    }
  }

  // For complicated matrix, we compare L and U with reference L and U
  auto h_L3 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), L3);
  auto h_U3 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), U3);
  for (int ib = 0; ib < Nb; ib++) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        EXPECT_NEAR_KK(h_L3(ib, i, j), h_L3_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_U3(ib, i, j), h_U3_ref(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched getrf test
///        LU factorization with partial pivoting
///
/// \param N [in] Batch size of matrices
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename AlgoTagType>
void impl_test_batched_getrf(const int N, const int BlkSize) {
  using ats            = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType       = typename ats::mag_type;
  using RealView2DType = Kokkos::View<RealType **, LayoutType, DeviceType>;
  using View3DType     = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using PivView2DType  = Kokkos::View<int **, LayoutType, DeviceType>;

  View3DType A("A", N, BlkSize, BlkSize), A_reconst("A_reconst", N, BlkSize, BlkSize), NL("NL", N, BlkSize, BlkSize),
      L("L", N, BlkSize, BlkSize), U("U", N, BlkSize, BlkSize), LU("LU", N, BlkSize, BlkSize),
      I("I", N, BlkSize, BlkSize);
  RealView2DType ones(Kokkos::view_alloc("ones", Kokkos::WithoutInitializing), N, BlkSize);
  PivView2DType ipiv("ipiv", N, BlkSize);

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  // Initialize A_reconst with random matrix
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::deep_copy(LU, A);

  // Unit matrix I
  Kokkos::deep_copy(ones, RealType(1.0));
  create_diagonal_matrix(ones, I);

  Kokkos::fence();

  // getrf to factorize matrix A = P * L * U
  auto info = Functor_BatchedSerialGetrf<DeviceType, View3DType, PivView2DType, AlgoTagType>(LU, ipiv).run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // Reconstruct L and U from Factorized matrix A
  // Copy non-diagonal lower triangular components to NL
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Lower, KokkosBatched::Diag::NonUnit>(LU, NL,
                                                                                                             -1);

  // Copy upper triangular components to U
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Upper, KokkosBatched::Diag::NonUnit>(LU, U);

  // Copy I to L
  Kokkos::deep_copy(L, I);

  // Matrix matrix addition by Gemm
  // NL + I by NL * I + L (==I) (result stored in L)
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType>(1.0, NL, I, 1.0, L).run();

  // LU = L * U
  Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType>(1.0, L, U, 0.0, LU).run();

  Kokkos::fence();

  // permute A by ipiv
  auto h_ipiv = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv);
  auto h_A    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  for (int ib = 0; ib < N; ib++) {
    // Permute A by pivot vector
    for (int i = 0; i < BlkSize; i++) {
      for (int j = 0; j < BlkSize; j++) {
        Kokkos::kokkos_swap(h_A(ib, h_ipiv(ib, i), j), h_A(ib, i, j));
      }
    }
  }

  RealType eps = 1.0e1 * ats::epsilon();

  auto h_LU = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), LU);
  // Check if LU = A (permuted)
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      for (int j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_LU(ib, i, j), h_A(ib, i, j), eps);
      }
    }
  }
}

}  // namespace Getrf
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename AlgoTagType>
int test_batched_getrf() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Getrf::impl_test_batched_getrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(1);
    Test::Getrf::impl_test_batched_getrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Getrf::impl_test_batched_getrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Getrf::impl_test_batched_getrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Getrf::impl_test_batched_getrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(1);
    Test::Getrf::impl_test_batched_getrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Getrf::impl_test_batched_getrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Getrf::impl_test_batched_getrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_getrf_float) {
  using algo_tag_type = typename KokkosBatched::Algo::Getrf::Unblocked;

  test_batched_getrf<TestDevice, float, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_getrf_double) {
  using algo_tag_type = typename KokkosBatched::Algo::Getrf::Unblocked;

  test_batched_getrf<TestDevice, double, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_getrf_fcomplex) {
  using algo_tag_type = typename KokkosBatched::Algo::Getrf::Unblocked;

  test_batched_getrf<TestDevice, Kokkos::complex<float>, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_getrf_dcomplex) {
  using algo_tag_type = typename KokkosBatched::Algo::Getrf::Unblocked;

  test_batched_getrf<TestDevice, Kokkos::complex<double>, algo_tag_type>();
}
#endif
