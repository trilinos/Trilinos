// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Pbtrf.hpp"
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Pbtrf {

template <typename U>
struct ParamTag {
  using uplo = U;
};

template <typename DeviceType, typename ABViewType, typename ParamTagType, typename AlgoTagType>
struct Functor_BatchedSerialPbtrf {
  using execution_space = typename DeviceType::execution_space;
  ABViewType m_ab;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialPbtrf(const ABViewType &ab) : m_ab(ab) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k, int &info) const {
    auto sub_ab = Kokkos::subview(m_ab, k, Kokkos::ALL(), Kokkos::ALL());

    info += KokkosBatched::SerialPbtrf<typename ParamTagType::uplo, AlgoTagType>::invoke(sub_ab);
  }

  inline int run() {
    using value_type = typename ABViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPbtrf");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, m_ab.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename ScalarType, typename AViewType, typename BViewType, typename CViewType,
          typename ArgTransA, typename ArgTransB>
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
    auto aa = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(m_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(m_c, k, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialGemm<ArgTransA, ArgTransB, Algo::Gemm::Unblocked>::invoke(m_alpha, aa, bb, m_beta, cc);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPbtrf");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

/// \brief Implementation details of batched pbtrf analytical test
///        Confirm A = U**H * U or L * L**H, where
///        For conventional storage,
///        A: [[4, 1],
///            [1, 4]]
///        L: [[sqrt(4), 0],
///            [1/sqrt(4), sqrt(4 - (1/sqrt(4))**2)]
///        U: [[sqrt(4), 1/sqrt(4)],
///            [0, sqrt(4 - (1/sqrt(4))**2)]
///
///        For lower banded storage, Ab = Lb * Lb**H
///        Ab: [[4, 4],
///             [1, 0]]
///        Lb: [[sqrt(4), sqrt(4 - (1/sqrt(4))**2)],
///             [1/sqrt(4), 0]]
///
///        For upper banded storage, Ab = Ub**H * Ub
///        Ab: [[0, 1],
///             [4, 4]]
///        Ub: [[0, 1/sqrt(4)],
///             [sqrt(4), sqrt(4 - (1/sqrt(4))**2)]]
/// \param N [in] Batch size of AB
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_pbtrf_analytical(const int N) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  const int BlkSize = 2, k = 1;
  View3DType A("A", N, BlkSize, BlkSize), Ab("Ab", N, k + 1, BlkSize),
      Ab_ref("Ab_ref", N, k + 1, BlkSize);  // Banded matrix

  auto h_A = Kokkos::create_mirror_view(A);
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      for (int j = 0; j < BlkSize; j++) {
        h_A(ib, i, j) = i == j ? 4.0 : 1.0;
      }
    }
  }

  Kokkos::deep_copy(A, h_A);

  // Create banded triangluar matrix in normal and banded storage
  using ArgUplo = typename ParamTagType::uplo;
  create_banded_triangular_matrix<View3DType, View3DType, ArgUplo>(A, Ab, k, true);

  // Make a reference using the naive Cholesky decomposition
  // Cholesky decomposition for conventional storage
  // l_kk = np.sqrt( a_kk - sum_{i=1}^{k-1}( l_ik^2 ) )
  // l_ik = 1/l_kk * ( a_ik - sum_{j=1}^{k-1}( l_ij * l_kj ) )
  auto h_Ab_ref = Kokkos::create_mirror_view(Ab_ref);
  auto h_Ab     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Ab);
  if (std::is_same_v<typename ParamTagType::uplo, KokkosBatched::Uplo::Upper>) {
    // A = U**H * U
    for (int ib = 0; ib < N; ib++) {
      h_Ab_ref(ib, 1, 0) = Kokkos::sqrt(h_Ab(ib, 1, 0));
      h_Ab_ref(ib, 0, 1) = 1.0 / h_Ab_ref(ib, 1, 0);
      h_Ab_ref(ib, 1, 1) = Kokkos::sqrt(h_Ab(ib, 1, 1) - h_Ab_ref(ib, 0, 1) * h_Ab_ref(ib, 0, 1));
    }
  } else {
    // A = L * L**H
    for (int ib = 0; ib < N; ib++) {
      h_Ab_ref(ib, 0, 0) = Kokkos::sqrt(h_Ab(ib, 0, 0));
      h_Ab_ref(ib, 1, 0) = 1.0 / h_Ab_ref(ib, 0, 0);
      h_Ab_ref(ib, 0, 1) = Kokkos::sqrt(h_Ab(ib, 0, 1) - h_Ab_ref(ib, 1, 0) * h_Ab_ref(ib, 1, 0));
    }
  }

  // Factorize with Pbtrf: A = U**H * U or A = L * L**H
  auto info = Functor_BatchedSerialPbtrf<DeviceType, View3DType, ParamTagType, AlgoTagType>(Ab).run();
  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  // Check if Ab == Ub or Lb
  Kokkos::deep_copy(h_Ab, Ab);
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < k + 1; i++) {
      for (int j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_Ab(ib, i, j), h_Ab_ref(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched pbtrf test
///
/// \param N [in] Batch size of RHS
/// \param k [in] Number of superdiagonals or subdiagonals of matrix A
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_pbtrf(const int N, const int k, const int BlkSize) {
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  View3DType A("A", N, BlkSize, BlkSize), A_reconst("A_reconst", N, BlkSize, BlkSize),
      Ab("Ab", N, k + 1, BlkSize);  // Banded matrix

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  // Initialize A_reconst with random matrix
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);

  // Make the matrix Positive Definite Symmetric and Diagonal dominant
  random_to_pds(A, A_reconst);
  Kokkos::deep_copy(A, 0.0);

  // Create banded triangluar matrix in normal and banded storage
  using ArgUplo = typename ParamTagType::uplo;
  create_banded_pds_matrix<View3DType, View3DType, ArgUplo>(A_reconst, A, k, false);

  create_banded_triangular_matrix<View3DType, View3DType, ArgUplo>(A_reconst, Ab, k, true);

  // Clear matrix
  Kokkos::deep_copy(A_reconst, 0.0);

  // Factorize with Pbtrf: A = U**H * U or A = L * L**H
  auto info = Functor_BatchedSerialPbtrf<DeviceType, View3DType, ParamTagType, AlgoTagType>(Ab).run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // Compute A = U**H * U or A = L * L**H
  if (std::is_same_v<typename ParamTagType::uplo, KokkosBatched::Uplo::Upper>) {
    // A = U**H * U
    View3DType U("U", N, BlkSize, BlkSize);
    banded_to_dense<View3DType, View3DType, ArgUplo>(Ab, U, k);
    Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::ConjTranspose,
                              Trans::NoTranspose>(1.0, U, U, 0.0, A_reconst)
        .run();
  } else {
    // A = L * L**H
    View3DType L("L", N, BlkSize, BlkSize);
    banded_to_dense<View3DType, View3DType, ArgUplo>(Ab, L, k);
    Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose,
                              Trans::ConjTranspose>(1.0, L, L, 0.0, A_reconst)
        .run();
  }

  // this eps is about 10^-14
  using ats      = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType = typename ats::mag_type;
  RealType eps   = 1.0e3 * ats::epsilon();

  auto h_A         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  auto h_A_reconst = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A_reconst);

  // Check if A = U**H * U or A = L * L**H
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      for (int j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_A_reconst(ib, i, j), h_A(ib, i, j), eps);
      }
    }
  }
}

}  // namespace Pbtrf
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_pbtrf() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Pbtrf::impl_test_batched_pbtrf_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Pbtrf::impl_test_batched_pbtrf_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      int k = 1;
      Test::Pbtrf::impl_test_batched_pbtrf<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, k, i);
      Test::Pbtrf::impl_test_batched_pbtrf<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, k, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Pbtrf::impl_test_batched_pbtrf_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Pbtrf::impl_test_batched_pbtrf_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      int k = 1;
      Test::Pbtrf::impl_test_batched_pbtrf<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, k, i);
      Test::Pbtrf::impl_test_batched_pbtrf<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, k, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_pbtrf_l_float) {
  using algo_tag_type  = typename Algo::Pbtrf::Unblocked;
  using param_tag_type = ::Test::Pbtrf::ParamTag<Uplo::Lower>;

  test_batched_pbtrf<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pbtrf_u_float) {
  using algo_tag_type  = typename Algo::Pbtrf::Unblocked;
  using param_tag_type = ::Test::Pbtrf::ParamTag<Uplo::Upper>;

  test_batched_pbtrf<TestDevice, float, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_pbtrf_l_double) {
  using algo_tag_type  = typename Algo::Pbtrf::Unblocked;
  using param_tag_type = ::Test::Pbtrf::ParamTag<Uplo::Lower>;

  test_batched_pbtrf<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pbtrf_u_double) {
  using algo_tag_type  = typename Algo::Pbtrf::Unblocked;
  using param_tag_type = ::Test::Pbtrf::ParamTag<Uplo::Upper>;

  test_batched_pbtrf<TestDevice, double, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_pbtrf_l_fcomplex) {
  using algo_tag_type  = typename Algo::Pbtrf::Unblocked;
  using param_tag_type = ::Test::Pbtrf::ParamTag<Uplo::Lower>;

  test_batched_pbtrf<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pbtrf_u_fcomplex) {
  using algo_tag_type  = typename Algo::Pbtrf::Unblocked;
  using param_tag_type = ::Test::Pbtrf::ParamTag<Uplo::Upper>;

  test_batched_pbtrf<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_pbtrf_l_dcomplex) {
  using algo_tag_type  = typename Algo::Pbtrf::Unblocked;
  using param_tag_type = ::Test::Pbtrf::ParamTag<Uplo::Lower>;

  test_batched_pbtrf<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pbtrf_u_dcomplex) {
  using algo_tag_type  = typename Algo::Pbtrf::Unblocked;
  using param_tag_type = ::Test::Pbtrf::ParamTag<Uplo::Upper>;

  test_batched_pbtrf<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
#endif
