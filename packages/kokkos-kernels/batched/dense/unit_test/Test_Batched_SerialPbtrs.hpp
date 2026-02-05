// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas2_gemv.hpp>
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Pbtrf.hpp"
#include "KokkosBatched_Pbtrs.hpp"
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Pbtrs {

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
  void operator()(const ParamTagType &, const int k) const {
    auto sub_ab = Kokkos::subview(m_ab, k, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialPbtrf<typename ParamTagType::uplo, AlgoTagType>::invoke(sub_ab);
  }

  inline void run() {
    using value_type = typename ABViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPbtrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, m_ab.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename ABViewType, typename BViewType, typename ParamTagType, typename AlgoTagType>
struct Functor_BatchedSerialPbtrs {
  using execution_space = typename DeviceType::execution_space;
  ABViewType m_ab;
  BViewType m_b;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialPbtrs(const ABViewType &ab, const BViewType &b) : m_ab(ab), m_b(b) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k, int &info) const {
    auto sub_ab = Kokkos::subview(m_ab, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb     = Kokkos::subview(m_b, k, Kokkos::ALL());

    info += KokkosBatched::SerialPbtrs<typename ParamTagType::uplo, AlgoTagType>::invoke(sub_ab, bb);
  }

  inline int run() {
    using value_type = typename ABViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPbtrs");
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
    auto aa = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto xx = Kokkos::subview(m_x, k, Kokkos::ALL());
    auto yy = Kokkos::subview(m_y, k, Kokkos::ALL());

    KokkosBlas::SerialGemv<Trans::NoTranspose, Algo::Gemv::Unblocked>::invoke(m_alpha, aa, xx, m_beta, yy);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialPbtrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_x.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

/// \brief Implementation details of batched pbtrs test
/// Confirm A * x = b, where
///        A: [[4, 1, 0],
///            [1, 4, 1],
///            [0, 1, 4]]
///        b: [1, 1, 1]
///        x: [3/14, 1/7, 3/14]
///
/// This corresponds to the following system of equations:
///        4 x0 +   x1        = 1
///          x0 + 4 x1 +   x2 = 1
///                 x1 + 4 x2 = 1
///
/// We confirm this with the factorized band matrix Ub or Lb.
///        For upper banded storage, Ab = Ub**H * Ub
///        Ub: [[0, 1/sqrt(4), 1/sqrt(4 - (1/sqrt(4))**2)],
///             [sqrt(4), sqrt(4 - (1/sqrt(4))**2), sqrt(4 - 1/sqrt(4 - (1/sqrt(4))**2))],]
///        For lower banded storage, Ab = Lb * Lb**H
///        Lb: [[sqrt(4), sqrt(4 - (1/sqrt(4))**2), sqrt(4 - 1/sqrt(4 - (1/sqrt(4))**2))],
///             [1/sqrt(4), 1/sqrt(4 - (1/sqrt(4))**2), 0],]
///
/// \param N [in] Batch size of RHS (banded matrix can also be batched matrix)
/// \param k [in] Number of superdiagonals or subdiagonals of matrix A
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_pbtrs_analytical(const int N) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View2DType = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  const int BlkSize = 3, k = 1;
  View3DType Ab("Ab", N, k + 1, BlkSize);                       // In band storage
  View2DType x0("x0", N, BlkSize), x_ref("x_ref", N, BlkSize);  // Solutions

  auto h_Ab    = Kokkos::create_mirror_view(Ab);
  auto h_x_ref = Kokkos::create_mirror_view(x_ref);

  for (int ib = 0; ib < N; ib++) {
    if (std::is_same_v<typename ParamTagType::uplo, KokkosBatched::Uplo::Upper>) {
      // Ub
      h_Ab(ib, 1, 0) = Kokkos::sqrt(4.0);
      h_Ab(ib, 0, 1) = 1.0 / h_Ab(ib, 1, 0);
      h_Ab(ib, 1, 1) = Kokkos::sqrt(4.0 - h_Ab(ib, 0, 1) * h_Ab(ib, 0, 1));
      h_Ab(ib, 0, 2) = 1.0 / h_Ab(ib, 1, 1);
      h_Ab(ib, 1, 2) = Kokkos::sqrt(4.0 - h_Ab(ib, 0, 2) * h_Ab(ib, 0, 2));
    } else {
      // Lb
      h_Ab(ib, 0, 0) = Kokkos::sqrt(4.0);
      h_Ab(ib, 1, 0) = 1.0 / h_Ab(ib, 0, 0);
      h_Ab(ib, 0, 1) = Kokkos::sqrt(4.0 - h_Ab(ib, 1, 0) * h_Ab(ib, 1, 0));
      h_Ab(ib, 1, 1) = 1.0 / h_Ab(ib, 0, 1);
      h_Ab(ib, 0, 2) = Kokkos::sqrt(4.0 - h_Ab(ib, 1, 1) * h_Ab(ib, 1, 1));
    }

    h_x_ref(ib, 0) = 3.0 / 14.0;
    h_x_ref(ib, 1) = 1.0 / 7.0;
    h_x_ref(ib, 2) = 3.0 / 14.0;
  }

  Kokkos::fence();

  Kokkos::deep_copy(x0, ScalarType(1.0));
  Kokkos::deep_copy(Ab, h_Ab);

  // pbtrs (Note, Ab is a factorized matrix of A)
  auto info = Functor_BatchedSerialPbtrs<DeviceType, View3DType, View2DType, ParamTagType, AlgoTagType>(Ab, x0).run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();
  auto h_x0    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x0);

  // Check x0 = x1
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      Test::EXPECT_NEAR_KK_REL(h_x0(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched pbtrs test
///        Confirm A * x = b, where A is a real symmetric positive definitie
///        or complex Hermitian band matrix. A is storead in a band storage.
///        A must be factorized as A=U**H*U or A=L*L**H (Cholesky factorization)
///        by pbtrf.
///
/// \param N [in] Batch size of RHS (banded matrix can also be batched matrix)
/// \param k [in] Number of superdiagonals or subdiagonals of matrix A
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_pbtrs(const int N, const int k, const int BlkSize) {
  using ats        = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType   = typename ats::mag_type;
  using View2DType = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;

  View3DType A("A", N, BlkSize, BlkSize), A_reconst("A_reconst", N, BlkSize, BlkSize);
  View3DType Ab("Ab", N, k + 1, BlkSize);                                             // Banded matrix
  View2DType x0("x0", N, BlkSize), x_ref("x_ref", N, BlkSize), y0("y0", N, BlkSize);  // Solutions

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  // Initialize A_reconst with random matrix
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x0, rand_pool, randStart, randEnd);
  Kokkos::deep_copy(x_ref, x0);

  // Make the matrix Positive Definite Symmetric and Diagonal dominant
  random_to_pds(A, A_reconst);
  Kokkos::deep_copy(A, ScalarType(0.0));

  // Create banded triangluar matrix in normal and banded storage
  using ArgUplo = typename ParamTagType::uplo;
  create_banded_pds_matrix<View3DType, View3DType, ArgUplo>(A_reconst, A, k, false);

  create_banded_triangular_matrix<View3DType, View3DType, ArgUplo>(A_reconst, Ab, k, true);

  Kokkos::fence();

  // Factorize with Pbtrf: A = U**H * U or A = L * L**H
  Functor_BatchedSerialPbtrf<DeviceType, View3DType, ParamTagType, AlgoTagType>(Ab).run();

  // pbtrs (Note, Ab is a factorized matrix of A)
  auto info = Functor_BatchedSerialPbtrs<DeviceType, View3DType, View2DType, ParamTagType, AlgoTagType>(Ab, x0).run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // Gemv to compute A*x0, this should be identical to x_ref
  Functor_BatchedSerialGemv<DeviceType, ScalarType, View3DType, View2DType, View2DType>(1.0, A, x0, 0.0, y0).run();

  Kokkos::fence();

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  auto h_y0    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y0);
  auto h_x_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x_ref);

  // Check A * x0 = x_ref
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      Test::EXPECT_NEAR_KK_REL(h_y0(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

}  // namespace Pbtrs
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_pbtrs() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Pbtrs::impl_test_batched_pbtrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Pbtrs::impl_test_batched_pbtrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      int k = 1;
      Test::Pbtrs::impl_test_batched_pbtrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, k, i);
      Test::Pbtrs::impl_test_batched_pbtrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, k, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Pbtrs::impl_test_batched_pbtrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Pbtrs::impl_test_batched_pbtrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      int k = 1;
      Test::Pbtrs::impl_test_batched_pbtrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, k, i);
      Test::Pbtrs::impl_test_batched_pbtrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, k, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_pbtrs_l_float) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<KokkosBatched::Uplo::Lower>;

  test_batched_pbtrs<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pbtrs_u_float) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<KokkosBatched::Uplo::Upper>;

  test_batched_pbtrs<TestDevice, float, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_pbtrs_l_double) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<KokkosBatched::Uplo::Lower>;

  test_batched_pbtrs<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pbtrs_u_double) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<KokkosBatched::Uplo::Upper>;

  test_batched_pbtrs<TestDevice, double, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_pbtrs_l_fcomplex) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<KokkosBatched::Uplo::Lower>;

  test_batched_pbtrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pbtrs_u_fcomplex) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<KokkosBatched::Uplo::Upper>;

  test_batched_pbtrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_pbtrs_l_dcomplex) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<KokkosBatched::Uplo::Lower>;

  test_batched_pbtrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_pbtrs_u_dcomplex) {
  using algo_tag_type  = typename Algo::Pbtrs::Unblocked;
  using param_tag_type = ::Test::Pbtrs::ParamTag<KokkosBatched::Uplo::Upper>;

  test_batched_pbtrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
#endif
