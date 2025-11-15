// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Gbtrf.hpp>
#include <KokkosBatched_Gbtrs.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Gbtrs {

template <typename T>
struct ParamTag {
  using trans = T;
};

template <typename DeviceType, typename ABViewType, typename PivViewType>
struct Functor_BatchedSerialGbtrf {
  using execution_space = typename DeviceType::execution_space;
  ABViewType m_ab;
  PivViewType m_ipiv;
  int m_kl, m_ku;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGbtrf(const ABViewType &ab, const PivViewType &ipiv, int kl, int ku)
      : m_ab(ab), m_ipiv(ipiv), m_kl(kl), m_ku(ku) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto ab   = Kokkos::subview(m_ab, k, Kokkos::ALL(), Kokkos::ALL());
    auto ipiv = Kokkos::subview(m_ipiv, k, Kokkos::ALL());

    KokkosBatched::SerialGbtrf<Algo::Gbtrf::Unblocked>::invoke(ab, ipiv, m_kl, m_ku);
  }

  inline void run() {
    using value_type = typename ABViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGbtrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_ab.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename AViewType, typename BViewType, typename PivViewType, typename ParamTagType,
          typename AlgoTagType>
struct Functor_BatchedSerialGbtrs {
  using execution_space = typename DeviceType::execution_space;
  AViewType m_a;
  BViewType m_b;
  PivViewType m_ipiv;
  int m_kl, m_ku;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGbtrs(const AViewType &a, const PivViewType &ipiv, const BViewType &b, int kl, int ku)
      : m_a(a), m_b(b), m_ipiv(ipiv), m_kl(kl), m_ku(ku) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k, int &info) const {
    auto aa   = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb   = Kokkos::subview(m_b, k, Kokkos::ALL());
    auto ipiv = Kokkos::subview(m_ipiv, k, Kokkos::ALL());

    info += KokkosBatched::SerialGbtrs<typename ParamTagType::trans, AlgoTagType>::invoke(aa, ipiv, bb, m_kl, m_ku);
  }

  inline int run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGbtrs");
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
    std::string name_region("KokkosBatched::Test::SerialGbtrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, m_x.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

/// \brief Implementation details of batched gbtrs test
///        Confirm A * x = b, where
///        A: [[1, -3, -2,  0],
///            [-1, 1, -3, -2],
///            [2, -1,  1, -3],
///            [0,  2, -1,  1]]
///        b: [1, 1, 1, 1]
///        x: [67/81, 22/81, -40/81, -1/27] or [-1/27, -40/81, 22/81, 67/81]
///
///        This corresponds to the following system of equations:
///          x0 - 3 x1 - 2 x2        = 1
///        - x0 +   x1 - 3 x2 - 2 x3 = 1
///        2 x0 -   x1 +   x3 - 3 x3 = 1
///               2 x1 -   x2 +   x3 = 1
///
/// We confirm this with the factorized matrix LUB and pivot given by
///       LUB: [[0,       0,    0,    0],
///             [0,       0,    0,   -3],
///             [0,       0,    1,  1.5],
///             [0,      -1, -2.5, -3.2],
///             [2,    -2.5,   -3,  5.4],
///             [-0.5, -0.2,    1,    0],
///             [0.5,  -0.8,    0,    0]]
///       piv: [2, 2, 2, 3]
///
/// \param N [in] Batch size of RHS
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_gbtrs_analytical(const int N) {
  using ats         = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType    = typename ats::mag_type;
  using View2DType  = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType  = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using PivViewType = Kokkos::View<int **, LayoutType, DeviceType>;

  const int BlkSize = 4, kl = 2, ku = 2;
  const int ldab = 2 * kl + ku + 1;
  View3DType AB("AB", N, ldab, BlkSize);                        // Banded matrix
  View2DType x0("x0", N, BlkSize), x_ref("x_ref", N, BlkSize);  // Solutions
  PivViewType ipiv("ipiv", N, BlkSize);

  using ArgTrans = typename ParamTagType::trans;
  auto h_AB      = Kokkos::create_mirror_view(AB);
  auto h_ipiv    = Kokkos::create_mirror_view(ipiv);
  auto h_x_ref   = Kokkos::create_mirror_view(x_ref);
  for (int ib = 0; ib < N; ib++) {
    h_AB(ib, 1, 3) = -3.0;
    h_AB(ib, 2, 2) = 1.0;
    h_AB(ib, 2, 3) = 1.5;
    h_AB(ib, 3, 1) = -1.0;
    h_AB(ib, 3, 2) = -2.5;
    h_AB(ib, 3, 3) = -3.2;
    h_AB(ib, 4, 0) = 2.0;
    h_AB(ib, 4, 1) = -2.5;
    h_AB(ib, 4, 2) = -3.0;
    h_AB(ib, 4, 3) = 5.4;
    h_AB(ib, 5, 0) = -0.5;
    h_AB(ib, 5, 1) = -0.2;
    h_AB(ib, 5, 2) = 1.0;
    h_AB(ib, 6, 0) = 0.5;
    h_AB(ib, 6, 1) = -0.8;

    h_ipiv(ib, 0) = 2;
    h_ipiv(ib, 1) = 2;
    h_ipiv(ib, 2) = 2;
    h_ipiv(ib, 3) = 3;
    if (std::is_same_v<ArgTrans, KokkosBatched::Trans::NoTranspose>) {
      h_x_ref(ib, 0) = 67.0 / 81.0;
      h_x_ref(ib, 1) = 22.0 / 81.0;
      h_x_ref(ib, 2) = -40.0 / 81.0;
      h_x_ref(ib, 3) = -1.0 / 27.0;
    } else {
      h_x_ref(ib, 0) = -1.0 / 27.0;
      h_x_ref(ib, 1) = -40.0 / 81.0;
      h_x_ref(ib, 2) = 22.0 / 81.0;
      h_x_ref(ib, 3) = 67.0 / 81.0;
    }
  }

  Kokkos::fence();

  Kokkos::deep_copy(AB, h_AB);
  Kokkos::deep_copy(ipiv, h_ipiv);
  Kokkos::deep_copy(x0, ScalarType(1.0));

  // gbtrs (Note, Ab is a factorized matrix of A)
  auto info = Functor_BatchedSerialGbtrs<DeviceType, View3DType, View2DType, PivViewType, ParamTagType, AlgoTagType>(
                  AB, ipiv, x0, kl, ku)
                  .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();
  auto h_x0    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x0);

  // Check x0 = [67/81, 22/81, -40/81, -1/27] or [-1/27, -40/81, 22/81, 67/81]
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x0(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched gbtrs test
///        Confirm A * x = b, or A**T * x = b, or A**H * x = b, where
/// \param N [in] Batch size of RHS
/// \param k [in] Number of superdiagonals and subdiagonals of matrix A
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_gbtrs(const int N, const int k, const int BlkSize) {
  using ats         = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType    = typename ats::mag_type;
  using View2DType  = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType  = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using PivViewType = Kokkos::View<int **, LayoutType, DeviceType>;

  const int kl = k, ku = k;
  const int ldab = 2 * kl + ku + 1;
  View3DType A("A", N, BlkSize, BlkSize), tmp_A("tmp_A", N, BlkSize, BlkSize),
      AB("AB", N, ldab, BlkSize);                         // Banded matrix
  View2DType x0("x0", N, BlkSize), y0("y0", N, BlkSize);  // Solutions
  PivViewType ipiv("ipiv", N, BlkSize);

  // Create a random matrix A and make it Positive Definite Symmetric
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  // Initialize tmp_A with random matrix
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(tmp_A, rand_pool, randStart, randEnd);

  dense_to_banded(tmp_A, AB, kl, ku);  // In banded storage
  banded_to_dense(AB, A, kl, ku);      // In conventional storage

  Kokkos::fence();

  // Create an initial solution vector x0 = [1, 1, 1, ...]
  Kokkos::deep_copy(x0, ScalarType(1.0));

  // gbtrf to factorize matrix A = P * L * U
  Functor_BatchedSerialGbtrf<DeviceType, View3DType, PivViewType>(AB, ipiv, kl, ku).run();

  // gbtrs (Note, Ab is a factorized matrix of A)
  auto info = Functor_BatchedSerialGbtrs<DeviceType, View3DType, View2DType, PivViewType, ParamTagType, AlgoTagType>(
                  AB, ipiv, x0, kl, ku)
                  .run();
  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // Gemv to compute A*x0, this should be identical to x_ref
  Functor_BatchedSerialGemv<DeviceType, ScalarType, View3DType, View2DType, View2DType, ParamTagType>(1.0, A, x0, 0.0,
                                                                                                      y0)
      .run();

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  // Check A * x0 = x_ref
  auto h_y0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y0);
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_y0(ib, i), ScalarType(1.0), eps);
    }
  }
}

}  // namespace Gbtrs
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_gbtrs() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Gbtrs::impl_test_batched_gbtrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(0);
    Test::Gbtrs::impl_test_batched_gbtrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Gbtrs::impl_test_batched_gbtrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      for (int k = 1; k < 4; k++) {
        Test::Gbtrs::impl_test_batched_gbtrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(0, k, i);
        Test::Gbtrs::impl_test_batched_gbtrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, k, i);
        Test::Gbtrs::impl_test_batched_gbtrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, k, i);
      }
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Gbtrs::impl_test_batched_gbtrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(0);
    Test::Gbtrs::impl_test_batched_gbtrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Gbtrs::impl_test_batched_gbtrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      for (int k = 1; k < 4; k++) {
        Test::Gbtrs::impl_test_batched_gbtrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(0, k, i);
        Test::Gbtrs::impl_test_batched_gbtrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, k, i);
        Test::Gbtrs::impl_test_batched_gbtrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, k, i);
      }
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_gbtrs_nt_float) {
  using param_tag_type = ::Test::Gbtrs::ParamTag<Trans::NoTranspose>;
  using algo_tag_type  = typename Algo::Gbtrs::Unblocked;

  test_batched_gbtrs<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_gbtrs_t_float) {
  using param_tag_type = ::Test::Gbtrs::ParamTag<Trans::Transpose>;
  using algo_tag_type  = typename Algo::Gbtrs::Unblocked;

  test_batched_gbtrs<TestDevice, float, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_gbtrs_nt_double) {
  using param_tag_type = ::Test::Gbtrs::ParamTag<Trans::NoTranspose>;
  using algo_tag_type  = typename Algo::Gbtrs::Unblocked;

  test_batched_gbtrs<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_gbtrs_t_double) {
  using param_tag_type = ::Test::Gbtrs::ParamTag<Trans::Transpose>;
  using algo_tag_type  = typename Algo::Gbtrs::Unblocked;

  test_batched_gbtrs<TestDevice, double, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_gbtrs_nt_fcomplex) {
  using param_tag_type = ::Test::Gbtrs::ParamTag<Trans::NoTranspose>;
  using algo_tag_type  = typename Algo::Gbtrs::Unblocked;

  test_batched_gbtrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_gbtrs_t_fcomplex) {
  using param_tag_type = ::Test::Gbtrs::ParamTag<Trans::Transpose>;
  using algo_tag_type  = typename Algo::Gbtrs::Unblocked;

  test_batched_gbtrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_gbtrs_c_fcomplex) {
  using param_tag_type = ::Test::Gbtrs::ParamTag<Trans::ConjTranspose>;
  using algo_tag_type  = typename Algo::Gbtrs::Unblocked;

  test_batched_gbtrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_gbtrs_nt_dcomplex) {
  using param_tag_type = ::Test::Gbtrs::ParamTag<Trans::NoTranspose>;
  using algo_tag_type  = typename Algo::Gbtrs::Unblocked;

  test_batched_gbtrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_gbtrs_t_dcomplex) {
  using param_tag_type = ::Test::Gbtrs::ParamTag<Trans::Transpose>;
  using algo_tag_type  = typename Algo::Gbtrs::Unblocked;

  test_batched_gbtrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_gbtrs_c_dcomplex) {
  using param_tag_type = ::Test::Gbtrs::ParamTag<Trans::ConjTranspose>;
  using algo_tag_type  = typename Algo::Gbtrs::Unblocked;

  test_batched_gbtrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
#endif
