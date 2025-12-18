// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Getrf.hpp>
#include <KokkosBatched_Getrs.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Getrs {

template <typename T>
struct ParamTag {
  using trans = T;
};

template <typename DeviceType, typename AViewType, typename PivViewType>
struct Functor_BatchedSerialGetrf {
  using execution_space = typename DeviceType::execution_space;
  AViewType m_a;
  PivViewType m_ipiv;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGetrf(const AViewType &a, const PivViewType &ipiv) : m_a(a), m_ipiv(ipiv) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto aa   = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto ipiv = Kokkos::subview(m_ipiv, k, Kokkos::ALL());

    KokkosBatched::SerialGetrf<Algo::Getrf::Unblocked>::invoke(aa, ipiv);
  }

  inline void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGetrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename AViewType, typename PivViewType, typename BViewType, typename ParamTagType,
          typename AlgoTagType>
struct Functor_BatchedSerialGetrs {
  using execution_space = typename DeviceType::execution_space;
  AViewType m_a;
  BViewType m_b;
  PivViewType m_ipiv;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGetrs(const AViewType &a, const PivViewType &ipiv, const BViewType &b)
      : m_a(a), m_b(b), m_ipiv(ipiv) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k, int &info) const {
    auto aa   = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto ipiv = Kokkos::subview(m_ipiv, k, Kokkos::ALL());
    auto bb   = Kokkos::subview(m_b, k, Kokkos::ALL());

    info += KokkosBatched::SerialGetrs<typename ParamTagType::trans, AlgoTagType>::invoke(aa, ipiv, bb);
  }

  inline int run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGetrs");
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
    std::string name_region("KokkosBatched::Test::SerialGetrs");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, m_x.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

/// \brief Implementation details of batched getrs test
/// Confirm A * x = b, where
///        A: [[1, 1],
///            [1, -1]]
///        b: [2, 0]
///        x: [1, 1]
/// This corresponds to the following system of equations:
///        x0 + x1 = 2
///        x0 - x1 = 0
/// We confirm this with the factorized matrix LU and pivot given by
///       LU: [[1, 1],
///            [1, -2]]
///       piv: [0, 1]
///
/// \param N [in] Batch size of RHS (banded matrix can also be batched matrix)
/// \param k [in] Number of superdiagonals or subdiagonals of matrix A
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_getrs_analytical(const int N) {
  using ats           = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType      = typename ats::mag_type;
  using View2DType    = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType    = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using PivView2DType = Kokkos::View<int **, LayoutType, DeviceType>;

  const int BlkSize = 2;
  View3DType LU("LU", N, BlkSize, BlkSize);                   // Factorized matrix of A
  View2DType x("x", N, BlkSize), x_ref("x_ref", N, BlkSize);  // Solutions
  PivView2DType ipiv("ipiv", N, BlkSize);

  auto h_LU    = Kokkos::create_mirror_view(LU);
  auto h_ipiv  = Kokkos::create_mirror_view(ipiv);
  auto h_x     = Kokkos::create_mirror_view(x);
  auto h_x_ref = Kokkos::create_mirror_view(x_ref);
  Kokkos::deep_copy(h_LU, 1.0);
  for (int ib = 0; ib < N; ib++) {
    h_LU(ib, 1, 1) = -2.0;
    h_ipiv(ib, 0)  = 0;
    h_ipiv(ib, 1)  = 1;

    h_x(ib, 0)     = 2;
    h_x(ib, 1)     = 0;
    h_x_ref(ib, 0) = 1;
    h_x_ref(ib, 1) = 1;
  }

  Kokkos::deep_copy(LU, h_LU);
  Kokkos::deep_copy(ipiv, h_ipiv);
  Kokkos::deep_copy(x, h_x);

  // getrs (Note, LU is a factorized matrix of A)
  auto info = Functor_BatchedSerialGetrs<DeviceType, View3DType, PivView2DType, View2DType, ParamTagType, AlgoTagType>(
                  LU, ipiv, x)
                  .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  // Check if x = [1, 1]
  Kokkos::deep_copy(h_x, x);
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_x(ib, i), h_x_ref(ib, i), eps);
    }
  }
}

/// \brief Implementation details of batched getrs test
///
/// \param N [in] Batch size of RHS (banded matrix can also be batched matrix)
/// \param k [in] Number of superdiagonals or subdiagonals of matrix A
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename ParamTagType, typename AlgoTagType>
void impl_test_batched_getrs(const int N, const int BlkSize) {
  using ats           = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType      = typename ats::mag_type;
  using View2DType    = Kokkos::View<ScalarType **, LayoutType, DeviceType>;
  using View3DType    = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using PivView2DType = Kokkos::View<int **, LayoutType, DeviceType>;

  View3DType A("A", N, BlkSize, BlkSize), ref("Ref", N, BlkSize, BlkSize);
  View3DType LU("LU", N, BlkSize, BlkSize);                               // Factorized
  View2DType x("x", N, BlkSize), y("y", N, BlkSize), b("b", N, BlkSize);  // Solutions
  PivView2DType ipiv("ipiv", N, BlkSize);

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  // Initialize A_reconst with random matrix
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(x, rand_pool, randStart, randEnd);
  Kokkos::deep_copy(LU, A);
  Kokkos::deep_copy(b, x);

  // getrf to factorize matrix A = P * L * U
  Functor_BatchedSerialGetrf<DeviceType, View3DType, PivView2DType>(LU, ipiv).run();

  // getrs (Note, LU is a factorized matrix of A)
  auto info = Functor_BatchedSerialGetrs<DeviceType, View3DType, PivView2DType, View2DType, ParamTagType, AlgoTagType>(
                  LU, ipiv, x)
                  .run();

  Kokkos::fence();
  EXPECT_EQ(info, 0);

  // Gemv to compute A*x, this should be identical to b
  Functor_BatchedSerialGemv<DeviceType, ScalarType, View3DType, View2DType, View2DType, ParamTagType>(1.0, A, x, 0.0, y)
      .run();

  Kokkos::fence();

  // this eps is about 10^-14
  RealType eps = 1.0e3 * ats::epsilon();

  // Check if A * x = b
  auto h_y = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  auto h_b = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b);
  for (int ib = 0; ib < N; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_NEAR_KK(h_y(ib, i), h_b(ib, i), eps);
    }
  }
}

}  // namespace Getrs
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_getrs() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Getrs::impl_test_batched_getrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Getrs::impl_test_batched_getrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Getrs::impl_test_batched_getrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, i);
      Test::Getrs::impl_test_batched_getrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Getrs::impl_test_batched_getrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1);
    Test::Getrs::impl_test_batched_getrs_analytical<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Getrs::impl_test_batched_getrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(1, i);
      Test::Getrs::impl_test_batched_getrs<DeviceType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_getrs_nt_float) {
  using param_tag_type = ::Test::Getrs::ParamTag<Trans::NoTranspose>;
  using algo_tag_type  = typename Algo::Getrs::Unblocked;

  test_batched_getrs<TestDevice, float, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_getrs_t_float) {
  using param_tag_type = ::Test::Getrs::ParamTag<Trans::Transpose>;
  using algo_tag_type  = typename Algo::Getrs::Unblocked;

  test_batched_getrs<TestDevice, float, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_getrs_nt_double) {
  using param_tag_type = ::Test::Getrs::ParamTag<Trans::NoTranspose>;
  using algo_tag_type  = typename Algo::Getrs::Unblocked;

  test_batched_getrs<TestDevice, double, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_getrs_t_double) {
  using param_tag_type = ::Test::Getrs::ParamTag<Trans::Transpose>;
  using algo_tag_type  = typename Algo::Getrs::Unblocked;

  test_batched_getrs<TestDevice, double, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_getrs_nt_fcomplex) {
  using param_tag_type = ::Test::Getrs::ParamTag<Trans::NoTranspose>;
  using algo_tag_type  = typename Algo::Getrs::Unblocked;

  test_batched_getrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_getrs_t_fcomplex) {
  using param_tag_type = ::Test::Getrs::ParamTag<Trans::Transpose>;
  using algo_tag_type  = typename Algo::Getrs::Unblocked;

  test_batched_getrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_getrs_c_fcomplex) {
  using param_tag_type = ::Test::Getrs::ParamTag<Trans::ConjTranspose>;
  using algo_tag_type  = typename Algo::Getrs::Unblocked;

  test_batched_getrs<TestDevice, Kokkos::complex<float>, param_tag_type, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_getrs_nt_dcomplex) {
  using param_tag_type = ::Test::Getrs::ParamTag<Trans::NoTranspose>;
  using algo_tag_type  = typename Algo::Getrs::Unblocked;

  test_batched_getrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_getrs_t_dcomplex) {
  using param_tag_type = ::Test::Getrs::ParamTag<Trans::Transpose>;
  using algo_tag_type  = typename Algo::Getrs::Unblocked;

  test_batched_getrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
TEST_F(TestCategory, test_batched_getrs_c_dcomplex) {
  using param_tag_type = ::Test::Getrs::ParamTag<Trans::ConjTranspose>;
  using algo_tag_type  = typename Algo::Getrs::Unblocked;

  test_batched_getrs<TestDevice, Kokkos::complex<double>, param_tag_type, algo_tag_type>();
}
#endif
