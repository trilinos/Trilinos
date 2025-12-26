// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_TestVanilla.hpp"

namespace Test {
namespace Gemm {

template <typename TA, typename TB>
struct ParamTag {
  using transA = TA;
  using transB = TB;
};

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
struct Functor_TestBatchedSerialGemm {
  using execution_space = typename DeviceType::execution_space;
  ViewType m_a, m_b, m_c;
  ScalarType m_alpha, m_beta;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedSerialGemm(const ScalarType alpha, const ViewType &a, const ViewType &b, const ScalarType beta,
                                const ViewType &c)
      : m_a(a), m_b(b), m_c(c), m_alpha(alpha), m_beta(beta) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k, int &info) const {
    auto aa = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(m_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(m_c, k, Kokkos::ALL(), Kokkos::ALL());

    info +=
        KokkosBatched::SerialGemm<typename ParamTagType::transA, typename ParamTagType::transB, AlgoTagType>::invoke(
            m_alpha, aa, bb, m_beta, cc);
  }

  inline int run() {
    using value_type = typename ViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGemm");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, m_c.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

/// \brief Implementation details of batched gemm test
/// \param N [in] Batch size of matrices
/// \param matAdim1 [in] Number of rows of matrix A
/// \param matAdim2 [in] Number of columns of matrix A
/// \param matBdim1 [in] Number of rows of matrix B
/// \param matBdim2 [in] Number of columns of matrix B
/// \param matCdim1 [in] Number of rows of matrix C
/// \param matCdim2 [in] Number of columns of matrix C
template <typename DeviceType, typename ValueType, typename ScalarType, typename LayoutType, typename ParamTagType,
          typename AlgoTagType>
void impl_test_batched_gemm(const int N, const int matAdim1, const int matAdim2, const int matBdim1, const int matBdim2,
                            const int matCdim1, const int matCdim2) {
  using execution_space = typename DeviceType::execution_space;
  using transA          = typename ParamTagType::transA;
  using transB          = typename ParamTagType::transB;
  using ats             = KokkosKernels::ArithTraits<ValueType>;
  using ViewType        = Kokkos::View<ValueType ***, LayoutType, DeviceType>;

  /// randomized input testing views
  ScalarType alpha = ScalarType(1.5);
  ScalarType beta  = ScalarType(3.0);

  ViewType A("A", N, matAdim1, matAdim2), B("B", N, matBdim1, matBdim2), C("C", N, matCdim1, matCdim2),
      C_ref("C_ref", N, matCdim1, matCdim2);

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  ScalarType randStart, randEnd;
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(B, rand_pool, randStart, randEnd);
  Kokkos::fill_random(C, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(C_ref, C);

  Functor_BatchedVanillaGEMM<ViewType, ViewType, ViewType, execution_space> vgemm;
  vgemm.A_t   = !std::is_same_v<transA, KokkosBatched::Trans::NoTranspose>;
  vgemm.B_t   = !std::is_same_v<transB, KokkosBatched::Trans::NoTranspose>;
  vgemm.A_c   = std::is_same_v<transA, KokkosBatched::Trans::ConjTranspose>;
  vgemm.B_c   = std::is_same_v<transB, KokkosBatched::Trans::ConjTranspose>;
  vgemm.A_    = A;
  vgemm.B_    = B;
  vgemm.C_    = C_ref;
  vgemm.alpha = alpha;
  vgemm.beta  = beta;
  vgemm.run();  // Compute C_ref

  // Compute using gemm API
  auto info =
      Functor_TestBatchedSerialGemm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(alpha, A, B, beta, C)
          .run();
  EXPECT_EQ(info, 0);

  auto h_C     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), C);
  auto h_C_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), C_ref);

  // check C = C_ref
  using mag_type = typename ats::mag_type;
  mag_type sum(1), diff(0);

  mag_type eps = ats::epsilon();
  eps *= std::is_same_v<ValueType, Kokkos::Experimental::half_t> ||
                 std::is_same_v<ValueType, Kokkos::Experimental::bhalf_t>
             ? 4
             : 1e3;

  for (int k = 0; k < N; ++k)
    for (int i = 0; i < matCdim1; ++i)
      for (int j = 0; j < matCdim2; ++j) {
        sum += ats::abs(h_C_ref(k, i, j));
        diff += ats::abs(h_C_ref(k, i, j) - h_C(k, i, j));
      }
  EXPECT_NEAR_KK(diff / sum, 0, eps);
}
}  // namespace Gemm
}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_gemm() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
        0, 10, 10, 10, 10, 10, 10);
    for (int i = 0; i < 10; ++i) {
      Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
          1024, i, i, i, i, i, i);
    }
    for (int i = 0; i < 10; ++i) {
      int dimM = i;
      int dimN = 2 * i;
      int dimK = 3 * i;
      if ((std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
            1024, dimM, dimK, dimK, dimN, dimM, dimN);
      }
      if ((std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              !std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
            1024, dimM, dimK, dimN, dimK, dimM, dimN);
      }
      if ((!std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
            1024, dimK, dimM, dimK, dimN, dimM, dimN);
      }
      if ((!std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              !std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
            1024, dimK, dimM, dimN, dimK, dimM, dimN);
      }
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
        0, 10, 10, 10, 10, 10, 10);
    for (int i = 0; i < 10; ++i) {
      Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
          1024, i, i, i, i, i, i);
    }
    for (int i = 0; i < 10; ++i) {
      int dimM = i;
      int dimN = 2 * i;
      int dimK = 3 * i;
      if ((std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
            1024, dimM, dimK, dimK, dimN, dimM, dimN);
      }
      if ((std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              !std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
            1024, dimM, dimK, dimN, dimK, dimM, dimN);
      }
      if ((!std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
            1024, dimK, dimM, dimK, dimN, dimM, dimN);
      }
      if ((!std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              !std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::Gemm::impl_test_batched_gemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType, AlgoTagType>(
            1024, dimK, dimM, dimN, dimK, dimM, dimN);
      }
    }
  }
#endif

  return 0;
}
