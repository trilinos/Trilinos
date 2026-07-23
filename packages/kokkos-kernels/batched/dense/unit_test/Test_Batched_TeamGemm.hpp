// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Gemm_Decl.hpp"

#include "KokkosBatched_Gemm_Team_Impl.hpp"
#include "KokkosBatched_Gemm_TeamVector_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace TeamGemm {

template <typename Mode, typename TA, typename TB>
struct ParamTag {
  using mode   = Mode;
  using transA = TA;
  using transB = TB;
};

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
struct Functor_TestBatchedTeamGemm {
  using execution_space = typename DeviceType::execution_space;
  using member_type     = typename Kokkos::TeamPolicy<execution_space>::member_type;

  ViewType m_a, m_b, m_c;
  ScalarType m_alpha, m_beta;

  Functor_TestBatchedTeamGemm(const ScalarType alpha, const ViewType &a, const ViewType &b, const ScalarType beta,
                              const ViewType &c)
      : m_a(a), m_b(b), m_c(c), m_alpha(alpha), m_beta(beta) {}

  KOKKOS_INLINE_FUNCTION void operator()(const ParamTagType &, const member_type &member) const {
    const int k = member.league_rank();

    auto aa = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(m_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(m_c, k, Kokkos::ALL(), Kokkos::ALL());

    if constexpr (std::is_same_v<typename ParamTagType::mode, KokkosBatched::Mode::Team>) {
      KokkosBatched::TeamGemm<member_type, typename ParamTagType::transA, typename ParamTagType::transB,
                              AlgoTagType>::invoke(member, m_alpha, aa, bb, m_beta, cc);
    } else if constexpr (std::is_same_v<typename ParamTagType::mode, KokkosBatched::Mode::TeamVector>) {
      KokkosBatched::TeamVectorGemm<member_type, typename ParamTagType::transA, typename ParamTagType::transB,
                                    AlgoTagType>::invoke(member, m_alpha, aa, bb, m_beta, cc);
    }
  }

  inline void run() {
    using value_type                  = typename ViewType::non_const_value_type;
    std::string name_region           = std::is_same_v<typename ParamTagType::mode, KokkosBatched::Mode::Team>
                                            ? "KokkosBatched::Test::TeamGemm"
                                            : "KokkosBatched::Test::TeamVectorGemm";
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    const int league_size = m_c.extent(0);
    Kokkos::TeamPolicy<execution_space, ParamTagType> policy(league_size, Kokkos::AUTO);
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
  }
};

/// \brief Compute the error bound for the GEMM operation
/// C = alpha * op(A) * op(B) + beta * C
/// The error bound with machine precison eps is computed as follows:
/// 1. Error analysis of Matrix Multiplication D = op(A) * op(B)
///    delta [D (i, j)] <= gamma_k * E (i, j)
///    where E (i, j) = \sum_{k=1}^{ka} (A' (i, k) * B' (k, j)),
///    A' = op(A), B' = op(B), gamma_k = ka * eps / (1 - ka * eps)
/// 2. Error from Scaling by alpha
///    delta [alpha * D (i, j)] <= gamma_k * |alpha| * E (i, j) + eps * |alpha * D (i, j)|
/// 3. Error from Scaling by beta
///    delta [beta * C (i, j)] <= eps * |beta * C (i, j)|
/// 4. Error from the addition operation
///    delta [beta * C (i, j) + alpha * D (i, j)]
///    <= eps * |beta * C (i, j) + alpha * D (i, j)|
///    <= eps * (|beta * C (i, j)| + |alpha * D (i, j)|)
/// 5. Final error bound
///    delta [C (i, j)] <= gamma_k * |alpha| * E (i, j) + 2 * eps * (|beta * C (i, j)| + |alpha * D (i, j)|)
///
/// \tparam ArgTransA Transpose option for matrix A
/// \tparam ArgTransB Transpose option for matrix B
/// \tparam AViewType View type for matrix A
/// \tparam BViewType View type for matrix B
/// \tparam CViewType View type for matrix C
/// \tparam ErrorBoundViewType [out] View type for the error bound
/// \tparam ScalarType Type for scalar alpha and beta

/// \param A [in] Input matrix A
/// \param B [in] Input matrix B
/// \param C [in] Input matrix C
/// \param error_bound [out] Output matrix for the error bound
/// \param alpha [in] Scalar multiplier for matrix A
/// \param beta [in] Scalar multiplier for matrix C
/// \param eps [in] Machine precision
template <typename ArgTransA, typename ArgTransB, typename AViewType, typename BViewType, typename CViewType,
          typename ErrorBoundViewType, typename ScalarType>
void gemm_error_bound(const AViewType &A, const BViewType &B, const CViewType &C, const ErrorBoundViewType &error_bound,
                      ScalarType alpha, ScalarType beta, double eps) {
  using value_type = typename AViewType::non_const_value_type;
  using ats        = KokkosKernels::ArithTraits<value_type>;
  using mag_type   = typename ats::mag_type;
  const int nb     = C.extent_int(0);
  const int m      = C.extent_int(1);
  const int n      = C.extent_int(2);
  const int ka     = std::is_same_v<ArgTransA, Trans::NoTranspose> ? A.extent(2) : A.extent(1);

  auto h_A           = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  auto h_B           = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B);
  auto h_C           = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), C);
  auto h_error_bound = Kokkos::create_mirror_view(error_bound);

  for (int ib = 0; ib < nb; ib++) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        mag_type sum = 0;
        if (std::is_same_v<ArgTransA, Trans::NoTranspose>) {
          if (std::is_same_v<ArgTransB, Trans::NoTranspose>) {
            for (int ii = 0; ii < A.extent_int(2); ii++) {
              sum += ats::abs(alpha * h_A(ib, i, ii) * h_B(ib, ii, j));
            }
          } else {
            for (int ii = 0; ii < A.extent_int(2); ii++) {
              sum += ats::abs(alpha * h_A(ib, i, ii) * h_B(ib, j, ii));
            }
          }
        } else {
          if (std::is_same_v<ArgTransB, Trans::NoTranspose>) {
            for (int ii = 0; ii < A.extent_int(1); ii++) {
              sum += ats::abs(alpha * h_A(ib, ii, i) * h_B(ib, ii, j));
            }
          } else {
            for (int ii = 0; ii < A.extent_int(1); ii++) {
              sum += ats::abs(alpha * h_A(ib, ii, i) * h_B(ib, j, ii));
            }
          }
        }
        mag_type gamma_k        = static_cast<mag_type>(ka) * eps / (1.0 - static_cast<mag_type>(ka) * eps);
        h_error_bound(ib, i, j) = gamma_k * ats::abs(alpha) * ats::abs(sum) +
                                  2 * eps * (ats::abs(beta * h_C(ib, i, j)) + ats::abs(alpha * sum));
      }
    }
  }
  Kokkos::deep_copy(error_bound, h_error_bound);
}

/// \brief Implementation details of batched team gemm test
/// \param[in] N: Batch size of matrices
/// \param[in] matAdim1: Number of rows of matrix A
/// \param[in] matAdim2: Number of columns of matrix A
/// \param[in] matBdim1: Number of rows of matrix B
/// \param[in] matBdim2: Number of columns of matrix B
/// \param[in] matCdim1: Number of rows of matrix C
/// \param[in] matCdim2: Number of columns of matrix C
template <typename DeviceType, typename ValueType, typename ScalarType, typename LayoutType, typename ParamTagType,
          typename AlgoTagType>
void impl_test_batched_teamgemm(const int N, const int matAdim1, const int matAdim2, const int matBdim1,
                                const int matBdim2, const int matCdim1, const int matCdim2) {
  using execution_space = typename DeviceType::execution_space;
  using transA          = typename ParamTagType::transA;
  using transB          = typename ParamTagType::transB;
  using execution_space = typename DeviceType::execution_space;
  using ats             = KokkosKernels::ArithTraits<ValueType>;
  using ViewType        = Kokkos::View<ValueType ***, LayoutType, DeviceType>;
  using FP64Type =
      std::conditional_t<KokkosKernels::ArithTraits<ValueType>::is_complex, Kokkos::complex<double>, double>;
  using ViewFP64Type      = Kokkos::View<FP64Type ***, LayoutType, DeviceType>;
  using ErrorViewFP64Type = Kokkos::View<double ***, LayoutType, DeviceType>;

  /// randomized input testing views
  ScalarType alpha = ScalarType(1.5), beta = ScalarType(3.0);
  FP64Type alpha_fp64 = FP64Type(1.5), beta_fp64 = FP64Type(3.0);

  ViewType A("A", N, matAdim1, matAdim2), B("B", N, matBdim1, matBdim2), C("C", N, matCdim1, matCdim2);
  ViewFP64Type A_fp64("A_fp64", N, matAdim1, matAdim2), B_fp64("B_fp64", N, matBdim1, matBdim2),
      C_fp64("C_fp64", N, matCdim1, matCdim2);
  ErrorViewFP64Type ErrorBound("ErrorBound", N, matCdim1, matCdim2);

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  // Data in given precision
  ScalarType randStart, randEnd;
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(B, rand_pool, randStart, randEnd);
  Kokkos::fill_random(C, rand_pool, randStart, randEnd);

  Kokkos::deep_copy(A_fp64, A);
  Kokkos::deep_copy(B_fp64, B);
  Kokkos::deep_copy(C_fp64, C);

  Functor_BatchedVanillaGEMM<ViewFP64Type, ViewFP64Type, ViewFP64Type, execution_space> vgemm;
  vgemm.A_t   = !std::is_same<transA, KokkosBatched::Trans::NoTranspose>::value;
  vgemm.B_t   = !std::is_same<transB, KokkosBatched::Trans::NoTranspose>::value;
  vgemm.A_c   = std::is_same_v<transA, KokkosBatched::Trans::ConjTranspose>;
  vgemm.B_c   = std::is_same_v<transB, KokkosBatched::Trans::ConjTranspose>;
  vgemm.A_    = A_fp64;
  vgemm.B_    = B_fp64;
  vgemm.C_    = C_fp64;
  vgemm.alpha = alpha_fp64;
  vgemm.beta  = beta_fp64;
  vgemm.run();  // Compute C_fp64 (reference)

  // Compute using gemm API in given precision
  Functor_TestBatchedTeamGemm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(alpha, A, B, beta, C).run();

  // Machine precision for the value type
  double eps = static_cast<double>(ats::epsilon());
  gemm_error_bound<transA, transB>(A_fp64, B_fp64, C_fp64, ErrorBound, alpha_fp64, beta_fp64, eps);

  auto h_C          = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), C);
  auto h_C_fp64     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), C_fp64);
  auto h_ErrorBound = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ErrorBound);
  for (int k = 0; k < N; ++k) {
    for (int i = 0; i < matCdim1; ++i) {
      for (int j = 0; j < matCdim2; ++j) {
        double actual_err  = Kokkos::abs(h_C(k, i, j) - h_C_fp64(k, i, j));
        double error_bound = h_ErrorBound(k, i, j);
        EXPECT_NEAR_KK(actual_err, 0, error_bound);
      }
    }
  }
}
}  // namespace TeamGemm
}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_teamgemm() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                               AlgoTagType>(0, 10, 10, 10, 10, 10, 10);
    for (int i = 0; i < 10; ++i) {
      Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                                 AlgoTagType>(128, i, i, i, i, i, i);
    }
    for (int i = 0; i < 10; ++i) {
      int dimM = 3 * i;
      int dimN = 2 * i;
      int dimK = i;
      if ((std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                                   AlgoTagType>(128, dimM, dimK, dimK, dimN, dimM, dimN);
      }
      if ((std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              !std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                                   AlgoTagType>(128, dimM, dimK, dimN, dimK, dimM, dimN);
      }
      if ((!std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                                   AlgoTagType>(128, dimK, dimM, dimK, dimN, dimM, dimN);
      }
      if ((!std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              !std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                                   AlgoTagType>(128, dimK, dimM, dimN, dimK, dimM, dimN);
      }
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                               AlgoTagType>(0, 10, 10, 10, 10, 10, 10);
    for (int i = 0; i < 10; ++i) {
      Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                                 AlgoTagType>(128, i, i, i, i, i, i);
    }
    for (int i = 0; i < 10; ++i) {
      int dimM = 3 * i;
      int dimN = 2 * i;
      int dimK = i;
      if ((std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                                   AlgoTagType>(128, dimM, dimK, dimK, dimN, dimM, dimN);
      }
      if ((std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              !std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                                   AlgoTagType>(128, dimM, dimK, dimN, dimK, dimM, dimN);
      }
      if ((!std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                                   AlgoTagType>(128, dimK, dimM, dimK, dimN, dimM, dimN);
      }
      if ((!std::is_same_v<typename ParamTagType::transA, KokkosBatched::Trans::NoTranspose>)&&(
              !std::is_same_v<typename ParamTagType::transB, KokkosBatched::Trans::NoTranspose>)) {
        Test::TeamGemm::impl_test_batched_teamgemm<DeviceType, ValueType, ScalarType, LayoutType, ParamTagType,
                                                   AlgoTagType>(128, dimK, dimM, dimN, dimK, dimM, dimN);
      }
    }
  }
#endif

  return 0;
}
