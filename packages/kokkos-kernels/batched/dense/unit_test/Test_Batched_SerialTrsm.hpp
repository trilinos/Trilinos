// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Trsm_Serial_Impl.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "Test_Batched_DenseUtils.hpp"

using namespace KokkosBatched;

namespace Test {
namespace Trsm {

template <typename S, typename U, typename T, typename D>
struct ParamTag {
  typedef S side;
  typedef U uplo;
  typedef T trans;
  typedef D diag;
};

template <typename DeviceType, typename ViewType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
struct Functor_TestBatchedSerialTrsm {
  using execution_space = typename DeviceType::execution_space;
  ViewType m_a, m_b;

  ScalarType m_alpha;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedSerialTrsm(const ScalarType alpha, const ViewType &a, const ViewType &b)
      : m_a(a), m_b(b), m_alpha(alpha) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const ParamTagType &, const int k) const {
    auto aa = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(m_b, k, Kokkos::ALL(), Kokkos::ALL());

    SerialTrsm<typename ParamTagType::side, typename ParamTagType::uplo, typename ParamTagType::trans,
               typename ParamTagType::diag, AlgoTagType>::invoke(m_alpha, aa, bb);
  }

  inline void run() {
    typedef typename ViewType::value_type value_type;
    std::string name_region("KokkosBatched::Test::SerialTrsm");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space, ParamTagType> policy(0, m_b.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
    Kokkos::Profiling::popRegion();
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
    std::string name_region("KokkosBatched::Test::SerialTrsm");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename ScalarType, typename ValueType, typename LayoutType, typename ParamTagType,
          typename AlgoTagType>
void impl_test_batched_trsm_blocking(const int N, const int BlkSize, const int NumCols) {
  using ats      = KokkosKernels::ArithTraits<ValueType>;
  using ViewType = Kokkos::View<ValueType ***, LayoutType, DeviceType>;

  /// randomized input testing views
  ScalarType alpha(1.0);

  const bool is_side_right = std::is_same_v<typename ParamTagType::side, Side::Right>;
  const int b_nrows        = is_side_right ? NumCols : BlkSize;
  const int b_ncols        = is_side_right ? BlkSize : NumCols;
  ViewType a0("a0", N, BlkSize, BlkSize), a1("a1", N, BlkSize, BlkSize), b0("b0", N, b_nrows, b_ncols),
      b1("b1", N, b_nrows, b_ncols);

  Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
  Kokkos::fill_random(a0, random, ValueType(1.0));
  Kokkos::fill_random(b0, random, ValueType(1.0));

  Kokkos::fence();

  Kokkos::deep_copy(a1, a0);
  Kokkos::deep_copy(b1, b0);

  Functor_TestBatchedSerialTrsm<DeviceType, ViewType, ScalarType, ParamTagType, Algo::Trsm::Blocked>(alpha, a0, b0)
      .run();
  Functor_TestBatchedSerialTrsm<DeviceType, ViewType, ScalarType, ParamTagType, AlgoTagType>(alpha, a1, b1).run();

  Kokkos::fence();

  /// for comparison send it to host
  auto b0_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b0);
  auto b1_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b1);

  /// check b0 = b1 ; this eps is about 10^-14
  using mag_type = typename ats::mag_type;
  mag_type sum(1), diff(0);
  const mag_type eps = 1.0e3 * ats::epsilon();

  for (int k = 0; k < N; ++k)
    for (int i = 0; i < b_nrows; ++i)
      for (int j = 0; j < b_ncols; ++j) {
        sum += ats::abs(b0_host(k, i, j));
        diff += ats::abs(b0_host(k, i, j) - b1_host(k, i, j));
      }
  EXPECT_NEAR_KK(diff / sum, 0.0, eps);
}

/// \brief Implementation details of batched trsm analytical test
/// Confirm A * x = b, where
/// A: [[1, 1],
///     [1.5, 2]]
/// X: [[1, 1],
///     [1, 1]]
///
/// Upper & Non-Transpose & Left
/// 1 x00 + 1 x10 = 1.5
/// 1 x01 + 1 x11 = 1.5
///         2 x10 = 1.5
///         2 x11 = 1.5
/// x = [[3/4, 3/4], [3/4, 3/4]]
///
/// Upper & Transpose & Left
/// 1 x00         = 1.5
/// 1 x01         = 1.5
/// 1 x00 + 2 x10 = 1.5
/// 1 x01 + 2 x11 = 1.5
/// x = [[3/2, 3/2], [0, 0]]
///
/// Upper & Non-Transpose & Right
/// 1 x00         = 1.5
/// 1 x10         = 1.5
/// 1 x00 + 2 x01 = 1.5
/// 1 x10 + 2 x11 = 1.5
/// x = [[3/2, 0], [3/2, 0]]
///
/// Upper & Transpose & Right
/// 1 x00 + 1 x01 = 1.5
/// 1 x10 + 1 x11 = 1.5
///         2 x01 = 1.5
///         2 x11 = 1.5
/// x = [[3/4, 3/4], [3/4, 3/4]]
///
/// Lower & Non-Transpose & Left
/// 1 x00           = 1.5
/// 3/2 x00 + 2 x10 = 1.5
/// 1 x01           = 1.5
/// 3/2 x01 + 2 x11 = 1.5
/// x = [[3/2, 3/2], [-3/8, -3/8]]
///
/// Lower & Transpose & Left
/// 1 x00 + 3/2 x10 = 1.5
///           2 x10 = 1.5
/// 1 x01 + 3/2 x11 = 1.5
///           2 x11 = 1.5
/// x = [[3/8, 3/8], [3/4, 3/4]]
///
/// Lower & Non-Transpose & Right
/// 1 x00 + 3/2 x10 = 1.5
///           2 x01 = 1.5
/// 1 x10 + 3/2 x11 = 1.5
///           2 x11 = 1.5
/// x = [[3/8, 3/4], [3/8, 3/4]]
///
/// Lower & Transpose & Right
/// 1 x00           = 1.5
/// 3/2 x00 + 2 x01 = 1.5
/// 1 x10           = 1.5
/// 3/2 x10 + 2 x11 = 1.5
/// x = [[3/2, -3/8], [3/2, -3/8]]
///
/// \param N [in] Batch size of matrices and RHS
template <typename DeviceType, typename ScalarType, typename ValueType, typename LayoutType, typename ParamTagType,
          typename AlgoTagType>
void impl_test_batched_trsm_analytical(const std::size_t N) {
  using ats        = typename KokkosKernels::ArithTraits<ValueType>;
  using RealType   = typename ats::mag_type;
  using View3DType = Kokkos::View<ValueType ***, LayoutType, DeviceType>;

  const std::size_t m = 2;
  View3DType A("A", N, m, m), B("B", N, m, m), B_ref("B_ref", N, m, m);

  auto h_A     = Kokkos::create_mirror_view(A);
  auto h_B_ref = Kokkos::create_mirror_view(B_ref);

  for (std::size_t ib = 0; ib < N; ib++) {
    h_A(ib, 0, 0) = 1.0;
    h_A(ib, 0, 1) = 1.0;
    h_A(ib, 1, 0) = 1.5;
    h_A(ib, 1, 1) = 2.0;
    if (std::is_same_v<typename ParamTagType::side, KokkosBatched::Side::Left>) {
      // Left
      if (std::is_same_v<typename ParamTagType::uplo, KokkosBatched::Uplo::Upper>) {
        // Upper
        if (std::is_same_v<typename ParamTagType::trans, Trans::NoTranspose>) {
          // No-Transpose
          if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
            // Non-Unit
            h_B_ref(ib, 0, 0) = 3.0 / 4.0;
            h_B_ref(ib, 1, 0) = 3.0 / 4.0;
            h_B_ref(ib, 0, 1) = 3.0 / 4.0;
            h_B_ref(ib, 1, 1) = 3.0 / 4.0;
          } else {
            // Unit
            h_B_ref(ib, 0, 0) = 0.0;
            h_B_ref(ib, 1, 0) = 3.0 / 2.0;
            h_B_ref(ib, 0, 1) = 0.0;
            h_B_ref(ib, 1, 1) = 3.0 / 2.0;
          }
        } else {
          // Transpose/ConjTrans
          if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
            // Non-Unit
            h_B_ref(ib, 0, 0) = 3.0 / 2.0;
            h_B_ref(ib, 1, 0) = 0.0;
            h_B_ref(ib, 0, 1) = 3.0 / 2.0;
            h_B_ref(ib, 1, 1) = 0.0;
          } else {
            // Unit
            h_B_ref(ib, 0, 0) = 3.0 / 2.0;
            h_B_ref(ib, 1, 0) = 0.0;
            h_B_ref(ib, 0, 1) = 3.0 / 2.0;
            h_B_ref(ib, 1, 1) = 0.0;
          }
        }
      } else {
        // Lower
        if (std::is_same_v<typename ParamTagType::trans, Trans::NoTranspose>) {
          // No-Transpose
          if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
            // Non-Unit
            h_B_ref(ib, 0, 0) = 3.0 / 2.0;
            h_B_ref(ib, 1, 0) = -3.0 / 8.0;
            h_B_ref(ib, 0, 1) = 3.0 / 2.0;
            h_B_ref(ib, 1, 1) = -3.0 / 8.0;
          } else {
            // Unit
            h_B_ref(ib, 0, 0) = 3.0 / 2.0;
            h_B_ref(ib, 1, 0) = -3.0 / 4.0;
            h_B_ref(ib, 0, 1) = 3.0 / 2.0;
            h_B_ref(ib, 1, 1) = -3.0 / 4.0;
          }
        } else {
          // Transpose/ConjTrans
          if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
            // Non-Unit
            h_B_ref(ib, 0, 0) = 3.0 / 8.0;
            h_B_ref(ib, 1, 0) = 3.0 / 4.0;
            h_B_ref(ib, 0, 1) = 3.0 / 8.0;
            h_B_ref(ib, 1, 1) = 3.0 / 4.0;
          } else {
            // Unit
            h_B_ref(ib, 0, 0) = -3.0 / 4.0;
            h_B_ref(ib, 1, 0) = 3.0 / 2.0;
            h_B_ref(ib, 0, 1) = -3.0 / 4.0;
            h_B_ref(ib, 1, 1) = 3.0 / 2.0;
          }
        }
      }
    } else {
      // Right
      if (std::is_same_v<typename ParamTagType::uplo, KokkosBatched::Uplo::Upper>) {
        // Upper
        if (std::is_same_v<typename ParamTagType::trans, Trans::NoTranspose>) {
          // No-Transpose
          if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
            // Non-Unit
            h_B_ref(ib, 0, 0) = 3.0 / 2.0;
            h_B_ref(ib, 1, 0) = 3.0 / 2.0;
            h_B_ref(ib, 0, 1) = 0.0;
            h_B_ref(ib, 1, 1) = 0.0;
          } else {
            // Unit
            h_B_ref(ib, 0, 0) = 3.0 / 2.0;
            h_B_ref(ib, 1, 0) = 3.0 / 2.0;
            h_B_ref(ib, 0, 1) = 0.0;
            h_B_ref(ib, 1, 1) = 0.0;
          }
        } else {
          // Transpose/ConjTrans
          if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
            // Non-Unit
            h_B_ref(ib, 0, 0) = 3.0 / 4.0;
            h_B_ref(ib, 1, 0) = 3.0 / 4.0;
            h_B_ref(ib, 0, 1) = 3.0 / 4.0;
            h_B_ref(ib, 1, 1) = 3.0 / 4.0;
          } else {
            // Unit
            h_B_ref(ib, 0, 0) = 0.0;
            h_B_ref(ib, 1, 0) = 0.0;
            h_B_ref(ib, 0, 1) = 3.0 / 2.0;
            h_B_ref(ib, 1, 1) = 3.0 / 2.0;
          }
        }
      } else {
        // Lower
        if (std::is_same_v<typename ParamTagType::trans, Trans::NoTranspose>) {
          // No-Transpose
          if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
            // Non-Unit
            h_B_ref(ib, 0, 0) = 3.0 / 8.0;
            h_B_ref(ib, 1, 0) = 3.0 / 8.0;
            h_B_ref(ib, 0, 1) = 3.0 / 4.0;
            h_B_ref(ib, 1, 1) = 3.0 / 4.0;
          } else {
            // Unit
            h_B_ref(ib, 0, 0) = -3.0 / 4.0;
            h_B_ref(ib, 1, 0) = -3.0 / 4.0;
            h_B_ref(ib, 0, 1) = 3.0 / 2.0;
            h_B_ref(ib, 1, 1) = 3.0 / 2.0;
          }
        } else {
          // Transpose/ConjTrans
          if (std::is_same_v<typename ParamTagType::diag, Diag::NonUnit>) {
            // Non-Unit
            h_B_ref(ib, 0, 0) = 3.0 / 2.0;
            h_B_ref(ib, 1, 0) = 3.0 / 2.0;
            h_B_ref(ib, 0, 1) = -3.0 / 8.0;
            h_B_ref(ib, 1, 1) = -3.0 / 8.0;
          } else {
            // Unit
            h_B_ref(ib, 0, 0) = 3.0 / 2.0;
            h_B_ref(ib, 1, 0) = 3.0 / 2.0;
            h_B_ref(ib, 0, 1) = -3.0 / 4.0;
            h_B_ref(ib, 1, 1) = -3.0 / 4.0;
          }
        }
      }
    }
  }

  Kokkos::deep_copy(A, h_A);

  // Initialize B with ones
  Kokkos::deep_copy(B, 1.0);

  // Solve using trsm
  ScalarType alpha(1.5);
  Functor_TestBatchedSerialTrsm<DeviceType, View3DType, ScalarType, ParamTagType, AlgoTagType>(alpha, A, B).run();

  Kokkos::fence();

  // Check B = B_ref
  RealType eps = 1.0e1 * ats::epsilon();
  auto h_B     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B);
  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t i = 0; i < m; i++) {
      for (std::size_t j = 0; j < m; j++) {
        EXPECT_NEAR_KK(h_B(ib, i, j), h_B_ref(ib, i, j), eps);
      }
    }
  }
}

/// \brief Implementation details of batched trsv test
/// \param N [in] Batch size of matrices and RHS
/// \param m [in] Number of rows of matrix A
/// \param n [in] Number of columns of matrix A
template <typename DeviceType, typename ScalarType, typename ValueType, typename LayoutType, typename ParamTagType,
          typename AlgoTagType>
void impl_test_batched_trsm(const std::size_t N, const std::size_t m, const std::size_t n) {
  using ats        = typename KokkosKernels::ArithTraits<ValueType>;
  using RealType   = typename ats::mag_type;
  using View3DType = Kokkos::View<ValueType ***, LayoutType, DeviceType>;

  const bool is_side_right  = std::is_same_v<typename ParamTagType::side, Side::Right>;
  const std::size_t b_nrows = is_side_right ? n : m;
  const std::size_t b_ncols = is_side_right ? m : n;

  View3DType A("A", N, m, m), B("B", N, b_nrows, b_ncols), C("C", N, b_nrows, b_ncols),
      B_ref("B_ref", N, b_nrows, b_ncols);
  View3DType Atri("Atri", N, m, m);  // Triangular components of A

  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::fill_random(B, rand_pool, randStart, randEnd);

  create_triangular_matrix<View3DType, View3DType, typename ParamTagType::uplo, typename ParamTagType::diag>(A, Atri);

  // Keep original B in B_ref
  Kokkos::deep_copy(B_ref, B);

  // Solve using trsm
  ScalarType alpha(1.5);
  Functor_TestBatchedSerialTrsm<DeviceType, View3DType, ScalarType, ParamTagType, AlgoTagType>(alpha, A, B).run();

  using trans_type = typename ParamTagType::trans;

  if (is_side_right) {
    if constexpr (std::is_same_v<trans_type, Trans::ConjTranspose>) {
      // ConjTrans is not implemented in Gemm, need to compute conj and use
      // Transpose Atri -> conj(Atri)
      auto h_Atri = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Atri);
      for (std::size_t ib = 0; ib < N; ib++) {
        for (std::size_t i = 0; i < m; i++) {
          for (std::size_t j = 0; j < m; j++) {
            h_Atri(ib, i, j) = KokkosKernels::ArithTraits<ScalarType>::conj(h_Atri(ib, i, j));
          }
        }
      }
      Kokkos::deep_copy(Atri, h_Atri);

      // Compute 1/alpha * B * Op(A) => C
      Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose,
                                Trans::Transpose>(1.0 / alpha, B, Atri, 0.0, C)
          .run();
    } else {
      // Compute 1/alpha * B * Op(A) => C
      Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::NoTranspose,
                                typename ParamTagType::trans>(1.0 / alpha, B, Atri, 0.0, C)
          .run();
    }
  } else {
    if constexpr (std::is_same_v<trans_type, Trans::ConjTranspose>) {
      // ConjTrans is not implemented in Gemm, need to compute conj and use
      // Transpose Atri -> conj(Atri)
      auto h_Atri = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Atri);
      for (std::size_t ib = 0; ib < N; ib++) {
        for (std::size_t i = 0; i < m; i++) {
          for (std::size_t j = 0; j < m; j++) {
            h_Atri(ib, i, j) = KokkosKernels::ArithTraits<ScalarType>::conj(h_Atri(ib, i, j));
          }
        }
      }
      Kokkos::deep_copy(Atri, h_Atri);

      // Compute 1/alpha * Op(A) * B => C
      Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType, Trans::Transpose,
                                Trans::NoTranspose>(1.0 / alpha, Atri, B, 0.0, C)
          .run();
    } else {
      // Compute 1/alpha * Op(A) * B => C
      Functor_BatchedSerialGemm<DeviceType, ScalarType, View3DType, View3DType, View3DType,
                                typename ParamTagType::trans, Trans::NoTranspose>(1.0 / alpha, Atri, B, 0.0, C)
          .run();
    }
  }

  Kokkos::fence();

  // Check 1/alpha * Op(A) * X == B_ref
  RealType eps = 1.0e3 * ats::epsilon();
  auto h_C     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), C);
  auto h_B_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), B_ref);
  for (std::size_t ib = 0; ib < N; ib++) {
    for (std::size_t i = 0; i < b_nrows; i++) {
      for (std::size_t j = 0; j < b_ncols; j++) {
        EXPECT_NEAR_KK(h_C(ib, i, j), h_B_ref(ib, i, j), eps);
      }
    }
  }
}
}  // namespace Trsm
}  // namespace Test

template <typename DeviceType, typename ValueType, typename ScalarType, typename ParamTagType, typename AlgoTagType>
int test_batched_trsm() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Trsm::impl_test_batched_trsm_analytical<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(1);
    Test::Trsm::impl_test_batched_trsm_analytical<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(2);
    // FIXME: ConjTranspose with blocking is not implemented yet
    if constexpr (!std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>) {
      Test::Trsm::impl_test_batched_trsm_blocking<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(0, 10, 4);
    }

    std::vector<std::size_t> sizes = {0, 1, 2, 3, 4, 5};
    for (auto m : sizes) {
      for (auto n : sizes) {
        Test::Trsm::impl_test_batched_trsm<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType, AlgoTagType>(
            1, m, n);
        Test::Trsm::impl_test_batched_trsm<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType, AlgoTagType>(
            2, m, n);
      }
    }
    for (int i = 0; i < 10; ++i) {
      // FIXME: ConjTranspose with blocking is not implemented yet
      if constexpr (!std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>) {
        Test::Trsm::impl_test_batched_trsm_blocking<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                    AlgoTagType>(1024, i, 4);
        Test::Trsm::impl_test_batched_trsm_blocking<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                    AlgoTagType>(1024, i, 1);
      }
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Trsm::impl_test_batched_trsm_analytical<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(1);
    Test::Trsm::impl_test_batched_trsm_analytical<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(2);
    // FIXME: ConjTranspose with blocking is not implemented yet
    if constexpr (!std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>) {
      Test::Trsm::impl_test_batched_trsm_blocking<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                  AlgoTagType>(0, 10, 4);
    }

    std::vector<std::size_t> sizes = {0, 1, 2, 3, 4, 5};
    for (auto m : sizes) {
      for (auto n : sizes) {
        Test::Trsm::impl_test_batched_trsm<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType, AlgoTagType>(
            1, m, n);
        Test::Trsm::impl_test_batched_trsm<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType, AlgoTagType>(
            2, m, n);
      }
    }
    for (int i = 0; i < 10; ++i) {
      // FIXME: ConjTranspose with blocking is not implemented yet
      if constexpr (!std::is_same_v<typename ParamTagType::trans, Trans::ConjTranspose>) {
        Test::Trsm::impl_test_batched_trsm_blocking<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                    AlgoTagType>(1024, i, 4);
        Test::Trsm::impl_test_batched_trsm_blocking<DeviceType, ScalarType, ValueType, LayoutType, ParamTagType,
                                                    AlgoTagType>(1024, i, 1);
      }
    }
  }
#endif

  return 0;
}
