// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Getrf.hpp>
#include <KokkosBatched_Gbtrf.hpp>
#include "Test_Batched_DenseUtils.hpp"

namespace Test {
namespace Gbtrf {

template <typename DeviceType, typename ABViewType, typename PivViewType, typename AlgoTagType>
struct Functor_BatchedSerialGbtrf {
  using execution_space = typename DeviceType::execution_space;
  ABViewType m_ab;
  PivViewType m_ipiv;
  int m_kl, m_ku, m_m;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGbtrf(const ABViewType &ab, const PivViewType &ipiv, int kl, int ku, int m)
      : m_ab(ab), m_ipiv(ipiv), m_kl(kl), m_ku(ku), m_m(m) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k, int &info) const {
    auto sub_ab   = Kokkos::subview(m_ab, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_ipiv = Kokkos::subview(m_ipiv, k, Kokkos::ALL());

    info += KokkosBatched::SerialGbtrf<AlgoTagType>::invoke(sub_ab, sub_ipiv, m_kl, m_ku, m_m);
  }

  inline int run() {
    using value_type = typename ABViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGbtrf");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    int info_sum                      = 0;
    Kokkos::Profiling::pushRegion(name.c_str());
    Kokkos::RangePolicy<execution_space> policy(0, m_ab.extent(0));
    Kokkos::parallel_reduce(name.c_str(), policy, *this, info_sum);
    Kokkos::Profiling::popRegion();
    return info_sum;
  }
};

template <typename DeviceType, typename AViewType, typename PivViewType, typename AlgoTagType>
struct Functor_BatchedSerialGetrf {
  using execution_space = typename DeviceType::execution_space;
  AViewType m_a;
  PivViewType m_ipiv;

  KOKKOS_INLINE_FUNCTION
  Functor_BatchedSerialGetrf(const AViewType &a, const PivViewType &ipiv) : m_a(a), m_ipiv(ipiv) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int k) const {
    auto sub_a    = Kokkos::subview(m_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto sub_ipiv = Kokkos::subview(m_ipiv, k, Kokkos::ALL());

    KokkosBatched::SerialGetrf<KokkosBatched::Algo::Getrf::Unblocked>::invoke(sub_a, sub_ipiv);
  }

  void run() {
    using value_type = typename AViewType::non_const_value_type;
    std::string name_region("KokkosBatched::Test::SerialGbtrf");
    const std::string name_value_type = Test::value_type_name<value_type>();
    std::string name                  = name_region + name_value_type;
    Kokkos::RangePolicy<execution_space> policy(0, m_a.extent(0));
    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

/// \brief Implementation details of batched gbtrf analytical test
///
/// \param Nb [in] Batch size of matrices
///        4x4 matrix
///        which satisfies PA = LU
///        P = [[0, 0, 1, 0],
///             [1, 0, 0, 0],
///             [0, 1, 0, 0],
///             [0, 0, 0, 1]]
///        A: [[1, -3, -2,  0],           AB: [[0,   0,  0,  0],
///            [-1, 1, -3, -2],                [0,   0,  0,  0],
///            [2, -1,  1, -3],                [0,   0, -2, -2],
///            [0,  2, -1,  1]]                [0,  -3, -3, -3],
///                                            [1,   1,  1,  1],
///                                            [-1, -1, -1,  0],
///                                            [2,   2,  0,  0]]
///
///        L = [[1,       0, 0, 0],       LUB: [[0,       0,    0,    0],
///             [0.5,     1, 0, 0],             [0,       0,    0,   -3],
///             [-0.5, -0.2, 1, 0],             [0,       0,    1,  1.5],
///             [0,    -0.8, 1, 1]]             [0,      -1, -2.5, -3.2],
///        U =  [[2, -1,      1,  -3. ],        [2,    -2.5,   -3,  5.4],
///              [0, -2.5, -2.5,  1.5 ],        [-0.5, -0.2,    1,    0],
///              [0,  0,     -3,  -3.2],        [0.5,  -0.8,    0,    0]]
///              [0,  0,      0,   5.4]]
///        Note P is obtained by piv = [2 2 2 3]
///
///        3x4 matrix
///        which satisfies PA = LU
///        P = [[0, 1, 0],
///             [0, 0, 1],
///             [1, 0, 0]]
///        A: [[1, -3, -2,  0],           AB: [[ 0,  0,  0,  0],
///            [-1, 1, -3, -2],                [ 0,  0,  0,  0],
///            [2, -1,  1, -3]]                [ 0,  0, -2, -2],
///                                            [ 0, -3, -3, -3],
///                                            [ 1,  1,  1,  0],
///                                            [-1, -1,  0,  0],
///                                            [ 2,  0,  0,  0]]
///
///        L = [[1,       0, 0, 0],       LUB: [[0,       0,    0,    0],
///             [0.5,     1, 0, 0],             [0,       0,    0,   -3],
///             [-0.5, -0.2, 1, 0]]             [0,       0,    1,  1.5],
///                                             [0,      -1, -2.5, -3.2],
///        U =  [[2, -1,      1,  -3. ],        [2,    -2.5,   -3,    0],
///              [0, -2.5, -2.5,  1.5 ],        [-0.5, -0.2,    0,    0],
///              [0,  0,     -3,  -3.2],        [0.5,     0,    0,    0]]
///              [0,  0,      0,     0]]
///        Note P is obtained by piv = [2 2 2]
///        For simplicity, we represent the matrix U in 3x4 form instead of 4x4
///
///        5x4 matrix
///        which satisfies PA = LU
///        P = [[0. 1. 0. 0. 0.],
///             [0. 0. 1. 0. 0.],
///             [1. 0. 0. 0. 0.],
///             [0. 0. 0. 1. 0.],
///             [0. 0. 0. 0. 1.]]
///        A: [[1, -3, -2,  0],           AB: [[ 0,  0,  0,  0],
///            [-1, 1, -3, -2],                [ 0,  0,  0,  0],
///            [2, -1,  1, -3],                [ 0,  0, -2, -2],
///            [0,  2, -1,  1],                [ 0, -3, -3, -3],
///            [0,  0,  0,  0]]                [ 1,  1,  1,  1],
///                                            [-1, -1, -1,  0],
///                                            [ 2,  2,  0,  0]]
///
///        L = [[1,       0,     0,     0],  LUB: [[0,       0,     0,     0],
///             [0.5,     1,     0,     0],        [0,       0,     0,    -3],
///             [-0.5, -0.2,     1,     0],        [0,       0,     1,   1.5],
///             [   0, -0.8,     1,     1],        [0,      -1,  -2.5,  -3.2],
///             [   0,    0,     0,     0]]        [2,    -2.5,    -3,   5.4],
///        U =  [[2, -1,      1,  -3. ],           [-0.5, -0.2,     1,     0],
///              [0, -2.5, -2.5,  1.5 ],           [0.5,  -0.8,     0,     0]]
///              [0,  0,     -3,  -3.2],
///              [0,  0,      0,   5.4]]
///        Note P is obtained by piv = [2 2 2 3]
///        For simplicity, we represent the matrix U in 5x4 form instead of 4x4
template <typename DeviceType, typename ScalarType, typename LayoutType, typename AlgoTagType>
void impl_test_batched_gbtrf_analytical(const int Nb) {
  using ats           = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType      = typename ats::mag_type;
  using View3DType    = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using PivView2DType = Kokkos::View<int **, LayoutType, DeviceType>;

  const int BlkSize = 4, kl = 2, ku = 2;
  const int ldab = 2 * kl + ku + 1;
  View3DType A0("A0", Nb, BlkSize, BlkSize), NL0("NL0", Nb, BlkSize, BlkSize), NL0_ref("NL0_ref", Nb, BlkSize, BlkSize),
      U0("U0", Nb, BlkSize, BlkSize), U0_ref("U0_ref", Nb, BlkSize, BlkSize);
  View3DType A1("A1", Nb, BlkSize - 1, BlkSize), NL1("NL1", Nb, BlkSize - 1, BlkSize),
      U1("U1", Nb, BlkSize - 1, BlkSize), NL1_ref("NL1_ref", Nb, BlkSize - 1, BlkSize),
      U1_ref("U1_ref", Nb, BlkSize - 1, BlkSize);
  View3DType A2("A2", Nb, BlkSize + 1, BlkSize), NL2("NL2", Nb, BlkSize + 1, BlkSize),
      U2("U2", Nb, BlkSize + 1, BlkSize), NL2_ref("NL2_ref", Nb, BlkSize + 1, BlkSize),
      U2_ref("U2_ref", Nb, BlkSize + 1, BlkSize);

  View3DType AB0("AB0", Nb, ldab, BlkSize), AB1("AB1", Nb, ldab, BlkSize), AB2("AB2", Nb, ldab, BlkSize);
  PivView2DType ipiv0("ipiv0", Nb, BlkSize), ipiv0_ref("ipiv0_ref", Nb, BlkSize), ipiv1("ipiv1", Nb, BlkSize - 1),
      ipiv1_ref("ipiv1_ref", Nb, BlkSize - 1), ipiv2("ipiv2", Nb, BlkSize), ipiv2_ref("ipiv2_ref", Nb, BlkSize);

  // Only filling A2 and deep_copy from its subview
  auto h_A2        = Kokkos::create_mirror_view(A2);
  auto h_NL2_ref   = Kokkos::create_mirror_view(NL2_ref);
  auto h_U2_ref    = Kokkos::create_mirror_view(U2_ref);
  auto h_ipiv2_ref = Kokkos::create_mirror_view(ipiv2_ref);

  for (int ib = 0; ib < Nb; ib++) {
    h_A2(ib, 0, 0) = 1;
    h_A2(ib, 0, 1) = -3;
    h_A2(ib, 0, 2) = -2;
    h_A2(ib, 0, 3) = 0;
    h_A2(ib, 1, 0) = -1;
    h_A2(ib, 1, 1) = 1;
    h_A2(ib, 1, 2) = -3;
    h_A2(ib, 1, 3) = -2;
    h_A2(ib, 2, 0) = 2;
    h_A2(ib, 2, 1) = -1;
    h_A2(ib, 2, 2) = 1;
    h_A2(ib, 2, 3) = -3;
    h_A2(ib, 3, 0) = 0;
    h_A2(ib, 3, 1) = 2;
    h_A2(ib, 3, 2) = -1;
    h_A2(ib, 3, 3) = 1;
    h_A2(ib, 4, 0) = 0;
    h_A2(ib, 4, 1) = 0;

    h_U2_ref(ib, 0, 0) = 2;
    h_U2_ref(ib, 0, 1) = -1;
    h_U2_ref(ib, 0, 2) = 1;
    h_U2_ref(ib, 0, 3) = -3;
    h_U2_ref(ib, 1, 1) = -2.5;
    h_U2_ref(ib, 1, 2) = -2.5;
    h_U2_ref(ib, 1, 3) = 1.5;
    h_U2_ref(ib, 2, 2) = -3;
    h_U2_ref(ib, 2, 3) = -3.2;
    h_U2_ref(ib, 3, 3) = 5.4;

    h_NL2_ref(ib, 1, 0) = 0.5;
    h_NL2_ref(ib, 2, 0) = -0.5;
    h_NL2_ref(ib, 2, 1) = -0.2;
    h_NL2_ref(ib, 3, 1) = -0.8;
    h_NL2_ref(ib, 3, 2) = 1.0;

    h_ipiv2_ref(ib, 0) = 2;
    h_ipiv2_ref(ib, 1) = 2;
    h_ipiv2_ref(ib, 2) = 2;
    h_ipiv2_ref(ib, 3) = 3;
  }

  // Copy matrix A2 to device
  Kokkos::deep_copy(A2, h_A2);

  // Copy submatrix of A2 to device
  auto A2_m3        = Kokkos::subview(A2, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize - 1), Kokkos::ALL);
  auto A2_m4        = Kokkos::subview(A2, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize), Kokkos::ALL);
  auto h_NL2_ref_m3 = Kokkos::subview(h_NL2_ref, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize - 1), Kokkos::ALL);
  auto h_NL2_ref_m4 = Kokkos::subview(h_NL2_ref, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize), Kokkos::ALL);
  auto h_U2_ref_m3  = Kokkos::subview(h_U2_ref, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize - 1), Kokkos::ALL);
  auto h_U2_ref_m4  = Kokkos::subview(h_U2_ref, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize), Kokkos::ALL);

  Kokkos::deep_copy(A0, A2_m4);  // Extract 4x4 matrix
  Kokkos::deep_copy(A1, A2_m3);  // Extract 3x4 matrix

  auto h_NL0_ref = Kokkos::create_mirror_view(NL0_ref);
  auto h_NL1_ref = Kokkos::create_mirror_view(NL1_ref);
  auto h_U0_ref  = Kokkos::create_mirror_view(U0_ref);
  auto h_U1_ref  = Kokkos::create_mirror_view(U1_ref);
  Kokkos::deep_copy(h_NL0_ref, h_NL2_ref_m4);  // Extract 4x4 matrix
  Kokkos::deep_copy(h_NL1_ref, h_NL2_ref_m3);  // Extract 3x4 matrix
  Kokkos::deep_copy(h_U0_ref, h_U2_ref_m4);    // Extract 4x4 matrix
  Kokkos::deep_copy(h_U1_ref, h_U2_ref_m3);    // Extract 3x4 matrix

  // Copy submatrix of ipiv2 to device
  auto h_ipiv0_ref = Kokkos::create_mirror_view(ipiv0_ref);
  auto h_ipiv1_ref = Kokkos::create_mirror_view(ipiv1_ref);
  auto h_ipiv2_m3  = Kokkos::subview(h_ipiv2_ref, Kokkos::ALL, Kokkos::pair<int, int>(0, BlkSize - 1));
  Kokkos::deep_copy(h_ipiv0_ref, h_ipiv2_ref);
  Kokkos::deep_copy(h_ipiv1_ref, h_ipiv2_m3);

  // Convert into banded storage
  dense_to_banded(A0, AB0, kl, ku);
  dense_to_banded(A1, AB1, kl, ku);
  dense_to_banded(A2, AB2, kl, ku);

  // gbtrf to factorize matrix A = P * L * U
  auto info0 =
      Functor_BatchedSerialGbtrf<DeviceType, View3DType, PivView2DType, AlgoTagType>(AB0, ipiv0, kl, ku, BlkSize).run();
  auto info1 =
      Functor_BatchedSerialGbtrf<DeviceType, View3DType, PivView2DType, AlgoTagType>(AB1, ipiv1, kl, ku, BlkSize - 1)
          .run();
  auto info2 =
      Functor_BatchedSerialGbtrf<DeviceType, View3DType, PivView2DType, AlgoTagType>(AB2, ipiv2, kl, ku, BlkSize + 1)
          .run();

  Kokkos::fence();
  EXPECT_EQ(info0, 0);
  EXPECT_EQ(info1, 0);
  EXPECT_EQ(info2, 0);

  // Extract matrix U and L from AB
  // first convert it to the dense matrix (stored in A)
  banded_to_dense<View3DType, View3DType>(AB0, A0, kl, ku);
  banded_to_dense<View3DType, View3DType>(AB1, A1, kl, ku);
  banded_to_dense<View3DType, View3DType>(AB2, A2, kl, ku);

  // Copy upper triangular components to U
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Upper, KokkosBatched::Diag::NonUnit>(A0, U0);
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Upper, KokkosBatched::Diag::NonUnit>(A1, U1);
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Upper, KokkosBatched::Diag::NonUnit>(A2, U2);

  // Extract L
  // Apply pivot at host
  auto h_ipiv0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv0);
  auto h_ipiv1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv1);
  auto h_ipiv2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv2);
  auto h_A0    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A0);
  auto h_A1    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A1);
  Kokkos::deep_copy(h_A2, A2);
  for (int ib = 0; ib < Nb; ib++) {
    for (int j = 0; j < BlkSize - 1; j++) {
      for (int i = j + 1; i < BlkSize - 1; i++) {
        Kokkos::kokkos_swap(h_A0(ib, i, j), h_A0(ib, h_ipiv0(ib, i), j));
      }
    }

    for (int j = 0; j < BlkSize - 2; j++) {
      for (int i = j + 1; i < BlkSize - 1; i++) {
        Kokkos::kokkos_swap(h_A1(ib, i, j), h_A1(ib, h_ipiv1(ib, i), j));
      }
    }

    for (int j = 0; j < BlkSize; j++) {
      for (int i = j + 1; i < BlkSize - 1; i++) {
        Kokkos::kokkos_swap(h_A2(ib, i, j), h_A2(ib, h_ipiv2(ib, i), j));
      }
    }
  }
  Kokkos::deep_copy(A0, h_A0);
  Kokkos::deep_copy(A1, h_A1);
  Kokkos::deep_copy(A2, h_A2);

  // Copy non-diagonal lower triangular components to NL
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Lower, KokkosBatched::Diag::NonUnit>(A0, NL0,
                                                                                                             -1);
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Lower, KokkosBatched::Diag::NonUnit>(A1, NL1,
                                                                                                             -1);
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Lower, KokkosBatched::Diag::NonUnit>(A2, NL2,
                                                                                                             -1);

  // Check if U, NL and ipiv have expected values
  auto h_U0  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), U0);
  auto h_U1  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), U1);
  auto h_U2  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), U2);
  auto h_NL0 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NL0);
  auto h_NL1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NL1);
  auto h_NL2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NL2);

  RealType eps = 1.0e1 * ats::epsilon();
  for (int ib = 0; ib < Nb; ib++) {
    for (int j = 0; j < BlkSize; j++) {
      EXPECT_EQ(h_ipiv0(ib, j), h_ipiv0_ref(ib, j));
      EXPECT_EQ(h_ipiv2(ib, j), h_ipiv2_ref(ib, j));
      for (int i = 0; i < BlkSize; i++) {
        EXPECT_NEAR_KK(h_U0(ib, i, j), h_U0_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_NL0(ib, i, j), h_NL0_ref(ib, i, j), eps);
      }
      for (int i = 0; i < BlkSize - 1; i++) {
        EXPECT_NEAR_KK(h_U1(ib, i, j), h_U1_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_NL1(ib, i, j), h_NL1_ref(ib, i, j), eps);
      }
      for (int i = 0; i < BlkSize + 1; i++) {
        EXPECT_NEAR_KK(h_U2(ib, i, j), h_U2_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_NL2(ib, i, j), h_NL2_ref(ib, i, j), eps);
      }
    }
    for (int j = 0; j < BlkSize - 1; j++) {
      EXPECT_EQ(h_ipiv1(ib, j), h_ipiv1_ref(ib, j));
    }
  }
}

/// \brief Implementation details of batched gbtrf test
///
/// \param N [in] Batch size of matrix A
/// \param BlkSize [in] Block size of matrix A
template <typename DeviceType, typename ScalarType, typename LayoutType, typename AlgoTagType>
void impl_test_batched_gbtrf(const int Nb, const int BlkSize) {
  using ats           = typename KokkosKernels::ArithTraits<ScalarType>;
  using RealType      = typename ats::mag_type;
  using View3DType    = Kokkos::View<ScalarType ***, LayoutType, DeviceType>;
  using PivView2DType = Kokkos::View<int **, LayoutType, DeviceType>;

  const int kl = 2, ku = 2;
  const int ldab = 2 * kl + ku + 1;
  View3DType A("A", Nb, BlkSize, BlkSize), LU("LU", Nb, BlkSize, BlkSize), NL("NL", Nb, BlkSize, BlkSize),
      NL_ref("NL_ref", Nb, BlkSize, BlkSize), U("U", Nb, BlkSize, BlkSize), U_ref("U_ref", Nb, BlkSize, BlkSize);

  View3DType AB("AB", Nb, ldab, BlkSize), AB_upper("AB_upper", Nb, kl + ku + 1, BlkSize);
  PivView2DType ipiv("ipiv", Nb, BlkSize), ipiv_ref("ipiv_ref", Nb, BlkSize);

  // Create a random matrix A and make it Positive Definite Symmetric
  using execution_space = typename DeviceType::execution_space;
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);
  ScalarType randStart, randEnd;

  // Initialize LU with random matrix
  KokkosKernels::Impl::getRandomBounds(1.0, randStart, randEnd);
  Kokkos::fill_random(LU, rand_pool, randStart, randEnd);

  dense_to_banded(LU, AB, kl, ku);  // In banded storage
  banded_to_dense(AB, A, kl, ku);   // In conventional storage

  Kokkos::deep_copy(LU, A);  // for getrf

  // gbtrf to factorize matrix A = P * L * U
  Functor_BatchedSerialGbtrf<DeviceType, View3DType, PivView2DType, AlgoTagType>(AB, ipiv, kl, ku, BlkSize).run();

  // Extract matrix U and L from AB
  // first convert it to the dense matrix (stored in A)
  banded_to_dense<View3DType, View3DType>(AB, A, kl, ku);

  // Copy upper triangular components to U
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Upper, KokkosBatched::Diag::NonUnit>(A, U);

  // Apply pivot at host
  auto h_ipiv = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv);
  auto h_A    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
  for (int ib = 0; ib < Nb; ib++) {
    for (int j = 0; j < BlkSize - 1; j++) {
      for (int i = j + 1; i < BlkSize - 1; i++) {
        Kokkos::kokkos_swap(h_A(ib, i, j), h_A(ib, h_ipiv(ib, i), j));
      }
    }
  }
  Kokkos::deep_copy(A, h_A);

  // Copy non-diagonal lower triangular components to NL
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Lower, KokkosBatched::Diag::NonUnit>(A, NL, -1);

  // Reference is made by getrf
  // getrf to factorize matrix A = P * L * U
  Functor_BatchedSerialGetrf<DeviceType, View3DType, PivView2DType, AlgoTagType>(LU, ipiv_ref).run();

  // Copy non-diagonal lower triangular components to NL
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Lower, KokkosBatched::Diag::NonUnit>(LU, NL_ref,
                                                                                                             -1);

  // Copy upper triangular components to U_ref
  create_triangular_matrix<View3DType, View3DType, KokkosBatched::Uplo::Upper, KokkosBatched::Diag::NonUnit>(LU, U_ref);

  // Check if U, NL and ipiv have expected values
  auto h_U        = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), U);
  auto h_U_ref    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), U_ref);
  auto h_NL       = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NL);
  auto h_NL_ref   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), NL_ref);
  auto h_ipiv_ref = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv_ref);

  RealType eps = 1.0e3 * ats::epsilon();
  for (int ib = 0; ib < Nb; ib++) {
    for (int i = 0; i < BlkSize; i++) {
      EXPECT_EQ(h_ipiv(ib, i), h_ipiv_ref(ib, i));
      for (int j = 0; j < BlkSize; j++) {
        EXPECT_NEAR_KK(h_U(ib, i, j), h_U_ref(ib, i, j), eps);
        EXPECT_NEAR_KK(h_NL(ib, i, j), h_NL_ref(ib, i, j), eps);
      }
    }
  }
}

}  // namespace Gbtrf
}  // namespace Test

template <typename DeviceType, typename ScalarType, typename AlgoTagType>
int test_batched_gbtrf() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  {
    using LayoutType = Kokkos::LayoutLeft;
    Test::Gbtrf::impl_test_batched_gbtrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(1);
    Test::Gbtrf::impl_test_batched_gbtrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Gbtrf::impl_test_batched_gbtrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Gbtrf::impl_test_batched_gbtrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  {
    using LayoutType = Kokkos::LayoutRight;
    Test::Gbtrf::impl_test_batched_gbtrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(1);
    Test::Gbtrf::impl_test_batched_gbtrf_analytical<DeviceType, ScalarType, LayoutType, AlgoTagType>(2);
    for (int i = 0; i < 10; i++) {
      Test::Gbtrf::impl_test_batched_gbtrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(1, i);
      Test::Gbtrf::impl_test_batched_gbtrf<DeviceType, ScalarType, LayoutType, AlgoTagType>(2, i);
    }
  }
#endif

  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, test_batched_gbtrf_float) {
  using algo_tag_type = typename KokkosBatched::Algo::Gbtrf::Unblocked;

  test_batched_gbtrf<TestDevice, float, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, test_batched_gbtrf_double) {
  using algo_tag_type = typename KokkosBatched::Algo::Gbtrf::Unblocked;

  test_batched_gbtrf<TestDevice, double, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, test_batched_gbtrf_fcomplex) {
  using algo_tag_type = typename KokkosBatched::Algo::Gbtrf::Unblocked;

  test_batched_gbtrf<TestDevice, Kokkos::complex<float>, algo_tag_type>();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, test_batched_gbtrf_dcomplex) {
  using algo_tag_type = typename KokkosBatched::Algo::Gbtrf::Unblocked;

  test_batched_gbtrf<TestDevice, Kokkos::complex<double>, algo_tag_type>();
}
#endif
