//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <chrono>

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_TestUtils.hpp>
#include <KokkosBlas3_gemm.hpp>

#include <KokkosLapack_svd.hpp>

namespace Test {

template <class AMatrix, class SVector, class UMatrix, class VMatrix>
void check_triple_product(const AMatrix& A, const SVector& S, const UMatrix& U, const VMatrix& Vt,
                          typename Kokkos::ArithTraits<typename AMatrix::non_const_value_type>::mag_type tol) {
  // After a successful SVD decomposition we have A=U*S*V
  // So using gemm we should be able to compare the above
  // triple product to the original matrix A.
  using execution_space = typename AMatrix::execution_space;

  AMatrix temp("intermediate U*S product", A.extent(0), A.extent(1));
  AMatrix M("U*S*V product", A.extent(0), A.extent(1));

  // First compute the left side of the product: temp = U*S
  Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space, int>(0, U.extent_int(0)), KOKKOS_LAMBDA(const int& rowIdx) {
        for (int colIdx = 0; colIdx < U.extent_int(1); ++colIdx) {
          if (colIdx < S.extent_int(0)) {
            temp(rowIdx, colIdx) = U(rowIdx, colIdx) * S(colIdx);
          }
        }
      });

  // Second compute the right side of the product: M = temp*V = U*S*V
  KokkosBlas::gemm("N", "N", 1, temp, Vt, 0, M);

  typename AMatrix::HostMirror A_h = Kokkos::create_mirror_view(A);
  typename AMatrix::HostMirror M_h = Kokkos::create_mirror_view(M);
  Kokkos::deep_copy(A_h, A);
  Kokkos::deep_copy(M_h, M);
  for (int rowIdx = 0; rowIdx < A.extent_int(0); ++rowIdx) {
    for (int colIdx = 0; colIdx < A.extent_int(1); ++colIdx) {
      if (tol < Kokkos::abs(A_h(rowIdx, colIdx))) {
        EXPECT_NEAR_KK_REL(A_h(rowIdx, colIdx), M_h(rowIdx, colIdx), tol);
      } else {
        EXPECT_NEAR_KK(A_h(rowIdx, colIdx), M_h(rowIdx, colIdx), tol);
      }
    }
  }
}

template <class Matrix>
void check_unitary_orthogonal_matrix(
    const Matrix& M, typename Kokkos::ArithTraits<typename Matrix::non_const_value_type>::mag_type tol) {
  // After a successful SVD decomposition the matrices
  // U and V are unitary matrices. Thus we can check
  // the property UUt=UtU=I and VVt=VtV=I using gemm.
  using scalar_type = typename Matrix::non_const_value_type;

  Matrix I0("M*Mt", M.extent(0), M.extent(0));
  KokkosBlas::gemm("N", "C", 1, M, M, 0, I0);
  typename Matrix::HostMirror I0_h = Kokkos::create_mirror_view(I0);
  Kokkos::deep_copy(I0_h, I0);
  for (int rowIdx = 0; rowIdx < M.extent_int(0); ++rowIdx) {
    for (int colIdx = 0; colIdx < M.extent_int(0); ++colIdx) {
      if (rowIdx == colIdx) {
        EXPECT_NEAR_KK_REL(I0_h(rowIdx, colIdx), Kokkos::ArithTraits<scalar_type>::one(), tol);
      } else {
        EXPECT_NEAR_KK(I0_h(rowIdx, colIdx), Kokkos::ArithTraits<scalar_type>::zero(), tol);
      }
    }
  }

  Matrix I1("Mt*M", M.extent(1), M.extent(1));
  KokkosBlas::gemm("C", "N", 1, M, M, 0, I1);
  typename Matrix::HostMirror I1_h = Kokkos::create_mirror_view(I1);
  Kokkos::deep_copy(I1_h, I1);
  for (int rowIdx = 0; rowIdx < M.extent_int(1); ++rowIdx) {
    for (int colIdx = 0; colIdx < M.extent_int(1); ++colIdx) {
      if (rowIdx == colIdx) {
        EXPECT_NEAR_KK_REL(I1_h(rowIdx, colIdx), Kokkos::ArithTraits<scalar_type>::one(), tol);
      } else {
        EXPECT_NEAR_KK(I1_h(rowIdx, colIdx), Kokkos::ArithTraits<scalar_type>::zero(), tol);
      }
    }
  }
}

template <class AMatrix, class Device>
int impl_analytic_2x2_svd() {
  using scalar_type = typename AMatrix::value_type;
  using mag_type    = typename Kokkos::ArithTraits<scalar_type>::mag_type;
  using vector_type = Kokkos::View<mag_type*, typename AMatrix::array_layout, Device>;
  using KAT_S       = Kokkos::ArithTraits<scalar_type>;

  const mag_type eps = KAT_S::eps();

  AMatrix A("A", 2, 2), U("U", 2, 2), Vt("Vt", 2, 2), Aref("A ref", 2, 2);
  vector_type S("S", 2);

  typename AMatrix::HostMirror A_h = Kokkos::create_mirror_view(A);

  // A = [3  0]
  //     [4  5]
  // USV = 1/sqrt(10) [1  -3] * sqrt(5) [3  0] *  1/sqrt(2) [ 1  1]
  //                  [3   1]           [0  1]              [-1  1]
  A_h(0, 0) = 3;
  A_h(1, 0) = 4;
  A_h(1, 1) = 5;

  Kokkos::deep_copy(A, A_h);
  Kokkos::deep_copy(Aref, A_h);

  KokkosLapack::svd("A", "A", A, S, U, Vt);
  // Don't really need to fence here as we deep_copy right after...

  typename vector_type::HostMirror S_h = Kokkos::create_mirror_view(S);
  Kokkos::deep_copy(S_h, S);
  typename AMatrix::HostMirror U_h = Kokkos::create_mirror_view(U);
  Kokkos::deep_copy(U_h, U);
  typename AMatrix::HostMirror Vt_h = Kokkos::create_mirror_view(Vt);
  Kokkos::deep_copy(Vt_h, Vt);

  // The singular values for this problem
  // are known: sqrt(45) and sqrt(5)
  EXPECT_NEAR_KK_REL(S_h(0), static_cast<mag_type>(Kokkos::sqrt(45)), 100 * eps);
  EXPECT_NEAR_KK_REL(S_h(1), static_cast<mag_type>(Kokkos::sqrt(5)), 100 * eps);

  // The singular vectors should be identical
  // or of oposite sign we check the first
  // component of the vectors to determine
  // the proper signed comparison.
  std::vector<scalar_type> Uref = {
      static_cast<scalar_type>(1 / Kokkos::sqrt(10)), static_cast<scalar_type>(3 / Kokkos::sqrt(10)),
      static_cast<scalar_type>(-3 / Kokkos::sqrt(10)), static_cast<scalar_type>(1 / Kokkos::sqrt(10))};
  std::vector<scalar_type> Vtref = {
      static_cast<scalar_type>(1 / Kokkos::sqrt(2)), static_cast<scalar_type>(-1 / Kokkos::sqrt(2)),
      static_cast<scalar_type>(1 / Kokkos::sqrt(2)), static_cast<scalar_type>(1 / Kokkos::sqrt(2))};

  // Both rotations and reflections are valid
  // vector basis so we need to check both signs
  // to confirm proper SVD was achieved.
  Kokkos::View<mag_type**, Kokkos::HostSpace> U_real("U real", 2, 2), Vt_real("Vt real", 2, 2);
  if constexpr (KAT_S::is_complex) {
    U_real(0, 0) = U_h(0, 0).real();
    U_real(0, 1) = U_h(0, 1).real();
    U_real(1, 0) = U_h(1, 0).real();
    U_real(1, 1) = U_h(1, 1).real();

    Vt_real(0, 0) = Vt_h(0, 0).real();
    Vt_real(0, 1) = Vt_h(0, 1).real();
    Vt_real(1, 0) = Vt_h(1, 0).real();
    Vt_real(1, 1) = Vt_h(1, 1).real();
  } else {
    U_real(0, 0) = U_h(0, 0);
    U_real(0, 1) = U_h(0, 1);
    U_real(1, 0) = U_h(1, 0);
    U_real(1, 1) = U_h(1, 1);

    Vt_real(0, 0) = Vt_h(0, 0);
    Vt_real(0, 1) = Vt_h(0, 1);
    Vt_real(1, 0) = Vt_h(1, 0);
    Vt_real(1, 1) = Vt_h(1, 1);
  }

  const mag_type tol        = 100 * KAT_S::eps();
  const mag_type one_sqrt10 = static_cast<mag_type>(1 / Kokkos::sqrt(10));
  const mag_type one_sqrt2  = static_cast<mag_type>(1 / Kokkos::sqrt(2));

  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(0, 0)), one_sqrt10, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(0, 1)), 3 * one_sqrt10, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(1, 0)), 3 * one_sqrt10, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(1, 1)), one_sqrt10, tol);

  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(0, 0)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(0, 1)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(1, 0)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(1, 1)), one_sqrt2, tol);

  check_unitary_orthogonal_matrix(U, tol);
  check_unitary_orthogonal_matrix(Vt, tol);

  check_triple_product(Aref, S, U, Vt, tol);

  return 0;
}

template <class AMatrix, class Device>
int impl_analytic_2x3_svd() {
  using scalar_type = typename AMatrix::value_type;
  using mag_type    = typename Kokkos::ArithTraits<scalar_type>::mag_type;
  using vector_type = Kokkos::View<mag_type*, typename AMatrix::array_layout, Device>;
  using KAT_S       = Kokkos::ArithTraits<scalar_type>;

  const mag_type tol = 100 * KAT_S::eps();

  AMatrix A("A", 2, 3), U("U", 2, 2), Vt("Vt", 3, 3), Aref("A ref", 2, 3);
  vector_type S("S", 2);

  typename AMatrix::HostMirror A_h = Kokkos::create_mirror_view(A);

  // A = [3  2   2]
  //     [2  3  -2]
  // USVt = 1/sqrt(2) [1   1] * [5  0  0] *  1/(3*sqrt(2)) [        3 3 0]
  //                  [1  -1]   [0  3  0]                  [        1 -1 4]
  //                                                       [2*sqrt(2) -2*sqrt(2)
  //                                                       -sqrt(2)]
  A_h(0, 0) = 3;
  A_h(0, 1) = 2;
  A_h(0, 2) = 2;
  A_h(1, 0) = 2;
  A_h(1, 1) = 3;
  A_h(1, 2) = -2;

  Kokkos::deep_copy(A, A_h);
  Kokkos::deep_copy(Aref, A_h);

  try {
    KokkosLapack::svd("A", "A", A, S, U, Vt);
  } catch (const std::runtime_error& e) {
    std::string test_string = e.what();
    std::string cusolver_m_less_than_n =
        "CUSOLVER does not support SVD for matrices with more columns "
        "than rows, you can transpose you matrix first then compute "
        "SVD of that transpose: At=VSUt, and swap the output U and Vt"
        " and transpose them to recover the desired SVD.";

    if (test_string == cusolver_m_less_than_n) {
      return 0;
    }
  }
  // Don't really need to fence here as we deep_copy right after...

  typename vector_type::HostMirror S_h = Kokkos::create_mirror_view(S);
  Kokkos::deep_copy(S_h, S);
  typename AMatrix::HostMirror U_h = Kokkos::create_mirror_view(U);
  Kokkos::deep_copy(U_h, U);
  typename AMatrix::HostMirror Vt_h = Kokkos::create_mirror_view(Vt);
  Kokkos::deep_copy(Vt_h, Vt);

  // The singular values for this problem
  // are known: sqrt(45) and sqrt(5)
  EXPECT_NEAR_KK_REL(S_h(0), static_cast<mag_type>(5), tol);
  EXPECT_NEAR_KK_REL(S_h(1), static_cast<mag_type>(3), tol);

  // Both rotations and reflections are valid
  // vector basis so we need to check both signs
  // to confirm proper SVD was achieved.
  Kokkos::View<mag_type**, Kokkos::HostSpace> U_real("U real", 2, 2), Vt_real("Vt real", 3, 3);
  if constexpr (KAT_S::is_complex) {
    U_real(0, 0) = U_h(0, 0).real();
    U_real(0, 1) = U_h(0, 1).real();
    U_real(1, 0) = U_h(1, 0).real();
    U_real(1, 1) = U_h(1, 1).real();

    Vt_real(0, 0) = Vt_h(0, 0).real();
    Vt_real(0, 1) = Vt_h(0, 1).real();
    Vt_real(0, 2) = Vt_h(0, 2).real();
    Vt_real(1, 0) = Vt_h(1, 0).real();
    Vt_real(1, 1) = Vt_h(1, 1).real();
    Vt_real(1, 2) = Vt_h(1, 2).real();
    Vt_real(2, 0) = Vt_h(2, 0).real();
    Vt_real(2, 1) = Vt_h(2, 1).real();
    Vt_real(2, 2) = Vt_h(2, 2).real();
  } else {
    U_real(0, 0) = U_h(0, 0);
    U_real(0, 1) = U_h(0, 1);
    U_real(1, 0) = U_h(1, 0);
    U_real(1, 1) = U_h(1, 1);

    Vt_real(0, 0) = Vt_h(0, 0);
    Vt_real(0, 1) = Vt_h(0, 1);
    Vt_real(0, 2) = Vt_h(0, 2);
    Vt_real(1, 0) = Vt_h(1, 0);
    Vt_real(1, 1) = Vt_h(1, 1);
    Vt_real(1, 2) = Vt_h(1, 2);
    Vt_real(2, 0) = Vt_h(2, 0);
    Vt_real(2, 1) = Vt_h(2, 1);
    Vt_real(2, 2) = Vt_h(2, 2);
  }

  const mag_type one_sqrt2  = static_cast<mag_type>(1 / Kokkos::sqrt(2));
  const mag_type one_sqrt18 = static_cast<mag_type>(1 / Kokkos::sqrt(18));
  const mag_type one_third  = static_cast<mag_type>(1. / 3.);

  // Check values of U
  // Don't worry about the sign
  // it will be check with the
  // triple product
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(0, 0)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(0, 1)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(1, 0)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(1, 1)), one_sqrt2, tol);

  // Check values of Vt
  // Don't worry about the sign
  // it will be check with the
  // triple product
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(0, 0)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(0, 1)), one_sqrt2, tol);
  EXPECT_NEAR_KK(Kokkos::abs(Vt_real(0, 2)), 0, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(1, 0)), one_sqrt18, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(1, 1)), one_sqrt18, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(1, 2)), 4 * one_sqrt18, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(2, 0)), 2 * one_third, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(2, 1)), 2 * one_third, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(2, 2)), one_third, tol);

  check_unitary_orthogonal_matrix(U, tol);
  check_unitary_orthogonal_matrix(Vt, tol);

  check_triple_product(Aref, S, U, Vt, tol);

  return 0;
}

template <class AMatrix, class Device>
int impl_analytic_3x2_svd() {
  using scalar_type = typename AMatrix::value_type;
  using mag_type    = typename Kokkos::ArithTraits<scalar_type>::mag_type;
  using vector_type = Kokkos::View<mag_type*, typename AMatrix::array_layout, Device>;
  using KAT_S       = Kokkos::ArithTraits<scalar_type>;

  const mag_type tol = 100 * KAT_S::eps();

  AMatrix A("A", 3, 2), U("U", 3, 3), Vt("Vt", 2, 2), Aref("A ref", 3, 2);
  vector_type S("S", 2);

  typename AMatrix::HostMirror A_h = Kokkos::create_mirror_view(A);

  // Note this is simply the transpose of the 2x3 matrix in the test above
  // A = [3  2]
  //     [2  3]
  //     [2 -2]
  // USVt = 1/(3*sqrt(2)) [3   1   2*sqrt(2)] * [5  0] * 1/sqrt(2) [1   1]
  //                      [3  -1  -2*sqrt(2)]   [0  3]             [1  -1]
  //                      [0   4     sqrt(2)]   [0  0]
  A_h(0, 0) = 3;
  A_h(0, 1) = 2;
  A_h(1, 0) = 2;
  A_h(1, 1) = 3;
  A_h(2, 0) = 2;
  A_h(2, 1) = -2;

  Kokkos::deep_copy(A, A_h);
  Kokkos::deep_copy(Aref, A_h);

  KokkosLapack::svd("A", "A", A, S, U, Vt);
  // Don't really need to fence here as we deep_copy right after...

  typename vector_type::HostMirror S_h = Kokkos::create_mirror_view(S);
  Kokkos::deep_copy(S_h, S);
  typename AMatrix::HostMirror U_h = Kokkos::create_mirror_view(U);
  Kokkos::deep_copy(U_h, U);
  typename AMatrix::HostMirror Vt_h = Kokkos::create_mirror_view(Vt);
  Kokkos::deep_copy(Vt_h, Vt);

  // The singular values for this problem
  // are known: sqrt(45) and sqrt(5)
  EXPECT_NEAR_KK_REL(S_h(0), static_cast<mag_type>(5), tol);
  EXPECT_NEAR_KK_REL(S_h(1), static_cast<mag_type>(3), tol);

  // Both rotations and reflections are valid
  // vector basis so we need to check both signs
  // to confirm proper SVD was achieved.
  Kokkos::View<mag_type**, Kokkos::HostSpace> U_real("U real", 3, 3), Vt_real("Vt real", 2, 2);
  if constexpr (KAT_S::is_complex) {
    U_real(0, 0) = U_h(0, 0).real();
    U_real(0, 1) = U_h(0, 1).real();
    U_real(0, 2) = U_h(0, 2).real();
    U_real(1, 0) = U_h(1, 0).real();
    U_real(1, 1) = U_h(1, 1).real();
    U_real(1, 2) = U_h(1, 2).real();
    U_real(2, 0) = U_h(2, 0).real();
    U_real(2, 1) = U_h(2, 1).real();
    U_real(2, 2) = U_h(2, 2).real();

    Vt_real(0, 0) = Vt_h(0, 0).real();
    Vt_real(0, 1) = Vt_h(0, 1).real();
    Vt_real(1, 0) = Vt_h(1, 0).real();
    Vt_real(1, 1) = Vt_h(1, 1).real();
  } else {
    U_real(0, 0) = U_h(0, 0);
    U_real(0, 1) = U_h(0, 1);
    U_real(0, 2) = U_h(0, 2);
    U_real(1, 0) = U_h(1, 0);
    U_real(1, 1) = U_h(1, 1);
    U_real(1, 2) = U_h(1, 2);
    U_real(2, 0) = U_h(2, 0);
    U_real(2, 1) = U_h(2, 1);
    U_real(2, 2) = U_h(2, 2);

    Vt_real(0, 0) = Vt_h(0, 0);
    Vt_real(0, 1) = Vt_h(0, 1);
    Vt_real(1, 0) = Vt_h(1, 0);
    Vt_real(1, 1) = Vt_h(1, 1);
  }

  const mag_type one_sqrt2  = static_cast<mag_type>(1 / Kokkos::sqrt(2));
  const mag_type one_sqrt18 = static_cast<mag_type>(1 / Kokkos::sqrt(18));
  const mag_type one_third  = static_cast<mag_type>(1. / 3.);

  // Check values of U
  // Don't worry about the sign
  // it will be check with the
  // triple product
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(0, 0)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(0, 1)), one_sqrt18, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(0, 2)), 2 * one_third, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(1, 0)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(1, 1)), one_sqrt18, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(1, 2)), 2 * one_third, tol);
  EXPECT_NEAR_KK(Kokkos::abs(U_real(2, 0)), 0, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(2, 1)), 4 * one_sqrt18, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(U_real(2, 2)), one_third, tol);

  // Check values of Vt
  // Don't worry about the sign
  // it will be check with the
  // triple product
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(0, 0)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(0, 1)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(1, 0)), one_sqrt2, tol);
  EXPECT_NEAR_KK_REL(Kokkos::abs(Vt_real(1, 1)), one_sqrt2, tol);

  check_unitary_orthogonal_matrix(U, tol);
  check_unitary_orthogonal_matrix(Vt, tol);

  check_triple_product(Aref, S, U, Vt, tol);

  return 0;
}

template <class AMatrix, class Device>
int impl_test_svd(const int m, const int n) {
  using execution_space = typename Device::execution_space;
  using scalar_type     = typename AMatrix::value_type;
  using KAT_S           = Kokkos::ArithTraits<scalar_type>;
  using mag_type        = typename KAT_S::mag_type;
  using vector_type     = Kokkos::View<mag_type*, typename AMatrix::array_layout, Device>;

  const mag_type max_val = 10;
  const mag_type tol     = 2000 * max_val * KAT_S::eps();

  AMatrix A("A", m, n), U("U", m, m), Vt("Vt", n, n), Aref("A ref", m, n);
  vector_type S("S", Kokkos::min(m, n));

  const uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);

  // Initialize A with random numbers
  scalar_type randStart = 0, randEnd = 0;
  Test::getRandomBounds(max_val, randStart, randEnd);
  Kokkos::fill_random(A, rand_pool, randStart, randEnd);
  Kokkos::deep_copy(Aref, A);

  // Working around CUSOLVER constraint for m >= n
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)
  if constexpr (std::is_same_v<typename Device::execution_space, Kokkos::Cuda>) {
    if (m >= n) {
      KokkosLapack::svd("A", "A", A, S, U, Vt);
    } else {
      return 0;
    }
  } else {
    KokkosLapack::svd("A", "A", A, S, U, Vt);
  }
#else
  KokkosLapack::svd("A", "A", A, S, U, Vt);
#endif

  check_unitary_orthogonal_matrix(U, tol);
  check_unitary_orthogonal_matrix(Vt, tol);

  // For larger sizes with the triple product
  // we accumulate a bit more error apparently?
  check_triple_product(Aref, S, U, Vt, 100 * Kokkos::max(m, n) * tol);

  return 0;
}

}  // namespace Test

template <class ScalarA, class Device>
int test_svd() {
  int ret;

#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_layout_left = Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device>;

  ret = Test::impl_analytic_2x2_svd<view_type_a_layout_left, Device>();
  EXPECT_EQ(ret, 0);

  ret = Test::impl_analytic_2x3_svd<view_type_a_layout_left, Device>();
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_left, Device>(0, 0);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_left, Device>(1, 1);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_left, Device>(15, 15);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_left, Device>(100, 100);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_left, Device>(100, 70);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_left, Device>(70, 100);
  EXPECT_EQ(ret, 0);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_layout_right = Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device>;

  ret = Test::impl_analytic_2x2_svd<view_type_a_layout_right, Device>();
  EXPECT_EQ(ret, 0);

  ret = Test::impl_analytic_2x3_svd<view_type_a_layout_right, Device>();
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_right, Device>(0, 0);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_right, Device>(1, 1);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_right, Device>(15, 15);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_right, Device>(100, 100);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_right, Device>(100, 70);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_svd<view_type_a_layout_right, Device>(70, 100);
  EXPECT_EQ(ret, 0);
#endif

  return 1;
}

template <class Scalar, class Device>
int test_svd_wrapper() {
#if defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
  if constexpr (std::is_same_v<typename Device::memory_space, Kokkos::HostSpace>) {
    // Using a device side space with LAPACK/MKL
    return test_svd<Scalar, Device>();
  }
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)
  if constexpr (std::is_same_v<typename Device::execution_space, Kokkos::Cuda>) {
    // Using a Cuda device with CUSOLVER
    return test_svd<Scalar, Device>();
  }
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER)
  if constexpr (std::is_same_v<typename Device::execution_space, Kokkos::HIP>) {
    // Using a HIP device with ROCSOLVER
    return test_svd<Scalar, Device>();
  }
#endif

  std::cout << "No TPL support enabled, svd is not tested" << std::endl;
  return 0;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, svd_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::svd_float");
  test_svd_wrapper<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, svd_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::svd_double");
  test_svd_wrapper<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, svd_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::svd_complex_float");
  test_svd_wrapper<Kokkos::complex<float>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, svd_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::svd_complex_double");
  test_svd_wrapper<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
