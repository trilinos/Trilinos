// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \author Luc Berger-Vergiat (lberg@sandia.gov)
/// \author Cameron Smith (smithc11@rpi.edu)

#include "gtest/gtest.h"
#include <random>

#include <KokkosBatched_QR_Decl.hpp>      // KokkosBatched::QR
#include "KokkosBatched_ApplyQ_Decl.hpp"  // KokkosBatched::ApplyQ
#include "KokkosBatched_QR_FormQ_Serial_Internal.hpp"
#include <KokkosBatched_Util.hpp>  // KokkosBlas::Algo
#include "KokkosBatched_Gemm_Decl.hpp"
#include <Kokkos_Core.hpp>

template <class MatricesType, class TauViewType, class TmpViewType, class ErrorViewType>
struct qrFunctor {
  using Scalar = typename MatricesType::value_type;
  using KAT    = KokkosKernels::ArithTraits<Scalar>;

  MatricesType As;
  TauViewType taus;
  TmpViewType ws;
  MatricesType Qs, Qts, Is, Bs;
  ErrorViewType error;

  Scalar max_val;

  qrFunctor(MatricesType As_, TauViewType taus_, TmpViewType ws_, MatricesType Qs_, MatricesType Qts_, MatricesType Is_,
            MatricesType Bs_, ErrorViewType error_, Scalar max_val_)
      : As(As_), taus(taus_), ws(ws_), Qs(Qs_), Qts(Qts_), Is(Is_), Bs(Bs_), error(error_), max_val(max_val_) {}

  KOKKOS_FUNCTION
  void operator()(const int matIdx) const {
    const int max_rows_cols = Kokkos::max(As.extent_int(1), As.extent_int(2));

    auto A   = Kokkos::subview(As, matIdx, Kokkos::ALL, Kokkos::ALL);
    auto tau = Kokkos::subview(taus, matIdx, Kokkos::ALL);
    auto w   = Kokkos::subview(ws, Kokkos::pair<int, int>(matIdx * max_rows_cols, (matIdx + 1) * max_rows_cols));
    auto Q   = Kokkos::subview(Qs, matIdx, Kokkos::ALL, Kokkos::ALL);
    auto Qt  = Kokkos::subview(Qts, matIdx, Kokkos::ALL, Kokkos::ALL);
    auto I   = Kokkos::subview(Is, matIdx, Kokkos::ALL, Kokkos::ALL);
    auto B   = Kokkos::subview(Bs, matIdx, Kokkos::ALL, Kokkos::ALL);

    const Scalar SC_one              = KAT::one();
    const typename KAT::mag_type tol = KAT::eps() * max_rows_cols * max_rows_cols * Kokkos::abs(max_val);

    int error_lcl = 0;

    for (int idx = 0; idx < w.extent_int(0); ++idx) {
      w(idx) = 0.0;
    }
    KokkosBatched::SerialQR<KokkosBlas::Algo::QR::Unblocked>::invoke(A, tau, w);

    // Store identity in Q and Qt
    for (int diagIdx = 0; diagIdx < Q.extent_int(0); ++diagIdx) {
      Q(diagIdx, diagIdx)  = SC_one;
      Qt(diagIdx, diagIdx) = SC_one;
    }

    // Call ApplyQ on Q
    for (int idx = 0; idx < w.extent_int(0); ++idx) {
      w(idx) = 0.0;
    }
    KokkosBatched::SerialApplyQ<Side::Left, Trans::NoTranspose, Algo::ApplyQ::Unblocked>::invoke(A, tau, Q, w);

    // Copy Q into I
    for (int rowIdx = 0; rowIdx < Q.extent_int(0); ++rowIdx) {
      for (int colIdx = 0; colIdx < Q.extent_int(1); ++colIdx) {
        I(rowIdx, colIdx) = Q(rowIdx, colIdx);
      }
    }

    // Call ApplyQ with transpose mode on Qt
    for (int idx = 0; idx < w.extent_int(0); ++idx) {
      w(idx) = 0.0;
    }
    KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, tau, Qt, w);

    // Call ApplyQ with transpose mode on I
    for (int idx = 0; idx < w.extent_int(0); ++idx) {
      w(idx) = 0.0;
    }
    KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, tau, I, w);

    // At this point I stores Q'Q
    // which should be the identity matrix
    for (int rowIdx = 0; rowIdx < I.extent_int(0); ++rowIdx) {
      for (int colIdx = 0; colIdx < I.extent_int(1); ++colIdx) {
        if (Kokkos::abs(Q(rowIdx, colIdx) - KAT::conj(Qt(colIdx, rowIdx))) > tol) {
          ++error_lcl;
        }
        if (rowIdx == colIdx) {
          if (Kokkos::abs(I(rowIdx, colIdx) - SC_one) > tol) {
            error_lcl += 1;
          }
        } else {
          if (Kokkos::abs(I(rowIdx, colIdx)) > tol) {
            error_lcl += 1;
          }
        }
      }
    }

    // Store identity in Q
    for (int rowIdx = 0; rowIdx < Q.extent_int(0); ++rowIdx) {
      for (int colIdx = 0; colIdx < Q.extent_int(1); ++colIdx) {
        Q(rowIdx, colIdx) = (rowIdx == colIdx) ? SC_one : KokkosKernels::ArithTraits<Scalar>::zero();
      }
    }

    // Call ApplyQ on Q from the right side
    for (int idx = 0; idx < w.extent_int(0); ++idx) {
      w(idx) = 0.0;
    }
    KokkosBatched::SerialApplyQ<Side::Right, Trans::NoTranspose, Algo::ApplyQ::Unblocked>::invoke(A, tau, Q, w);
    for (int rowIdx = 0; rowIdx < I.extent_int(0); ++rowIdx) {
      for (int colIdx = 0; colIdx < I.extent_int(1); ++colIdx) {
        if (Kokkos::abs(Q(rowIdx, colIdx) - KAT::conj(Qt(colIdx, rowIdx))) > tol) {
          ++error_lcl;
        }
      }
    }

    // Apply Q' to B which holds a copy of the orginal A
    // Afterwards B should hold a copy of R and be zero below its diagonal
    for (int idx = 0; idx < w.extent_int(0); ++idx) {
      w(idx) = 0.0;
    }
    KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, tau, B, w);

    for (int rowIdx = 0; rowIdx < B.extent_int(0); ++rowIdx) {
      for (int colIdx = 0; colIdx < B.extent_int(1); ++colIdx) {
        if (rowIdx <= colIdx) {
          if (Kokkos::abs(B(rowIdx, colIdx) - A(rowIdx, colIdx)) > tol * Kokkos::abs(A(rowIdx, colIdx))) {
            error_lcl += 1;
          }
        } else {
          if (Kokkos::abs(B(rowIdx, colIdx)) > tol) {
            error_lcl += 1;
          }
        }
      }
    }
    error(matIdx) = error_lcl;
  }
};

template <class Device, class Scalar, class AlgoTagType>
void test_QR_square() {
  // Analytical test with a square matrix
  //     [12, -51,   4]              [150, -69,  58]        [14,  21, -14]
  // A = [ 6, 167, -68]    Q = 1/175 [ 75, 158,  -6]    R = [ 0, 175, -70]
  //     [-4,  24, -41]              [-50,  30, 165]        [ 0,   0, -35]
  //
  // Expected outputs:
  //                          [   14,  21, -14]
  // tau = [7/13, , 1/2]  A = [-6/26, 175, -70]
  //                          [ 4/26,    , -35]
  //

  using MatrixViewType    = Kokkos::View<Scalar**>;
  using ColVectorViewType = Kokkos::View<Scalar*>;
  using ColWorkViewType   = Kokkos::View<Scalar*>;

  const Scalar tol = 20 * KokkosKernels::ArithTraits<Scalar>::eps();
  constexpr int m = 3, n = 3;

  MatrixViewType A("A", m, n), B("B", m, n), Q("Q", m, m);
  ColVectorViewType t("t", n);
  ColWorkViewType w("w", n);

  auto A_h  = Kokkos::create_mirror_view(A);
  A_h(0, 0) = 12;
  A_h(0, 1) = -51;
  A_h(0, 2) = 4;
  A_h(1, 0) = 6;
  A_h(1, 1) = 167;
  A_h(1, 2) = -68;
  A_h(2, 0) = -4;
  A_h(2, 1) = 24;
  A_h(2, 2) = -41;

  Kokkos::deep_copy(A, A_h);
  Kokkos::deep_copy(B, A_h);

  Kokkos::parallel_for(
      "serialQR", 1, KOKKOS_LAMBDA(int) {
        // compute the QR factorization of A and store the results in A and t
        // (tau) - see the lapack dgeqp3(...) documentation:
        // www.netlib.org/lapack/explore-html-3.6.1/dd/d9a/group__double_g_ecomputational_ga1b0500f49e03d2771b797c6e88adabbb.html
        KokkosBatched::SerialQR<AlgoTagType>::invoke(A, t, w);
      });

  Kokkos::fence();
  Kokkos::deep_copy(A_h, A);
  auto tau_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, t);

  Test::EXPECT_NEAR_KK_REL(A_h(0, 0), static_cast<Scalar>(-14), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(0, 1), static_cast<Scalar>(-21), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(0, 2), static_cast<Scalar>(14), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(1, 0), static_cast<Scalar>(6. / 26.), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(1, 1), static_cast<Scalar>(-175), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(1, 2), static_cast<Scalar>(70), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(2, 0), static_cast<Scalar>(-4.0 / 26.0), tol);
  // Test::EXPECT_NEAR_KK_REL(A_h(2, 1),   35.0, tol);      // Analytical expression too painful to compute...
  Test::EXPECT_NEAR_KK_REL(A_h(2, 2), static_cast<Scalar>(35), tol);

  Test::EXPECT_NEAR_KK_REL(tau_h(0), static_cast<Scalar>(7. / 13.), tol);
  // Test::EXPECT_NEAR_KK_REL(tau_h(1), 25. / 32., tol);      // Analytical expression too painful to compute...
  Test::EXPECT_NEAR_KK_REL(tau_h(2), static_cast<Scalar>(1. / 2.), tol);

  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, t, B, w);
      });
  auto B_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, B);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-14), B_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-21), B_h(0, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(14), B_h(0, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0), B_h(1, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-175), B_h(1, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(70), B_h(1, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0), B_h(2, 0), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0), B_h(2, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(35), B_h(2, 2), tol);

  Kokkos::parallel_for(
      "serialFormQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialQR_FormQ_Internal::invoke(m, n, A.data(), A.stride(0), A.stride(1), t.data(), t.stride(0),
                                                       Q.data(), Q.stride(0), Q.stride(1), w.data(), w.stride(0));
      });
  auto Q_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Q);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-6. / 7.), Q_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(69. / 175.), Q_h(0, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-58. / 175.), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-3. / 7.), Q_h(1, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-158. / 175.), Q_h(1, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(6. / 175.), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(2. / 7.), Q_h(2, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-6. / 35.), Q_h(2, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-33. / 35.), Q_h(2, 2), tol);

  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, t, Q, w);
      });
  Kokkos::deep_copy(Q_h, Q);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(1.), Q_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(1.), Q_h(1, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(1.), Q_h(2, 2), tol);

  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(0, 1), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(1, 0), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(2, 0), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(2, 1), tol);

  Kokkos::deep_copy(Q_h, 0);
  Q_h(0, 0) = 1.0;
  Q_h(1, 1) = 1.0;
  Q_h(2, 2) = 1.0;
  Kokkos::deep_copy(Q, Q_h);
  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Right, Trans::NoTranspose, Algo::ApplyQ::Unblocked>::invoke(A, t, Q, w);
      });
  Kokkos::deep_copy(Q_h, Q);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-6. / 7.), Q_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(69. / 175.), Q_h(0, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-58. / 175.), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-3. / 7.), Q_h(1, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-158. / 175.), Q_h(1, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(6. / 175.), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(2. / 7.), Q_h(2, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-6. / 35.), Q_h(2, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-33. / 35.), Q_h(2, 2), tol);
}

template <class Device, class Scalar, class AlgoTagType>
void test_QR_rectangular() {
  // Analytical test with a rectangular matrix
  //     [3,  5]        [-0.60, -0.64,  0.48]        [-5, -3]
  // A = [4,  0]    Q = [-0.80, -0.48, -0.36]    R = [ 0,  5]
  //     [0, -3]        [ 0.00, -0.60, -0.80]        [ 0,  0]
  //
  // Expected outputs:
  //                         [ -5,  -3]
  // tau = [5/8, 10/18]  A = [1/2,   5]
  //                         [  0, 1/3]
  //

  using MatrixViewType    = Kokkos::View<Scalar**>;
  using ColVectorViewType = Kokkos::View<Scalar*>;
  using ColWorkViewType   = Kokkos::View<Scalar*>;

  const Scalar tol = 10 * KokkosKernels::ArithTraits<Scalar>::eps();
  constexpr int m = 3, n = 2;

  MatrixViewType A("A", m, n), B("B", m, n), Q("Q", m, m);
  ColVectorViewType t("t", n);
  ColWorkViewType w("w", Kokkos::max(m, n));

  auto A_h  = Kokkos::create_mirror_view(A);
  A_h(0, 0) = 3;
  A_h(0, 1) = 5;
  A_h(1, 0) = 4;
  A_h(1, 1) = 0;
  A_h(2, 0) = 0;
  A_h(2, 1) = -3;

  Kokkos::deep_copy(A, A_h);
  Kokkos::deep_copy(B, A_h);

  Kokkos::parallel_for(
      "serialQR", 1, KOKKOS_LAMBDA(int) {
        // compute the QR factorization of A and store the results in A and t
        // (tau) - see the lapack dgeqp3(...) documentation:
        // www.netlib.org/lapack/explore-html-3.6.1/dd/d9a/group__double_g_ecomputational_ga1b0500f49e03d2771b797c6e88adabbb.html
        KokkosBatched::SerialQR<AlgoTagType>::invoke(A, t, w);
      });

  Kokkos::fence();
  Kokkos::deep_copy(A_h, A);
  auto tau_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, t);

  Test::EXPECT_NEAR_KK_REL(A_h(0, 0), static_cast<Scalar>(-5.0), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(0, 1), static_cast<Scalar>(-3.0), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(1, 0), static_cast<Scalar>(0.5), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(1, 1), static_cast<Scalar>(5.0), tol);
  Test::EXPECT_NEAR_KK(A_h(2, 0), static_cast<Scalar>(0.), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(2, 1), static_cast<Scalar>(1. / 3.), tol);

  Test::EXPECT_NEAR_KK_REL(tau_h(0), static_cast<Scalar>(5. / 8.), tol);
  Test::EXPECT_NEAR_KK_REL(tau_h(1), static_cast<Scalar>(10. / 18.), tol);

  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, t, B, w);
      });
  auto B_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, B);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-5.0), B_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-3.0), B_h(0, 1), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), B_h(1, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(5.0), B_h(1, 1), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), B_h(2, 0), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), B_h(2, 1), tol);

  Kokkos::parallel_for(
      "serialFormQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialQR_FormQ_Internal::invoke(m, n, A.data(), A.stride(0), A.stride(1), t.data(), t.stride(0),
                                                       Q.data(), Q.stride(0), Q.stride(1), w.data(), w.stride(0));
      });
  auto Q_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Q);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.64), Q_h(0, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.48), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.80), Q_h(1, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.48), Q_h(1, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.36), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(2, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(2, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.80), Q_h(2, 2), tol);

  Kokkos::deep_copy(Q_h, 0.0);
  Q_h(0, 0) = 1.0;
  Q_h(1, 1) = 1.0;
  Q_h(2, 2) = 1.0;
  Kokkos::deep_copy(Q, Q_h);
  Kokkos::deep_copy(w, 0.0);
  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, t, Q, w);
      });
  Kokkos::deep_copy(Q_h, Q);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.64), Q_h(1, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.48), Q_h(2, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.80), Q_h(0, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.48), Q_h(1, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.36), Q_h(2, 1), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.80), Q_h(2, 2), tol);

  Kokkos::deep_copy(Q_h, 0.0);
  Q_h(0, 0) = 1.0;
  Q_h(1, 1) = 1.0;
  Q_h(2, 2) = 1.0;
  Kokkos::deep_copy(Q, Q_h);
  Kokkos::deep_copy(w, 0.0);
  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Left, Trans::NoTranspose, Algo::ApplyQ::Unblocked>::invoke(A, t, Q, w);
      });
  Kokkos::deep_copy(Q_h, Q);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.64), Q_h(0, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.48), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.80), Q_h(1, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.48), Q_h(1, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.36), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(2, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(2, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.80), Q_h(2, 2), tol);

  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, t, Q, w);
      });
  Kokkos::deep_copy(Q_h, Q);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(1.0), Q_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(1.0), Q_h(1, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(1.0), Q_h(2, 2), tol);

  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(0, 1), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(1, 0), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(2, 0), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(2, 1), tol);

  Kokkos::deep_copy(Q_h, 0);
  Q_h(0, 0) = 1.0;
  Q_h(1, 1) = 1.0;
  Q_h(2, 2) = 1.0;
  Kokkos::deep_copy(Q, Q_h);
  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Right, Trans::NoTranspose, Algo::ApplyQ::Unblocked>::invoke(A, t, Q, w);
      });
  Kokkos::deep_copy(Q_h, Q);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.64), Q_h(0, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.48), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.80), Q_h(1, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.48), Q_h(1, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.36), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(2, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(2, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.80), Q_h(2, 2), tol);
}

template <class Device, class Scalar, class AlgoTagType>
void test_QR_rectangular_cplx() {
  // Analytical test with a complex rectangular matrix
  //     [3+3i,  15]        [-0.60, -0.64,  0.48]        [-6, (-22 + 21i) / 3]
  // A = [4+i,    0]    Q = [-0.80, -0.48, -0.36]    R = [ 0, sqrt(170**2 + 68**2 + 80**2 + 54**2) / 18]
  //     [0+i, -3-i]        [ 0.00, -0.60, -0.80]        [ 0,           0]
  //
  // Expected outputs:
  //      [           -6,  (-22 + 21i) / 3]
  //  A = [(13 - i) / 13,  norm_y]
  //      [(1 + 3i) / 30,  (-80 - 54i) / (-170 - 18*norm_y +68i)]
  //
  // tau = [(3 + i)/2, norm_y / (norm_y - (-170 + 68i) / 18)]
  //
  // with norm_y = sqrt(170**2 + 68**2 + 80**2 + 54**2) / 18

  using KAT      = KokkosKernels::ArithTraits<Scalar>;
  using mag_type = typename KAT::mag_type;

  using MatrixViewType    = Kokkos::View<Scalar**>;
  using ColVectorViewType = Kokkos::View<Scalar*>;
  using ColWorkViewType   = Kokkos::View<Scalar*>;

  const mag_type tol = 10 * KAT::eps();
  constexpr int m = 3, n = 2;

  MatrixViewType A("A", m, n), B("B", m, n), Q("Q", m, m);
  ColVectorViewType t("t", n);
  ColWorkViewType w("w", Kokkos::max(m, n));

  auto A_h  = Kokkos::create_mirror_view(A);
  A_h(0, 0) = Scalar(3, 3);
  A_h(0, 1) = Scalar(15, 0);
  A_h(1, 0) = Scalar(4, 1);
  A_h(1, 1) = Scalar(0, 0);
  A_h(2, 0) = Scalar(0, 1);
  A_h(2, 1) = Scalar(-3, -1);

  Kokkos::deep_copy(A, A_h);
  Kokkos::deep_copy(B, A_h);

  Kokkos::parallel_for(
      "serialQR", 1, KOKKOS_LAMBDA(int) {
        // compute the QR factorization of A and store the results in A and t
        // (tau) - see the lapack dgeqp3(...) documentation:
        // www.netlib.org/lapack/explore-html-3.6.1/dd/d9a/group__double_g_ecomputational_ga1b0500f49e03d2771b797c6e88adabbb.html
        KokkosBatched::SerialQR<AlgoTagType>::invoke(A, t, w);
      });

  Kokkos::fence();
  Kokkos::deep_copy(A_h, A);
  auto tau_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, t);

  const mag_type norm_y = Kokkos::sqrt(170 * 170 + 68 * 68 + 80 * 80 + 54 * 54) / 18;

  Test::EXPECT_NEAR_KK_REL(A_h(0, 0), Scalar(-6, 0), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(0, 1), Scalar(-22. / 3, 21. / 3), tol);
  Test::EXPECT_NEAR_KK_REL(A_h(1, 0), Scalar(13, -1) / 30, tol);
  Test::EXPECT_NEAR_KK_REL(A_h(1, 1), norm_y, tol);
  Test::EXPECT_NEAR_KK_REL(A_h(2, 0), Scalar(1, 3) / 30, tol);
  Test::EXPECT_NEAR_KK_REL(A_h(2, 1), Scalar(-80, -54) / Scalar(-170 - 18 * norm_y, 68), tol);

  Test::EXPECT_NEAR_KK_REL(tau_h(0), KAT::one() / Scalar(3. / 2., 1. / 2.), tol);
  Test::EXPECT_NEAR_KK_REL(tau_h(1), norm_y / (norm_y - Scalar(-170, 68) / 18), tol);

  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, t, B, w);
      });
  auto B_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, B);

  Test::EXPECT_NEAR_KK_REL(Scalar(-6, 0), B_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(Scalar(-22, 21) / 3, B_h(0, 1), tol);
  Test::EXPECT_NEAR_KK(Scalar(0, 0), B_h(1, 0), tol);
  Test::EXPECT_NEAR_KK_REL(Scalar(norm_y, 0), B_h(1, 1), tol);
  Test::EXPECT_NEAR_KK(Scalar(0, 0), B_h(2, 0), tol);
  Test::EXPECT_NEAR_KK(Scalar(0, 0), B_h(2, 1), tol);

  Kokkos::parallel_for(
      "serialFormQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialQR_FormQ_Internal::invoke(m, n, A.data(), A.stride(0), A.stride(1), t.data(), t.stride(0),
                                                       Q.data(), Q.stride(0), Q.stride(1), w.data(), w.stride(0));
      });
  auto Q_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Q);

  Test::EXPECT_NEAR_KK_REL(Scalar(-0.5, -0.5), Q_h(0, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.64), Q_h(0, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.48), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK_REL(Scalar(-2.0 / 3.0, -1.0 / 6.0), Q_h(1, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.48), Q_h(1, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.36), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK_REL(Scalar(0.0, -1.0 / 6.0), Q_h(2, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(2, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.80), Q_h(2, 2), tol);

  // Kokkos::deep_copy(Q_h, 0.0);
  // Q_h(0, 0) = 1.0;
  // Q_h(1, 1) = 1.0;
  // Q_h(2, 2) = 1.0;
  // Kokkos::deep_copy(Q, Q_h);
  // Kokkos::deep_copy(w, 0.0);
  // Kokkos::parallel_for(
  //     "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
  //       KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, t, Q, w);
  //     });
  // Kokkos::deep_copy(Q_h, Q);

  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(0, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.64), Q_h(1, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.48), Q_h(2, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.80), Q_h(0, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.48), Q_h(1, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.36), Q_h(2, 1), tol);
  // Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(0, 2), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(1, 2), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.80), Q_h(2, 2), tol);

  Kokkos::deep_copy(Q_h, 0.0);
  Q_h(0, 0) = 1.0;
  Q_h(1, 1) = 1.0;
  Q_h(2, 2) = 1.0;
  Kokkos::deep_copy(Q, Q_h);
  Kokkos::deep_copy(w, 0.0);
  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Left, Trans::NoTranspose, Algo::ApplyQ::Unblocked>::invoke(A, t, Q, w);
      });
  Kokkos::deep_copy(Q_h, Q);

  Test::EXPECT_NEAR_KK_REL(Scalar(-0.5, -0.5), Q_h(0, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.64), Q_h(0, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.48), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK_REL(Scalar(-2.0 / 3.0, -1.0 / 6.0), Q_h(1, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.48), Q_h(1, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.36), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK_REL(Scalar(0.0, -1.0 / 6.0), Q_h(2, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(2, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.80), Q_h(2, 2), tol);

  Kokkos::parallel_for(
      "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
        KokkosBatched::SerialApplyQ<Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked>::invoke(A, t, Q, w);
      });
  Kokkos::deep_copy(Q_h, Q);

  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(1.0), Q_h(0, 0), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(1.0), Q_h(1, 1), tol);
  Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(1.0), Q_h(2, 2), tol);

  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(0, 1), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(0, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(1, 0), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(1, 2), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(2, 0), tol);
  Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(2, 1), tol);

  // Kokkos::deep_copy(Q_h, 0);
  // Q_h(0, 0) = 1.0;
  // Q_h(1, 1) = 1.0;
  // Q_h(2, 2) = 1.0;
  // Kokkos::deep_copy(Q, Q_h);
  // Kokkos::parallel_for(
  //     "serialApplyQ", 1, KOKKOS_LAMBDA(int) {
  //       KokkosBatched::SerialApplyQ<Side::Right, Trans::NoTranspose, Algo::ApplyQ::Unblocked>::invoke(A, t, Q, w);
  //     });
  // Kokkos::deep_copy(Q_h, Q);

  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(0, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.64), Q_h(0, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.48), Q_h(0, 2), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.80), Q_h(1, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.48), Q_h(1, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.36), Q_h(1, 2), tol);
  // Test::EXPECT_NEAR_KK(static_cast<Scalar>(0.), Q_h(2, 0), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(-0.60), Q_h(2, 1), tol);
  // Test::EXPECT_NEAR_KK_REL(static_cast<Scalar>(0.80), Q_h(2, 2), tol);
}

template <class Device, class Scalar, class AlgoTagType>
void test_QR_batch(const int numMat, const int numRows, const int numCols) {
  // Generate a batch of matrices
  // Compute QR factorization
  // Verify that R is triangular
  // Verify that Q is unitary
  // Check that Q*R = A

  using ExecutionSpace = typename Device::execution_space;
  using size_type      = typename Kokkos::View<Scalar*, ExecutionSpace>::size_type;

  {
    Kokkos::View<Scalar**, ExecutionSpace> tau("tau", numMat, numCols);
    Kokkos::View<Scalar*, ExecutionSpace> tmp("work buffer",
                                              static_cast<size_type>(numMat) * Kokkos::max(numRows, numCols));
    Kokkos::View<Scalar***, ExecutionSpace> As("A matrices", numMat, numRows, numCols);
    Kokkos::View<Scalar***, ExecutionSpace> Bs("B matrices", numMat, numRows, numCols);
    Kokkos::View<Scalar***, ExecutionSpace> Qs("Q matrices", numMat, numRows, numRows);
    Kokkos::View<Scalar***, ExecutionSpace> Qts("Q transpose matrices", numMat, numRows, numRows);
    Kokkos::View<Scalar***, ExecutionSpace> Is("Identity matrices", numMat, numRows, numRows);
    Kokkos::View<int*, ExecutionSpace> error("global number of error", numMat);

    Kokkos::Random_XorShift64_Pool<ExecutionSpace> rand_pool(2718);
    constexpr double max_val = 1000;
    {
      Scalar randStart, randEnd;
      Test::getRandomBounds(max_val, randStart, randEnd);
      Kokkos::fill_random(ExecutionSpace{}, As, rand_pool, randStart, randEnd);
    }
    Kokkos::fence();
    Kokkos::deep_copy(Bs, As);

    auto As_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, As);

    qrFunctor myFunc(As, tau, tmp, Qs, Qts, Is, Bs, error, max_val);
    Kokkos::parallel_for("KokkosBatched::test_QR_batch", Kokkos::RangePolicy<ExecutionSpace>(0, numMat), myFunc);
    Kokkos::fence();

    auto error_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, error);

    int global_error = 0;
    for (int matIdx = 0; matIdx < numMat; ++matIdx) {
      global_error += error_h(matIdx);
    }
    EXPECT_EQ(global_error, 0);
  }
}

template <class Device, class Scalar, class AlgoTagType>
void test_QR_batch(const int numMat, const int numRows) {
  test_QR_batch<Device, Scalar, AlgoTagType>(numMat, numRows, numRows);
}

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, serial_qr_square_analytic_float) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_square<TestDevice, float, AlgoTagType>();
}
TEST_F(TestCategory, serial_qr_rectangular_analytic_float) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_rectangular<TestDevice, float, AlgoTagType>();
}
TEST_F(TestCategory, serial_qr_batch_float) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_batch<TestDevice, float, AlgoTagType>(314, 36);
  test_QR_batch<TestDevice, float, AlgoTagType>(10, 42, 36);
  test_QR_batch<TestDevice, float, AlgoTagType>(100, 42, 36);
  test_QR_batch<TestDevice, float, AlgoTagType>(200, 42, 36);
  test_QR_batch<TestDevice, float, AlgoTagType>(250, 42, 36);
  test_QR_batch<TestDevice, float, AlgoTagType>(300, 42, 36);
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, serial_qr_square_analytic_double) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_square<TestDevice, double, AlgoTagType>();
}
TEST_F(TestCategory, serial_qr_rectangular_analytic_double) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_rectangular<TestDevice, double, AlgoTagType>();
}
TEST_F(TestCategory, serial_qr_batch_double) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_batch<TestDevice, double, AlgoTagType>(314, 36);
  test_QR_batch<TestDevice, double, AlgoTagType>(10, 42, 36);
  test_QR_batch<TestDevice, double, AlgoTagType>(100, 42, 36);
  test_QR_batch<TestDevice, double, AlgoTagType>(200, 42, 36);
  test_QR_batch<TestDevice, double, AlgoTagType>(250, 42, 36);
  test_QR_batch<TestDevice, double, AlgoTagType>(300, 42, 36);
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT)
TEST_F(TestCategory, serial_qr_square_analytic_scomplex) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_square<TestDevice, Kokkos::complex<float>, AlgoTagType>();
}
TEST_F(TestCategory, serial_qr_rectangular_analytic_scomplex) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_rectangular<TestDevice, Kokkos::complex<float>, AlgoTagType>();
}
TEST_F(TestCategory, serial_qr_batch_scomplex) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_batch<TestDevice, Kokkos::complex<float>, AlgoTagType>(314, 36);
  test_QR_batch<TestDevice, Kokkos::complex<float>, AlgoTagType>(10, 42, 36);
  test_QR_batch<TestDevice, Kokkos::complex<float>, AlgoTagType>(100, 42, 36);
  test_QR_batch<TestDevice, Kokkos::complex<float>, AlgoTagType>(200, 42, 36);
  test_QR_batch<TestDevice, Kokkos::complex<float>, AlgoTagType>(250, 42, 36);
  test_QR_batch<TestDevice, Kokkos::complex<float>, AlgoTagType>(300, 42, 36);
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE)
TEST_F(TestCategory, serial_qr_square_analytic_dcomplex) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_square<TestDevice, Kokkos::complex<double>, AlgoTagType>();
}
TEST_F(TestCategory, serial_qr_rectangular_analytic_dcomplex) {
  using AlgoTagType = KokkosBlas::Algo::QR::Unblocked;
  test_QR_rectangular_cplx<TestDevice, Kokkos::complex<double>, AlgoTagType>();
}
TEST_F(TestCategory, serial_qr_batch_dcomplex) {
  typedef KokkosBlas::Algo::QR::Unblocked AlgoTagType;
  test_QR_batch<TestDevice, Kokkos::complex<double>, AlgoTagType>(314, 36);
  test_QR_batch<TestDevice, Kokkos::complex<double>, AlgoTagType>(1, 4, 3);
  test_QR_batch<TestDevice, Kokkos::complex<double>, AlgoTagType>(100, 42, 36);
  test_QR_batch<TestDevice, Kokkos::complex<double>, AlgoTagType>(200, 42, 36);
  test_QR_batch<TestDevice, Kokkos::complex<double>, AlgoTagType>(250, 42, 36);
  test_QR_batch<TestDevice, Kokkos::complex<double>, AlgoTagType>(300, 42, 36);
}
#endif
