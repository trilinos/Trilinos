// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_TEST_DENSE_LINEAR_ALGEBRA_HPP__
#define __TACHO_TEST_DENSE_LINEAR_ALGEBRA_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <Kokkos_Random.hpp>
#include <Kokkos_Timer.hpp>

#include "Tacho_Blas_External.hpp"
#include "Tacho_Blas_Team.hpp"
#include "Tacho_Util.hpp"

#include "Tacho_Lapack_External.hpp"
#include "Tacho_Lapack_Team.hpp"

using namespace Tacho;

typedef Kokkos::DualView<value_type **, Kokkos::LayoutLeft, device_type> matrix_type;

namespace Test {
struct Functor_TeamGemm {
  char _transa, _transb;
  int _m, _n, _k;
  matrix_type _A, _B, _C;
  value_type _alpha, _beta;

  Functor_TeamGemm(const char transa, const char transb, const int m, const int n, const int k, const value_type alpha,
                   const matrix_type &A, const matrix_type &B, const value_type beta, const matrix_type &C)
      : _transa(transa), _transb(transb), _m(m), _n(n), _k(k), _A(A), _B(B), _C(C), _alpha(alpha), _beta(beta) {}

  template <typename MemberType> KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    ::BlasTeam<value_type>::gemm(member, _transa, _transb, _m, _n, _k, _alpha, (const value_type *)_A.d_view.data(),
                                 (int)_A.d_view.stride_1(), (const value_type *)_B.d_view.data(),
                                 (int)_B.d_view.stride_1(), _beta, (value_type *)_C.d_view.data(),
                                 (int)_C.d_view.stride_1());
  }

  inline void run() {
    _A.sync_device();
    _B.sync_device();

    _C.sync_device();
    _C.modify_device();

    Kokkos::parallel_for(Kokkos::TeamPolicy<typename device_type::execution_space>(1, Kokkos::AUTO), *this);

    _C.sync_host();
  }
};
} // namespace Test

TEST(DenseLinearAlgebra, team_gemm_nn) {

  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'N', transb = 'N';
  const value_type alpha = 1.3, beta = 2.5;

  // test problem setup
  matrix_type A1("A1", m, k), B1("B1", k, n), C1("C1", m, n);
  matrix_type A2("A2", m, k), B2("B2", k, n), C2("C2", m, n);

  A1.modify_device();
  B1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(B1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  // tacho test
  ::Test::Functor_TeamGemm test(transa, transb, m, n, k, alpha, A1, B1, beta, C1);
  test.run();

  // reference test
  Blas<value_type>::gemm(transa, transb, m, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), B2.h_view.data(),
                         B2.h_view.stride_1(), beta, C2.h_view.data(), C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_gemm_nt) {

  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'N', transb = 'T';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", m, k), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", m, k), B2("B2", n, k), C2("C2", m, n);

  // test problem setup
  A1.modify_device();
  B1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(B1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  // tacho test
  ::Test::Functor_TeamGemm test(transa, transb, m, n, k, alpha, A1, B1, beta, C1);
  test.run();

  // reference test
  Blas<value_type>::gemm(transa, transb, m, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), B2.h_view.data(),
                         B2.h_view.stride_1(), beta, C2.h_view.data(), C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_gemm_nc) {

  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'N', transb = 'C';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", m, k), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", m, k), B2("B2", n, k), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(B1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamGemm test(transa, transb, m, n, k, alpha, A1, B1, beta, C1);
  test.run();

  Blas<value_type>::gemm(transa, transb, m, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), B2.h_view.data(),
                         B2.h_view.stride_1(), beta, C2.h_view.data(), C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_gemm_tn) {

  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'T', transb = 'N';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", k, n), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", k, n), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(B1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamGemm test(transa, transb, m, n, k, alpha, A1, B1, beta, C1);
  test.run();

  Blas<value_type>::gemm(transa, transb, m, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), B2.h_view.data(),
                         B2.h_view.stride_1(), beta, C2.h_view.data(), C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_gemm_tt) {

  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'T', transb = 'T';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  A1.modify_device();
  B1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(B1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamGemm test(transa, transb, m, n, k, alpha, A1, B1, beta, C1);
  test.run();

  Blas<value_type>::gemm(transa, transb, m, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), B2.h_view.data(),
                         B2.h_view.stride_1(), beta, C2.h_view.data(), C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_gemm_tc) {

  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'T', transb = 'C';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  A1.modify_device();
  B1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(B1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamGemm test(transa, transb, m, n, k, alpha, A1, B1, beta, C1);
  test.run();

  Blas<value_type>::gemm(transa, transb, m, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), B2.h_view.data(),
                         B2.h_view.stride_1(), beta, C2.h_view.data(), C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_gemm_cn) {

  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'C', transb = 'N';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", k, n), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", k, n), C2("C2", m, n);

  A1.modify_device();
  B1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(B1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamGemm test(transa, transb, m, n, k, alpha, A1, B1, beta, C1);
  test.run();

  Blas<value_type>::gemm(transa, transb, m, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), B2.h_view.data(),
                         B2.h_view.stride_1(), beta, C2.h_view.data(), C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_gemm_ct) {

  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'C', transb = 'T';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  A1.modify_device();
  B1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(B1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamGemm test(transa, transb, m, n, k, alpha, A1, B1, beta, C1);
  test.run();

  Blas<value_type>::gemm(transa, transb, m, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), B2.h_view.data(),
                         B2.h_view.stride_1(), beta, C2.h_view.data(), C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_gemm_cc) {

  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'C', transb = 'C';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  A1.modify_device();
  B1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(B1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamGemm test(transa, transb, m, n, k, alpha, A1, B1, beta, C1);
  test.run();

  Blas<value_type>::gemm(transa, transb, m, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), B2.h_view.data(),
                         B2.h_view.stride_1(), beta, C2.h_view.data(), C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

namespace Test {
struct Functor_TeamGemv {
  char _trans;
  int _m, _n;
  matrix_type _A, _x, _y;
  value_type _alpha, _beta;

  Functor_TeamGemv(const char trans, const int m, const int n, const value_type alpha, const matrix_type &A,
                   const matrix_type &x, const value_type beta, const matrix_type &y)
      : _trans(trans), _m(m), _n(n), _A(A), _x(x), _y(y), _alpha(alpha), _beta(beta) {}

  template <typename MemberType> KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    ::BlasTeam<value_type>::gemv(member, _trans, _m, _n, _alpha, (const value_type *)_A.d_view.data(),
                                 (int)_A.d_view.stride_1(), (const value_type *)_x.d_view.data(),
                                 (int)_x.d_view.stride_0(), _beta, (value_type *)_y.d_view.data(),
                                 (int)_y.d_view.stride_0());
  }

  inline void run() {
    _A.sync_device();
    _x.sync_device();

    _y.sync_device();
    _y.modify_device();

    Kokkos::parallel_for(Kokkos::TeamPolicy<typename device_type::execution_space>(1, Kokkos::AUTO), *this);

    _y.sync_host();
  }
};
} // namespace Test

TEST(DenseLinearAlgebra, team_gemv_n) {

  const ordinal_type m = 20, n = 10;
  const char trans = 'N';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", m, n), x1("B1", n, 1), y1("C1", m, 1);
  matrix_type A2("A2", m, n), x2("B2", n, 1), y2("C2", m, 1);

  A1.modify_device();
  x1.modify_device();
  y1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(x1.d_view, random, value_type(1));
  Kokkos::fill_random(y1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(x2.h_view, x1.d_view);
  Kokkos::deep_copy(y2.h_view, y1.d_view);

  ::Test::Functor_TeamGemv test(trans, m, n, alpha, A1, x1, beta, y1);
  test.run();

  Blas<value_type>::gemv(trans, m, n, alpha, A2.h_view.data(), A2.h_view.stride_1(), x2.h_view.data(),
                         x2.h_view.stride_0(), beta, y2.h_view.data(), y2.h_view.stride_0());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    EXPECT_NEAR(ats::abs(y1.h_view(i, 0)), ats::abs(y2.h_view(i, 0)), eps);
}

TEST(DenseLinearAlgebra, team_gemv_t) {

  const ordinal_type m = 20, n = 10;
  const char trans = 'T';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", m, n), x1("x1", m, 1), y1("y1", n, 1);
  matrix_type A2("A2", m, n), x2("x2", m, 1), y2("y2", n, 1);

  A1.modify_device();
  x1.modify_device();
  y1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(x1.d_view, random, value_type(1));
  Kokkos::fill_random(y1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(x2.h_view, x1.d_view);
  Kokkos::deep_copy(y2.h_view, y1.d_view);

  ::Test::Functor_TeamGemv test(trans, m, n, alpha, A1, x1, beta, y1);
  test.run();

  Blas<value_type>::gemv(trans, m, n, alpha, A2.h_view.data(), A2.h_view.stride_1(), x2.h_view.data(),
                         x2.h_view.stride_0(), beta, y2.h_view.data(), y2.h_view.stride_0());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < n; ++i)
    EXPECT_NEAR(ats::abs(y1.h_view(i, 0)), ats::abs(y2.h_view(i, 0)), eps);
}

TEST(DenseLinearAlgebra, team_gemv_c) {

  const ordinal_type m = 20, n = 10;
  const char trans = 'C';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", m, n), x1("x1", m, 1), y1("y1", n, 1);
  matrix_type A2("A2", m, n), x2("x2", m, 1), y2("y2", n, 1);

  A1.modify_device();
  x1.modify_device();
  y1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(x1.d_view, random, value_type(1));
  Kokkos::fill_random(y1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(x2.h_view, x1.d_view);
  Kokkos::deep_copy(y2.h_view, y1.d_view);

  ::Test::Functor_TeamGemv test(trans, m, n, alpha, A1, x1, beta, y1);
  test.run();

  Blas<value_type>::gemv(trans, m, n, alpha, A2.h_view.data(), A2.h_view.stride_1(), x2.h_view.data(),
                         x2.h_view.stride_0(), beta, y2.h_view.data(), y2.h_view.stride_0());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < n; ++i)
    EXPECT_NEAR(ats::abs(y1.h_view(i, 0)), ats::abs(y2.h_view(i, 0)), eps);
}

namespace Test {
struct Functor_TeamHerk {
  char _uplo, _trans;
  int _n, _k;
  matrix_type _A, _C;
  value_type _alpha, _beta;

  Functor_TeamHerk(const char uplo, const char trans, const int n, const int k, const value_type alpha,
                   const matrix_type &A, const value_type beta, const matrix_type &C)
      : _uplo(uplo), _trans(trans), _n(n), _k(k), _A(A), _C(C), _alpha(alpha), _beta(beta) {}

  template <typename MemberType> KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    ::BlasTeam<value_type>::herk(member, _uplo, _trans, _n, _k, _alpha, (const value_type *)_A.d_view.data(),
                                 (int)_A.d_view.stride_1(), _beta, (value_type *)_C.d_view.data(),
                                 (int)_C.d_view.stride_1());
  }

  inline void run() {
    _A.sync_device();

    _C.sync_device();
    _C.modify_device();

    Kokkos::parallel_for(Kokkos::TeamPolicy<typename device_type::execution_space>(1, Kokkos::AUTO), *this);

    _C.sync_host();
  }
};
} // namespace Test

TEST(DenseLinearAlgebra, team_herk_un) {

  const ordinal_type n = 20, k = 10;
  const char uplo = 'U';
  const char trans = 'N';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", n, k), C1("C1", n, n);
  matrix_type A2("A2", n, k), C2("C2", n, n);

  A1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamHerk test(uplo, trans, n, k, alpha, A1, beta, C1);
  test.run();

  Blas<value_type>::herk(uplo, trans, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), beta, C2.h_view.data(),
                         C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_herk_uc) {

  const ordinal_type n = 20, k = 10;
  const char uplo = 'U';
  const char trans = 'C';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, n), C1("C1", n, n);
  matrix_type A2("A2", k, n), C2("C2", n, n);

  A1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamHerk test(uplo, trans, n, k, alpha, A1, beta, C1);
  test.run();

  Blas<value_type>::herk(uplo, trans, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), beta, C2.h_view.data(),
                         C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_herk_ln) {

  const ordinal_type n = 20, k = 10;
  const char uplo = 'L';
  const char trans = 'N';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", n, k), C1("C1", n, n);
  matrix_type A2("A2", n, k), C2("C2", n, n);

  A1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamHerk test(uplo, trans, n, k, alpha, A1, beta, C1);
  test.run();

  Blas<value_type>::herk(uplo, trans, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), beta, C2.h_view.data(),
                         C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

TEST(DenseLinearAlgebra, team_herk_lc) {

  const ordinal_type n = 20, k = 10;
  const char uplo = 'L';
  const char trans = 'C';
  const value_type alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, n), C1("C1", n, n);
  matrix_type A2("A2", k, n), C2("C2", n, n);

  A1.modify_device();
  C1.modify_device();

  Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);

  Kokkos::fill_random(A1.d_view, random, value_type(1));
  Kokkos::fill_random(C1.d_view, random, value_type(1));

  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamHerk test(uplo, trans, n, k, alpha, A1, beta, C1);
  test.run();

  Blas<value_type>::herk(uplo, trans, n, k, alpha, A2.h_view.data(), A2.h_view.stride_1(), beta, C2.h_view.data(),
                         C2.h_view.stride_1());

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      EXPECT_NEAR(ats::abs(C1.h_view(i, j)), ats::abs(C2.h_view(i, j)), eps);
}

namespace Test {
struct Functor_TeamTrsv {
  char _uplo, _trans, _diag;
  int _m;
  matrix_type _A, _b;

  Functor_TeamTrsv(const char uplo, const char trans, const char diag, const int m, const matrix_type &A,
                   const matrix_type &b)
      : _uplo(uplo), _trans(trans), _diag(diag), _m(m), _A(A), _b(b) {}

  template <typename MemberType> KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    ::BlasTeam<value_type>::trsv(member, _uplo, _trans, _diag, _m, (const value_type *)_A.d_view.data(),
                                 (int)_A.d_view.stride_1(), (value_type *)_b.d_view.data(), (int)_b.d_view.stride_0());
  }

  inline void run() {
    _A.sync_device();

    _b.sync_device();
    _b.modify_device();

    Kokkos::parallel_for(Kokkos::TeamPolicy<typename device_type::execution_space>(1, Kokkos::AUTO), *this);

    _b.sync_host();
  }
};
} // namespace Test

#define TEAM_TRSV_TEST_BODY                                                                                            \
  {                                                                                                                    \
    matrix_type A1("A1", m, m), b1("b1", m, 1);                                                                        \
    matrix_type A2("A2", m, m), b2("b2", m, 1);                                                                        \
                                                                                                                       \
    A1.modify_device();                                                                                                \
    b1.modify_device();                                                                                                \
                                                                                                                       \
    Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);                               \
                                                                                                                       \
    Kokkos::fill_random(A1.d_view, random, value_type(1));                                                             \
    Kokkos::fill_random(b1.d_view, random, value_type(1));                                                             \
                                                                                                                       \
    Kokkos::deep_copy(A2.h_view, A1.d_view);                                                                           \
    Kokkos::deep_copy(b2.h_view, b1.d_view);                                                                           \
                                                                                                                       \
    ::Test::Functor_TeamTrsv test(uplo, trans, diag, m, A1, b1);                                                       \
    test.run();                                                                                                        \
                                                                                                                       \
    Blas<value_type>::trsv(uplo, trans, diag, m, A2.h_view.data(), A2.h_view.stride_1(), b2.h_view.data(),             \
                           b2.h_view.stride_0());                                                                      \
                                                                                                                       \
    const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 10000;                                 \
    for (int i = 0; i < m; ++i)                                                                                        \
      EXPECT_NEAR(ats::abs(b1.h_view(i, 0)), ats::abs(b2.h_view(i, 0)), eps *ats::abs(b2.h_view(i, 0)));               \
  }

TEST(DenseLinearAlgebra, team_trsv_unu) {

  const ordinal_type m = 4;
  const char uplo = 'U';
  const char trans = 'N';
  const char diag = 'U';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_unn) {

  const ordinal_type m = 20;
  const char uplo = 'U';
  const char trans = 'N';
  const char diag = 'N';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_utu) {

  const ordinal_type m = 20;
  const char uplo = 'U';
  const char trans = 'T';
  const char diag = 'U';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_utn) {

  const ordinal_type m = 20;
  const char uplo = 'U';
  const char trans = 'T';
  const char diag = 'N';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_ucu) {

  const ordinal_type m = 20;
  const char uplo = 'U';
  const char trans = 'C';
  const char diag = 'U';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_ucn) {

  const ordinal_type m = 20;
  const char uplo = 'U';
  const char trans = 'C';
  const char diag = 'N';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_lnu) {

  const ordinal_type m = 20;
  const char uplo = 'L';
  const char trans = 'N';
  const char diag = 'U';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_lnn) {

  const ordinal_type m = 20;
  const char uplo = 'L';
  const char trans = 'N';
  const char diag = 'N';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_ltu) {

  const ordinal_type m = 20;
  const char uplo = 'L';
  const char trans = 'T';
  const char diag = 'U';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_ltn) {

  const ordinal_type m = 20;
  const char uplo = 'L';
  const char trans = 'T';
  const char diag = 'N';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_lcu) {

  const ordinal_type m = 20;
  const char uplo = 'L';
  const char trans = 'C';
  const char diag = 'U';

  TEAM_TRSV_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsv_lcn) {

  const ordinal_type m = 20;
  const char uplo = 'L';
  const char trans = 'C';
  const char diag = 'N';

  TEAM_TRSV_TEST_BODY;
}
#undef TEAM_TRSV_TEST_BODY

namespace Test {
struct Functor_TeamTrsm {
  char _side, _uplo, _trans, _diag;
  int _m, _n;
  matrix_type _A, _B;
  value_type _alpha;

  Functor_TeamTrsm(const char side, const char uplo, const char trans, const char diag, const int m, const int n,
                   const value_type alpha, const matrix_type &A, const matrix_type &B)
      : _side(side), _uplo(uplo), _trans(trans), _diag(diag), _m(m), _n(n), _A(A), _B(B), _alpha(alpha) {}

  template <typename MemberType> KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    ::BlasTeam<value_type>::trsm(member, _side, _uplo, _trans, _diag, _m, _n, _alpha,
                                 (const value_type *)_A.d_view.data(), (int)_A.d_view.stride_1(),
                                 (value_type *)_B.d_view.data(), (int)_B.d_view.stride_1());
  }

  inline void run() {
    _A.sync_device();

    _B.sync_device();
    _B.modify_device();

    Kokkos::parallel_for(Kokkos::TeamPolicy<typename device_type::execution_space>(1, Kokkos::AUTO), *this);

    _B.sync_host();
  }
};
} // namespace Test

#define TEAM_TRSM_TEST_BODY                                                                                            \
  {                                                                                                                    \
    matrix_type A1("A1", m, m), B1("b1", m, n);                                                                        \
    matrix_type A2("A2", m, m), B2("b2", m, n);                                                                        \
                                                                                                                       \
    A1.modify_device();                                                                                                \
    B1.modify_device();                                                                                                \
                                                                                                                       \
    Kokkos::Random_XorShift64_Pool<typename device_type::execution_space> random(13718);                               \
                                                                                                                       \
    Kokkos::fill_random(A1.d_view, random, value_type(1));                                                             \
    Kokkos::fill_random(B1.d_view, random, value_type(1));                                                             \
                                                                                                                       \
    Kokkos::deep_copy(A2.h_view, A1.d_view);                                                                           \
    Kokkos::deep_copy(B2.h_view, B1.d_view);                                                                           \
                                                                                                                       \
    ::Test::Functor_TeamTrsm test(side, uplo, trans, diag, m, n, alpha, A1, B1);                                       \
    test.run();                                                                                                        \
                                                                                                                       \
    Blas<value_type>::trsm(side, uplo, trans, diag, m, n, alpha, A2.h_view.data(), A2.h_view.stride_1(),               \
                           B2.h_view.data(), B2.h_view.stride_1());                                                    \
                                                                                                                       \
    const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 100000;                                \
    for (int i = 0; i < m; ++i)                                                                                        \
      for (int j = 0; j < n; ++j)                                                                                      \
        EXPECT_NEAR(ats::abs(B1.h_view(i, j)), ats::abs(B2.h_view(i, j)), eps *ats::abs(B2.h_view(i, j)));             \
  }

TEST(DenseLinearAlgebra, team_trsm_lunu) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'U';
  const char trans = 'N';
  const char diag = 'U';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_lunn) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'U';
  const char trans = 'N';
  const char diag = 'N';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_lutu) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'U';
  const char trans = 'T';
  const char diag = 'U';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_lutn) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'U';
  const char trans = 'T';
  const char diag = 'N';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_lucu) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'U';
  const char trans = 'C';
  const char diag = 'U';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_lucn) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'U';
  const char trans = 'C';
  const char diag = 'N';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_llnu) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'L';
  const char trans = 'N';
  const char diag = 'U';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_llnn) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'L';
  const char trans = 'N';
  const char diag = 'N';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_lltu) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'L';
  const char trans = 'T';
  const char diag = 'U';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_lltn) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'L';
  const char trans = 'T';
  const char diag = 'N';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_llcu) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'L';
  const char trans = 'C';
  const char diag = 'U';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}

TEST(DenseLinearAlgebra, team_trsm_llcn) {

  const ordinal_type m = 20, n = 10;
  const char side = 'L';
  const char uplo = 'L';
  const char trans = 'C';
  const char diag = 'N';
  const value_type alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}
#undef TEAM_TRSM_TEST_BODY

namespace Test {
struct Functor_TeamChol {
  char _uplo;
  int _m;
  matrix_type _A;

  Functor_TeamChol(const char uplo, const int m, const matrix_type &A) : _uplo(uplo), _m(m), _A(A) {}

  template <typename MemberType> KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    int r_val = 0;
    ::LapackTeam<value_type>::potrf(member, _uplo, _m, (value_type *)_A.d_view.data(), (int)_A.d_view.stride_1(),
                                    &r_val);
  }

  inline void run() {
    _A.sync_device();
    _A.modify_device();

    Kokkos::parallel_for(Kokkos::TeamPolicy<typename device_type::execution_space>(1, Kokkos::AUTO), *this);

    _A.sync_host();
  }
};
} // namespace Test

TEST(DenseLinearAlgebra, team_chol_u) {

  const ordinal_type m = 20;
  const char uplo = 'U';

  matrix_type A1("A1", m, m);
  matrix_type A2("A2", m, m);

  A2.modify_host();

  for (int i = 0; i < m; ++i)
    A2.h_view(i, i) = 4.0;
  for (int i = 0; i < (m - 1); ++i) {
    A2.h_view(i, i + 1) = -1.0;
    A2.h_view(i + 1, i) = -1.0;
  }

  A1.modify_device();
  Kokkos::deep_copy(A1.d_view, A2.h_view);

  ::Test::Functor_TeamChol test(uplo, m, A1);
  test.run();

  int r_val = 0;
  Lapack<value_type>::potrf(uplo, m, (value_type *)A2.h_view.data(), (int)A2.h_view.stride_1(), &r_val);

  const magnitude_type eps = std::numeric_limits<magnitude_type>::epsilon() * 1000;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < m; ++j)
      EXPECT_NEAR(ats::abs(A1.h_view(i, j)), ats::abs(A2.h_view(i, j)), eps);
}

#endif
