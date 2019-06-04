#ifndef __TACHO_TEST_DENSE_LINEAR_ALGEBRA_HPP__
#define __TACHO_TEST_DENSE_LINEAR_ALGEBRA_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <Kokkos_Random.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_Blas_Team.hpp"
#include "Tacho_Blas_External.hpp"

#include "Tacho_Lapack_Team.hpp"
#include "Tacho_Lapack_External.hpp"

using namespace Tacho;
using std::abs;
using Kokkos::abs;

typedef Kokkos::DualView<ValueType**,Kokkos::LayoutLeft,DeviceSpaceType> matrix_type;

namespace Test {
  struct Functor_TeamGemm {
    char _transa, _transb;
    int _m, _n, _k;
    matrix_type _A, _B, _C;
    ValueType _alpha, _beta;

    Functor_TeamGemm(const char transa, const char transb,
                     const int m, const int n, const int k, 
                     const ValueType alpha, 
                     const matrix_type &A,
                     const matrix_type &B,
                     const ValueType beta,
                     const matrix_type &C) 
      : _transa(transa), _transb(transb), 
        _m(m), _n(n), _k(k), 
        _A(A), _B(B), _C(C), 
        _alpha(alpha), _beta(beta) {}
    
    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &member) const {
      ::BlasTeam<ValueType>::gemm(member,
                                  _transa, _transb,
                                  _m, _n, _k,
                                  _alpha, 
                                  (const ValueType*)_A.d_view.data(), (int)_A.d_view.stride_1(),
                                  (const ValueType*)_B.d_view.data(), (int)_B.d_view.stride_1(),
                                  _beta,
                                  (      ValueType*)_C.d_view.data(), (int)_C.d_view.stride_1());
    }
    
    inline
    void run() {
      _A.sync<DeviceSpaceType>();
      _B.sync<DeviceSpaceType>();

      _C.sync<DeviceSpaceType>();
      _C.modify<DeviceSpaceType>();

      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);

      _C.sync<HostSpaceType>();
    }
  };
}


TEST( DenseLinearAlgebra, team_gemm_nn ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'N', transb = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  // test problem setup
  matrix_type A1("A1", m, k), B1("B1", k, n), C1("C1", m, n);
  matrix_type A2("A2", m, k), B2("B2", k, n), C2("C2", m, n);

  A1.modify<DeviceSpaceType>();
  B1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);

  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(B1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);
  
  // tacho test
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();

  // reference test 
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        B2.h_view.data(), B2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());
  
  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_gemm_nt ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'N', transb = 'T';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", m, k), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", m, k), B2("B2", n, k), C2("C2", m, n);

  // test problem setup
  A1.modify<DeviceSpaceType>();
  B1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(B1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  // tacho test 
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  // reference test 
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        B2.h_view.data(), B2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_gemm_nc ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'N', transb = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", m, k), B1("B1", n, k), C1("C1", m, n);
  matrix_type   A2("A2", m, k), B2("B2", n, k), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(B1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        B2.h_view.data(), B2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());
  
  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_gemm_tn ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'T', transb = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", k, n), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", k, n), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(B1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        B2.h_view.data(), B2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_gemm_tt ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'T', transb = 'T';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  A1.modify<DeviceSpaceType>();
  B1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(B1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        B2.h_view.data(), B2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_gemm_tc ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'T', transb = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  A1.modify<DeviceSpaceType>();
  B1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();
  
  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(B1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        B2.h_view.data(), B2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_gemm_cn ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'C', transb = 'N';
  const ValueType alpha = 1.3, beta = 2.5;
  
  matrix_type A1("A1", k, m), B1("B1", k, n), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", k, n), C2("C2", m, n);

  A1.modify<DeviceSpaceType>();
  B1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(B1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        B2.h_view.data(), B2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_gemm_ct ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'C', transb = 'T';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  A1.modify<DeviceSpaceType>();
  B1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(B1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        B2.h_view.data(), B2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_gemm_cc ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'C', transb = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  A1.modify<DeviceSpaceType>();
  B1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(B1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(B2.h_view, B1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        B2.h_view.data(), B2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}

namespace Test {
  struct Functor_TeamGemv {
    char _trans;
    int _m, _n;
    matrix_type _A, _x, _y;
    ValueType _alpha, _beta;

    Functor_TeamGemv(const char trans,
                     const int m, const int n, 
                     const ValueType alpha, 
                     const matrix_type &A,
                     const matrix_type &x,
                     const ValueType beta,
                     const matrix_type &y) 
      : _trans(trans), 
        _m(m), _n(n),
        _A(A), _x(x), _y(y), 
        _alpha(alpha), _beta(beta) {}
    
    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &member) const {
      ::BlasTeam<ValueType>::gemv(member,
                                  _trans,
                                  _m, _n,
                                  _alpha, 
                                  (const ValueType*)_A.d_view.data(), (int)_A.d_view.stride_1(),
                                  (const ValueType*)_x.d_view.data(), (int)_x.d_view.stride_0(),
                                  _beta,
                                  (      ValueType*)_y.d_view.data(), (int)_y.d_view.stride_0());
    }
    
    inline
    void run() {
      _A.sync<DeviceSpaceType>();
      _x.sync<DeviceSpaceType>();
      
      _y.sync<DeviceSpaceType>();
      _y.modify<DeviceSpaceType>();

      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);

      _y.sync<HostSpaceType>();
    }
  };
}


TEST( DenseLinearAlgebra, team_gemv_n ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char trans = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", m, n), x1("B1", n, 1), y1("C1", m, 1);
  matrix_type A2("A2", m, n), x2("B2", n, 1), y2("C2", m, 1);

  A1.modify<DeviceSpaceType>();
  x1.modify<DeviceSpaceType>();
  y1.modify<DeviceSpaceType>();

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(x1.d_view, random, ValueType(1));
  Kokkos::fill_random(y1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(x2.h_view, x1.d_view);
  Kokkos::deep_copy(y2.h_view, y1.d_view);
  
  ::Test::Functor_TeamGemv test(trans,
                                m, n, 
                                alpha, A1, x1, beta, y1);
  test.run();
  
  Blas<ValueType>::gemv(trans,
                        m, n, 
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        x2.h_view.data(), x2.h_view.stride_0(),
                        beta,
                        y2.h_view.data(), y2.h_view.stride_0());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    EXPECT_NEAR(::abs(y1.h_view(i,0)), ::abs(y2.h_view(i,0)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_gemv_t ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char trans = 'T';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", m, n), x1("x1", m, 1), y1("y1", n, 1);
  matrix_type A2("A2", m, n), x2("x2", m, 1), y2("y2", n, 1);

  A1.modify<DeviceSpaceType>();
  x1.modify<DeviceSpaceType>();
  y1.modify<DeviceSpaceType>();

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(x1.d_view, random, ValueType(1));
  Kokkos::fill_random(y1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(x2.h_view, x1.d_view);
  Kokkos::deep_copy(y2.h_view, y1.d_view);

  ::Test::Functor_TeamGemv test(trans,
                                m, n, 
                                alpha, A1, x1, beta, y1);
  test.run();
  
  Blas<ValueType>::gemv(trans,
                        m, n, 
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        x2.h_view.data(), x2.h_view.stride_0(),
                        beta,
                        y2.h_view.data(), y2.h_view.stride_0());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    EXPECT_NEAR(::abs(y1.h_view(i,0)), ::abs(y2.h_view(i,0)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_gemv_c ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char trans = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", m, n), x1("x1", m, 1), y1("y1", n, 1);
  matrix_type A2("A2", m, n), x2("x2", m, 1), y2("y2", n, 1);

  A1.modify<DeviceSpaceType>();
  x1.modify<DeviceSpaceType>();
  y1.modify<DeviceSpaceType>();

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(x1.d_view, random, ValueType(1));
  Kokkos::fill_random(y1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(x2.h_view, x1.d_view);
  Kokkos::deep_copy(y2.h_view, y1.d_view);

  ::Test::Functor_TeamGemv test(trans,
                                m, n, 
                                alpha, A1, x1, beta, y1);
  test.run();

  Blas<ValueType>::gemv(trans,
                        m, n, 
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        x2.h_view.data(), x2.h_view.stride_0(),
                        beta,
                        y2.h_view.data(), y2.h_view.stride_0());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    EXPECT_NEAR(::abs(y1.h_view(i,0)), ::abs(y2.h_view(i,0)), eps);
  TEST_END;
}

namespace Test {
  struct Functor_TeamHerk {
    char _uplo, _trans;
    int _n, _k;
    matrix_type _A, _C;
    ValueType _alpha, _beta;

    Functor_TeamHerk(const char uplo, const char trans,
                     const int n, const int k, 
                     const ValueType alpha, 
                     const matrix_type &A,
                     const ValueType beta,
                     const matrix_type &C) 
      : _uplo(uplo), _trans(trans), 
        _n(n), _k(k),
        _A(A), _C(C),
        _alpha(alpha), _beta(beta) {}
    
    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &member) const {
      ::BlasTeam<ValueType>::herk(member,
                                  _uplo, _trans,
                                  _n, _k,
                                  _alpha, 
                                  (const ValueType*)_A.d_view.data(), (int)_A.d_view.stride_1(),
                                  _beta,
                                  (      ValueType*)_C.d_view.data(), (int)_C.d_view.stride_1());
    }
    
    inline
    void run() {
      _A.sync<DeviceSpaceType>();

      _C.sync<DeviceSpaceType>();
      _C.modify<DeviceSpaceType>();

      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);

      _C.sync<HostSpaceType>();
    }
  };
}


TEST( DenseLinearAlgebra, team_herk_un ) {
  TEST_BEGIN;
  const ordinal_type n = 20, k = 10;
  const char uplo  = 'U';
  const char trans = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", n, k), C1("C1", n, n);
  matrix_type A2("A2", n, k), C2("C2", n, n);

  A1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();
  
  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamHerk test(uplo, trans,
                                n, k, 
                                alpha, A1, beta, C1);
  test.run();

  Blas<ValueType>::herk(uplo, trans,
                        n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());
  
  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_herk_uc ) {
  TEST_BEGIN;
  const ordinal_type n = 20, k = 10;
  const char uplo  = 'U';
  const char trans = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, n), C1("C1", n, n);
  matrix_type A2("A2", k, n), C2("C2", n, n);

  A1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();
  
  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamHerk test(uplo, trans,
                                n, k, 
                                alpha, A1, beta, C1);
  test.run();

  Blas<ValueType>::herk(uplo, trans,
                        n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());
  
  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_herk_ln ) {
  TEST_BEGIN;
  const ordinal_type n = 20, k = 10;
  const char uplo  = 'L';
  const char trans = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", n, k), C1("C1", n, n);
  matrix_type A2("A2", n, k), C2("C2", n, n);

  A1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();
  
  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamHerk test(uplo, trans,
                                n, k, 
                                alpha, A1, beta, C1);
  test.run();

  Blas<ValueType>::herk(uplo, trans,
                        n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


TEST( DenseLinearAlgebra, team_herk_lc ) {
  TEST_BEGIN;
  const ordinal_type n = 20, k = 10;
  const char uplo  = 'L';
  const char trans = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type A1("A1", k, n), C1("C1", n, n);
  matrix_type A2("A2", k, n), C2("C2", n, n);

  A1.modify<DeviceSpaceType>();
  C1.modify<DeviceSpaceType>();
  
  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1.d_view, random, ValueType(1));
  Kokkos::fill_random(C1.d_view, random, ValueType(1));
  
  Kokkos::deep_copy(A2.h_view, A1.d_view);
  Kokkos::deep_copy(C2.h_view, C1.d_view);

  ::Test::Functor_TeamHerk test(uplo, trans,
                                n, k, 
                                alpha, A1, beta, C1);
  test.run();

  Blas<ValueType>::herk(uplo, trans,
                        n, k,
                        alpha, 
                        A2.h_view.data(), A2.h_view.stride_1(),
                        beta,
                        C2.h_view.data(), C2.h_view.stride_1());

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1.h_view(i,j)), ::abs(C2.h_view(i,j)), eps);
  TEST_END;
}


namespace Test {
  struct Functor_TeamTrsv {
    char _uplo, _trans, _diag;
    int _m;
    matrix_type _A, _b;

    Functor_TeamTrsv(const char uplo, const char trans, const char diag,
                     const int m,
                     const matrix_type &A,
                     const matrix_type &b) 
      : _uplo(uplo), _trans(trans), _diag(diag),
        _m(m),
        _A(A), _b(b) {}

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &member) const {
      ::BlasTeam<ValueType>::trsv(member,
                                  _uplo, _trans, _diag,
                                  _m,
                                  (const ValueType*)_A.d_view.data(), (int)_A.d_view.stride_1(),
                                  (      ValueType*)_b.d_view.data(), (int)_b.d_view.stride_0());
    }
    
    inline
    void run() {
      _A.sync<DeviceSpaceType>();

      _b.sync<DeviceSpaceType>();
      _b.modify<DeviceSpaceType>();
      
      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);
      
      _b.sync<HostSpaceType>();
    }
  };
}


#define TEAM_TRSV_TEST_BODY {                                           \
    matrix_type A1("A1", m, m), b1("b1", m, 1);                         \
    matrix_type A2("A2", m, m), b2("b2", m, 1);                         \
                                                                        \
    A1.modify<DeviceSpaceType>();                                       \
    b1.modify<DeviceSpaceType>();                                       \
                                                                        \
    Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);      \
                                                                        \
    Kokkos::fill_random(A1.d_view, random, ValueType(1));               \
    Kokkos::fill_random(b1.d_view, random, ValueType(1));               \
                                                                        \
    Kokkos::deep_copy(A2.h_view, A1.d_view);                            \
    Kokkos::deep_copy(b2.h_view, b1.d_view);                            \
                                                                        \
    ::Test::Functor_TeamTrsv test(uplo, trans, diag,                    \
                                  m,                                    \
                                  A1, b1);                              \
    test.run();                                                         \
                                                                        \
    Blas<ValueType>::trsv(uplo, trans, diag,                            \
                          m,                                            \
                          A2.h_view.data(), A2.h_view.stride_1(),       \
                          b2.h_view.data(), b2.h_view.stride_0());      \
                                                                        \
    const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 10000; \
    for (int i=0;i<m;++i)                                               \
      EXPECT_NEAR(::abs(b1.h_view(i,0)), ::abs(b2.h_view(i,0)), eps*::abs(b2.h_view(i,0))); \
  }

TEST( DenseLinearAlgebra, team_trsv_unu ) {
  TEST_BEGIN;
  const ordinal_type m = 4;
  const char uplo  = 'U';
  const char trans = 'N';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsv_unn ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'N';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsv_utu ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'T';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}

TEST( DenseLinearAlgebra, team_trsv_utn ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'T';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}

TEST( DenseLinearAlgebra, team_trsv_ucu ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'C';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}

TEST( DenseLinearAlgebra, team_trsv_ucn ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'C';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}





TEST( DenseLinearAlgebra, team_trsv_lnu ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'N';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsv_lnn ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'N';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsv_ltu ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'T';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}

TEST( DenseLinearAlgebra, team_trsv_ltn ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'T';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}

TEST( DenseLinearAlgebra, team_trsv_lcu ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'C';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}

TEST( DenseLinearAlgebra, team_trsv_lcn ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'C';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
  TEST_END;
}
#undef TEAM_TRSV_TEST_BODY


namespace Test {
  struct Functor_TeamTrsm {
    char _side, _uplo, _trans, _diag;
    int _m, _n;
    matrix_type _A, _B;
    ValueType _alpha;

    Functor_TeamTrsm(const char side, const char uplo, const char trans, const char diag,
                     const int m, const int n,
                     const ValueType alpha,
                     const matrix_type &A,
                     const matrix_type &B) 
      : _side(side), _uplo(uplo), _trans(trans), _diag(diag),
        _m(m), _n(n),
        _A(A), _B(B),
        _alpha(alpha) {}
    
    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &member) const {
      ::BlasTeam<ValueType>::trsm(member,
                                  _side, _uplo, _trans, _diag,
                                  _m, _n,
                                  _alpha,
                                  (const ValueType*)_A.d_view.data(), (int)_A.d_view.stride_1(),
                                  (      ValueType*)_B.d_view.data(), (int)_B.d_view.stride_1());
    }
    
    inline
    void run() {
      _A.sync<DeviceSpaceType>();

      _B.sync<DeviceSpaceType>();
      _B.modify<DeviceSpaceType>();
      
      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);

      _B.sync<HostSpaceType>();
    }
  };
}

#define TEAM_TRSM_TEST_BODY {                                           \
    matrix_type A1("A1", m, m), B1("b1", m, n);                         \
    matrix_type A2("A2", m, m), B2("b2", m, n);                         \
                                                                        \
    A1.modify<DeviceSpaceType>();                                       \
    B1.modify<DeviceSpaceType>();                                       \
                                                                        \
    Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);      \
                                                                        \
    Kokkos::fill_random(A1.d_view, random, ValueType(1));               \
    Kokkos::fill_random(B1.d_view, random, ValueType(1));               \
                                                                        \
    Kokkos::deep_copy(A2.h_view, A1.d_view);                            \
    Kokkos::deep_copy(B2.h_view, B1.d_view);                            \
                                                                        \
    ::Test::Functor_TeamTrsm test(side, uplo, trans, diag,              \
                                  m, n,                                 \
                                  alpha,                                \
                                  A1, B1);                              \
    test.run();                                                         \
                                                                        \
    Blas<ValueType>::trsm(side, uplo, trans, diag,                      \
                          m, n,                                         \
                          alpha,                                        \
                          A2.h_view.data(), A2.h_view.stride_1(),       \
                          B2.h_view.data(), B2.h_view.stride_1());      \
                                                                        \
    const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 100000; \
    for (int i=0;i<m;++i)                                               \
      for (int j=0;j<n;++j)                                             \
        EXPECT_NEAR(::abs(B1.h_view(i,j)), ::abs(B2.h_view(i,j)), eps*::abs(B2.h_view(i,j))); \
  }

TEST( DenseLinearAlgebra, team_trsm_lunu ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'N';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_lunn ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'N';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_lutu ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'T';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_lutn ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'T';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_lucu ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'C';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_lucn ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'C';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_llnu ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'N';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_llnn ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'N';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_lltu ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'T';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_lltn ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'T';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_llcu ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'C';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}


TEST( DenseLinearAlgebra, team_trsm_llcn ) {
  TEST_BEGIN;
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'C';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
  TEST_END;
}
#undef TEAM_TRSM_TEST_BODY

namespace Test {
  struct Functor_TeamChol {
    char _uplo;
    int _m;
    matrix_type _A;

    Functor_TeamChol(const char uplo,
                     const int m, 
                     const matrix_type &A) 
      : _uplo(uplo),
        _m(m),
        _A(A) {}
    
    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &member) const {
      int r_val = 0;
      ::LapackTeam<ValueType>::potrf(member,
                                     _uplo,
                                     _m,
                                     (ValueType*)_A.d_view.data(), (int)_A.d_view.stride_1(),
                                     &r_val);
    }
    
    inline
    void run() {
      _A.sync<DeviceSpaceType>();
      _A.modify<DeviceSpaceType>();
      
      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);

      _A.sync<HostSpaceType>();
    }
  };
}


TEST( DenseLinearAlgebra, team_chol_u ) {
  TEST_BEGIN;
  const ordinal_type m = 20;
  const char uplo = 'U';
  
  matrix_type A1("A1", m, m);
  matrix_type A2("A2", m, m);

  A2.modify<HostSpaceType>();
  
  for (int i=0;i<m;++i) 
    A2.h_view(i,i) = 4.0;
  for (int i=0;i<(m-1);++i) {
    A2.h_view(i,i+1) = -1.0;
    A2.h_view(i+1,i) = -1.0;
  }

  A1.modify<DeviceSpaceType>();  
  Kokkos::deep_copy(A1.d_view, A2.h_view);
  
  ::Test::Functor_TeamChol test(uplo,
                                m, 
                                A1);
  test.run();
  
  int r_val = 0;
  Lapack<ValueType>::potrf(uplo,
                           m, 
                           (ValueType*)A2.h_view.data(), (int)A2.h_view.stride_1(),
                           &r_val);
  
  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<m;++j) 
      EXPECT_NEAR(::abs(A1.h_view(i,j)), ::abs(A2.h_view(i,j)), eps);
  TEST_END;
}




#endif
