#ifndef __TACHO_TEST_DENSE_LINEAR_ALGEBRA_HPP__
#define __TACHO_TEST_DENSE_LINEAR_ALGEBRA_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_Blas_Team.hpp"
#include "TachoExp_Blas_External.hpp"

using namespace Tacho::Experimental;
using std::abs;
using Kokkos::abs;

typedef Kokkos::View<ValueType**,Kokkos::LayoutLeft,DeviceSpaceType> matrix_type_device;
typedef typename matrix_type_device::HostMirror matrix_type_host;

namespace Test {
  struct Functor_TeamGemm {
    char _transa, _transb;
    int _m, _n, _k;
    matrix_type_device _A, _B, _C;
    ValueType _alpha, _beta;

    Functor_TeamGemm(const char transa, const char transb,
                     const int m, const int n, const int k, 
                     const ValueType alpha, 
                     const matrix_type_device &A,
                     const matrix_type_device &B,
                     const ValueType beta,
                     const matrix_type_device &C) 
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
                                  (const ValueType*)_A.data(), (int)_A.stride_1(),
                                  (const ValueType*)_B.data(), (int)_B.stride_1(),
                                  _beta,
                                  (      ValueType*)_C.data(), (int)_C.stride_1());
    }
    
    inline
    void run() {
      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);
    }
  };
}


TEST( DenseLinearAlgebra, team_gemm_nn ) {
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'N', transb = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", m, k), B1("B1", k, n), C1("C1", m, n);
  matrix_type_host   A2("A2", m, k), B2("B2", k, n), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(B1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(B2, B1);
  Kokkos::deep_copy(C2, C1);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        B2.data(), B2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_gemm_nt ) {
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'N', transb = 'T';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", m, k), B1("B1", n, k), C1("C1", m, n);
  matrix_type_host   A2("A2", m, k), B2("B2", n, k), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(B1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(B2, B1);
  Kokkos::deep_copy(C2, C1);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        B2.data(), B2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_gemm_nc ) {
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'N', transb = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", m, k), B1("B1", n, k), C1("C1", m, n);
  matrix_type_host   A2("A2", m, k), B2("B2", n, k), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(B1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(B2, B1);
  Kokkos::deep_copy(C2, C1);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        B2.data(), B2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_gemm_tn ) {
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'T', transb = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", k, m), B1("B1", k, n), C1("C1", m, n);
  matrix_type_host   A2("A2", k, m), B2("B2", k, n), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(B1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(B2, B1);
  Kokkos::deep_copy(C2, C1);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        B2.data(), B2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_gemm_tt ) {
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'T', transb = 'T';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type_host   A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(B1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(B2, B1);
  Kokkos::deep_copy(C2, C1);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        B2.data(), B2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_gemm_tc ) {
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'T', transb = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type_host   A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(B1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(B2, B1);
  Kokkos::deep_copy(C2, C1);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        B2.data(), B2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_gemm_cn ) {
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'C', transb = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", k, m), B1("B1", k, n), C1("C1", m, n);
  matrix_type_host   A2("A2", k, m), B2("B2", k, n), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(B1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(B2, B1);
  Kokkos::deep_copy(C2, C1);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        B2.data(), B2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_gemm_ct ) {
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'C', transb = 'T';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type_host   A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(B1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(B2, B1);
  Kokkos::deep_copy(C2, C1);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        B2.data(), B2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_gemm_cc ) {
  const ordinal_type m = 20, n = 10, k = 15;
  const char transa = 'C', transb = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", k, m), B1("B1", n, k), C1("C1", m, n);
  matrix_type_host   A2("A2", k, m), B2("B2", n, k), C2("C2", m, n);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(B1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(B2, B1);
  Kokkos::deep_copy(C2, C1);
  
  ::Test::Functor_TeamGemm test(transa, transb, 
                                m, n, k,
                                alpha, A1, B1, beta, C1);
  test.run();
  
  Blas<ValueType>::gemm(transa, transb, 
                        m, n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        B2.data(), B2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}

namespace Test {
  struct Functor_TeamGemv {
    char _trans;
    int _m, _n;
    matrix_type_device _A, _x, _y;
    ValueType _alpha, _beta;

    Functor_TeamGemv(const char trans,
                     const int m, const int n, 
                     const ValueType alpha, 
                     const matrix_type_device &A,
                     const matrix_type_device &x,
                     const ValueType beta,
                     const matrix_type_device &y) 
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
                                  (const ValueType*)_A.data(), (int)_A.stride_1(),
                                  (const ValueType*)_x.data(), (int)_x.stride_0(),
                                  _beta,
                                  (      ValueType*)_y.data(), (int)_y.stride_0());
    }
    
    inline
    void run() {
      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);
    }
  };
}


TEST( DenseLinearAlgebra, team_gemv_n ) {
  const ordinal_type m = 20, n = 10;
  const char trans = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", m, n), x1("B1", n, 1), y1("C1", m, 1);
  matrix_type_host   A2("A2", m, n), x2("B2", n, 1), y2("C2", m, 1);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(x1, random, ValueType(1));
  Kokkos::fill_random(y1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(x2, x1);
  Kokkos::deep_copy(y2, y1);
  
  ::Test::Functor_TeamGemv test(trans,
                                m, n, 
                                alpha, A1, x1, beta, y1);
  test.run();
  
  Blas<ValueType>::gemv(trans,
                        m, n, 
                        alpha, 
                        A2.data(), A2.stride_1(),
                        x2.data(), x2.stride_0(),
                        beta,
                        y2.data(), y2.stride_0());
  
  auto y1_host = Kokkos::create_mirror_view(y1);
  Kokkos::deep_copy(y1_host, y1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<m;++i)
    EXPECT_NEAR(::abs(y1_host(i,0)), ::abs(y2(i,0)), eps);
}


TEST( DenseLinearAlgebra, team_gemv_t ) {
  const ordinal_type m = 20, n = 10;
  const char trans = 'T';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", m, n), x1("x1", m, 1), y1("y1", n, 1);
  matrix_type_host   A2("A2", m, n), x2("x2", m, 1), y2("y2", n, 1);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(x1, random, ValueType(1));
  Kokkos::fill_random(y1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(x2, x1);
  Kokkos::deep_copy(y2, y1);

  ::Test::Functor_TeamGemv test(trans,
                                m, n, 
                                alpha, A1, x1, beta, y1);
  test.run();
  
  Blas<ValueType>::gemv(trans,
                        m, n, 
                        alpha, 
                        A2.data(), A2.stride_1(),
                        x2.data(), x2.stride_0(),
                        beta,
                        y2.data(), y2.stride_0());
  
  auto y1_host = Kokkos::create_mirror_view(y1);
  Kokkos::deep_copy(y1_host, y1);
  
  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    EXPECT_NEAR(::abs(y1_host(i,0)), ::abs(y2(i,0)), eps);
}


TEST( DenseLinearAlgebra, team_gemv_c ) {
  const ordinal_type m = 20, n = 10;
  const char trans = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", m, n), x1("x1", m, 1), y1("y1", n, 1);
  matrix_type_host   A2("A2", m, n), x2("x2", m, 1), y2("y2", n, 1);

  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(x1, random, ValueType(1));
  Kokkos::fill_random(y1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(x2, x1);
  Kokkos::deep_copy(y2, y1);

  ::Test::Functor_TeamGemv test(trans,
                                m, n, 
                                alpha, A1, x1, beta, y1);
  test.run();

  Blas<ValueType>::gemv(trans,
                        m, n, 
                        alpha, 
                        A2.data(), A2.stride_1(),
                        x2.data(), x2.stride_0(),
                        beta,
                        y2.data(), y2.stride_0());
  
  auto y1_host = Kokkos::create_mirror_view(y1);
  Kokkos::deep_copy(y1_host, y1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    EXPECT_NEAR(::abs(y1_host(i,0)), ::abs(y2(i,0)), eps);
}

namespace Test {
  struct Functor_TeamHerk {
    char _uplo, _trans;
    int _n, _k;
    matrix_type_device _A, _C;
    ValueType _alpha, _beta;

    Functor_TeamHerk(const char uplo, const char trans,
                     const int n, const int k, 
                     const ValueType alpha, 
                     const matrix_type_device &A,
                     const ValueType beta,
                     const matrix_type_device &C) 
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
                                  (const ValueType*)_A.data(), (int)_A.stride_1(),
                                  _beta,
                                  (      ValueType*)_C.data(), (int)_C.stride_1());
    }
    
    inline
    void run() {
      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);
    }
  };
}


TEST( DenseLinearAlgebra, team_herk_un ) {
  const ordinal_type n = 20, k = 10;
  const char uplo  = 'U';
  const char trans = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", n, k), C1("C1", n, n);
  matrix_type_host   A2("A2", n, k), C2("C2", n, n);
  
  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(C2, C1);

  ::Test::Functor_TeamHerk test(uplo, trans,
                                n, k, 
                                alpha, A1, beta, C1);
  test.run();

  Blas<ValueType>::herk(uplo, trans,
                        n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_herk_uc ) {
  const ordinal_type n = 20, k = 10;
  const char uplo  = 'U';
  const char trans = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", k, n), C1("C1", n, n);
  matrix_type_host   A2("A2", k, n), C2("C2", n, n);
  
  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(C2, C1);

  ::Test::Functor_TeamHerk test(uplo, trans,
                                n, k, 
                                alpha, A1, beta, C1);
  test.run();

  Blas<ValueType>::herk(uplo, trans,
                        n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_herk_ln ) {
  const ordinal_type n = 20, k = 10;
  const char uplo  = 'L';
  const char trans = 'N';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", n, k), C1("C1", n, n);
  matrix_type_host   A2("A2", n, k), C2("C2", n, n);
  
  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(C2, C1);

  ::Test::Functor_TeamHerk test(uplo, trans,
                                n, k, 
                                alpha, A1, beta, C1);
  test.run();

  Blas<ValueType>::herk(uplo, trans,
                        n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());

  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


TEST( DenseLinearAlgebra, team_herk_lc ) {
  const ordinal_type n = 20, k = 10;
  const char uplo  = 'L';
  const char trans = 'C';
  const ValueType alpha = 1.3, beta = 2.5;

  matrix_type_device A1("A1", k, n), C1("C1", n, n);
  matrix_type_host   A2("A2", k, n), C2("C2", n, n);
  
  Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);
  
  Kokkos::fill_random(A1, random, ValueType(1));
  Kokkos::fill_random(C1, random, ValueType(1));
  
  Kokkos::deep_copy(A2, A1);
  Kokkos::deep_copy(C2, C1);

  ::Test::Functor_TeamHerk test(uplo, trans,
                                n, k, 
                                alpha, A1, beta, C1);
  test.run();

  Blas<ValueType>::herk(uplo, trans,
                        n, k,
                        alpha, 
                        A2.data(), A2.stride_1(),
                        beta,
                        C2.data(), C2.stride_1());
  
  auto C1_host = Kokkos::create_mirror_view(C1);
  Kokkos::deep_copy(C1_host, C1);

  const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 1000;
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j) 
      EXPECT_NEAR(::abs(C1_host(i,j)), ::abs(C2(i,j)), eps);
}


namespace Test {
  struct Functor_TeamTrsv {
    char _uplo, _trans, _diag;
    int _m;
    matrix_type_device _A, _b;

    Functor_TeamTrsv(const char uplo, const char trans, const char diag,
                     const int m,
                     const matrix_type_device &A,
                     const matrix_type_device &b) 
      : _uplo(uplo), _trans(trans), _diag(diag),
        _m(m),
        _A(A), _b(b) {}
    
    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &member) const {
      ::BlasTeam<ValueType>::trsv(member,
                                  _uplo, _trans, _diag,
                                  _m,
                                  (const ValueType*)_A.data(), (int)_A.stride_1(),
                                  (      ValueType*)_b.data(), (int)_b.stride_0());
    }
    
    inline
    void run() {
      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);
    }
  };
}


#define TEAM_TRSV_TEST_BODY {                                           \
    matrix_type_device A1("A1", m, m), b1("b1", m, 1);                  \
    matrix_type_host   A2("A2", m, m), b2("b2", m, 1);                  \
                                                                        \
    Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);      \
                                                                        \
    Kokkos::fill_random(A1, random, ValueType(1));                      \
    Kokkos::fill_random(b1, random, ValueType(1));                      \
                                                                        \
    Kokkos::deep_copy(A2, A1);                                          \
    Kokkos::deep_copy(b2, b1);                                          \
                                                                        \
    ::Test::Functor_TeamTrsv test(uplo, trans, diag,                    \
                                  m,                                    \
                                  A1, b1);                              \
    test.run();                                                         \
                                                                        \
    Blas<ValueType>::trsv(uplo, trans, diag,                            \
                          m,                                            \
                          A2.data(), A2.stride_1(),                     \
                          b2.data(), b2.stride_0());                    \
                                                                        \
    auto b1_host = Kokkos::create_mirror_view(b1);                      \
    Kokkos::deep_copy(b1_host, b1);                                     \
                                                                        \
    const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 10000; \
    for (int i=0;i<m;++i)                                               \
      EXPECT_NEAR(::abs(b1_host(i,0)), ::abs(b2(i,0)), eps*::abs(b2(i,0))); \
  }

TEST( DenseLinearAlgebra, team_trsv_unu ) {
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'N';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsv_unn ) {
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'N';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsv_utu ) {
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'T';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
}

TEST( DenseLinearAlgebra, team_trsv_utn ) {
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'T';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
}

TEST( DenseLinearAlgebra, team_trsv_ucu ) {
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'C';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
}

TEST( DenseLinearAlgebra, team_trsv_ucn ) {
  const ordinal_type m = 20;
  const char uplo  = 'U';
  const char trans = 'C';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
}





TEST( DenseLinearAlgebra, team_trsv_lnu ) {
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'N';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsv_lnn ) {
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'N';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsv_ltu ) {
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'T';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
}

TEST( DenseLinearAlgebra, team_trsv_ltn ) {
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'T';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
}

TEST( DenseLinearAlgebra, team_trsv_lcu ) {
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'C';
  const char diag  = 'U';

  TEAM_TRSV_TEST_BODY;
}

TEST( DenseLinearAlgebra, team_trsv_lcn ) {
  const ordinal_type m = 20;
  const char uplo  = 'L';
  const char trans = 'C';
  const char diag  = 'N';

  TEAM_TRSV_TEST_BODY;
}
#undef TEAM_TRSV_TEST_BODY


namespace Test {
  struct Functor_TeamTrsm {
    char _side, _uplo, _trans, _diag;
    int _m, _n;
    matrix_type_device _A, _B;
    ValueType _alpha;

    Functor_TeamTrsm(const char side, const char uplo, const char trans, const char diag,
                     const int m, const int n,
                     const ValueType alpha,
                     const matrix_type_device &A,
                     const matrix_type_device &B) 
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
                                  (const ValueType*)_A.data(), (int)_A.stride_1(),
                                  (      ValueType*)_B.data(), (int)_B.stride_1());
    }
    
    inline
    void run() {
      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceSpaceType>(1, Kokkos::AUTO), *this);
    }
  };
}

#define TEAM_TRSM_TEST_BODY {                                           \
    matrix_type_device A1("A1", m, m), B1("b1", m, n);                  \
    matrix_type_host   A2("A2", m, m), B2("b2", m, n);                  \
                                                                        \
    Kokkos::Random_XorShift64_Pool<DeviceSpaceType> random(13718);      \
                                                                        \
    Kokkos::fill_random(A1, random, ValueType(1));                      \
    Kokkos::fill_random(B1, random, ValueType(1));                      \
                                                                        \
    Kokkos::deep_copy(A2, A1);                                          \
    Kokkos::deep_copy(B2, B1);                                          \
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
                          A2.data(), A2.stride_1(),                     \
                          B2.data(), B2.stride_1());                    \
                                                                        \
    auto B1_host = Kokkos::create_mirror_view(B1);                      \
    Kokkos::deep_copy(B1_host, B1);                                     \
                                                                        \
    const MagnitudeType eps = std::numeric_limits<MagnitudeType>::epsilon() * 100000; \
    for (int i=0;i<m;++i)                                               \
      for (int j=0;j<n;++j)                                             \
        EXPECT_NEAR(::abs(B1_host(i,j)), ::abs(B2(i,j)), eps*::abs(B2(i,j))); \
  }

TEST( DenseLinearAlgebra, team_trsm_lunu ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'N';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_lunn ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'N';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_lutu ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'T';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_lutn ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'T';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_lucu ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'C';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_lucn ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'U';
  const char trans = 'C';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_llnu ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'N';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_llnn ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'N';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_lltu ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'T';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_lltn ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'T';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_llcu ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'C';
  const char diag  = 'U';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}


TEST( DenseLinearAlgebra, team_trsm_llcn ) {
  const ordinal_type m = 20, n = 10;
  const char side  = 'L';
  const char uplo  = 'L';
  const char trans = 'C';
  const char diag  = 'N';
  const ValueType alpha = 1.2;

  TEAM_TRSM_TEST_BODY;
}



#endif
