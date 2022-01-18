#ifndef __KOKKOSBATCHED_INNER_MULTIPLE_DOT_PRODUCT_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_INNER_MULTIPLE_DOT_PRODUCT_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_InnerMultipleDotProduct_Decl.hpp"

namespace KokkosBatched {

  ///
  /// Dot Product for GEMV
  /// ====================

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int 
  InnerMultipleDotProduct<5>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ x,
                const int n, 
                /**/  ValueType *__restrict__ y) {
    if (n <= 0) return 0;

    const int 
      i0 = 0*_as0, i1 = 1*_as0, i2 = 2*_as0, i3 = 3*_as0, i4 = 4*_as0;

    // unroll by rows
    ValueType
      y_0 = 0, y_1 = 0, y_2 = 0, y_3 = 0, y_4 = 0;

                    
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int j=0;j<n;++j) {
      const int jj = j*_as1;
      const ValueType x_j = x[j*_xs0];

      y_0 += A[i0+jj]*x_j;
      y_1 += A[i1+jj]*x_j;
      y_2 += A[i2+jj]*x_j;
      y_3 += A[i3+jj]*x_j;
      y_4 += A[i4+jj]*x_j;
    }

    y[0*_ys0] += alpha*y_0;
    y[1*_ys0] += alpha*y_1;
    y[2*_ys0] += alpha*y_2;
    y[3*_ys0] += alpha*y_3;
    y[4*_ys0] += alpha*y_4;

    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int 
  InnerMultipleDotProduct<4>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ x,
                const int n, 
                /**/  ValueType *__restrict__ y) {
    if (!n) return 0;

    const int 
      i0 = 0*_as0, i1 = 1*_as0, i2 = 2*_as0, i3 = 3*_as0;

    // unroll by rows
    ValueType
      y_0 = 0, y_1 = 0, y_2 = 0, y_3 = 0;
            
                    
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int j=0;j<n;++j) {
      const int jj = j*_as1;
      const ValueType x_j = x[j*_xs0];

      y_0 += A[i0+jj]*x_j;
      y_1 += A[i1+jj]*x_j;
      y_2 += A[i2+jj]*x_j;
      y_3 += A[i3+jj]*x_j;
    }

    y[0*_ys0] += alpha*y_0;
    y[1*_ys0] += alpha*y_1;
    y[2*_ys0] += alpha*y_2;
    y[3*_ys0] += alpha*y_3;

    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int 
  InnerMultipleDotProduct<3>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ x,
                const int n, 
                /**/  ValueType *__restrict__ y) {
    if (n <= 0) return 0;

    const int 
      i0 = 0*_as0, i1 = 1*_as0, i2 = 2*_as0;

    // unroll by rows
    ValueType
      y_0 = 0, y_1 = 0, y_2 = 0;
            
                    
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int j=0;j<n;++j) {
      const int jj = j*_as1;
      const ValueType x_j = x[j*_xs0];

      y_0 += A[i0+jj]*x_j;
      y_1 += A[i1+jj]*x_j;
      y_2 += A[i2+jj]*x_j;
    }

    y[0*_ys0] += alpha*y_0;
    y[1*_ys0] += alpha*y_1;
    y[2*_ys0] += alpha*y_2;

    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int 
  InnerMultipleDotProduct<2>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ x,
                const int n, 
                /**/  ValueType *__restrict__ y) {
    if (n <= 0) return 0;

    const int 
      i0 = 0*_as0, i1 = 1*_as0;

    // unroll by rows
    ValueType
      y_0 = 0, y_1 = 0;
            
                    
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int j=0;j<n;++j) {
      const int jj = j*_as1;
      const ValueType x_j = x[j*_xs0];

      y_0 += A[i0+jj]*x_j;
      y_1 += A[i1+jj]*x_j;
    }

    y[0*_ys0] += alpha*y_0;
    y[1*_ys0] += alpha*y_1;

    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int 
  InnerMultipleDotProduct<1>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ x,
                const int n, 
                /**/  ValueType *__restrict__ y) {
    if (n <= 0) return 0;

    // unroll by rows
    ValueType
      y_0 = 0;
            
                    
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int j=0;j<n;++j) 
      y_0 += A[j*_as1]*x[j*_xs0];

    y[0] += alpha*y_0;

    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int 
  InnerMultipleDotProduct<5>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ x,
                const int m, const int n, 
                /**/  ValueType *__restrict__ y) {
    if (m <= 0 || n <= 0) return 0;
    switch (m) {
    case 5: { InnerMultipleDotProduct<5> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    case 4: { InnerMultipleDotProduct<4> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    case 3: { InnerMultipleDotProduct<3> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    case 2: { InnerMultipleDotProduct<2> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    case 1: { InnerMultipleDotProduct<1> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    }
    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int 
  InnerMultipleDotProduct<4>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ x,
                const int m, const int n, 
                /**/  ValueType *__restrict__ y) {
    if (m <= 0 || n <= 0) return 0;
    switch (m) {
    case 4: { InnerMultipleDotProduct<4> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    case 3: { InnerMultipleDotProduct<3> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    case 2: { InnerMultipleDotProduct<2> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    case 1: { InnerMultipleDotProduct<1> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    }
    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int 
  InnerMultipleDotProduct<3>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ x,
                const int m, const int n, 
                /**/  ValueType *__restrict__ y) {
    if (m <= 0 || n <= 0) return 0;
    switch (m) {
    case 3: { InnerMultipleDotProduct<3> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    case 2: { InnerMultipleDotProduct<2> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    case 1: { InnerMultipleDotProduct<1> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    }
    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int 
  InnerMultipleDotProduct<2>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ x,
                const int m, const int n, 
                /**/  ValueType *__restrict__ y) {
    if (m <= 0 || n <= 0) return 0;
    switch (m) {
    case 2: { InnerMultipleDotProduct<2> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    case 1: { InnerMultipleDotProduct<1> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    }
    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int 
  InnerMultipleDotProduct<1>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ x,
                const int m, const int n, 
                /**/  ValueType *__restrict__ y) {
    if (m <= 0 || n <= 0) return 0;
    switch (m) {
    case 1: { InnerMultipleDotProduct<1> inner(_as0, _as1, _xs0, _ys0); inner.serial_invoke(alpha, A, x, n, y); break; }
    }
    return 0;
  }
}



#endif
