#ifndef __KOKKOSKERNELS_GEMV_SERIAL_IMPL_HPP__
#define __KOKKOSKERNELS_GEMV_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_Util.hpp"

namespace KokkosKernels {

  ///
  /// Serial Impl
  /// ===========
  namespace Serial {

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    int
    Gemv<Trans::NoTranspose,Algo::Gemv::CompactMKL>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y) {
      typedef typename yViewType::value_type vector_type;
      typedef typename vector_type::value_type value_type;

      const int
        m = A.dimension(0),
        n = 1,
        k = A.dimension(1),
        vl = vector_type::vector_length;

      // no error check
      cblas_dgemm_compact(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                          m, n, k, 
                          alpha, 
                          (const double*)A.data(), A.stride_0(), 
                          (const double*)x.data(), x.stride_0(), 
                          beta,
                          (double*)y.data(), y.stride_0(),
                          (MKL_INT)vl, (MKL_INT)1);
    }
    
    template<>
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    int
    Gemv<Trans::NoTranspose,Algo::Gemv::Unblocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y) {
      // y = beta y + alpha A x
      // y (m), A(m x n), B(n)
      
      typedef typename yViewType::value_type value_type;
      
      if      (beta == 0) Util::set  (y, value_type(0)   );
      else if (beta != 1) Util::scale(y, value_type(beta));
      
      if (alpha != 0) {
        const int
          m = A.dimension_0(),
          n = A.dimension_1();

        if (!(m && n)) return 0;
        
        const int
          as0 = A.stride_0(),
          as1 = A.stride_1(),
          xs0 = x.stride_0(),
          ys0 = y.stride_0();
        
        value_type
          *__restrict__ pY = &y(0);
        const value_type
          *__restrict__ pA = &A(0,0),
          *__restrict__ pX = &x(0);

        for (int i=0;i<m;++i) {
          value_type t(0);
          const value_type 
            *__restrict__ tA = (pA + i*as0);
          for (int j=0;j<n;++j)
            t += tA[j*as1]*pX[j*xs0];
          pY[i*ys0] += alpha*t;
        }
      }
      return 0;
    }

    template<int mb>
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int 
    InnerMultipleDotProduct<mb>::
    invoke(const ScalarType alpha,
           const ValueType *__restrict__ A,
           const ValueType *__restrict__ x,
           const int m, const int n, 
           /**/  ValueType *__restrict__ y) {
      if (!(m && n)) return 0;

      for (int i=0;i<m;++i) {
        ValueType t(0);
        const ValueType
          *__restrict__ tA = A + i*_as0;
        for (int j=0;j<n;++j)
          t += tA[j*_as1]*x[j*_xs0];
        y[i*_ys0] += alpha*t;
      }
    }

    template<>
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int 
    InnerMultipleDotProduct<4>::
    invoke(const ScalarType alpha,
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

      // // unroll by rows and cols

      // ValueType
      //   y_00 = 0, y_10 = 0, y_20 = 0, y_30 = 0,
      //   y_01 = 0, y_11 = 0, y_21 = 0, y_31 = 0;

      // const int nn = (n/2)*2;
      // for (int j=0;j<nn;j+=2) {
      //   const int j0 = j*_as1, j1 = (j+1)*_as1;
      //   const ValueType x_j0 = x[j*_xs0], x_j1 = x[(j+1)*_xs0];
        
      //   y_00 += A[i0+j0]*x_j0;
      //   y_10 += A[i1+j0]*x_j0;
      //   y_20 += A[i2+j0]*x_j0;
      //   y_30 += A[i3+j0]*x_j0;

      //   y_01 += A[i0+j1]*x_j1;
      //   y_11 += A[i1+j1]*x_j1;
      //   y_21 += A[i2+j1]*x_j1;
      //   y_31 += A[i3+j1]*x_j1;
      // }
      // if (n%2) {
      //   const int j0 = nn*_as1;
      //   const ValueType x_j0 = x[nn*_xs0];
      //   y_00 += A[i0+j0]*x_j0;
      //   y_10 += A[i1+j0]*x_j0;
      //   y_20 += A[i2+j0]*x_j0;
      //   y_30 += A[i3+j0]*x_j0;
      // }
        
      // y[0*_ys0] += alpha*(y_00+=y_01);
      // y[1*_ys0] += alpha*(y_10+=y_11);
      // y[2*_ys0] += alpha*(y_20+=y_21);
      // y[3*_ys0] += alpha*(y_30+=y_31);

      return 0;
    }

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    int
    Gemv<Trans::NoTranspose,Algo::Gemv::Blocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y) {
      // y = beta y + alpha A x
      // y (m), A(m x n), B(n)
      
      typedef typename yViewType::value_type value_type;
      
      if      (beta == 0) Util::set  (y, value_type(0)   );
      else if (beta != 1) Util::scale(y, value_type(beta));
      
      if (alpha != 0) {
        const int
          m = A.dimension_0(),
          n = A.dimension_1();

        if (!(m && n)) return 0;
        
        const int
          as0 = A.stride_0(),
          as1 = A.stride_1(),
          xs0 = x.stride_0(),
          ys0 = y.stride_0();

        enum : int {
          mb = Algo::Gemm::Blocked::mb };

        InnerMultipleDotProduct<mb> inner(as0, as1,
                                          xs0, 
                                          ys0);

        const value_type
          *__restrict__ pX = &x(0);
        
        const int mm = (m/mb)*mb;
        for (int i=0;i<mm;i+=mb) 
          inner.invoke(alpha, &A(i,0), pX, n, &y(i));
        
        const int mp = (m%mb);
        if (mp) inner.invoke(alpha, &A(mm,0), pX, mp, n, &y(mm));
      }
      return 0;
    }
  }
} // end namespace KokkosKernels

#endif
