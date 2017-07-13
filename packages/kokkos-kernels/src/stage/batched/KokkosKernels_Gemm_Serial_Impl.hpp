#ifndef __KOKKOSKERNELS_GEMM_SERIAL_IMPL_HPP__
#define __KOKKOSKERNELS_GEMM_SERIAL_IMPL_HPP__


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
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::CompactMKL>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      typedef typename CViewType::value_type vector_type;
      typedef typename vector_type::value_type value_type;

      const int
        m = C.dimension(0),
        n = C.dimension(1),
        k = A.dimension(1),
        vl = vector_type::vector_length;

      // no error check
      cblas_dgemm_compact(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                          m, n, k, 
                          alpha, 
                          (const double*)A.data(), A.stride_0(), 
                          (const double*)B.data(), B.stride_0(), 
                          beta,
                          (double*)C.data(), C.stride_0(),
                          (MKL_INT)vl, (MKL_INT)1);
    }

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Unblocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      
      typedef typename CViewType::value_type value_type;
      
      if      (beta == 0) Util::set  (C, value_type(0)   );
      else if (beta != 1) Util::scale(C, value_type(beta));
      
      if (alpha != 0) {
        const int
          m = C.dimension(0),
          n = C.dimension(1),
          k = A.dimension(1);

        if (!(m && n && k)) return 0;
        
        const int
          as0 = A.stride_0(),
          bs1 = B.stride_1(),
          cs0 = C.stride_0(),
          cs1 = C.stride_1();
        
        value_type
          *__restrict__ pC = &C(0,0);
        for (int p=0;p<k;++p) {
          const value_type
            *__restrict__ pA = &A(0,p),
            *__restrict__ pB = &B(p,0);
          for (int i=0;i<m;++i) {
            const value_type tA(alpha*pA[i*as0]);
            for (int j=0;j<n;++j)
              pC[i*cs0+j*cs1] += tA*pB[j*bs1];
          }
        }
      }
      return 0;
    }

    template<>
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int
    Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Blocked>::
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      
      typedef typename CViewType::value_type value_type;
      
      if      (beta == 0) Util::set  (C, value_type(0)   );
      else if (beta != 1) Util::scale(C, value_type(beta));
      
      if (alpha != 0) {
        const int
          m = C.dimension(0),
          n = C.dimension(1),
          k = A.dimension(1);
        
        if (!(m && n && k)) return 0;

        const int
          as0 = A.stride_0(),
          as1 = A.stride_1(),
          bs0 = B.stride_0(),
          bs1 = B.stride_1(),
          cs0 = C.stride_0(),
          cs1 = C.stride_1();
        
	enum : int {
          mb = Algo::Gemm::Blocked::mb,
          nb = Algo::Gemm::Blocked::nb };
        
        InnerRankUpdate<mb,nb> inner(as0, as1,
                                     bs0, bs1,
                                     cs0, bs1);
        const int mm = (m/mb)*mb, nn = (n/nb)*nb;
        for (int i=0;i<mm;i+=mb)
          for (int j=0;j<nn;j+=nb)
            inner.serial_invoke(alpha, &A(i,0), &B(0,j), k, &C(i,j));
        
        const int mp = (m%mb), np = (n%nb);
        if (mp      ) inner.serial_invoke(alpha, &A(mm, 0), &B( 0, 0), mp, nn, k, &C(mm, 0));
        if (      np) inner.serial_invoke(alpha, &A( 0, 0), &B( 0,nn), mm, np, k, &C( 0,nn));
        if (mp && np) inner.serial_invoke(alpha, &A(mm, 0), &B( 0,nn), mp, np, k, &C(mm,nn));
      }
      return 0;
    }

  } // end namespace Serial


  ///
  /// Inner kernel for remainders
  /// ===========================

  template<int mb, int nb>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  InnerRankUpdate<mb,nb>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ B,
                const int m, const int n, const int k,
                /**/  ValueType *__restrict__ C) {
    if (!(m && n && k)) return 0;

    for (int p=0;p<k;++p) {
      const ValueType
        *__restrict__ pA = A + p*_as1,
        *__restrict__ pB = B + p*_bs0;
      for (int i=0;i<m;++i) {
        const ValueType tA(alpha*pA[i*_as0]);
        for (int j=0;j<n;++j)
          C[i*_cs0+j*_cs1] += tA*pB[j*_bs1];
      }
    }
    return 0;
  }

  ///
  /// Inner kernel (unrolled version, 5x5, 4x4, 3x3, 2x2)
  /// ===================================================
  
  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  InnerRankUpdate<5,5>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ B,
                const int k,
                /**/  ValueType *__restrict__ C) {
    if (!k) return 0;

    ValueType
      a_0p, b_p0, c_00 = 0, c_01 = 0, c_02 = 0, c_03 = 0, c_04 = 0,
      a_1p, b_p1, c_10 = 0, c_11 = 0, c_12 = 0, c_13 = 0, c_14 = 0,
      a_2p, b_p2, c_20 = 0, c_21 = 0, c_22 = 0, c_23 = 0, c_24 = 0,
      a_3p, b_p3, c_30 = 0, c_31 = 0, c_32 = 0, c_33 = 0, c_34 = 0,
      a_4p, b_p4, c_40 = 0, c_41 = 0, c_42 = 0, c_43 = 0, c_44 = 0;
    
    const int
      i0 = 0*_as0, i1 = 1*_as0, i2 = 2*_as0, i3 = 3*_as0, i4 = 4*_as0,
      j0 = 0*_bs1, j1 = 1*_bs1, j2 = 2*_bs1, j3 = 3*_bs1, j4 = 4*_bs1;
    
    for (int p=0;p<k;++p) {
      a_0p = A[i0+p*_as1]; b_p0 = B[p*_bs0+j0];
      a_1p = A[i1+p*_as1]; b_p1 = B[p*_bs0+j1];
      a_2p = A[i2+p*_as1]; b_p2 = B[p*_bs0+j2];
      a_3p = A[i3+p*_as1]; b_p3 = B[p*_bs0+j3];
      a_4p = A[i4+p*_as1]; b_p4 = B[p*_bs0+j4];
      
      c_00 += a_0p * b_p0; c_01 += a_0p * b_p1; c_02 += a_0p * b_p2; c_03 += a_0p * b_p3; c_04 += a_0p * b_p4;
      c_10 += a_1p * b_p0; c_11 += a_1p * b_p1; c_12 += a_1p * b_p2; c_13 += a_1p * b_p3; c_14 += a_1p * b_p4;
      c_20 += a_2p * b_p0; c_21 += a_2p * b_p1; c_22 += a_2p * b_p2; c_23 += a_2p * b_p3; c_24 += a_2p * b_p4;
      c_30 += a_3p * b_p0; c_31 += a_3p * b_p1; c_32 += a_3p * b_p2; c_33 += a_3p * b_p3; c_34 += a_3p * b_p4;
      c_40 += a_4p * b_p0; c_41 += a_4p * b_p1; c_42 += a_4p * b_p2; c_43 += a_4p * b_p3; c_44 += a_4p * b_p4;
    }
    
    C[0*_cs0+0*_cs1] += alpha * c_00; C[0*_cs0+1*_cs1] += alpha * c_01; C[0*_cs0+2*_cs1] += alpha * c_02; C[0*_cs0+3*_cs1] += alpha * c_03; C[0*_cs0+4*_cs1] += alpha * c_04;
    C[1*_cs0+0*_cs1] += alpha * c_10; C[1*_cs0+1*_cs1] += alpha * c_11; C[1*_cs0+2*_cs1] += alpha * c_12; C[1*_cs0+3*_cs1] += alpha * c_13; C[1*_cs0+4*_cs1] += alpha * c_14;
    C[2*_cs0+0*_cs1] += alpha * c_20; C[2*_cs0+1*_cs1] += alpha * c_21; C[2*_cs0+2*_cs1] += alpha * c_22; C[2*_cs0+3*_cs1] += alpha * c_23; C[2*_cs0+4*_cs1] += alpha * c_24;
    C[3*_cs0+0*_cs1] += alpha * c_30; C[3*_cs0+1*_cs1] += alpha * c_31; C[3*_cs0+2*_cs1] += alpha * c_32; C[3*_cs0+3*_cs1] += alpha * c_33; C[3*_cs0+4*_cs1] += alpha * c_34;
    C[4*_cs0+0*_cs1] += alpha * c_40; C[4*_cs0+1*_cs1] += alpha * c_41; C[4*_cs0+2*_cs1] += alpha * c_42; C[4*_cs0+3*_cs1] += alpha * c_43; C[4*_cs0+4*_cs1] += alpha * c_44;
    
    return 0;
  }
  
  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  InnerRankUpdate<4,4>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ B,
                const int k,
                /**/  ValueType *__restrict__ C) {
    if (!k) return 0;

    ValueType
      a_0p, b_p0, c_00 = 0, c_01 = 0, c_02 = 0, c_03 = 0,
      a_1p, b_p1, c_10 = 0, c_11 = 0, c_12 = 0, c_13 = 0,
      a_2p, b_p2, c_20 = 0, c_21 = 0, c_22 = 0, c_23 = 0,
      a_3p, b_p3, c_30 = 0, c_31 = 0, c_32 = 0, c_33 = 0;

    const int
      i0 = 0*_as0, i1 = 1*_as0, i2 = 2*_as0, i3 = 3*_as0,
      j0 = 0*_bs1, j1 = 1*_bs1, j2 = 2*_bs1, j3 = 3*_bs1;

    for (int p=0;p<k;++p) {
      a_0p = A[i0+p*_as1]; b_p0 = B[p*_bs0+j0];
      a_1p = A[i1+p*_as1]; b_p1 = B[p*_bs0+j1];
      a_2p = A[i2+p*_as1]; b_p2 = B[p*_bs0+j2];
      a_3p = A[i3+p*_as1]; b_p3 = B[p*_bs0+j3];

      c_00 += a_0p * b_p0; c_01 += a_0p * b_p1; c_02 += a_0p * b_p2; c_03 += a_0p * b_p3;
      c_10 += a_1p * b_p0; c_11 += a_1p * b_p1; c_12 += a_1p * b_p2; c_13 += a_1p * b_p3;
      c_20 += a_2p * b_p0; c_21 += a_2p * b_p1; c_22 += a_2p * b_p2; c_23 += a_2p * b_p3;
      c_30 += a_3p * b_p0; c_31 += a_3p * b_p1; c_32 += a_3p * b_p2; c_33 += a_3p * b_p3;
    }

    C[0*_cs0+0*_cs1] += alpha * c_00; C[0*_cs0+1*_cs1] += alpha * c_01; C[0*_cs0+2*_cs1] += alpha * c_02; C[0*_cs0+3*_cs1] += alpha * c_03;
    C[1*_cs0+0*_cs1] += alpha * c_10; C[1*_cs0+1*_cs1] += alpha * c_11; C[1*_cs0+2*_cs1] += alpha * c_12; C[1*_cs0+3*_cs1] += alpha * c_13;
    C[2*_cs0+0*_cs1] += alpha * c_20; C[2*_cs0+1*_cs1] += alpha * c_21; C[2*_cs0+2*_cs1] += alpha * c_22; C[2*_cs0+3*_cs1] += alpha * c_23;
    C[3*_cs0+0*_cs1] += alpha * c_30; C[3*_cs0+1*_cs1] += alpha * c_31; C[3*_cs0+2*_cs1] += alpha * c_32; C[3*_cs0+3*_cs1] += alpha * c_33;

    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  InnerRankUpdate<3,3>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ B,
                const int k,
                /**/  ValueType *__restrict__ C) {
    if (!k) return 0;

    ValueType
      a_0p, b_p0, c_00 = 0, c_01 = 0, c_02 = 0,
      a_1p, b_p1, c_10 = 0, c_11 = 0, c_12 = 0,
      a_2p, b_p2, c_20 = 0, c_21 = 0, c_22 = 0;

    const int
      i0 = 0*_as0, i1 = 1*_as0, i2 = 2*_as0,
      j0 = 0*_bs1, j1 = 1*_bs1, j2 = 2*_bs1;

    for (int p=0;p<k;++p) {
      a_0p = A[i0+p*_as1]; b_p0 = B[p*_bs0+j0];
      a_1p = A[i1+p*_as1]; b_p1 = B[p*_bs0+j1];
      a_2p = A[i2+p*_as1]; b_p2 = B[p*_bs0+j2];

      c_00 += a_0p * b_p0; c_01 += a_0p * b_p1; c_02 += a_0p * b_p2;
      c_10 += a_1p * b_p0; c_11 += a_1p * b_p1; c_12 += a_1p * b_p2;
      c_20 += a_2p * b_p0; c_21 += a_2p * b_p1; c_22 += a_2p * b_p2;
    }

    C[0*_cs0+0*_cs1] += alpha * c_00; C[0*_cs0+1*_cs1] += alpha * c_01; C[0*_cs0+2*_cs1] += alpha * c_02;
    C[1*_cs0+0*_cs1] += alpha * c_10; C[1*_cs0+1*_cs1] += alpha * c_11; C[1*_cs0+2*_cs1] += alpha * c_12;
    C[2*_cs0+0*_cs1] += alpha * c_20; C[2*_cs0+1*_cs1] += alpha * c_21; C[2*_cs0+2*_cs1] += alpha * c_22;

    return 0;
  }

  template<>
  template<typename ScalarType,
           typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  InnerRankUpdate<2,2>::
  serial_invoke(const ScalarType alpha,
                const ValueType *__restrict__ A,
                const ValueType *__restrict__ B,
                const int k,
                /**/  ValueType *__restrict__ C) {
    if (!k) return 0;

    ValueType
      a_0p, b_p0, c_00 = 0, c_01 = 0,
      a_1p, b_p1, c_10 = 0, c_11 = 0,
      a_2p, b_p2, c_20 = 0, c_21 = 0;

    const int
      i0 = 0*_as0, i1 = 1*_as0,
      j0 = 0*_bs1, j1 = 1*_bs1;

    for (int p=0;p<k;++p) {
      a_0p = A[i0+p*_as1]; b_p0 = B[p*_bs0+j0];
      a_1p = A[i1+p*_as1]; b_p1 = B[p*_bs0+j1];

      c_00 += a_0p * b_p0; c_01 += a_0p * b_p1;
      c_10 += a_1p * b_p0; c_11 += a_1p * b_p1;
    }

    C[0*_cs0+0*_cs1] += alpha * c_00; C[0*_cs0+1*_cs1] += alpha * c_01;
    C[1*_cs0+0*_cs1] += alpha * c_10; C[1*_cs0+1*_cs1] += alpha * c_11;

    return 0;
  }

} // end namespace KokkosKernels

#endif
