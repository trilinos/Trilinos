#ifndef __KOKKOSKERNELS_TRSV_SERIAL_IMPL_HPP__
#define __KOKKOSKERNELS_TRSV_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_Util.hpp"

#include "KokkosKernels_Gemv_Decl.hpp"
#include "KokkosKernels_Gemv_Serial_Impl.hpp"

namespace KokkosKernels {

  namespace Serial {

    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsvLowerUnitDiag<bmn>::
    invoke(const ValueType *__restrict__ A,
           const int m, 
           /**/  ValueType *__restrict__ B) {
      // L(m x m) X(m) = B(m)
      
      if (m == 1 || m == 0) return 0;

      ValueType 
        *__restrict__ a21 = (ValueType*)A + _as0,
        *__restrict__ beta1 = B,
        *__restrict__ b2 = beta1 + _bs0;

      const int asd = _as0 + _as1;
      for (int p=0;p<m;++p,beta1+=_bs0,b2+=_bs0,a21+=asd) {
        const int iend = m-p-1;
        for (int i=0;i<iend;++i)
          b2[i*_bs0] -= a21[i*_as0] * *beta1;
      }
      return 0;
    }
    
    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsvLowerUnitDiag<4>::
    invoke(const ValueType *__restrict__ A,
           /**/  ValueType *__restrict__ B) {
      // L(m x m) X(m) = B(m)

      const ValueType 
        a_10 = A[1*_as0+0*_as1], 
        a_20 = A[2*_as0+0*_as1], a_21 = A[2*_as0+1*_as1], 
        a_30 = A[3*_as0+0*_as1], a_31 = A[3*_as0+1*_as1], a_32 = A[3*_as0+2*_as1];

      ValueType b_0, b_1, b_2, b_3;

      // load
      b_0 = B[0*_bs0];
      b_1 = B[1*_bs0];
      b_2 = B[2*_bs0];
      b_3 = B[3*_bs0];
      
      // 0 iteration
      b_1 -= a_10 * b_0; 
      b_2 -= a_20 * b_0; 
      b_3 -= a_30 * b_0; 
      
      // 1 iteration
      b_2 -= a_21 * b_1;
      b_3 -= a_31 * b_1;
      
      // 2 iteration
      b_3 -= a_32 * b_2; 
      
      // store
      B[1*_bs0] = b_1;
      B[2*_bs0] = b_2;
      B[3*_bs0] = b_3;
      
      return 0;
    }

    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsvLowerNonUnitDiag<bmn>::
    invoke(const ValueType *__restrict__ A,
           const int m, 
           /**/  ValueType *__restrict__ B) {
      // L(m x m) X(m) = B (m)

      if (m == 0) return 0;
      
      ValueType 
        *__restrict__ alpha11 = (ValueType*)A,
        *__restrict__ a21 = (ValueType*)A + _as0,
        *__restrict__ beta1 = B,
        *__restrict__ b2 = beta1 + _bs0;

      const int asd = _as0 + _as1;
      for (int p=0;p<m;++p,beta1+=_bs0,b2+=_bs0,alpha11+=asd,a21+=asd) {
        *beta1 /= *alpha11;

        const int iend = m-p-1;
        for (int i=0;i<iend;++i)
          b2[i*_bs0] -= a21[i*_as0] * *beta1;
      }
      return 0;
    }

    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsvLowerNonUnitDiag<4>::
    invoke(const ValueType *__restrict__ A,
           /**/  ValueType *__restrict__ B) {
      // L(m x m) X(m) = B (m)

      const ValueType 
        a_10 = A[1*_as0+0*_as1], 
        a_20 = A[2*_as0+0*_as1], a_21 = A[2*_as0+1*_as1], 
        a_30 = A[3*_as0+0*_as1], a_31 = A[3*_as0+1*_as1], a_32 = A[3*_as0+2*_as1];

      const ValueType 
        a_00 = A[0*_as0+0*_as1], 
        a_11 = A[1*_as0+1*_as1], 
        a_22 = A[2*_as0+2*_as1], 
        a_33 = A[3*_as0+3*_as1];

      ValueType b_0, b_1, b_2, b_3;

      // load
      b_0 = B[0*_bs0]; 
      b_1 = B[1*_bs0]; 
      b_2 = B[2*_bs0]; 
      b_3 = B[3*_bs0]; 
      
      // 0 iteration
      b_0 /= a_00;
      b_1 -= a_10 * b_0;                  
      b_2 -= a_20 * b_0;                  
      b_3 -= a_30 * b_0;                    
      
      // 1 iteration                         
      b_1 /= a_11;
      b_2 -= a_21 * b_1;                   
      b_3 -= a_31 * b_1;                   
      
      // 2 iteration                         
      b_2 /= a_22; 
      b_3 -= a_32 * b_2;                   
      
      // 3 iteration                         
      b_3 /= a_33;
      
      // store
      B[0*_bs0] = b_0; 
      B[1*_bs0] = b_1; 
      B[2*_bs0] = b_2; 
      B[3*_bs0] = b_3; 

      return 0;
    }

    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsvUpperUnitDiag<bmn>::
    invoke(const ValueType *__restrict__ A,
           const int m, 
           /**/  ValueType *__restrict__ B) {
      // U(m x m) X(m) = B (m)

      if (m == 1 || m == 0) return 0;
      
      ValueType 
        *__restrict__ a01 = (ValueType*)A + (m-1)*_as1,
        *__restrict__ b0 = B,
        *__restrict__ beta1 = B + (m-1)*_bs0;
      
      for (int p=(m-1);p>=0;--p,beta1-=_bs0,a01-=_as1) {
        const int iend = p;
        for (int i=0;i<iend;++i)
          b0[i*_bs0] -= a01[i*_as0] * (*beta1);
      }
      return 0;
    }

    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsvUpperUnitDiag<4>::
    invoke(const ValueType *__restrict__ A,
           /**/  ValueType *__restrict__ B) {
      // U(m x m) X(m x n) = B (m x n)

      const ValueType 
        a_01 = A[0*_as0+1*_as1], a_02 = A[0*_as0+2*_as1], a_03 = A[0*_as0+3*_as1], 
        /**/                     a_12 = A[1*_as0+2*_as1], a_13 = A[1*_as0+3*_as1], 
        /**/                                              a_23 = A[2*_as0+3*_as1];

      ValueType b_0, b_1, b_2, b_3;

      // load
      b_0 = B[0*_bs0];
      b_1 = B[1*_bs0];
      b_2 = B[2*_bs0];
      b_3 = B[3*_bs0];
      
      // 0 iteration
      b_0 -= a_03 * b_3;
      b_1 -= a_13 * b_3;
      b_2 -= a_23 * b_3;
      
      // 1 iteration
      b_0 -= a_02 * b_2; 
      b_1 -= a_12 * b_2;
      
      // 2 iteration
      b_0 -= a_01 * b_1;
      
      // store
      B[0*_bs0] = b_0;
      B[1*_bs0] = b_1;
      B[2*_bs0] = b_2;

      return 0;
    }

    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsvUpperNonUnitDiag<bmn>::
    invoke(const ValueType *__restrict__ A,
           const int m, 
           /**/  ValueType *__restrict__ B) {
      // U(m x m) X(m) = B (m)

      if (m == 0) return 0;

      const int asd = _as0 + _as1;

      ValueType 
        *__restrict__ alpha11 = (ValueType*)A + (m-1)*asd,
        *__restrict__ a01 = (ValueType*)A + (m-1)*_as1,
        *__restrict__ b0 = B,
        *__restrict__ beta1 = B + (m-1)*_bs0;
      
      for (int p=(m-1);p>=0;--p,beta1-=_bs0,alpha11-=asd,a01-=_as1) {
        *beta1 /= *alpha11;

        const int iend = p;
        for (int i=0;i<iend;++i)
          b0[i*_bs0] -= a01[i*_as0] * *beta1;
      }
      return 0;
    }

    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsvUpperNonUnitDiag<4>::
    invoke(const ValueType *__restrict__ A,
           /**/  ValueType *__restrict__ B) {
      // U(m x m) X(m) = B(m)
      const ValueType 
        a_01 = A[0*_as0+1*_as1], a_02 = A[0*_as0+2*_as1], a_03 = A[0*_as0+3*_as1], 
        /**/                     a_12 = A[1*_as0+2*_as1], a_13 = A[1*_as0+3*_as1], 
        /**/                                              a_23 = A[2*_as0+3*_as1];

      const ValueType 
        a_00 = A[0*_as0+0*_as1], 
        a_11 = A[1*_as0+1*_as1], 
        a_22 = A[2*_as0+2*_as1], 
        a_33 = A[3*_as0+3*_as1];

      ValueType b_0, b_1, b_2, b_3;

      // load
      b_0 = B[0*_bs0];
      b_1 = B[1*_bs0];
      b_2 = B[2*_bs0];
      b_3 = B[3*_bs0];
      
      // 0 iteration
      b_3 /= a_33;
      b_2 -= a_23 * b_3;
      b_1 -= a_13 * b_3;
      b_0 -= a_03 * b_3;
      
      // 1 iteration
      b_2 /= a_22;
      b_1 -= a_12 * b_2;
      b_0 -= a_02 * b_2; 
      
      // 2 iteration
      b_1 /= a_11;
      b_0 -= a_01 * b_1;
      
      // 3 iteration                         
      b_0 /= a_00;
      
      // store
      B[0*_bs0] = b_0;
      B[1*_bs0] = b_1;
      B[2*_bs0] = b_2;
      B[3*_bs0] = b_3; 

      return 0;
    }

    template<typename ArgDiag>
    struct Trsv<Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsv::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        typedef typename bViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = b.dimension_0(),
          n = 1,
          vl = vector_type::vector_length;

        // no error check
        cblas_dtrsm_compact(CblasRowMajor, 
                            CblasLeft, CblasLower, CblasNoTrans, 
                            ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                            m, n, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (double*)b.data(), b.stride_0(), 
                            (MKL_INT)vl, (MKL_INT)1);

      }
    };

    template<typename ArgDiag>
    struct Trsv<Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsv::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {

        typedef typename bViewType::value_type value_type;
        
        if (alpha == 0) {
          Util::set(b, 0);
        } else {
          if (alpha != 1)
            Util::scale(b, alpha);
          
          // B (m x 1), A(m x m)
          const int
            m = b.dimension(0);

          if (m == 0 || (m == 1 && ArgDiag::use_unit_diag)) return 0;
          
          const int
            as0 = A.stride_0(),
            as1 = A.stride_1(),
            asd = as0 + as1,
            bs0 = b.stride_0();

          value_type
            *__restrict__ alpha11 = &A(0,0),
            *__restrict__ beta1 = &b(0);

          for (int p=0;p<m;++p,alpha11+=asd,beta1+=bs0) {
            if (!ArgDiag::use_unit_diag) 
              *beta1 /= *alpha11;
            
            value_type
              *__restrict__ a21 = alpha11 + as0,
              *__restrict__ b2 = beta1 + bs0;
            
            const int iend = m-p-1;
            for (int i=0;i<iend;++i)
              b2[i*bs0] -= a21[i*as0] * *beta1;
          }
        }
        
        return 0;
      }
    };

    template<typename ArgDiag>
    struct Trsv<Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsv::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {

        typedef typename bViewType::value_type value_type;
        
        if (alpha == 0) {
          Util::set(b, 0);
        } else {
          if (alpha != 1)
            Util::scale(b, alpha);
          
          // Lower(A) (m x m), B(m x 1)
          const int
            m = b.dimension(0);

          if (m == 0 || (m == 1 && ArgDiag::use_unit_diag)) return 0;

          const int
            as0 = A.stride_0(),
            as1 = A.stride_1(),
            bs0 = b.stride_0();

	  enum : int {
            mb = Algo::Trsv::Blocked::mb };
          
          InnerTrsvLowerUnitDiag<mb> trsv_lu(as0, as1,
                                             bs0);
          
          InnerTrsvLowerNonUnitDiag<mb> trsv_ln(as0, as1,
                                                bs0);

          InnerMultipleDotProduct<mb> gemv_n(as0, as1,
                                             bs0, 
                                             bs0);
          
          const int mm = (m/mb)*mb;
          for (int p=0;p<mm;p+=mb) {
            // trsv update
            if (ArgDiag::use_unit_diag) 
              trsv_lu.invoke(&A(p,p), &b(p));
            else 
              trsv_ln.invoke(&A(p,p), &b(p));
            
            // gemv update
            for (int i=p+mb;i<mm;i+=mb)
              gemv_n.invoke(-1, &A(i,p), &b(p), mb, mb, &b(i));
          }

          // remainder
          const int mp = (m%mb);
          if (mp) {
            gemv_n.invoke(-1, &A(mm,0), &b(0), mp, mm, &b(mm)); 
            if (ArgDiag::use_unit_diag) 
              trsv_lu.invoke(&A(mm,mm), mp, &b(mm));
            else 
              trsv_ln.invoke(&A(mm,mm), mp, &b(mm));
          }
        }
        return 0;
      }
    };

    template<typename ArgDiag>
    struct Trsv<Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsv::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        typedef typename bViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = b.dimension_0(),
          n = 1,
          vl = vector_type::vector_length;

        // no error check
        cblas_dtrsm_compact(CblasRowMajor, 
                            CblasLeft, CblasUpper, CblasNoTrans, 
                            ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                            m, n, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (double*)b.data(), b.stride_0(), 
                            (MKL_INT)vl, (MKL_INT)1);

      }
    };

    template<typename ArgDiag>
    struct Trsv<Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsv::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {

        typedef typename bViewType::value_type value_type;        

        if (alpha == 0) {
          Util::set(b, 0);
        } else {
          if (alpha != 1)
            Util::scale(b, alpha);
          
          // B (m x 1), A(m x m)
          const int
            m = b.dimension(0);

          if (m == 0 || (m == 1 && ArgDiag::use_unit_diag)) return 0;
          
          const int
            as0 = A.stride_0(),
            as1 = A.stride_1(),
            asd = as0 + as1,
            bs0 = b.stride_0();

          value_type
            *__restrict__ a01 = &A(0,m-1),
            *__restrict__ alpha11 = &A(m-1,m-1),
            *__restrict__ beta1 = &b(m-1),
            *__restrict__ b0 = &b(0);
          
          for (int p=(m-1);p>=0;--p,alpha11-=asd,a01-=as1,beta1-=bs0) {
            if (!ArgDiag::use_unit_diag) 
              *beta1 /= *alpha11;
            
            const int iend = p;
            for (int i=0;i<iend;++i)
              b0[i*bs0] -= a01[i*as0] * *beta1;
          }
        }        
        return 0;
      }
    };

    template<typename ArgDiag>
    struct Trsv<Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsv::Blocked> {
      
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {

        typedef typename bViewType::value_type value_type;        

        if (alpha == 0) {
          Util::set(b, 0);
        } else {
          if (alpha != 1)
            Util::scale(b, alpha);
          
          // B (m x n), A(m x m)
          const int
            m = b.dimension(0);

          if (m == 0 || (m == 1 && ArgDiag::use_unit_diag)) return 0;
          
          const int
            as0 = A.stride_0(),
            as1 = A.stride_1(),
            bs0 = b.stride_0();

	  enum : int {
            mb = Algo::Trsv::Blocked::mb };
          
          InnerTrsvUpperUnitDiag<mb> trsv_uu(as0, as1,
                                             bs0);
          
          InnerTrsvUpperNonUnitDiag<mb> trsv_un(as0, as1,
                                                bs0);

          InnerMultipleDotProduct<mb> gemv_n(as0, as1,
                                             bs0, 
                                             bs0);
          
          const int 
            mm = (m/mb)*mb;

          for (int pp=0;pp<mm;pp+=mb) {
            const int p = m - pp - mb;

            // trsm update
            if (ArgDiag::use_unit_diag) 
              trsv_uu.invoke(&A(p,p), &b(p));
            else
              trsv_un.invoke(&A(p,p), &b(p));
            
            // gemm update
            for (int i=(p-mb);i>=0;i-=mb)
              gemv_n.invoke(-1, &A(i,p), &b(p), mb, mb, &b(i));
          }
          
          // remainder
          const int 
            mp = (m%mb);

          if (mp) {
            if (mm) 
              gemv_n.invoke(-1, &A(0,mp), &b(mp), mp, mm, &b(0)); 
            if (ArgDiag::use_unit_diag) 
              trsv_uu.invoke(&A(0,0), mp, &b(0));
            else 
              trsv_un.invoke(&A(0,0), mp, &b(0));
          }
        }
        
        return 0;
      }
      
    };

  }


}

#endif
