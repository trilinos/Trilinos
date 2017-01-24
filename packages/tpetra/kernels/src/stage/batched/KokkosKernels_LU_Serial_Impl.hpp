#ifndef __KOKKOSKERNELS_LU_SERIAL_IMPL_HPP__
#define __KOKKOSKERNELS_LU_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_Util.hpp"

#include "KokkosKernels_Gemm_Decl.hpp"
#include "KokkosKernels_Trsm_Serial_Decl.hpp"

#include "KokkosKernels_Gemm_Serial_Impl.hpp"
#include "KokkosKernels_Trsm_Serial_Impl.hpp"

namespace KokkosKernels {

  namespace Serial {

    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerLU<bmn>::
    invoke(const int m, const int n,
           ValueType *__restrict__ A) {
      // A(m x n)
      const int
        k = (m < n ? m : n);

      for (int p=0;p<k;++p) {
        const ValueType 
          // inv_alpha11 = 1.0/A[p*_as0+p*_as1],  
          alpha11 = A[p*_as0+p*_as1],
          *__restrict__ a12t = A + (p  )*_as0 + (p+1)*_as1;

        ValueType
          *__restrict__ a21  = A + (p+1)*_as0 + (p  )*_as1,
          *__restrict__ A22  = A + (p+1)*_as0 + (p+1)*_as1;

        const int
          iend = m-p-1,
          jend = n-p-1;
        
        for (int i=0;i<iend;++i) {
          // a21[i*_as0] *= inv_alpha11; 
          a21[i*_as0] /= alpha11;
          for (int j=0;j<jend;++j)
            A22[i*_as0+j*_as1] -= a21[i*_as0] * a12t[j*_as1];
        }
      }
      return 0;
    }
    
    // let's do something stupid...
    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerLU<4>::
    invoke(ValueType *__restrict__ A) {
      // load
      ValueType
        a_00 = A[0*_as0+0*_as1], a_01 = A[0*_as0+1*_as1], a_02 = A[0*_as0+2*_as1], a_03 = A[0*_as0+3*_as1],
        a_10 = A[1*_as0+0*_as1], a_11 = A[1*_as0+1*_as1], a_12 = A[1*_as0+2*_as1], a_13 = A[1*_as0+3*_as1],
        a_20 = A[2*_as0+0*_as1], a_21 = A[2*_as0+1*_as1], a_22 = A[2*_as0+2*_as1], a_23 = A[2*_as0+3*_as1],
        a_30 = A[3*_as0+0*_as1], a_31 = A[3*_as0+1*_as1], a_32 = A[3*_as0+2*_as1], a_33 = A[3*_as0+3*_as1];

      // ValueType inv_alpha;

      // 0 iteration
      //inv_alpha = 1.0/a_00;
      /* a_10 *= inv_alpha; */ a_10 /= a_00;  a_11 -= a_10*a_01; a_12 -= a_10*a_02; a_13 -= a_10*a_03;
      /* a_20 *= inv_alpha; */ a_20 /= a_00;  a_21 -= a_20*a_01; a_22 -= a_20*a_02; a_23 -= a_20*a_03;
      /* a_30 *= inv_alpha; */ a_30 /= a_00;  a_31 -= a_30*a_01; a_32 -= a_30*a_02; a_33 -= a_30*a_03;
      
      // 1 iteration
      //inv_alpha = 1.0/a_11;
      /* a_21 *= inv_alpha; */ a_21 /= a_11;  a_22 -= a_21*a_12; a_23 -= a_21*a_13;
      /* a_31 *= inv_alpha; */ a_31 /= a_11;  a_32 -= a_31*a_12; a_33 -= a_31*a_13;

      // 2 iteration
      a_32 /= a_22; a_33 -= a_32*a_23;

      // store
      A[1*_as0+0*_as1] = a_10; A[1*_as0+1*_as1] = a_11; A[1*_as0+2*_as1] = a_12; A[1*_as0+3*_as1] = a_13;
      A[2*_as0+0*_as1] = a_20; A[2*_as0+1*_as1] = a_21; A[2*_as0+2*_as1] = a_22; A[2*_as0+3*_as1] = a_23;
      A[3*_as0+0*_as1] = a_30; A[3*_as0+1*_as1] = a_31; A[3*_as0+2*_as1] = a_32; A[3*_as0+3*_as1] = a_33;

      return 0;
    }
    
    template<>
    template<typename AViewType>
    KOKKOS_INLINE_FUNCTION
    int
    LU<Algo::LU::Unblocked>::
    invoke(AViewType A) {
      // A = LU (A)
      // A (m x n)

      typedef typename AViewType::value_type value_type;

      const int
        m = A.dimension(0),
        n = A.dimension(1),
        k = (m < n ? m : n);

      const int
        as0 = A.stride_0(),
        as1 = A.stride_1();

      for (int p=0;p<k;++p) {
        const value_type 
          // inv_alpha11 = 1.0/A(p,p),
          alpha11 = A(p,p),
          *__restrict__ a12t = &A(p,  p+1);
        
        value_type
          *__restrict__ a21  = &A(p+1,p  ),
          *__restrict__ A22  = &A(p+1,p+1);
        
        const int
          iend = m-p-1,
          jend = n-p-1;
        
        for (int i=0;i<iend;++i) {
          // a21[i*as0] *= inv_alpha11; 
          a21[i*as0] /= alpha11;
          for (int j=0;j<jend;++j)
            A22[i*as0+j*as1] -= a21[i*as0] * a12t[j*as1];
        }
      }
      return 0;
    }
    
    template<>
    template<typename AViewType>
    KOKKOS_INLINE_FUNCTION
    int
    LU<Algo::LU::Blocked>::
    invoke(AViewType A) {
      // A = LU (A)
      // A (m x n)

      typedef typename AViewType::value_type value_type;

      // should be square matrix
      const int
        m = A.dimension(0),
        n = A.dimension(1),
        k = (m < n ? m : n);
      
      const int
        as0 = A.stride_0(),
        as1 = A.stride_1();

      enum : int {
        mb = Algo::LU::Blocked::mb };

      InnerLU<mb> lu(as0, as1);

      InnerTriSolveLeftLowerUnitDiag<mb> trsm_llu(as0, as1,
                                                  as0, as1);

      InnerTriSolveRightUpperNonUnitDiag<mb> trsm_run(as0, as1,
                                                      as0, as1);

      InnerRankUpdate<mb,mb> gemm_nn(as0, as1,
                                     as0, as1,
                                     as0, as1);

      const int
        mm = (m/mb)*mb,
        nn = (n/mb)*mb,
        kk = (k/mb)*mb;

      for (int p=0;p<kk;p+=mb) {
        // block lu (for now use remainder version)
        //lu.invoke(mb, mb, &A(p,p));
        lu.invoke(&A(p,p));

        // trsm update
        trsm_llu.invoke(&A(p,p), nn-p-mb, &A(p,   p+mb));
        trsm_run.invoke(&A(p,p), mm-p-mb, &A(p+mb,p   ));

        // gemm update
        for (int i=p+mb;i<mm;i+=mb)
          for (int j=p+mb;j<nn;j+=mb)
            gemm_nn.serial_invoke(-1, &A(i,p), &A(p,j), mb, &A(i,j));
      }

      const int
        mp = (m%mb),
        np = (n%mb),
        kp = (k%mb),
        mk = m-kk,
        nk = n-kk;

      if (      np) trsm_llu.invoke(&A(0,0), kk, np, &A( 0,nn));
      if (mp      ) trsm_run.invoke(&A(0,0), mp, kk, &A(mm, 0));
      if (mp && np) {
        gemm_nn.serial_invoke(-1, &A(kk,0), &A(0,kk), mk, nk, kk, &A(kk,kk));
        lu.invoke(mk, nk, &A(kk,kk));
      }

      return 0;
    }

  }
}

#endif
