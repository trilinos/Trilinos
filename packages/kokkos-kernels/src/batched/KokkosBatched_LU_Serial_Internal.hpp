#ifndef __KOKKOSBATCHED_LU_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_LU_SERIAL_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_InnerLU_Serial_Impl.hpp"
#include "KokkosBatched_InnerTrsm_Serial_Impl.hpp"
#include "KokkosBatched_Gemm_Serial_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {
  
    ///
    /// Serial Internal Impl
    /// ====================

    template<typename AlgoType>
    struct SerialLU_Internal {
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static int 
      invoke(const int m, const int n,
             ValueType *__restrict__ A, const int as0, const int as1);
    };

    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialLU_Internal<Algo::LU::Unblocked>::
    invoke(const int m, const int n,
           ValueType *__restrict__ A, const int as0, const int as1) {
      const int k = (m < n ? m : n);
      if (k <= 0) return 0;

      for (int p=0;p<k;++p) {
        const int iend = m-p-1, jend = n-p-1;

        const ValueType 
          // inv_alpha11 = 1.0/A(p,p),
          alpha11 = A[p*as0+p*as1],
          *__restrict__ a12t = A+(p  )*as0+(p+1)*as1;
        
        ValueType
          *__restrict__ a21  = A+(p+1)*as0+(p  )*as1,
          *__restrict__ A22  = A+(p+1)*as0+(p+1)*as1;
        
        for (int i=0;i<iend;++i) {
          // a21[i*as0] *= inv_alpha11; 
          a21[i*as0] /= alpha11;

              
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
          for (int j=0;j<jend;++j)
            A22[i*as0+j*as1] -= a21[i*as0] * a12t[j*as1];
        }
      }
      return 0;
    }
    
    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    SerialLU_Internal<Algo::LU::Blocked>::
    invoke(const int m, const int n,
           ValueType *__restrict__ A, const int as0, const int as1) {
      enum : int {
        mbAlgo = Algo::LU::Blocked::mb<Kokkos::Impl::ActiveExecutionMemorySpace>()
      };
      const typename Kokkos::Details::ArithTraits<ValueType>::mag_type one(1.0), minus_one(-1.0);          
      const int k = (m < n ? m : n);
      if (k <= 0) return 0;

      InnerLU<mbAlgo> lu(as0, as1);
          
      InnerTrsmLeftLowerUnitDiag<mbAlgo>    trsm_llu(as0, as1, as0, as1);
      InnerTrsmLeftLowerNonUnitDiag<mbAlgo> trsm_run(as1, as0, as1, as0);

      auto lu_factorize = [&](const int ib,
                              const int jb,
                              ValueType *__restrict__ AA) {
        const int mb = mbAlgo;
        const int kb = ib < jb ? ib : jb; 
        for (int p=0;p<kb;p+=mb) {
          const int pb = (p+mb) > kb ? (kb-p) : mb;

          // diagonal block
          ValueType *__restrict__ Ap = AA+p*as0+p*as1;

          // lu on a block             
          lu.serial_invoke(pb, Ap);

          // dimension ABR
          const int m_abr = ib-p-mb, n_abr = jb-p-mb;

          // trsm update
          trsm_llu.serial_invoke(Ap, pb, n_abr, Ap+mb*as1);
          trsm_run.serial_invoke(Ap, pb, m_abr, Ap+mb*as0);
            
          // gemm update
          SerialGemmInternal<Algo::Gemm::Blocked>::
            invoke(m_abr, n_abr, pb,
                   minus_one,
                   Ap+mb*as0, as0, as1,
                   Ap+mb*as1, as0, as1,
                   one,
                   Ap+mb*as0+mb*as1, as0, as1);
        }
      };

      const bool is_small = true; //(m*n <= 64*64);
      if (is_small) {
        lu_factorize(m, n, A);
      } else {
        // // some cache blocking may need (not priority yet);
        // lu_factorize(m, n, A);
      }

      return 0;
    }

  }
}

#endif
