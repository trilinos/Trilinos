#ifndef __KOKKOSBATCHED_SHIFTED_TRSV_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_SHIFTED_TRSV_SERIAL_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Set_Internal.hpp"
#include "KokkosBatched_Scale_Internal.hpp"

#include "KokkosBatched_InnerTrsm_Serial_Impl.hpp"
#include "KokkosBatched_Gemv_Serial_Internal.hpp"

namespace KokkosBatched {

  ///
  /// Serial Internal Impl
  /// ====================

  ///
  /// Lower
  ///

  struct SerialShiftedTrsvInternalLower {
    template<typename ScalarType,
             typename ValueTypeA,
             typename ValueTypeB>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const int m,
           const ScalarType lambda,
           const ValueTypeA *__restrict__ A, const int as0, const int as1,
           /* */ ValueTypeB *__restrict__ b, const int bs0,
           const int *__restrict__ blks) {
      const int as = as0+as1;

      int p = 0;
      for (;p<m;) {
        const int blk = blks[p], iend = m-p-blk;
        assert( ((blk == 1) || (blk == 2)) && "ShiftedTrsvLower: blocks are not 1x1 or 2x2");
        if (blk == 1) {
          const auto alpha11 = A[p*as]-lambda;
          ValueTypeB *__restrict__ beta1 = b+p*bs0;

          // with __restrict__ a compiler assumes that the pointer is not accessed by others
          // op(/=) uses this pointer and changes the associated values, which brings a compiler problem
          *beta1 = *beta1 / alpha11;

          if (iend) {
            const ValueTypeA *__restrict__ a21 = A+p*as+as0; 
            ValueTypeB *__restrict__ b2 = beta1+bs0;
            for (int i=0;i<iend;++i)
              b2[i*bs0] -= a21[i*as0] * (*beta1);              
          }
        } else {
          const int p_plus_one = p+1;
          const auto alpha11 = A[p*as]-lambda;
          const auto alpha12 = A[p*as+as1];
          const auto alpha21 = A[p*as+as0];
          const auto alpha22 = A[p_plus_one*as]-lambda;
          const auto det = alpha11*alpha22-alpha12*alpha21;

          ValueTypeB *__restrict__ beta1 = b+p*bs0;            
          ValueTypeB *__restrict__ beta2 = b+p_plus_one*bs0;            

          const ValueTypeB beta_one = *beta1;
          const ValueTypeB beta_two = *beta2;           

          *beta1 = ( alpha22*(beta_one)-alpha12*(beta_two))/det;
          *beta2 = (-alpha21*(beta_one)+alpha11*(beta_two))/det;
            
          if (iend) {
            const ValueTypeA *__restrict__ A21 = A+p*as+2*as0; 
            ValueTypeB *__restrict__ b2 = beta1+2*bs0;
              
            for (int i=0;i<iend;++i)
              b2[i*bs0] -= (A21[i*as0]*(*beta1) + A21[i*as0+as1]*(*beta2));              
          }
        }
        p += blk;
      }
      return 0;
    }
  };

  ///
  /// Upper
  ///
    
  struct SerialShiftedTrsvInternalUpper {
    template<typename ScalarType,
             typename ValueTypeA,
             typename ValueTypeB>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const int m,
           const ScalarType lambda,
           const ValueTypeA *__restrict__ A, const int as0, const int as1,
           /**/  ValueTypeB *__restrict__ b, const int bs0,
           const int *__restrict__ blks) {
      const int as = as0+as1;

      ValueTypeB *__restrict__ b0 = b;
        
      int p = m-1;
      for (;p>=0;) {
        const int blk = blks[p], iend = p+1-blk;
        assert( ((blk == 1) || (blk == 2)) && "ShiftedTrsvUpper: blocks are not 1x1 or 2x2");          
        if (blk == 1) {
          const auto alpha11 = A[p*as]-lambda;
          /**/  ValueTypeB *__restrict__ beta1 = b+p*bs0;

          // with __restrict__ a compiler assumes that the pointer is not accessed by others
          // op(/=) uses this pointer and changes the associated values, which brings a compiler problem
          *beta1 = *beta1 / alpha11;

          if (iend) {
            const ValueTypeA *__restrict__ a01   = A+p*as1;
            for (int i=0;i<iend;++i)
              b0[i*bs0] -= a01[i*as0] * (*beta1);
          }
        } else {

          const int p_minus_one = p-1;
          const auto alpha11 = A[p_minus_one*as]-lambda;
          const auto alpha12 = A[p_minus_one*as+as1];
          const auto alpha21 = A[p_minus_one*as+as0];
          const auto alpha22 = A[p*as]-lambda;
          const auto det = alpha11*alpha22-alpha12*alpha21;

          ValueTypeB *__restrict__ beta1 = b+p_minus_one*bs0;            
          ValueTypeB *__restrict__ beta2 = b+p*bs0;            

          const ValueTypeB beta_one = *beta1;
          const ValueTypeB beta_two = *beta2;           

          *beta1 = ( alpha22*(beta_one)-alpha12*(beta_two))/det;
          *beta2 = (-alpha21*(beta_one)+alpha11*(beta_two))/det;
            
          if (iend) {
            const ValueTypeA *__restrict__ A01 = A+p_minus_one*as1; 
            for (int i=0;i<iend;++i)
              b0[i*bs0] -= (A01[i*as0]*(*beta1) + A01[i*as0+as1]*(*beta2));              
          }
        }
        p -= blk;
      }
      return 0;
    }
  };
    
}


#endif
