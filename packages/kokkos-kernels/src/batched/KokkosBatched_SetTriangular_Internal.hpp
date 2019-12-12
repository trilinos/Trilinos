#ifndef __KOKKOSBATCHED_SET_TRIANGULAR_INTERNAL_HPP__
#define __KOKKOSBATCHED_SET_TRIANGULAR_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

  ///
  /// Serial Internal Impl
  /// ==================== 
  struct SerialSetLowerTriangularInternal {
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const int m, const int n,
           const int dist,
           const ScalarType alpha,
           /* */ ValueType *__restrict__ A, const int as0, const int as1) {
      for (int j=0;j<n;++j) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
        for (int i=(j+dist);i<m;++i) {
          A[i*as0+j*as1] = alpha;
        }
      }
        
      return 0;
    }
  };

} // end namespace KokkosBatched


#endif
