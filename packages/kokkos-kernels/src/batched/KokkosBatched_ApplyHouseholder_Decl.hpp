#ifndef __KOKKOSBATCHED_APPLY_HOUSEHOLDER_DECL_HPP__
#define __KOKKOSBATCHED_APPLY_HOUSEHOLDER_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

  ///
  /// Serial Householder 
  ///

  // level 1 operation (no blocking algorithm info avail)
  template<typename ArgSide>
  struct SerialApplyHouseholder {
    template<typename uViewType,
             typename tauViewType,
             typename AViewType,
             typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const uViewType &u2,
           const tauViewType &tau,
           const AViewType
           const wViewType &w);
  };

}


#endif
