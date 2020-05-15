#ifndef __KOKKOSBATCHED_HOUSEHOLDER_DECL_HPP__
#define __KOKKOSBATCHED_HOUSEHOLDER_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

  ///
  /// Serial Householder 
  ///

  // level 1 operation (no blocking algorithm info avail)
  template<typename ArgSide>
  struct SerialHouseholder {
    template<typename aViewType,
             typename tauViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const aViewType &a,
           const tauViewType &tau);
  };

}

#endif
