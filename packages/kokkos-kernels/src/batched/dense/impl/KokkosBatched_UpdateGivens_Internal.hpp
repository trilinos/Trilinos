#ifndef __KOKKOSBATCHED_UPDATE_GIVENS_INTERNAL_HPP__
#define __KOKKOSBATCHED_UPDATE_GIVENS_INTERNAL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

  ///
  /// Serial Internal Impl
  /// ==================== 
  ///
  /// this impl follows the flame interface of householder transformation
  ///
  struct SerialUpdateGivensInternal {
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const Kokkos::pair<ValueType,ValueType> &S,
           /* */ Kokkos::pair<ValueType,ValueType> &G) {
      const ValueType
        tmp    = S.first*G.first -S.second*G.second;
      G.second = S.first*G.second+S.second*G.first;
      G.first  = tmp;

      return 0;
    }
  };

} // end namespace KokkosBatched


#endif
