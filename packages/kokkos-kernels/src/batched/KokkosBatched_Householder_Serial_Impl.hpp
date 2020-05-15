#ifndef __KOKKOSBATCHED_HOUSEHOLDER_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_HOUSEHOLDER_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Householder_Serial_Internal.hpp"

namespace KokkosBatched {

  ///
  /// Serial Impl
  /// ===========

  template<>
  template<typename aViewType,
           typename tauViewType>
  KOKKOS_INLINE_FUNCTION
  int
  SerialHouseholder<Side::Left>::
  invoke(const aViewType &a,
         const tauViewType &tau) {
    return SerialLeftHouseholderInternal::
      invoke(a.extent(0)-1,
             a.data(),
             a.data()+a.stride(0), a.stride(0),
             tau.data());
  }
        
}


#endif
