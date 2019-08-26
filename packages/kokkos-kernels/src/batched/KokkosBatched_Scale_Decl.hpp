#ifndef __KOKKOSBATCHED_SCALE_DECL_HPP__
#define __KOKKOSBATCHED_SCALE_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

  ///
  /// Serial Scale
  ///

  struct SerialScale {
    template<typename ScalarType,
             typename AViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const ScalarType alpha,
           const AViewType &A);
  };

  ///
  /// Team Scale
  ///

  template<typename MemberType>
  struct TeamScale {
    template<typename ScalarType,
             typename AViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const ScalarType alpha,
           const AViewType &A);
  };

}


#endif
