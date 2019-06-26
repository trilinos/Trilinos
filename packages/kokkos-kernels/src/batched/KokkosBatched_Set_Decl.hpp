#ifndef __KOKKOSBATCHED_SET_DECL_HPP__
#define __KOKKOSBATCHED_SET_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {
  ///
  /// Serial Set
  ///

  struct SerialSet {
    template<typename ScalarType,
             typename AViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const ScalarType alpha,
           const AViewType &A);
  };

  ///
  /// Team Set
  ///

  template<typename MemberType>
  struct TeamSet {
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
