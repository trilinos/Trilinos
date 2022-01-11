#ifndef __KOKKOSBATCHED_ADD_RADIAL_DECL_HPP__
#define __KOKKOSBATCHED_ADD_RADIAL_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

  ///
  /// This add tiny values on diagonals so the absolute values of diagonals become larger
  ///

  ///
  /// Serial AddRadial
  ///

  struct SerialAddRadial {
    template<typename ScalarType,
             typename AViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const ScalarType tiny,
           const AViewType &A);
  };

  ///
  /// Team Set
  ///

  template<typename MemberType>
  struct TeamAddRadial {
    template<typename ScalarType,
             typename AViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const ScalarType tiny,
           const AViewType &A);
  };

}


#endif
