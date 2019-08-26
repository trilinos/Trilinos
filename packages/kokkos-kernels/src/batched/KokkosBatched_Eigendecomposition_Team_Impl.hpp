#ifndef __KOKKOSBATCHED_EIGENDECOMPOSITION_TEAM_IMPL_HPP__
#define __KOKKOSBATCHED_EIGENDECOMPOSITION_TEAM_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Eigendecomposition_Team_Internal.hpp"

namespace KokkosBatched {

  ///
  /// Team Impl
  /// =========

  template<typename MemberType>
  struct TeamEigendecomposition {
    template<typename AViewType,
             typename EViewType,
             typename UViewType,
             typename WViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const AViewType &A,
           const EViewType &er, const EViewType &ei,
           const UViewType &UL, const UViewType &UR,
           const WViewType &W) {
      /// view checking
      const int m = A.extent(0);
      assert(m == A.extent(1)  && "Eigendecomposition: A is not square");
      assert(m == er.extent(0) && "Eigendecomposition: Length of er does not match to A's dimension");
      assert(m == ei.extent(0) && "Eigendecomposition: Length of ei does not match to A's dimension");
      assert(m == UL.extent(0) && "Eigendecomposition: Length of UL does not match to A's dimension");
      assert(m == UL.extent(1) && "Eigendecomposition: Width of UL does not match to A's dimension");
      assert(m == UR.extent(0) && "Eigendecomposition: Length of UR does not match to A's dimension");
      assert(m == UR.extent(1) && "Eigendecomposition: Width of UR does not match to A's dimension");
      assert(W.extent(0) >= (2*m*m+5*m) && "Eigendecomposition: workspace size is too small");
      assert(W.stride(0) == 1  && "Eigendecomposition: Provided workspace is not contiguous");
    
      return 0;
    }
  };

} /// end namespace KokkosBatched


#endif
