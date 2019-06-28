#ifndef __KOKKOSBATCHED_LU_DECL_HPP__
#define __KOKKOSBATCHED_LU_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {
      
  template<typename ArgAlgo>
  struct SerialLU {
    // no piv version
    template<typename AViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const AViewType &A,
           const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny = 0);
  };       

  template<typename MemberType,
           typename ArgAlgo>
  struct TeamLU {
    // no piv version
    template<typename AViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const AViewType &A,
           const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny = 0);
  };       

  ///
  /// Selective Interface
  ///
  template<typename MemberType,
           typename ArgMode, typename ArgAlgo>
  struct LU {
    // no piv version
    template<typename AViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const AViewType &A,
           const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny = 0) {
      int r_val = 0;
      if (std::is_same<ArgMode,Mode::Serial>::value) {
        r_val = SerialLU<ArgAlgo>::invoke(A, tiny);
      } else if (std::is_same<ArgMode,Mode::Team>::value) {
        r_val = TeamLU<MemberType,ArgAlgo>::invoke(member, A, tiny);
      } 
      return r_val;
    }
  };           
      
}

#endif
