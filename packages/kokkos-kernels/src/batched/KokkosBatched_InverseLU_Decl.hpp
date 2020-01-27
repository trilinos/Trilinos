#ifndef __KOKKOSBATCHED_INVERSELU_DECL_HPP__
#define __KOKKOSBATCHED_INVERSELU_DECL_HPP__


/// \author Vinh Dang (vqdang@sandia.gov)

#include "KokkosBatched_Vector.hpp"
#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_Copy_Impl.hpp"
#include "KokkosBatched_SetIdentity_Decl.hpp"
#include "KokkosBatched_SetIdentity_Impl.hpp"
#include "KokkosBatched_SolveLU_Decl.hpp"

namespace KokkosBatched {
      
  template<typename ArgAlgo>
  struct SerialInverseLU {
    template<typename AViewType,
             typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const AViewType &A,
           const wViewType &w) {
      typedef typename wViewType::value_type value_type;
      // workspace w is always 1D view; reinterpret it 
      Kokkos::View<value_type**,Kokkos::LayoutRight,Kokkos::AnonymousSpace>
        W(w.data(), A.extent(0), A.extent(1));

      int r_val[3] = {};
      r_val[0] = SerialCopy<Trans::NoTranspose>::invoke(A, W);
      r_val[1] = SerialSetIdentity::invoke(A);
      r_val[2] = SerialSolveLU<Trans::NoTranspose,ArgAlgo>::invoke(W, A);        
      return r_val[0]+r_val[1]+r_val[2];
    }
  };       

  template<typename MemberType,
           typename ArgAlgo>
  struct TeamInverseLU {
    template<typename AViewType,
             typename wViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const AViewType &A,
           const wViewType &w) {
      typedef typename wViewType::value_type value_type;
      // workspace w is always 1D view; reinterpret it 
      Kokkos::View<value_type**,Kokkos::LayoutRight,Kokkos::AnonymousSpace>
        W(w.data(), A.extent(0), A.extent(1));

      int r_val[3] = {};
      r_val[0] = TeamCopy<MemberType,Trans::NoTranspose>::invoke(member, A, W);
      r_val[1] = TeamSetIdentity<MemberType>::invoke(member, A);
      r_val[2] = TeamSolveLU<MemberType,Trans::NoTranspose,ArgAlgo>::invoke(member, W, A);        
      return r_val[0]+r_val[1]+r_val[2];
    }
  };       
      
}


#endif
