#ifndef __KOKKOSBATCHED_GEMM_DECL_HPP__
#define __KOKKOSBATCHED_GEMM_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

  ///
  /// Serial Gemm 
  ///

  template<typename ArgTransA,
           typename ArgTransB,
           typename ArgAlgo>
  struct SerialGemm {
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C);
  };

  ///
  /// Team Gemm
  ///

  template<typename MemberType,
           typename ArgTransA,
           typename ArgTransB,
           typename ArgAlgo>
  struct TeamGemm {
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C);
  };


  ///
  /// Selective Interface
  ///
  template<typename MemberType,
           typename ArgTransA,
           typename ArgTransB,
           typename ArgMode,
           typename ArgAlgo>
  struct Gemm {
    template<typename ScalarType,
             typename AViewType,
             typename BViewType,
             typename CViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const ScalarType alpha,
           const AViewType &A,
           const BViewType &B,
           const ScalarType beta,
           const CViewType &C) {
      int r_val = 0;
      if (std::is_same<ArgMode,Mode::Serial>::value) {
        r_val = SerialGemm<ArgTransA,ArgTransB,ArgAlgo>::invoke(alpha, A, B, beta, C);
      } else if (std::is_same<ArgMode,Mode::Team>::value) {
        r_val = TeamGemm<MemberType,ArgTransA,ArgTransB,ArgAlgo>::invoke(member, alpha, A, B, beta, C);
      } 
      return r_val;
    }
  };      

}


#endif
