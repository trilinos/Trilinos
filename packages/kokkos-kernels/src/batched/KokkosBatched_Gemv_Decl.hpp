#ifndef __KOKKOSBATCHED_GEMV_DECL_HPP__
#define __KOKKOSBATCHED_GEMV_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

  ///
  /// Serial Gemv 
  ///

  template<typename ArgTrans,
           typename ArgAlgo>
  struct SerialGemv {
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y);
  };
    
  ///
  /// Team Gemv 
  ///

  template<typename MemberType,
           typename ArgTrans,
           typename ArgAlgo>
  struct TeamGemv {
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y);
  };

  ///
  /// Selective Interface
  ///
  template<typename MemberType,
           typename ArgTrans,
           typename ArgMode, typename ArgAlgo>
  struct Gemv {
    template<typename ScalarType,
             typename AViewType,
             typename xViewType,
             typename yViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const ScalarType alpha,
           const AViewType &A,
           const xViewType &x,
           const ScalarType beta,
           const yViewType &y) {
      int r_val = 0;
      if (std::is_same<ArgMode,Mode::Serial>::value) {
        r_val = SerialGemv<ArgTrans,ArgAlgo>::invoke(alpha, A, x, beta, y);
      } else if (std::is_same<ArgMode,Mode::Team>::value) {
        r_val = TeamGemv<MemberType,ArgTrans,ArgAlgo>::invoke(member, alpha, A, x, beta, y);
      } 
      return r_val;
    }      
  };

}


#define KOKKOSBATCHED_SERIAL_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS) \
  KokkosBatched::SerialGemvInternal<ALGOTYPE>                           \
  ::invoke(M, N, ALPHA, A, AS0, AS1, X, XS, BETA, Y, YS)

#define KOKKOSBATCHED_SERIAL_GEMV_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS) \
  KokkosBatched::SerialGemvInternal<ALGOTYPE>                           \
  ::invoke(N, M, ALPHA, A, AS1, AS0, X, XS, BETA, Y, YS)

#define KOKKOSBATCHED_TEAM_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS) \
  KokkosBatched::TeamGemvInternal<ALGOTYPE>                             \
  ::invoke(MEMBER, M, N, ALPHA, A, AS0, AS1, X, XS, BETA, Y, YS)

#define KOKKOSBATCHED_TEAM_GEMV_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS) \
  KokkosBatched::TeamGemvInternal<ALGOTYPE>                             \
  ::invoke(MEMBER, N, M, ALPHA, A, AS1, AS0, X, XS, BETA, Y, YS)
  
#define KOKKOSBATCHED_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE(MODETYPE,ALGOTYPE,MEMBER,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS) \
  if (std::is_same<MODETYPE,KokkosBatched::Mode::Serial>::value) {      \
    KOKKOSBATCHED_SERIAL_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS); \
  } else if (std::is_same<MODETYPE,KokkosBatched::Mode::Team>::value) { \
    KOKKOSBATCHED_TEAM_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS); \
  }                                                                   

#define KOKKOSBATCHED_GEMV_TRANSPOSE_INTERNAL_INVOKE(MODETYPE,ALGOTYPE,MEMBER,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS) \
  if (std::is_same<MODETYPE,KokkosBatched::Mode::Serial>::value) {      \
    KOKKOSBATCHED_SERIAL_GEMV_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS); \
  } else if (std::is_same<MODETYPE,KokkosBatched::Mode::Team>::value) { \
    KOKKOSBATCHED_TEAM_GEMV_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS); \
  }

#endif
