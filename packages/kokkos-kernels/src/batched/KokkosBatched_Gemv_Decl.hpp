#ifndef __KOKKOSBATCHED_GEMV_DECL_HPP__
#define __KOKKOSBATCHED_GEMV_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)


namespace KokkosBatched {
  namespace Experimental {

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
    
#define KOKKOSBATCHED_SERIAL_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS) \
    KokkosBatched::Experimental::SerialGemvInternal<ALGOTYPE>           \
    ::invoke(M, N, ALPHA, A, AS0, AS1, X, XS, BETA, Y, YS)
    
#define KOKKOSBATCHED_SERIAL_GEMV_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS) \
    KokkosBatched::Experimental::SerialGemvInternal<ALGOTYPE>           \
    ::invoke(N, M, ALPHA, A, AS1, AS0, X, XS, BETA, Y, YS)
    
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

#define KOKKOSBATCHED_TEAM_GEMV_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS) \
    KokkosBatched::Experimental::TeamGemvInternal<ALGOTYPE>           \
    ::invoke(MEMBER, M, N, ALPHA, A, AS0, AS1, X, XS, BETA, Y, YS)
    
#define KOKKOSBATCHED_TEAM_GEMV_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,M,N,ALPHA,A,AS0,AS1,X,XS,BETA,Y,YS) \
    KokkosBatched::Experimental::TeamGemvInternal<ALGOTYPE>           \
    ::invoke(MEMBER, N, M, ALPHA, A, AS1, AS0, X, XS, BETA, Y, YS)

  }
}

#endif
