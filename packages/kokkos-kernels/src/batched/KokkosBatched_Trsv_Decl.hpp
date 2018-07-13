#ifndef __KOKKOSBATCHED_TRSV_DECL_HPP__
#define __KOKKOSBATCHED_TRSV_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)


namespace KokkosBatched {
  namespace Experimental {  
    ///
    /// Serial Trsv
    ///

    template<typename ArgUplo,
             typename ArgTrans,
             typename ArgDiag,
             typename ArgAlgo>
    struct SerialTrsv {

      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const bViewType &b);
    };

#define KOKKOSBATCHED_SERIAL_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
    KokkosBatched::Experimental::SerialTrsvInternalLower<ALGOTYPE>::invoke(DIAG::use_unit_diag, M, ALPHA, A, AS0, AS1, B, BS)
    
#define KOKKOSBATCHED_SERIAL_TRSV_LOWER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
    KokkosBatched::Experimental::SerialTrsvInternalUpper<ALGOTYPE>::invoke(DIAG::use_unit_diag, N, ALPHA, A, AS1, AS0, B, BS)
    
#define KOKKOSBATCHED_SERIAL_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
    KokkosBatched::Experimental::SerialTrsvInternalUpper<ALGOTYPE>::invoke(DIAG::use_unit_diag, M, ALPHA, A, AS0, AS1, B, BS)
    
#define KOKKOSBATCHED_SERIAL_TRSV_UPPER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
    KokkosBatched::Experimental::SerialTrsvInternalLower<ALGOTYPE>::invoke(DIAG::use_unit_diag, N, ALPHA, A, AS1, AS0, B, BS) 
    
    ///
    /// Team Trsv
    ///

    template<typename MemberType, 
             typename ArgUplo,
             typename ArgTrans,
             typename ArgDiag,
             typename ArgAlgo>
    struct TeamTrsv {

      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member,
             const ScalarType alpha,
             const AViewType &A,
             const bViewType &b);
    };

#define KOKKOSBATCHED_TEAM_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
    KokkosBatched::Experimental::TeamTrsvInternalLower<ALGOTYPE>::invoke(MEMBER,DIAG::use_unit_diag, M, ALPHA, A, AS0, AS1, B, BS)
    
#define KOKKOSBATCHED_TEAM_TRSV_LOWER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
    KokkosBatched::Experimental::TeamTrsvInternalUpper<ALGOTYPE>::invoke(MEMBER,DIAG::use_unit_diag, N, ALPHA, A, AS1, AS0, B, BS)
    
#define KOKKOSBATCHED_TEAM_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
    KokkosBatched::Experimental::TeamTrsvInternalUpper<ALGOTYPE>::invoke(MEMBER,DIAG::use_unit_diag, M, ALPHA, A, AS0, AS1, B, BS)
    
#define KOKKOSBATCHED_TEAM_TRSV_UPPER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
    KokkosBatched::Experimental::TeamTrsvInternalLower<ALGOTYPE>::invoke(MEMBER,DIAG::use_unit_diag, N, ALPHA, A, AS1, AS0, B, BS) 

  }
}

#endif
