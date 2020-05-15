#ifndef __KOKKOSBATCHED_TRSV_DECL_HPP__
#define __KOKKOSBATCHED_TRSV_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

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


  ///
  /// Selective Interface
  ///
  template<typename MemberType, 
           typename ArgUplo,
           typename ArgTrans,
           typename ArgDiag,
           typename ArgMode, typename ArgAlgo>
  struct Trsv {      
    template<typename ScalarType,
             typename AViewType,
             typename bViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const ScalarType alpha,
           const AViewType &A,
           const bViewType &b) {
      int r_val = 0;
      if (std::is_same<ArgMode,Mode::Serial>::value) {
        r_val = SerialTrsv<ArgUplo,ArgTrans,ArgDiag,ArgAlgo>::invoke(alpha, A, b);
      } else if (std::is_same<ArgMode,Mode::Team>::value) {
        r_val = TeamTrsv<MemberType,ArgUplo,ArgTrans,ArgDiag,ArgAlgo>::invoke(member, alpha, A, b);
      } 
      return r_val;
    }      
  };

}

#define KOKKOSBATCHED_SERIAL_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  KokkosBatched::SerialTrsvInternalLower<ALGOTYPE>::invoke(DIAG::use_unit_diag, M, ALPHA, A, AS0, AS1, B, BS)

#define KOKKOSBATCHED_SERIAL_TRSV_LOWER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  KokkosBatched::SerialTrsvInternalUpper<ALGOTYPE>::invoke(DIAG::use_unit_diag, N, ALPHA, A, AS1, AS0, B, BS)

#define KOKKOSBATCHED_SERIAL_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  KokkosBatched::SerialTrsvInternalUpper<ALGOTYPE>::invoke(DIAG::use_unit_diag, M, ALPHA, A, AS0, AS1, B, BS)

#define KOKKOSBATCHED_SERIAL_TRSV_UPPER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  KokkosBatched::SerialTrsvInternalLower<ALGOTYPE>::invoke(DIAG::use_unit_diag, N, ALPHA, A, AS1, AS0, B, BS) 

#define KOKKOSBATCHED_TEAM_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  KokkosBatched::TeamTrsvInternalLower<ALGOTYPE>::invoke(MEMBER,DIAG::use_unit_diag, M, ALPHA, A, AS0, AS1, B, BS)

#define KOKKOSBATCHED_TEAM_TRSV_LOWER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  KokkosBatched::TeamTrsvInternalUpper<ALGOTYPE>::invoke(MEMBER,DIAG::use_unit_diag, N, ALPHA, A, AS1, AS0, B, BS)

#define KOKKOSBATCHED_TEAM_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  KokkosBatched::TeamTrsvInternalUpper<ALGOTYPE>::invoke(MEMBER,DIAG::use_unit_diag, M, ALPHA, A, AS0, AS1, B, BS)

#define KOKKOSBATCHED_TEAM_TRSV_UPPER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  KokkosBatched::TeamTrsvInternalLower<ALGOTYPE>::invoke(MEMBER,DIAG::use_unit_diag, N, ALPHA, A, AS1, AS0, B, BS) 

#define KOKKOSBATCHED_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE(MODETYPE,ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  if (std::is_same<MODETYPE,KokkosBatched::Mode::Serial>::value) {      \
    KOKKOSBATCHED_SERIAL_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS); \
  } else if (std::is_same<MODETYPE,KokkosBatched::Mode::Team>::value) { \
    KOKKOSBATCHED_TEAM_TRSV_LOWER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS); \
  }                                                                   

#define KOKKOSBATCHED_TRSV_LOWER_TRANSPOSE_INTERNAL_INVOKE(MODETYPE,ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  if (std::is_same<MODETYPE,KokkosBatched::Mode::Serial>::value) {      \
    KOKKOSBATCHED_SERIAL_TRSV_LOWER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS); \
  } else if (std::is_same<MODETYPE,KokkosBatched::Mode::Team>::value) { \
    KOKKOSBATCHED_TEAM_TRSV_LOWER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS); \
  }

#define KOKKOSBATCHED_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE(MODETYPE,ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  if (std::is_same<MODETYPE,KokkosBatched::Mode::Serial>::value) {      \
    KOKKOSBATCHED_SERIAL_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS); \
  } else if (std::is_same<MODETYPE,KokkosBatched::Mode::Team>::value) { \
    KOKKOSBATCHED_TEAM_TRSV_UPPER_NO_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS); \
  }

#define KOKKOSBATCHED_TRSV_UPPER_TRANSPOSE_INTERNAL_INVOKE(MODETYPE,ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS) \
  if (std::is_same<MODETYPE,KokkosBatched::Mode::Serial>::value) {      \
    KOKKOSBATCHED_SERIAL_TRSV_UPPER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS); \
  } else if (std::is_same<MODETYPE,KokkosBatched::Mode::Team>::value) { \
    KOKKOSBATCHED_TEAM_TRSV_UPPER_TRANSPOSE_INTERNAL_INVOKE(ALGOTYPE,MEMBER,DIAG,M,N,ALPHA,A,AS0,AS1,B,BS); \
  }

#endif
