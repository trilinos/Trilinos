#ifndef __KOKKOSBATCHED_COPY_DECL_HPP__
#define __KOKKOSBATCHED_COPY_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

  ///
  /// Serial Copy
  ///

  template<typename ArgTrans>
  struct SerialCopy {
    template<typename AViewType,
             typename BViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const AViewType &A,
           /* */ BViewType &B);
  };

  ///
  /// Team Copy
  ///

  template<typename MemberType, typename ArgTrans>
  struct TeamCopy {
    template<typename AViewType,
             typename BViewType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const AViewType &A,
           const BViewType &B);
  };


  ///
  /// Selective Interface
  ///
  template<typename MemberType,
           typename ArgTrans,
           typename ArgMode>
  struct Copy {
    template<typename AViewType,
             typename BViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    static int
    invoke(const MemberType &member,
           const AViewType &A,
           const BViewType &B) {
      int r_val = 0;
      if (std::is_same<ArgMode,Mode::Serial>::value) {
        r_val = SerialCopy<ArgTrans>::invoke(A, B);
      } else if (std::is_same<ArgMode,Mode::Team>::value) {
        r_val = TeamCopy<MemberType,ArgTrans>::invoke(member, A, B);
      } 
      return r_val;
    }
  };      

}


#define KOKKOSBATCHED_SERIAL_COPY_MATRIX_NO_TRANSPOSE_INTERNAL_INVOKE(M,N,A,AS0,AS1,B,BS0,BS1) \
  KokkosBatched::SerialCopyInternal                                     \
  ::invoke(M, N, A, AS0, AS1, B, BS0, BS1)

#define KOKKOSBATCHED_TEAM_COPY_MATRIX_NO_TRANSPOSE_INTERNAL_INVOKE(MEMBER,M,N,A,AS0,AS1,B,BS0,BS1) \
  KokkosBatched::TeamCopyInternal                                       \
  ::invoke(MEMBER, M, N, A, AS0, AS1, B, BS0, BS1)

#define KOKKOSBATCHED_SERIAL_COPY_VECTOR_INTERNAL_INVOKE(M,A,AS,B,BS)	\
  KokkosBatched::SerialCopyInternal                                     \
  ::invoke(M, A, AS, B, BS)

#define KOKKOSBATCHED_TEAM_COPY_VECTOR_NO_TRANSPOSE_INTERNAL_INVOKE(MEMBER,M,A,AS,B,BS) \
  KokkosBatched::TeamCopyInternal                                       \
  ::invoke(MEMBER, M, A, AS, B, BS)

#define KOKKOSBATCHED_COPY_VECTOR_NO_TRANSPOSE_INTERNAL_INVOKE(MODETYPE,MEMBER,M,A,AS,B,BS) \
  if (std::is_same<MODETYPE,KokkosBatched::Mode::Serial>::value) {      \
    KOKKOSBATCHED_SERIAL_COPY_VECTOR_INTERNAL_INVOKE(M,A,AS,B,BS);	\
  } else if (std::is_same<MODETYPE,KokkosBatched::Mode::Team>::value) { \
    KOKKOSBATCHED_TEAM_COPY_VECTOR_NO_TRANSPOSE_INTERNAL_INVOKE(MEMBER,M,A,AS,B,BS); \
  }

#endif
