#ifndef __KOKKOSBATCHED_TRSV_TEAM_IMPL_HPP__
#define __KOKKOSBATCHED_TRSV_TEAM_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trsv_Team_Internal.hpp"


namespace KokkosBatched {
  namespace Experimental {

    ///
    /// Team Impl
    /// ===========

    ///
    /// Implemented:
    /// L/NT, U/NT, L/T, U/T
    /// 
    /// Not yet implemented
    /// L/CT, U/CT 

    ///
    /// L/NT
    ///
    
    template<typename MemberType, typename ArgDiag>
    struct TeamTrsv<MemberType,Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsv::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return TeamTrsvInternalLower<Algo::Trsv::Unblocked>::
          invoke(member,
                 ArgDiag::use_unit_diag,
                 A.dimension_0(), 
                 alpha,
                 A.data(), A.stride_0(), A.stride_1(),
                 b.data(), b.stride_0());
      }
    };

    template<typename MemberType, typename ArgDiag>
    struct TeamTrsv<MemberType,Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsv::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return TeamTrsvInternalLower<Algo::Trsv::Blocked>::
          invoke(member,
                 ArgDiag::use_unit_diag,
                 A.dimension_0(), 
                 alpha,
                 A.data(), A.stride_0(), A.stride_1(),
                 b.data(), b.stride_0());
      }
    };

    ///
    /// L/T
    ///

    template<typename MemberType, typename ArgDiag>
    struct TeamTrsv<MemberType,Uplo::Lower,Trans::Transpose,ArgDiag,Algo::Trsv::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return TeamTrsvInternalUpper<Algo::Trsv::Unblocked>::
          invoke(member,
                 ArgDiag::use_unit_diag,
                 A.dimension_1(), 
                 alpha,
                 A.data(), A.stride_1(), A.stride_0(),
                 b.data(), b.stride_0());
      }
    };

    template<typename MemberType, typename ArgDiag>
    struct TeamTrsv<MemberType,Uplo::Lower,Trans::Transpose,ArgDiag,Algo::Trsv::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member,
             const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return TeamTrsvInternalUpper<Algo::Trsv::Blocked>::
          invoke(ArgDiag::use_unit_diag,
                 A.dimension_1(), 
                 alpha,
                 A.data(), A.stride_1(), A.stride_0(),
                 b.data(), b.stride_0());
      }
    };

    ///
    /// U/NT
    ///

    template<typename MemberType, typename ArgDiag>
    struct TeamTrsv<MemberType,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsv::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return TeamTrsvInternalUpper<Algo::Trsv::Unblocked>::
          invoke(member,
                 ArgDiag::use_unit_diag,
                 A.dimension_0(), 
                 alpha,
                 A.data(), A.stride_0(), A.stride_1(),
                 b.data(), b.stride_0());
      }
    };

    template<typename MemberType, typename ArgDiag>
    struct TeamTrsv<MemberType,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsv::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member, 
             const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return TeamTrsvInternalUpper<Algo::Trsv::Blocked>::
          invoke(member,
                 ArgDiag::use_unit_diag,
                 A.dimension_0(), 
                 alpha,
                 A.data(), A.stride_0(), A.stride_1(),
                 b.data(), b.stride_0());
      }
    };

    ///
    /// U/T
    ///

    template<typename MemberType, typename ArgDiag>
    struct TeamTrsv<MemberType,Uplo::Upper,Trans::Transpose,ArgDiag,Algo::Trsv::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member,
             const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return TeamTrsvInternalLower<Algo::Trsv::Unblocked>::
          invoke(member,
                 ArgDiag::use_unit_diag,
                 A.dimension_1(), 
                 alpha,
                 A.data(), A.stride_1(), A.stride_0(),
                 b.data(), b.stride_0());
      }
    };

    template<typename MemberType, typename ArgDiag>
    struct TeamTrsv<MemberType,Uplo::Upper,Trans::Transpose,ArgDiag,Algo::Trsv::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename bViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member,
             const ScalarType alpha,
             const AViewType &A,
             const bViewType &b) {
        return TeamTrsvInternalLower<Algo::Trsv::Blocked>::
          invoke(member,
                 ArgDiag::use_unit_diag,
                 A.dimension_1(), 
                 alpha,
                 A.data(), A.stride_1(), A.stride_0(),
                 b.data(), b.stride_0());
      }
    };
  }
}


#endif
