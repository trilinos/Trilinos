#ifndef __KOKKOSBATCHED_GEMV_TEAMVECTOR_IMPL_HPP__
#define __KOKKOSBATCHED_GEMV_TEAMVECTOR_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Gemv_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// Team Impl
/// =========

///
/// Implemented:
/// NT, T
///
/// Not yet implemented
/// CT

///
/// NT
///

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::NoTranspose, Algo::Gemv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename xViewType,
            typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const ScalarType alpha, const AViewType &A,
      const xViewType &x, const ScalarType beta, const yViewType &y) {
    return TeamVectorGemvInternal<Algo::Gemv::Unblocked>::invoke(
        member, A.extent(0), A.extent(1), alpha, A.data(), A.stride_0(),
        A.stride_1(), x.data(), x.stride_0(), beta, y.data(), y.stride_0());
  }
};

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::NoTranspose, Algo::Gemv::Blocked> {
  template <typename ScalarType, typename AViewType, typename xViewType,
            typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const ScalarType alpha, const AViewType &A,
      const xViewType &x, const ScalarType beta, const yViewType &y) {
    return TeamVectorGemvInternal<Algo::Gemv::Blocked>::invoke(
        member, A.extent(0), A.extent(1), alpha, A.data(), A.stride_0(),
        A.stride_1(), x.data(), x.stride_0(), beta, y.data(), y.stride_0());
  }
};

///
/// T
///

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::Transpose, Algo::Gemv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename xViewType,
            typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const ScalarType alpha, const AViewType &A,
      const xViewType &x, const ScalarType beta, const yViewType &y) {
    return TeamVectorGemvInternal<Algo::Gemv::Unblocked>::invoke(
        member, A.extent(1), A.extent(0), alpha, A.data(), A.stride_1(),
        A.stride_0(), x.data(), x.stride_0(), beta, y.data(), y.stride_0());
  }
};

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::Transpose, Algo::Gemv::Blocked> {
  template <typename ScalarType, typename AViewType, typename xViewType,
            typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const ScalarType alpha, const AViewType &A,
      const xViewType &x, const ScalarType beta, const yViewType &y) {
    return TeamVectorGemvInternal<Algo::Gemv::Blocked>::invoke(
        member, A.extent(1), A.extent(0), alpha, A.data(), A.stride_1(),
        A.stride_0(), x.data(), x.stride_0(), beta, y.data(), y.stride_0());
  }
};

}  // namespace KokkosBatched

#endif
