//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSBLAS2_TEAM_GEMV_SPEC_HPP_
#define KOKKOSBLAS2_TEAM_GEMV_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>
#include <KokkosBlas2_team_gemv_impl.hpp>

namespace KokkosBlas {

template <typename MemberType, typename ArgTrans, typename ArgAlgo = Algo::Gemv::Default>
struct TeamGemv {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& /*member*/, const ScalarType /*alpha*/,
                                           const AViewType& /*A*/, const xViewType& /*x*/, const ScalarType /*beta*/,
                                           const yViewType& /*y*/);
};

template <typename MemberType, typename ArgTrans, typename ArgAlgo>
struct TeamVectorGemv {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& /*member*/, const ScalarType /*alpha*/,
                                           const AViewType& /*A*/, const xViewType& /*x*/, const ScalarType /*beta*/,
                                           const yViewType& /*y*/);
};

///
/// NT
///

template <typename MemberType>
struct TeamGemv<MemberType, Trans::NoTranspose, Algo::Gemv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A,
                                           const xViewType& x, const ScalarType beta, const yViewType& y) {
    static_assert(AViewType::rank == 2, "KokkosBlas::TeamGemv requires rank-2 A matrix");
    return Impl::TeamGemvInternal<Algo::Gemv::Unblocked>::invoke(member, A.extent(0), A.extent(1), alpha, A.data(),
                                                                 A.stride_0(), A.stride_1(), x.data(), x.stride_0(),
                                                                 beta, y.data(), y.stride_0());
  }
};

template <typename MemberType>
struct TeamGemv<MemberType, Trans::NoTranspose, Algo::Gemv::Blocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A,
                                           const xViewType& x, const ScalarType beta, const yViewType& y) {
    static_assert(AViewType::rank == 2, "KokkosBlas::TeamGemv requires rank-2 A matrix");
    return Impl::TeamGemvInternal<Algo::Gemv::Blocked>::invoke(member, A.extent(0), A.extent(1), alpha, A.data(),
                                                               A.stride_0(), A.stride_1(), x.data(), x.stride_0(), beta,
                                                               y.data(), y.stride_0());
  }
};

///
/// T
///

template <typename MemberType>
struct TeamGemv<MemberType, Trans::Transpose, Algo::Gemv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A,
                                           const xViewType& x, const ScalarType beta, const yViewType& y) {
    static_assert(AViewType::rank == 2, "BLAS TeamGemv requires rank-2 A matrix");
    return Impl::TeamGemvInternal<Algo::Gemv::Unblocked>::invoke(member, A.extent(1), A.extent(0), alpha, A.data(),
                                                                 A.stride_1(), A.stride_0(), x.data(), x.stride_0(),
                                                                 beta, y.data(), y.stride_0());
  }
};

template <typename MemberType>
struct TeamGemv<MemberType, Trans::Transpose, Algo::Gemv::Blocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A,
                                           const xViewType& x, const ScalarType beta, const yViewType& y) {
    static_assert(AViewType::rank == 2, "BLAS TeamGemv requires rank-2 A matrix");
    return Impl::TeamGemvInternal<Algo::Gemv::Blocked>::invoke(member, A.extent(1), A.extent(0), alpha, A.data(),
                                                               A.stride_1(), A.stride_0(), x.data(), x.stride_0(), beta,
                                                               y.data(), y.stride_0());
  }
};

///
/// CT
///

template <typename MemberType>
struct TeamGemv<MemberType, Trans::ConjTranspose, Algo::Gemv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A,
                                           const xViewType& x, const ScalarType beta, const yViewType& y) {
    static_assert(AViewType::rank == 2, "BLAS TeamGemv requires rank-2 A matrix");
    return Impl::TeamGemvInternal<Algo::Gemv::Unblocked>::invoke(member, Impl::OpConj{}, A.extent(1), A.extent(0),
                                                                 alpha, A.data(), A.stride_1(), A.stride_0(), x.data(),
                                                                 x.stride_0(), beta, y.data(), y.stride_0());
  }
};

template <typename MemberType>
struct TeamGemv<MemberType, Trans::ConjTranspose, Algo::Gemv::Blocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A,
                                           const xViewType& x, const ScalarType beta, const yViewType& y) {
    static_assert(AViewType::rank == 2, "BLAS TeamGemv requires rank-2 A matrix");
    return Impl::TeamGemvInternal<Algo::Gemv::Blocked>::invoke(member, Impl::OpConj{}, A.extent(1), A.extent(0), alpha,
                                                               A.data(), A.stride_1(), A.stride_0(), x.data(),
                                                               x.stride_0(), beta, y.data(), y.stride_0());
  }
};

///
/// NT
///

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::NoTranspose, Algo::Gemv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A,
                                           const xViewType& x, const ScalarType beta, const yViewType& y) {
    static_assert(AViewType::rank == 2, "Batched TeamVectorGemv requires rank-2 A matrix");
    return Impl::TeamVectorGemvInternal<Algo::Gemv::Unblocked>::invoke(member, A.extent(0), A.extent(1), alpha,
                                                                       A.data(), A.stride_0(), A.stride_1(), x.data(),
                                                                       x.stride_0(), beta, y.data(), y.stride_0());
  }
};

///
/// T
///

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::Transpose, Algo::Gemv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A,
                                           const xViewType& x, const ScalarType beta, const yViewType& y) {
    static_assert(AViewType::rank == 2, "Batched TeamVectorGemv requires rank-2 A matrix");
    return Impl::TeamVectorGemvInternal<Algo::Gemv::Unblocked>::invoke(member, A.extent(1), A.extent(0), alpha,
                                                                       A.data(), A.stride_1(), A.stride_0(), x.data(),
                                                                       x.stride_0(), beta, y.data(), y.stride_0());
  }
};

///
/// CT
///

template <typename MemberType>
struct TeamVectorGemv<MemberType, Trans::ConjTranspose, Algo::Gemv::Unblocked> {
  template <typename ScalarType, typename AViewType, typename xViewType, typename yViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A,
                                           const xViewType& x, const ScalarType beta, const yViewType& y) {
    static_assert(AViewType::rank == 2, "Batched TeamVectorGemv requires rank-2 A matrix");
    return Impl::TeamVectorGemvInternal<Algo::Gemv::Unblocked>::invoke(
        member, Impl::OpConj{}, A.extent(1), A.extent(0), alpha, A.data(), A.stride_1(), A.stride_0(), x.data(),
        x.stride_0(), beta, y.data(), y.stride_0());
  }
};

}  // namespace KokkosBlas

#endif
