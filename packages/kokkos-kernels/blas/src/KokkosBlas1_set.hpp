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

#ifndef KOKKOSBLAS1_SET_HPP_
#define KOKKOSBLAS1_SET_HPP_

#include <KokkosBlas1_set_impl.hpp>

namespace KokkosBlas {

///
/// Serial Set
///

struct SerialSet {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A) {
    return Impl::SerialSetInternal::invoke(A.extent(0), A.extent(1), alpha, A.data(), A.stride_0(), A.stride_1());
  }
};

///
/// Team Set
///

template <typename MemberType>
struct TeamSet {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A) {
    return Impl::TeamSetInternal::invoke(member, A.extent(0), A.extent(1), alpha, A.data(), A.stride_0(), A.stride_1());
  }
};

///
/// TeamVector Set
///

template <typename MemberType>
struct TeamVectorSet {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A) {
    return Impl::TeamVectorSetInternal::invoke(member, A.extent(0), A.extent(1), alpha, A.data(), A.stride_0(),
                                               A.stride_1());
  }
};

}  // namespace KokkosBlas

#endif
