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
#ifndef __KOKKOSBATCHED_TRSM_DECL_HPP__
#define __KOKKOSBATCHED_TRSM_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

template <typename ArgSide, typename ArgUplo, typename ArgTrans, typename ArgDiag, typename ArgAlgo>
struct SerialTrsm {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B);
};

template <typename MemberType, typename ArgSide, typename ArgUplo, typename ArgTrans, typename ArgDiag,
          typename ArgAlgo>
struct TeamTrsm {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B);
};

template <typename MemberType, typename ArgSide, typename ArgUplo, typename ArgTrans, typename ArgDiag,
          typename ArgAlgo>
struct TeamVectorTrsm {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B);
};

///
/// Selective Interface
///
template <typename MemberType, typename ArgSide, typename ArgUplo, typename ArgTrans, typename ArgDiag,
          typename ArgMode, typename ArgAlgo>
struct Trsm {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                                const BViewType &B) {
    int r_val = 0;
    if (std::is_same<ArgMode, Mode::Serial>::value) {
      r_val = SerialTrsm<ArgSide, ArgUplo, ArgTrans, ArgDiag, ArgAlgo>::invoke(alpha, A, B);
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      r_val = TeamTrsm<MemberType, ArgSide, ArgUplo, ArgTrans, ArgDiag, ArgAlgo>::invoke(member, alpha, A, B);
    }
    return r_val;
  }
};

}  // namespace KokkosBatched

#include "KokkosBatched_Trsm_Serial_Impl.hpp"
#include "KokkosBatched_Trsm_Team_Impl.hpp"
#include "KokkosBatched_Trsm_TeamVector_Impl.hpp"

#endif
