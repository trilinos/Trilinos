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
#ifndef __KOKKOSBATCHED_QR_WITH_COLUMNPIVOTING_DECL_HPP__
#define __KOKKOSBATCHED_QR_WITH_COLUMNPIVOTING_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// TeamVector QR
///

template <typename MemberType, typename ArgAlgo>
struct TeamVectorQR_WithColumnPivoting {
  template <typename AViewType, typename tViewType, typename pViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const pViewType &p, const wViewType &w,
                                           /* */ int &matrix_rank);
};

}  // namespace KokkosBatched

#include "KokkosBatched_QR_WithColumnPivoting_TeamVector_Impl.hpp"

#endif
