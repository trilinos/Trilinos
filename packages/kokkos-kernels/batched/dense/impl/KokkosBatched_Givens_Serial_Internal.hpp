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
#ifndef __KOKKOSBATCHED_GIVENS_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_GIVENS_SERIAL_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
///
/// this impl follows the flame interface of householder transformation
///
struct SerialGivensInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ValueType chi1, const ValueType chi2,
                                           /* */ Kokkos::pair<ValueType, ValueType>* G,
                                           /* */ ValueType* chi1_new) {
    typedef ValueType value_type;
    const value_type zero(0), one(1);
    /// compute G = [ gamma -sigma;
    ///               sigma  gamma ];
    /// G.first = gamma and G.second = sigma
    /// this rotation satisfy the following
    ///   G' [chi1; = [ alpha;
    ///       chi2]     zero ];
    value_type cs, sn, r;
    if (chi2 == zero) {
      r  = chi1;
      cs = one;
      sn = zero;
    } else if (chi1 == zero) {
      r  = chi2;
      cs = zero;
      sn = one;
    } else {
      // here we do not care overflow caused by the division although it is
      // probable....
      r  = Kokkos::ArithTraits<value_type>::sqrt(chi1 * chi1 + chi2 * chi2);
      cs = chi1 / r;
      sn = chi2 / r;

      if (Kokkos::ArithTraits<value_type>::abs(chi1) > Kokkos::ArithTraits<value_type>::abs(chi2) && cs < zero) {
        cs = -cs;
        sn = -sn;
        r  = -r;
      }
    }

    G->first  = cs;
    G->second = sn;
    *chi1_new = r;

    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
