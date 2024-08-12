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
#ifndef __KOKKOSBATCHED_NORMALIZE_INTERNAL_HPP__
#define __KOKKOSBATCHED_NORMALIZE_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
struct SerialNormalizeInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m,
                                           /* */ ValueType *KOKKOS_RESTRICT v, const int vs) {
    typedef ValueType value_type;
    typedef Kokkos::ArithTraits<value_type> ats;
    typedef typename ats::mag_type mag_type;

    mag_type norm(0);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) {
      const auto v_at_i = v[i * vs];
      norm += ats::real(v_at_i * ats::conj(v_at_i));
    }
    norm = Kokkos::ArithTraits<mag_type>::sqrt(norm);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) v[i * vs] /= norm;

    return 0;
  }

  template <typename RealType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m,
                                           /* */ RealType *KOKKOS_RESTRICT vr, const int vrs,
                                           /* */ RealType *KOKKOS_RESTRICT vi, const int vis) {
    typedef RealType real_type;
    typedef Kokkos::ArithTraits<real_type> ats;
    typedef typename ats::mag_type mag_type;

    mag_type norm(0);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) {
      const auto vr_at_i = vr[i * vrs];
      const auto vi_at_i = vi[i * vis];
      norm += vr_at_i * vr_at_i + vi_at_i * vi_at_i;
    }
    norm = Kokkos::ArithTraits<mag_type>::sqrt(norm);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) {
      vr[i * vrs] /= norm;
      vi[i * vis] /= norm;
    }

    return 0;
  }
};

}  // end namespace KokkosBatched

#endif
