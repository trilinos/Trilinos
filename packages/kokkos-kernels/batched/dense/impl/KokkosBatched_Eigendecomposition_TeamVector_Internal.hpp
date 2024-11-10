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
#ifndef __KOKKOSBATCHED_EIGENDECOMPOSITION_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_EIGENDECOMPOSITION_TEAMVECTOR_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_SetIdentity_Internal.hpp"
#include "KokkosBatched_SetTriangular_Internal.hpp"
#include "KokkosBatched_Normalize_Internal.hpp"

#include "KokkosBatched_Eigendecomposition_Serial_Internal.hpp"

// #include "KokkosBatched_Hessenberg_TeamVector_Internal.hpp"
// #include "KokkosBatched_ApplyQ_TeamVector_Internal.hpp"
// #include "KokkosBatched_Schur_TeamVector_Internal.hpp"
// #include "KokkosBatched_RightEigenvectorFromSchur_TeamVector_Internal.hpp"
// #include "KokkosBatched_LeftEigenvectorFromSchur_TeamVector_Internal.hpp"
// #include "KokkosBatched_Gemm_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Internal Impl
/// ====================

struct TeamVectorEigendecompositionInternal {
  template <typename MemberType, typename RealType>
  KOKKOS_INLINE_FUNCTION static int device_invoke(const MemberType &member, const int m, RealType *A, const int as0,
                                                  const int as1, RealType *er, const int ers, RealType *ei,
                                                  const int eis, RealType *UL, const int uls0, const int uls1,
                                                  RealType *UR, const int urs0, const int urs1, RealType *w,
                                                  const int wlen) {
    /// not yet implemented
    return 0;
  }

  /// Given a general nonsymmetric matrix A (m x m), it performs
  /// eigendecomposition of the matrix.
  ///
  /// Parameters:
  ///   [in]m
  ///     A dimension of the square matrix H.
  ///   [in/out]A, [in]as0, [in]as1
  ///     Real general nonsymmetric matrix A(m x m) with strides as0 and as1.
  ///     A is first condensed to a upper Hessenberg form. Then, the Francis
  ///     double shift QR algorithm is applied to compute its Schur form.
  ///     On exit, A stores a quasi upper triangular matrix of the Schur
  ///     decomposition.
  ///   [out]er, [in]ers, [out]ei, [in]eis
  ///     A complex vector er(m)+ei(m)i with a stride ers and eis to store
  ///     computed eigenvalues. For a complex eigen pair, it stores a+bi and
  ///     a-bi consecutively.
  ///   [out]UL, [in]uls0, [in] uls1
  ///     Left eigenvectors UL(m x m) with strides uls0 and uls1. When UL is
  ///     NULL, the left eigenvectors are not computed.
  ///   [out]UR, [in]urs0, [in] urs1
  ///     Right eigenvectors UR(m x m) with strides urs0 and urs1. When UR is
  ///     NULL, the right eigenvectors are not computed.
  ///   [out]w, [in]wlen
  ///     Workspace
  template <typename MemberType, typename RealType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, RealType *A, const int as0,
                                           const int as1, RealType *er, const int ers, RealType *ei, const int eis,
                                           RealType *UL, const int uls0, const int uls1, RealType *UR, const int urs0,
                                           const int urs1, RealType *w, const int wlen) {
    static_assert(false, "TeamVector eigendecomposition is not implemented yet.");
    /*
    // DO NOT USE
    //
    //     #ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    //     <host code>
    //     #else
    //     <device code>
    //     #endif
    //
    // DO THIS INSTEAD
    //
    //     KOKKOS_IF_ON_HOST((<host code>))
    //     KOKKOS_IF_ON_DEVICE((<device code>))
    //
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    if (as0 == 1 || as1 == 1) {
      /// column major or row major and it runs on host
      /// potentially it can run tpls internally
      Kokkos::single
        (Kokkos::PerTeam(member),
         [&]() {
           SerialEigendecompositionInternal::host_invoke(m,
                                                         A, as0, as1,
                                                         er, ers,
                                                         ei, eis,
                                                         UL, uls0, uls1,
                                                         UR, urs0, urs1,
                                                         w, wlen);
         });
    } else {
      /// arbitrary strides should be handled by native implementation
      device_invoke(member, m,
                    A, as0, as1,
                    er, ers,
                    ei, eis,
                    UL, uls0, uls1,
                    UR, urs0, urs1,
                    w, wlen);
      throw std::runtime_error("TeamVector eigendecomposition is not implemented
yet.");
    }
#else
    /// device code runs
    device_invoke(member, m,
                  A, as0, as1,
                  er, ers,
                  ei, eis,
                  UL, uls0, uls1,
                  UR, urs0, urs1,
                  w, wlen);
#endif
*/
    return 0;
  }
};

}  // namespace KokkosBatched

#endif
