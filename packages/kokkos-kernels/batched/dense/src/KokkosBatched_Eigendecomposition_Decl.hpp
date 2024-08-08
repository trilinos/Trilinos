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
#ifndef __KOKKOSBATCHED_EIGENDECOMPOSITION_DECL_HPP__
#define __KOKKOSBATCHED_EIGENDECOMPOSITION_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

/// Given a general nonsymmetric matrix A (m x m), it performs
/// eigendecomposition of the matrix.
///
/// Parameters:
///   [in] member
///     Team interface only has this argument. Partial specialization can be
///     applied for a different type of team member.
///   [in/out]A
///     Real general nonsymmetric rank 2 view A(m x m).
///     A is first condensed to a upper Hessenberg form. Then, the Francis
///     double shift QR algorithm is applied to compute its Schur form.
///     On exit, A stores a quasi upper triangular matrix of the Schur
///     decomposition.
///   [out]er, [out]ei
///     A real and imaginary eigenvalues, which forms er(m)+ei(m)i
///     For a complex eigen pair, it stores a+bi and a-bi consecutively.
///   [out]UL, [out]UR
///     Left/right eigenvectors are stored in (m x m) matrices. If zero span
///     view is provided, it does not compute the corresponding eigenvectors.
///     However, both UL and UR cannot have zero span. If eigenvalues are only
///     requested, use the Eigenvalue interface which simplifies computations
///   [out]W
///     1D contiguous workspace. The minimum size is (2*m*m+5*m) where m is the
///     dimension of matrix A.

struct SerialEigendecomposition {
  template <typename AViewType, typename EViewType, typename UViewType, typename WViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const EViewType &er, const EViewType &ei,
                                           const UViewType &UL, const UViewType &UR, const WViewType &W);
};

template <typename MemberType>
struct TeamVectorEigendecomposition {
  template <typename AViewType, typename EViewType, typename UViewType, typename WViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const EViewType &er,
                                           const EViewType &ei, const UViewType &UL, const UViewType &UR,
                                           const WViewType &W);
};

}  // namespace KokkosBatched

#endif
