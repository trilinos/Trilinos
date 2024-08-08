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
#ifndef __KOKKOSBATCHED_GESV_HPP__
#define __KOKKOSBATCHED_GESV_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

struct Gesv {
  struct StaticPivoting {};
  struct NoPivoting {};

  using Default = StaticPivoting;
};

/// \brief Serial Batched GESV:
///
/// Solve A_l x_l = b_l for all l = 0, ..., N
/// using a batched LU decomposition, 2 batched triangular solves, and a batched
/// static pivoting.
///
/// \tparam MatrixType: Input type for the matrix, needs to be a 2D view
/// \tparam VectorType: Input type for the right-hand side and the solution,
/// needs to be a 1D view
///
/// \param A [in]: matrix, a rank 2 view
/// \param X [out]: solution, a rank 1 view
/// \param B [in]: right-hand side, a rank 1 view
/// \param tmp [in]: a rank 2 view used to store temporary variable; dimension
/// must be n x (n+4) where n is the number of rows.
///
///
/// Two versions are available (those are chosen based on ArgAlgo):
///
///   1. NoPivoting: the solver does not use a pivoting strategy,
///   2. StaticPivoting: the solver uses a static pivoting strategy that relies
///   on using
///      maximal absolute value of row and column to choose pivots and apply
///      them before calling the LU decomposition. Known limitation: the
///      currently implemented strategy would not work with some matrices such
///      as [[2, 1], [1, 0]], when this is the case, the Gesv (if used with
///      pivoting), will return 1 and print an error message.
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgAlgo>
struct SerialGesv {
  template <typename MatrixType, typename XVectorType, typename YVectorType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MatrixType A, const XVectorType X, const YVectorType Y,
                                           const MatrixType tmp);

  template <typename MatrixType, typename VectorType>
  [[deprecated]] KOKKOS_INLINE_FUNCTION static int invoke(const MatrixType A, const VectorType X, const VectorType Y,
                                                          const MatrixType tmp) {
    return invoke(A, X, Y, tmp);
  }
};

/// \brief Team Batched GESV:
///
/// Solve A_l x_l = b_l for all l = 0, ..., N
/// using a batched LU decomposition, 2 batched triangular solves, and a batched
/// static pivoting.
///
/// \tparam MatrixType: Input type for the matrix, needs to be a 2D view
/// \tparam VectorType: Input type for the right-hand side and the solution,
/// needs to be a 1D view
///
/// \param member [in]: TeamPolicy member
/// \param A [in]: matrix, a rank 2 view
/// \param X [out]: solution, a rank 1 view
/// \param B [in]: right-hand side, a rank 1 view
///
/// Two versions are available (those are chosen based on ArgAlgo):
///
///   1. NoPivoting: the solver does not use a pivoting strategy,
///   2. StaticPivoting: the solver uses a static pivoting strategy that relies
///   on using
///      maximal absolute value of row and column to choose pivots and apply
///      them before calling the LU decomposition. Known limitation: the
///      currently implemented strategy would not work with some matrices such
///      as [[2, 1], [1, 0]], when this is the case, the Gesv (if used with
///      pivoting), will return 1 and print an error message.
///
/// A nested parallel_for with TeamThreadRange is used.
///

template <typename MemberType, typename ArgAlgo>
struct TeamGesv {
  template <typename MatrixType, typename VectorType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const MatrixType A, const VectorType X,
                                           const VectorType Y);
};

/// \brief Team Vector Batched GESV:
///
/// Solve A_l x_l = b_l for all l = 0, ..., N
/// using a batched LU decomposition, 2 batched triangular solves, and a batched
/// static pivoting.
///
/// \tparam MatrixType: Input type for the matrix, needs to be a 2D view
/// \tparam VectorType: Input type for the right-hand side and the solution,
/// needs to be a 1D view
///
/// \param member [in]: TeamPolicy member
/// \param A [in]: matrix, a rank 2 view
/// \param X [out]: solution, a rank 1 view
/// \param B [in]: right-hand side, a rank 1 view
///
/// Two versions are available (those are chosen based on ArgAlgo):
///
///   1. NoPivoting: the solver does not use a pivoting strategy,
///   2. StaticPivoting: the solver uses a static pivoting strategy that relies
///   on using
///      maximal absolute value of row and column to choose pivots and apply
///      them before calling the LU decomposition. Known limitation: the
///      currently implemented strategy would not work with some matrices such
///      as [[2, 1], [1, 0]], when this is the case, the Gesv (if used with
///      pivoting), will return 1 and print an error message.
///
///   Two nested parallel_for with both TeamVectorRange and ThreadVectorRange
///   (or one with TeamVectorRange) are used inside.
///

template <typename MemberType, typename ArgAlgo>
struct TeamVectorGesv {
  template <typename MatrixType, typename VectorType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const MatrixType A, const VectorType X,
                                           const VectorType Y);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Gesv_Impl.hpp"

#endif
