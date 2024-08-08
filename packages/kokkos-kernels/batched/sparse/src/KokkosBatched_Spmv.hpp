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
#ifndef __KOKKOSBATCHED_SPMV_HPP__
#define __KOKKOSBATCHED_SPMV_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

/// \brief Serial Batched SPMV:
///   y_l <- alpha_l * A_l * x_l + beta_l * y_l for all l = 1, ..., N
/// where:
///   * N is the number of matrices,
///   * A_1, ..., A_N are N sparse matrices which share the same sparsity
///   pattern,
///   * x_1, ..., x_N are the N input vectors,
///   * y_1, ..., y_N are the N output vectors,
///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N,
///   * beta_1, ..., beta_N are N scaling factors for y_1, ..., y_N.
///
/// \tparam ValuesViewType: Input type for the values of the batched crs matrix,
/// needs to be a 2D view \tparam IntView: Input type for row offset array and
/// column-index array, needs to be a 1D view \tparam xViewType: Input type for
/// X, needs to be a 2D view \tparam yViewType: Input type for Y, needs to be a
/// 2D view \tparam alphaViewType: Input type for alpha, needs to be a 1D view
/// \tparam betaViewType: Input type for beta, needs to be a 1D view
/// \tparam dobeta: Int which sepcifies if beta_l * y_l is used or not (if
/// dobeta == 0, beta_l * y_l is not added to the result of alpha_l * A_l * x_l)
///
/// \param alpha [in]: input coefficient for X, a rank 1 view
/// \param values [in]: values of the batched crs matrix, a rank 2 view
/// \param row_ptr [in]: row offset array of the batched crs matrix, a rank 1
/// view \param colIndices [in]: column-index array of the batched crs matrix, a
/// rank 1 view \param X [in]: Input vector X, a rank 2 view \param beta [in]:
/// input coefficient for Y (if dobeta != 0), a rank 1 view \param Y [in/out]:
/// Output vector Y, a rank 2 view
///
/// The matrices are represented using a Compressed Row Storage (CRS) format and
/// the shared sparsity pattern is reused from one matrix to the others.
///
/// Concretely, instead of providing an array of N matrices to the batched SPMV
/// kernel, the user provides one row offset array (1D view), one column-index
/// array (1D view), and one value array (2D view, one dimension for the
/// non-zero indices and one for the matrix indices).
///
/// No nested parallel_for is used inside of the function.
///

template <typename ArgTrans = Trans::NoTranspose>
struct SerialSpmv {
  template <typename ValuesViewType, typename IntView, typename xViewType, typename yViewType, typename alphaViewType,
            typename betaViewType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(const alphaViewType &alpha, const ValuesViewType &values,
                                           const IntView &row_ptr, const IntView &colIndices, const xViewType &x,
                                           const betaViewType &beta, const yViewType &Y);

  template <typename ValuesViewType, typename IntView, typename xViewType, typename yViewType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type &alpha,
      const ValuesViewType &values, const IntView &row_ptr, const IntView &colIndices, const xViewType &X,
      const typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type &beta,
      const yViewType &Y);
};

/// \brief Team Batched SPMV:
///   y_l <- alpha_l * A_l * x_l + beta_l * y_l for all l = 1, ..., N
/// where:
///   * N is the number of matrices,
///   * A_1, ..., A_N are N sparse matrices which share the same sparsity
///   pattern,
///   * x_1, ..., x_N are the N input vectors,
///   * y_1, ..., y_N are the N output vectors,
///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N,
///   * beta_1, ..., beta_N are N scaling factors for y_1, ..., y_N.
///
/// \tparam ValuesViewType: Input type for the values of the batched crs matrix,
/// needs to be a 2D view \tparam IntView: Input type for row offset array and
/// column-index array, needs to be a 1D view \tparam xViewType: Input type for
/// X, needs to be a 2D view \tparam yViewType: Input type for Y, needs to be a
/// 2D view \tparam alphaViewType: Input type for alpha, needs to be a 1D view
/// \tparam betaViewType: Input type for beta, needs to be a 1D view
/// \tparam dobeta: Int which sepcifies if beta_l * y_l is used or not (if
/// dobeta == 0, beta_l * y_l is not added to the result of alpha_l * A_l * x_l)
///
/// \param member [in]: TeamPolicy member
/// \param alpha [in]: input coefficient for X, a rank 1 view
/// \param values [in]: values of the batched crs matrix, a rank 2 view
/// \param row_ptr [in]: row offset array of the batched crs matrix, a rank 1
/// view \param colIndices [in]: column-index array of the batched crs matrix, a
/// rank 1 view \param X [in]: Input vector X, a rank 2 view \param beta [in]:
/// input coefficient for Y (if dobeta != 0), a rank 1 view \param Y [in/out]:
/// Output vector Y, a rank 2 view
///
/// The matrices are represented using a Compressed Row Storage (CRS) format and
/// the shared sparsity pattern is reused from one matrix to the others.
///
/// Concretely, instead of providing an array of N matrices to the batched SPMV
/// kernel, the user provides one row offset array (1D view), one column-index
/// array (1D view), and one value array (2D view, one dimension for the
/// non-zero indices and one for the matrix indices).
///
/// A nested parallel_for with TeamThreadRange is used.
///

template <typename MemberType, typename ArgTrans = Trans::NoTranspose>
struct TeamSpmv {
  template <typename ValuesViewType, typename IntView, typename xViewType, typename yViewType, typename alphaViewType,
            typename betaViewType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const alphaViewType &alpha,
                                           const ValuesViewType &values, const IntView &row_ptr,
                                           const IntView &colIndices, const xViewType &x, const betaViewType &beta,
                                           const yViewType &y);

  template <typename ValuesViewType, typename IntView, typename xViewType, typename yViewType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member,
      const typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type &alpha,
      const ValuesViewType &values, const IntView &row_ptr, const IntView &colIndices, const xViewType &x,
      const typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type &beta,
      const yViewType &y);
};

/// \brief TeamVector Batched SPMV:
///   y_l <- alpha_l * A_l * x_l + beta_l * y_l for all l = 1, ..., N
/// where:
///   * N is the number of matrices,
///   * A_1, ..., A_N are N sparse matrices which share the same sparsity
///   pattern,
///   * x_1, ..., x_N are the N input vectors,
///   * y_1, ..., y_N are the N output vectors,
///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N,
///   * beta_1, ..., beta_N are N scaling factors for y_1, ..., y_N.
///
/// \tparam ValuesViewType: Input type for the values of the batched crs matrix,
/// needs to be a 2D view \tparam IntView: Input type for row offset array and
/// column-index array, needs to be a 1D view \tparam xViewType: Input type for
/// X, needs to be a 2D view \tparam yViewType: Input type for Y, needs to be a
/// 2D view \tparam alphaViewType: Input type for alpha, needs to be a 1D view
/// \tparam betaViewType: Input type for beta, needs to be a 1D view
/// \tparam dobeta: Int which sepcifies if beta_l * y_l is used or not (if
/// dobeta == 0, beta_l * y_l is not added to the result of alpha_l * A_l * x_l)
///
/// \param member [in]: TeamPolicy member
/// \param alpha [in]: input coefficient for X, a rank 1 view
/// \param values [in]: values of the batched crs matrix, a rank 2 view
/// \param row_ptr [in]: row offset array of the batched crs matrix, a rank 1
/// view \param colIndices [in]: column-index array of the batched crs matrix, a
/// rank 1 view \param X [in]: Input vector X, a rank 2 view \param beta [in]:
/// input coefficient for Y (if dobeta != 0), a rank 1 view \param Y [in/out]:
/// Output vector Y, a rank 2 view
///
/// The matrices are represented using a Compressed Row Storage (CRS) format and
/// the shared sparsity pattern is reused from one matrix to the others.
///
/// Concretely, instead of providing an array of N matrices to the batched SPMV
/// kernel, the user provides one row offset array (1D view), one column-index
/// array (1D view), and one value array (2D view, one dimension for the
/// non-zero indices and one for the matrix indices).
///
/// Two nested parallel_for with both TeamThreadRange and ThreadVectorRange
/// (or one with TeamVectorRange) are used inside.
///

template <typename MemberType, typename ArgTrans = Trans::NoTranspose, unsigned N_team = 1>
struct TeamVectorSpmv {
  template <typename ValuesViewType, typename IntView, typename xViewType, typename yViewType, typename alphaViewType,
            typename betaViewType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const alphaViewType &alpha,
                                           const ValuesViewType &values, const IntView &row_ptr,
                                           const IntView &colIndices, const xViewType &x, const betaViewType &beta,
                                           const yViewType &y);

  template <typename ValuesViewType, typename IntView, typename xViewType, typename yViewType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member,
      const typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type &alpha,
      const ValuesViewType &values, const IntView &row_ptr, const IntView &colIndices, const xViewType &x,
      const typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type &beta,
      const yViewType &y);
};

/// \brief Batched SPMV: Selective Interface
///   y_l <- alpha_l * A_l * x_l + beta_l * y_l for all l = 1, ..., N
/// where:
///   * N is the number of matrices,
///   * A_1, ..., A_N are N sparse matrices which share the same sparsity
///   pattern,
///   * x_1, ..., x_N are the N input vectors,
///   * y_1, ..., y_N are the N output vectors,
///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N,
///   * beta_1, ..., beta_N are N scaling factors for y_1, ..., y_N.
///
/// \tparam ValuesViewType: Input type for the values of the batched crs matrix,
/// needs to be a 2D view \tparam IntView: Input type for row offset array and
/// column-index array, needs to be a 1D view \tparam xViewType: Input type for
/// X, needs to be a 2D view \tparam yViewType: Input type for Y, needs to be a
/// 2D view \tparam alphaViewType: Input type for alpha, needs to be a 1D view
/// \tparam betaViewType: Input type for beta, needs to be a 1D view
/// \tparam dobeta: Int which sepcifies if beta_l * y_l is used or not (if
/// dobeta == 0, beta_l * y_l is not added to the result of alpha_l * A_l * x_l)
///
/// \param member [in]: TeamPolicy member
/// \param alpha [in]: input coefficient for X, a rank 1 view
/// \param values [in]: values of the batched crs matrix, a rank 2 view
/// \param row_ptr [in]: row offset array of the batched crs matrix, a rank 1
/// view \param colIndices [in]: column-index array of the batched crs matrix, a
/// rank 1 view \param X [in]: Input vector X, a rank 2 view \param beta [in]:
/// input coefficient for Y (if dobeta != 0), a rank 1 view \param Y [in/out]:
/// Output vector Y, a rank 2 view

template <typename MemberType, typename ArgTrans, typename ArgMode>
struct Spmv {
  template <typename ValuesViewType, typename IntView, typename xViewType, typename yViewType, typename alphaViewType,
            typename betaViewType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const alphaViewType &alpha,
                                           const ValuesViewType &values, const IntView &row_ptr,
                                           const IntView &colIndices, const xViewType &x, const betaViewType &beta,
                                           const yViewType &y) {
    int r_val = 0;
    if (std::is_same<ArgMode, Mode::Serial>::value) {
      r_val =
          SerialSpmv<ArgTrans>::template invoke<ValuesViewType, IntView, xViewType, yViewType, alphaViewType,
                                                betaViewType, dobeta>(alpha, values, row_ptr, colIndices, x, beta, y);
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      r_val = TeamSpmv<MemberType, ArgTrans>::template invoke<ValuesViewType, IntView, xViewType, yViewType,
                                                              alphaViewType, betaViewType, dobeta>(
          member, alpha, values, row_ptr, colIndices, x, beta, y);
    } else if (std::is_same<ArgMode, Mode::TeamVector>::value) {
      r_val = TeamVectorSpmv<MemberType, ArgTrans>::template invoke<ValuesViewType, IntView, xViewType, yViewType,
                                                                    alphaViewType, betaViewType, dobeta>(
          member, alpha, values, row_ptr, colIndices, x, beta, y);
    }
    return r_val;
  }

  template <typename ValuesViewType, typename IntView, typename xViewType, typename yViewType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member,
      const typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type &alpha,
      const ValuesViewType &values, const IntView &row_ptr, const IntView &colIndices, const xViewType &x,
      const typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type &beta,
      const yViewType &y) {
    int r_val = 0;
    if (std::is_same<ArgMode, Mode::Serial>::value) {
      r_val = SerialSpmv<ArgTrans>::template invoke<ValuesViewType, IntView, xViewType, yViewType, dobeta>(
          alpha, values, row_ptr, colIndices, x, beta, y);
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      r_val = TeamSpmv<MemberType, ArgTrans>::template invoke<ValuesViewType, IntView, xViewType, yViewType, dobeta>(
          member, alpha, values, row_ptr, colIndices, x, beta, y);
    } else if (std::is_same<ArgMode, Mode::TeamVector>::value) {
      r_val =
          TeamVectorSpmv<MemberType, ArgTrans>::template invoke<ValuesViewType, IntView, xViewType, yViewType, dobeta>(
              member, alpha, values, row_ptr, colIndices, x, beta, y);
    }
    return r_val;
  }
};
}  // namespace KokkosBatched

#include "KokkosBatched_Spmv_Serial_Impl.hpp"
#include "KokkosBatched_Spmv_Team_Impl.hpp"
#include "KokkosBatched_Spmv_TeamVector_Impl.hpp"
#endif
