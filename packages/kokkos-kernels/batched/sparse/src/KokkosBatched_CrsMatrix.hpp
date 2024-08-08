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
#ifndef __KOKKOSBATCHED_CRSMATRIX_HPP__
#define __KOKKOSBATCHED_CRSMATRIX_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

namespace KokkosBatched {

/// \brief Batched CrsMatrix:
///
/// \tparam ValuesViewType: Input type for the values of the batched crs matrix,
/// needs to be a 2D view \tparam IntView: Input type for row offset array and
/// column-index array, needs to be a 1D view

template <class ValuesViewType, class IntViewType>
class CrsMatrix {
 public:
  using ScalarType    = typename ValuesViewType::non_const_value_type;
  using MagnitudeType = typename Kokkos::ArithTraits<ScalarType>::mag_type;

 private:
  ValuesViewType values;
  IntViewType row_ptr;
  IntViewType colIndices;
  int n_operators;
  int n_rows;
  int n_colums;

 public:
  KOKKOS_INLINE_FUNCTION
  CrsMatrix(const ValuesViewType &_values, const IntViewType &_row_ptr, const IntViewType &_colIndices)
      : values(_values), row_ptr(_row_ptr), colIndices(_colIndices) {
    n_operators = _values.extent(0);
    n_rows      = _row_ptr.extent(0) - 1;
    n_colums    = n_rows;
  }

  KOKKOS_INLINE_FUNCTION
  ~CrsMatrix() {}

  /// \brief apply version that uses constant coefficients alpha and beta
  ///
  ///   y_l <- alpha * A_l * x_l + beta * y_l for all l = 1, ..., N
  /// where:
  ///   * N is the number of matrices,
  ///   * A_1, ..., A_N are N sparse matrices which share the same sparsity
  ///   pattern,
  ///   * x_1, ..., x_N are the N input vectors,
  ///   * y_1, ..., y_N are the N output vectors,
  ///   * alpha is a scaling factor for x_1, ..., x_N,
  ///   * beta is a scaling factor for y_1, ..., y_N.
  ///
  /// \tparam MemberType: Input type for the TeamPolicy member
  /// \tparam XViewType: Input type for X, needs to be a 2D view
  /// \tparam YViewType: Input type for Y, needs to be a 2D view
  /// \tparam ArgTrans: Argument for transpose or notranspose
  /// \tparam ArgMode: Argument for the parallelism used in the apply
  ///
  /// \param member [in]: TeamPolicy member
  /// \param alpha [in]: input coefficient for X (default value 1.)
  /// \param X [in]: Input vector X, a rank 2 view
  /// \param beta [in]: input coefficient for Y (default value 0.)
  /// \param Y [in/out]: Output vector Y, a rank 2 view

  template <typename ArgTrans, typename ArgMode, typename MemberType, typename XViewType, typename YViewType>
  KOKKOS_INLINE_FUNCTION void apply(const MemberType &member, const XViewType &X, const YViewType &Y,
                                    MagnitudeType alpha = Kokkos::ArithTraits<MagnitudeType>::one(),
                                    MagnitudeType beta  = Kokkos::ArithTraits<MagnitudeType>::zero()) const {
    if (beta == Kokkos::ArithTraits<MagnitudeType>::zero()) {
      if (member.team_size() == 1 && n_operators == 8)
        KokkosBatched::TeamVectorSpmv<MemberType, ArgTrans, 8>::template invoke<ValuesViewType, IntViewType, XViewType,
                                                                                YViewType, 0>(
            member, alpha, values, row_ptr, colIndices, X, beta, Y);
      else
        KokkosBatched::TeamVectorSpmv<MemberType, ArgTrans>::template invoke<ValuesViewType, IntViewType, XViewType,
                                                                             YViewType, 0>(
            member, alpha, values, row_ptr, colIndices, X, beta, Y);
    } else {
      if (member.team_size() == 1 && n_operators == 8)
        KokkosBatched::TeamVectorSpmv<MemberType, ArgTrans, 8>::template invoke<ValuesViewType, IntViewType, XViewType,
                                                                                YViewType, 1>(
            member, alpha, values, row_ptr, colIndices, X, beta, Y);
      else
        KokkosBatched::TeamVectorSpmv<MemberType, ArgTrans>::template invoke<ValuesViewType, IntViewType, XViewType,
                                                                             YViewType, 1>(
            member, alpha, values, row_ptr, colIndices, X, beta, Y);
    }
  }

  template <typename ArgTrans, typename XViewType, typename YViewType>
  KOKKOS_INLINE_FUNCTION void apply(const XViewType &X, const YViewType &Y,
                                    MagnitudeType alpha = Kokkos::ArithTraits<MagnitudeType>::one(),
                                    MagnitudeType beta  = Kokkos::ArithTraits<MagnitudeType>::zero()) const {
    if (beta == Kokkos::ArithTraits<MagnitudeType>::zero())
      KokkosBatched::SerialSpmv<ArgTrans>::template invoke<ValuesViewType, IntViewType, XViewType, YViewType, 0>(
          alpha, values, row_ptr, colIndices, X, beta, Y);
    else
      KokkosBatched::SerialSpmv<ArgTrans>::template invoke<ValuesViewType, IntViewType, XViewType, YViewType, 1>(
          alpha, values, row_ptr, colIndices, X, beta, Y);
  }
};

}  // namespace KokkosBatched

#endif