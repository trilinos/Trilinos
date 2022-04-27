//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.4
//       Copyright (2021) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
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
  using ScalarType = typename ValuesViewType::non_const_value_type;
  using MagnitudeType =
      typename Kokkos::Details::ArithTraits<ScalarType>::mag_type;

 private:
  ValuesViewType values;
  IntViewType row_ptr;
  IntViewType colIndices;
  int n_operators;
  int n_rows;
  int n_colums;

 public:
  KOKKOS_INLINE_FUNCTION
  CrsMatrix(const ValuesViewType &_values, const IntViewType &_row_ptr,
            const IntViewType &_colIndices)
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

  template <typename MemberType, typename XViewType, typename YViewType,
            typename ArgTrans, typename ArgMode>
  KOKKOS_INLINE_FUNCTION void apply(
      const MemberType &member, const XViewType &X, const YViewType &Y,
      MagnitudeType alpha = Kokkos::Details::ArithTraits<MagnitudeType>::one(),
      MagnitudeType beta =
          Kokkos::Details::ArithTraits<MagnitudeType>::zero()) const {
    if (beta == 0)
      KokkosBatched::Spmv<MemberType, ArgTrans, ArgMode>::template invoke<
          ValuesViewType, IntViewType, XViewType, YViewType, 0>(
          member, alpha, values, row_ptr, colIndices, X, beta, Y);
    else
      KokkosBatched::Spmv<MemberType, ArgTrans, ArgMode>::template invoke<
          ValuesViewType, IntViewType, XViewType, YViewType, 1>(
          member, alpha, values, row_ptr, colIndices, X, beta, Y);
  }

  /// \brief apply version that uses variable coefficient alpha and no beta
  ///   y_l <- alpha_l * A_l * x_l  for all l = 1, ..., N
  /// where:
  ///   * N is the number of matrices,
  ///   * A_1, ..., A_N are N sparse matrices which share the same sparsity
  ///   pattern,
  ///   * x_1, ..., x_N are the N input vectors,
  ///   * y_1, ..., y_N are the N output vectors,
  ///   * alpha_1, ..., alpha_N are N scaling factors for x_1, ..., x_N.
  ///
  /// \tparam MemberType: Input type for the TeamPolicy member
  /// \tparam XViewType: Input type for X, needs to be a 2D view
  /// \tparam YViewType: Input type for Y, needs to be a 2D view
  /// \tparam ArgTrans: Argument for transpose or notranspose
  /// \tparam ArgMode: Argument for the parallelism used in the apply
  ///
  /// \param member [in]: TeamPolicy member
  /// \param alpha [in]: input coefficient for X, a rank 1 view
  /// \param X [in]: Input vector X, a rank 2 view
  /// \param Y [out]: Output vector Y, a rank 2 view

  template <typename MemberType, typename XViewType, typename YViewType,
            typename NormViewType, typename ArgTrans, typename ArgMode>
  KOKKOS_INLINE_FUNCTION void apply(const MemberType &member,
                                    const XViewType &X, const YViewType &Y,
                                    NormViewType alpha) const {
    KokkosBatched::Spmv<MemberType, ArgTrans, ArgMode>::template invoke<
        ValuesViewType, IntViewType, XViewType, YViewType, NormViewType,
        NormViewType, 0>(member, alpha, values, row_ptr, colIndices, X, alpha,
                         Y);
  }

  /// \brief apply version that uses variable coefficients alpha and beta
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
  /// \tparam MemberType: Input type for the TeamPolicy member
  /// \tparam XViewType: Input type for X, needs to be a 2D view
  /// \tparam YViewType: Input type for Y, needs to be a 2D view
  /// \tparam NormViewType: Input type for alpha and beta, needs to be a 1D view
  /// \tparam ArgTrans: Argument for transpose or notranspose
  /// \tparam ArgMode: Argument for the parallelism used in the apply
  ///
  /// \param member [in]: TeamPolicy member
  /// \param alpha [in]: input coefficient for X, a rank 1 view
  /// \param X [in]: Input vector X, a rank 2 view
  /// \param beta [in]: input coefficient for Y, a rank 1 view
  /// \param Y [in/out]: Output vector Y, a rank 2 view

  template <typename MemberType, typename XViewType, typename YViewType,
            typename NormViewType, typename ArgTrans, typename ArgMode>
  KOKKOS_INLINE_FUNCTION void apply(const MemberType &member,
                                    const XViewType &X, const YViewType &Y,
                                    const NormViewType &alpha,
                                    const NormViewType &beta) const {
    KokkosBatched::Spmv<MemberType, ArgTrans, ArgMode>::template invoke<
        ValuesViewType, IntViewType, XViewType, YViewType, NormViewType,
        NormViewType, 1>(member, alpha, values, row_ptr, colIndices, X, beta,
                         Y);
  }
};

}  // namespace KokkosBatched

#endif